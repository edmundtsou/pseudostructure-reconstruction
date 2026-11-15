"""Section A of 00_bin2cell_annotation_aggregation.ipynb -- Bin2Cell reconstruction of one or more sections.

Written by the notebook's %%writefile cell and executed in the dedicated `bin2cell` conda env
(TensorFlow 2.15 + StarDist 0.9.1 + bin2cell 0.3.4, numpy<2), which cannot share a kernel with the
rest of the notebook. The code is the reference notebook 01_bin_to_cell_segmentation (calibrated on
Ctrl_1A2 / Ctrl_1A4) with the parameters unchanged; the only additions are the per-nucleus
centroid (obs['nucleus_x_fullres_px' / 'nucleus_y_fullres_px']) used to place cells in tubules,
and one registration check that validates its coordinate transform.

Copied into this repository on 2026-11-15 from the mouse bin2cell_results bundle, unchanged except for the two human
entries in SAMPLE_INPUTS. Here the env has opencv 4.10.0.84 instead of 4.13.0: the 4.13 Linux wheel needs numpy>=2,
which TensorFlow 2.15 cannot load.

Usage: PYTHONNOUSERSITE=1 <bin2cell python> bin2cell_run.py --visium-dir DIR --out DIR --samples A,B [--force]
"""
import argparse
import json
import os
import platform
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "2")
if os.environ.get("PYTHONNOUSERSITE") != "1":
    raise RuntimeError("run with PYTHONNOUSERSITE=1 (a user-site scanpy can shadow the env's copy)")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
from PIL import Image

import bin2cell as b2c

Image.MAX_IMAGE_PIXELS = None     # the full-resolution H&E exceeds PIL's default bomb limit
pd.set_option("display.max_colwidth", 120); pd.set_option("display.width", 220)

_ARGS = argparse.ArgumentParser()
_ARGS.add_argument("--visium-dir", required=True, type=Path)
_ARGS.add_argument("--out", required=True, type=Path)
_ARGS.add_argument("--samples", required=True)
_ARGS.add_argument("--force", action="store_true")
ARGS = _ARGS.parse_args()
B2C_DIR = ARGS.out
FORCE = ARGS.force

VISIUM_DIR = ARGS.visium_dir
# Per-sample INPUTS only. Every Bin2Cell parameter below is shared by all four samples. All four
# sections: Space Ranger 3.1.1, Visium HD v1, Mouse Transcriptome Probe Set v2.0, mm10-2020-A. The
# hires H&E sits in square_002um/spatial for the controls and in outs/spatial for the IR sections.
_COMMON = dict(chemistry="Visium HD v1", spaceranger="3.1.1",
               assay="probe-based (Visium Mouse Transcriptome Probe Set v2.0)", genome="mm10 (refdata-gex-mm10-2020-A)")
SAMPLE_INPUTS: dict[str, dict] = {
    "Ctrl_1A2": dict(condition="control", outs=VISIUM_DIR / "Ctrl_1A2" / "outs",
                     source_image=VISIUM_DIR / "Ctrl_1A2" / "outs" / "spatial" / "Ctrl_1A2.tif",
                     image_dir=VISIUM_DIR / "Ctrl_1A2" / "outs" / "binned_outputs" / "square_002um" / "spatial", **_COMMON),
    "Ctrl_1A4": dict(condition="control", outs=VISIUM_DIR / "Ctrl_1A4" / "outs",
                     source_image=VISIUM_DIR / "Ctrl_1A4" / "outs" / "spatial" / "Ctrl_1A4.tif",
                     image_dir=VISIUM_DIR / "Ctrl_1A4" / "outs" / "binned_outputs" / "square_002um" / "spatial", **_COMMON),
    "IR_2A2": dict(condition="ischemia-reperfusion", outs=VISIUM_DIR / "IR_2A2" / "outs",
                   source_image=VISIUM_DIR / "IR_2A2" / "outs" / "spatial" / "IR_2A2.tif",
                   image_dir=VISIUM_DIR / "IR_2A2" / "outs" / "spatial", **_COMMON),
    "IR_2A4": dict(condition="ischemia-reperfusion", outs=VISIUM_DIR / "IR_2A4" / "outs",
                   source_image=VISIUM_DIR / "IR_2A4" / "outs" / "spatial" / "IR_2A4.tif",
                   image_dir=VISIUM_DIR / "IR_2A4" / "outs" / "spatial", **_COMMON),
}
# Human healthy cortex, one donor (both slices are cortex; "MED" is the slice name only). Space Ranger 4.0.1,
# Human Transcriptome Probe Set v2.1.0, GRCh38-2024-A. The full-res H&E sits beside outs/ as <sample>/<sample>.tif.
_HUMAN = dict(chemistry="Visium HD v1", spaceranger="4.0.1",
              assay="probe-based (Visium Human Transcriptome Probe Set v2.1.0)", genome="GRCh38-2024-A")
for _s in ("HUK1_COR_1", "HUK1_MED_1"):
    SAMPLE_INPUTS[_s] = dict(condition="human control", outs=VISIUM_DIR / _s / "outs",
                             source_image=VISIUM_DIR / _s / f"{_s}.tif",
                             image_dir=VISIUM_DIR / _s / "outs" / "spatial", **_HUMAN)
SAMPLES = [s for s in ARGS.samples.split(",") if s]

# ---- Bin2Cell parameters (identical for every sample; recorded to <sample>/params.json) ----
MPP = 0.5                          # microns per pixel of the H&E render StarDist segments (bin2cell demo value)
CROP_BUFFER = 150                  # full-res px of margin around the capture area in the render
DESTRIPE_QUANTILE = 0.99           # b2c.destripe default
STARDIST_HE_MODEL = "2D_versatile_he"
STARDIST_HE_PROB_THRESH = 0.1     # bin2cell demo value; the model default 0.69 under-calls on H&E (ch increased from 0.01 to .5) and back to .25

STARDIST_BLOCK_SIZE, STARDIST_MIN_OVERLAP, STARDIST_CONTEXT = 4096, 128, 128
MAX_BIN_DISTANCE = 4               # expand each nucleus by <= 2 bins = 4 um (b2c demo value; ch increased from 2 to 4)
EXPAND_K = 5                       # candidate labelled bins per unlabelled bin (b2c default) (ch increased from 4 to 5)
GEX_SIGMA = 6                      # gaussian sigma (px at MPP) on the total-count image (b2c demo value; ch increased from 4 to 6)
STARDIST_GEX_MODEL = "2D_versatile_fluo"
STARDIST_GEX_PROB_THRESH, STARDIST_GEX_NMS_THRESH = 0.5, 0.5 # bin2cell demo values; the model defaults 0.5 and 0.4 under-call on GEX (ch increased from 0.01 to .5)
BIN_SIDE_UM = 2.0
RANDOM_SEED = 0
MAX_FRAC_BINS_OUTSIDE_IMAGE = 0.001   # the capture area can overhang the image by a few rows of bins
TISSUE_PIXEL_MAX_MEAN_RGB = 225.0     # H&E background is near white
# Preview of the notebook-02 gates (REVIEW ONLY; pipeline-default floors, not the V10 150/100 floors)
PREVIEW_MIN_COUNTS, PREVIEW_MIN_GENES, PREVIEW_MAD_NMADS = 50, 25, 3.0
PREVIEW_MIN_AREA_UM2, PREVIEW_MAX_AREA_UM2 = 10.0, 600.0
N_ZOOM_REGIONS, ZOOM_SIDE_UM, OVERLAY_DPI = 6, 300.0, 250
# Added for the tubule assignment: the 2 um bin nearest a nucleus centroid must lie inside that
# nucleus (labels_he == its id) for at least this fraction of nucleus-anchored cells.
MIN_FRAC_NUCLEUS_CENTROID_IN_OWN_NUCLEUS = 0.90

PARAMS = {k: v for k, v in dict(globals()).items()
          if k.isupper() and not k.startswith("_") and isinstance(v, (int, float, str, tuple))}
CHECK_NAME_ALIASES = {   # registration checks renamed 2026-09-11, after Ctrl_1A2 was first run
    "all bins inside the full-res image": "bins inside the full-res image (capture area may overhang the image edge)",
    "all bins inside the cropped 0.5 mpp render": "bins inside the cropped render",
}
print("B2C_DIR:", B2C_DIR, "| samples:", SAMPLES, "| FORCE:", FORCE, flush=True)

# ---- shared utilities (duplicated in every notebook on purpose: no local imports) ----
_T0 = time.time()


def log(msg: str = "") -> None:
    """Print a progress line prefixed with the seconds elapsed since the kernel started.

    Parameters
    ----------
    msg : str
        Text to print. Flushed immediately so long cells show progress.
    """
    print(f"[{time.time() - _T0:7.1f}s] {msg}", flush=True)


def banner(msg: str) -> None:
    """Print a section banner (a line of '=' above and below ``msg``)."""
    print("\n" + "=" * 86 + f"\n{msg}\n" + "=" * 86, flush=True)


def warn_banner(msg: str) -> None:
    """Print a warning that cannot be missed in a long cell output."""
    print("\n" + "!" * 86 + f"\nWARNING: {msg}\n" + "!" * 86 + "\n", flush=True)


def show_png(path: Path, width: int = 1000) -> None:
    """Display a PNG inline when running under Jupyter; print its path otherwise.

    Parameters
    ----------
    path : Path
        Image file to show.
    width : int
        Display width in pixels (notebook only).
    """
    try:
        from IPython import get_ipython
        from IPython.display import Image as _IPImage, display
        if get_ipython() is None:
            raise RuntimeError
        display(_IPImage(filename=str(path), width=width))
    except Exception:  # noqa: BLE001 - plain-python execution (smoke test)
        print(f"  [figure] {path}", flush=True)


@dataclass
class SampleIO:
    """Every path one Bin2Cell run reads or writes.

    Attributes
    ----------
    sample : str
        Sample key in ``SAMPLE_INPUTS``.
    condition : str
        Human-readable condition, copied into ``obs['condition']``.
    bin_dir : Path
        Space Ranger ``binned_outputs/square_002um`` directory (input, read-only).
    source_image : Path
        Full-resolution H&E given to Space Ranger (input, read-only).
    image_dir : Path
        Directory holding ``tissue_hires_image.png`` (input, read-only).
    out : Path
        Output directory for this sample; ``out / 'qc'`` holds tables and figures.
    """
    sample: str
    condition: str
    bin_dir: Path
    source_image: Path
    image_dir: Path
    out: Path

    @property
    def he_render(self) -> Path:
        return self.out / "he_0.5mpp_crop.tiff"

    @property
    def gex_render(self) -> Path:
        return self.out / "gex_0.5mpp.tiff"

    @property
    def npz_he(self) -> Path:
        return self.out / "labels_he.npz"

    @property
    def npz_gex(self) -> Path:
        return self.out / "labels_gex.npz"

    @property
    def cells_h5ad(self) -> Path:
        return self.out / f"{self.sample}_bin2cell_cells_unfiltered.h5ad"


def sample_io(sample: str, out_base: Path) -> SampleIO:
    """Build the ``SampleIO`` of one sample from ``SAMPLE_INPUTS``.

    Parameters
    ----------
    sample : str
        Key of ``SAMPLE_INPUTS``.
    out_base : Path
        Root output directory; the sample writes into ``out_base / sample``.
    """
    cfg = SAMPLE_INPUTS[sample]
    return SampleIO(sample=sample, condition=cfg["condition"],
                    bin_dir=cfg["outs"] / "binned_outputs" / "square_002um",
                    source_image=cfg["source_image"], image_dir=cfg["image_dir"],
                    out=Path(out_base) / sample)


def record_environment(out: Path) -> dict:
    """Write ``environment_versions.json`` and ``pip_freeze.txt`` for this run.

    Parameters
    ----------
    out : Path
        Sample output directory.

    Returns
    -------
    dict
        Package versions and the GPUs TensorFlow can see.
    """
    import anndata, bin2cell, cv2, scanpy, skimage, stardist, tensorflow as tf  # noqa: E401
    gpus = tf.config.list_physical_devices("GPU")
    info = dict(
        python=platform.python_version(), executable=sys.executable,
        bin2cell=bin2cell.__version__, stardist=stardist.__version__,
        tensorflow=tf.__version__, numpy=np.__version__, scanpy=scanpy.__version__,
        anndata=anndata.__version__, opencv=cv2.__version__, scikit_image=skimage.__version__,
        pandas=pd.__version__, scipy=__import__("scipy").__version__,
        gpus=[g.name for g in gpus], gpu_used=bool(gpus),
        stardist_models=dict(he=STARDIST_HE_MODEL, gex=STARDIST_GEX_MODEL,
                             weights_dir=str(Path.home() / ".keras" / "models" / "StarDist2D"),
                             source="StarDist2D.from_pretrained() -> github.com/stardist/stardist-models"),
        PYTHONNOUSERSITE=os.environ.get("PYTHONNOUSERSITE", ""),
        hostname=platform.node(), date=time.strftime("%Y-%m-%d %H:%M:%S"),
    )
    (out / "environment_versions.json").write_text(json.dumps(info, indent=2))
    try:
        req = subprocess.run([sys.executable, "-m", "pip", "freeze"], capture_output=True,
                             text=True, env={**os.environ, "PYTHONNOUSERSITE": "1"}).stdout
        (out / "pip_freeze.txt").write_text(req)
    except Exception as exc:  # noqa: BLE001
        log(f"  (pip freeze failed: {exc})")
    return info


def read_bins(io: SampleIO):
    """Read the 2 um filtered bin matrix with coordinates in FULL-RES pixels (x = col, y = row).

    Parameters
    ----------
    io : SampleIO
        Paths of the sample.

    Returns
    -------
    AnnData
        One row per 2 um bin; ``obs['n_counts']`` / ``obs['n_genes']`` added. Zero-count bins
        are kept: they still occupy area and receive labels.
    """
    adata = b2c.read_visium(str(io.bin_dir), source_image_path=str(io.source_image),
                            spaceranger_image_path=str(io.image_dir))
    adata.var_names_make_unique()
    lib = list(adata.uns["spatial"])[0]
    sf = adata.uns["spatial"][lib]["scalefactors"]
    log(f"  {adata.n_obs:,} bins x {adata.n_vars:,} genes; library '{lib}'; "
        f"scalefactors {json.dumps({k: round(v, 6) for k, v in sf.items()})}")
    n_in = int(adata.obs["in_tissue"].sum())
    assert n_in == adata.n_obs, f"{adata.n_obs - n_in} bins in the filtered matrix are not in_tissue"
    adata.obs["n_counts"] = np.asarray(adata.X.sum(axis=1)).ravel()
    adata.obs["n_genes"] = np.asarray((adata.X > 0).sum(axis=1)).ravel()
    log(f"  bins with zero counts: {int((adata.obs.n_counts == 0).sum()):,} (kept)")
    return adata


def safe_destripe(adata, quantile: float = DESTRIPE_QUANTILE) -> None:
    """``b2c.destripe()`` with a NaN-safe count adjustment for zero-count bins.

    The library's ``destripe_counts()`` divides adjusted by raw counts per bin, which is 0/0
    for a bin with no UMIs. Those bins are kept (they carry area), so the factor is computed
    by the library and applied here with 0 in place of NaN.

    Parameters
    ----------
    adata : AnnData
        Bin object from ``read_bins``; modified in place. Raw counts move to
        ``layers['counts_raw']`` and ``X`` becomes the destriped matrix.
    quantile : float
        Per-row / per-column quantile used by ``b2c.destripe``.
    """
    b2c.destripe(adata, quantile=quantile, adjust_counts=False)
    raw = adata.obs["n_counts"].to_numpy(float)
    adj = adata.obs["n_counts_adjusted"].to_numpy(float)
    factor = np.where(raw > 0, adj / np.where(raw > 0, raw, 1.0), 0.0)
    adata.obs["n_counts_adjusted"] = np.where(raw > 0, adj, 0.0)
    adata.layers["counts_raw"] = adata.X.copy()
    adata.X = sp.diags(factor).dot(adata.X).tocsr()
    log(f"  destriped: quantile {quantile}; median per-bin factor {np.median(factor[raw > 0]):.3f}; "
        f"bins with factor==0 (no counts) {int((factor == 0).sum()):,}")


def run_stardist(image_path: Path, npz_path: Path, force: bool, **kw) -> None:
    """Run ``b2c.stardist`` unless its label file already exists.

    Parameters
    ----------
    image_path : Path
        Image to segment (H&E render or smoothed total-count image).
    npz_path : Path
        Sparse label matrix to write.
    force : bool
        Recompute even if ``npz_path`` exists.
    **kw
        Passed to ``b2c.stardist`` (``stardist_model``, ``prob_thresh``, ``nms_thresh``).
    """
    if npz_path.exists() and not force:
        log(f"  reusing {npz_path.name}")
        return
    t0 = time.time()
    b2c.stardist(image_path=str(image_path), labels_npz_path=str(npz_path),
                 block_size=STARDIST_BLOCK_SIZE, min_overlap=STARDIST_MIN_OVERLAP,
                 context=STARDIST_CONTEXT, **kw)
    log(f"  StarDist {kw.get('stardist_model')} -> {npz_path.name} in {time.time() - t0:.0f}s")


def label_image_2um(adata, key: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Dense label image on the 2 um array grid (rows x cols), plus the row/col index arrays."""
    r = adata.obs["array_row"].to_numpy(int); c = adata.obs["array_col"].to_numpy(int)
    img = np.zeros((r.max() + 1, c.max() + 1), dtype=np.int64)
    img[r, c] = adata.obs[key].to_numpy(int)
    return img, r, c


def footprint_morphometrics(adata, key: str) -> pd.DataFrame:
    """Area, perimeter and circularity of each label's bin footprint on the 2 um grid.

    The perimeter is the exact CRACK perimeter: the number of bin edges whose 4-neighbour
    carries a different label (or lies outside the grid), times the bin side. skimage's
    ``perimeter`` under-estimates small footprints (circularity > 1 for most cells here), so
    it is not used. A full n x n square of bins has circularity pi/4 = 0.785, the ceiling.

    Parameters
    ----------
    adata : AnnData
        Bin object carrying the integer label column ``obs[key]`` (0 = unassigned).
    key : str
        Label column, e.g. ``'labels_joint'``.

    Returns
    -------
    DataFrame
        One row per object: ``object_id``, ``cell_area_um2``, ``cell_perimeter_um``,
        ``cell_circularity``, ``eccentricity`` and the footprint bounding-box size.
    """
    from skimage.measure import regionprops_table
    img, r, c = label_image_2um(adata, key)
    props = regionprops_table(img, properties=("label", "area", "bbox", "eccentricity"))
    df = pd.DataFrame(props).rename(columns={"label": "object_id"})
    lab = img[r, c]
    exposed = np.zeros(len(r), dtype=np.int64)
    for dr, dc in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        rr, cc = r + dr, c + dc
        inside = (rr >= 0) & (rr < img.shape[0]) & (cc >= 0) & (cc < img.shape[1])
        nbr = np.full(len(r), -1, dtype=np.int64)
        nbr[inside] = img[rr[inside], cc[inside]]
        exposed += (nbr != lab)
    edges = pd.Series(exposed[lab > 0]).groupby(lab[lab > 0]).sum()
    df["cell_area_um2"] = df["area"] * BIN_SIDE_UM ** 2
    df["cell_perimeter_um"] = edges.reindex(df["object_id"]).to_numpy() * BIN_SIDE_UM
    with np.errstate(divide="ignore", invalid="ignore"):
        circ = 4.0 * np.pi * df["cell_area_um2"] / np.square(df["cell_perimeter_um"])
    df["cell_circularity"] = np.clip(circ.fillna(0.0), 0, 1)
    df["footprint_bbox_rows"] = df["bbox-2"] - df["bbox-0"]
    df["footprint_bbox_cols"] = df["bbox-3"] - df["bbox-1"]
    return df.drop(columns=["area", "bbox-0", "bbox-1", "bbox-2", "bbox-3"])


def nucleus_pixel_areas(npz_path: Path, mpp: float) -> pd.Series:
    """Nucleus area in um^2 from the StarDist mask pixel counts at ``mpp`` microns per pixel."""
    lab = sp.load_npz(npz_path)
    counts = np.bincount(lab.data.astype(np.int64))
    ids = np.flatnonzero(counts)
    return pd.Series(counts[ids].astype(float) * mpp ** 2, index=ids, name="nucleus_area_um2")


def nucleus_centroids_fullres(npz_path: Path, scalef: float, crop_origin) -> pd.DataFrame:
    """Centroid of every StarDist nucleus in FULL-RES H&E pixels (the frame of obsm['spatial'] and the GeoJSONs).

    The label matrix indexes the cropped render (row = y, col = x at ``MPP``); a render pixel maps
    back as ``render / scalef + crop_origin`` -- the inverse of the transform
    ``registration_checks`` uses to put bins on the render.

    Parameters
    ----------
    npz_path : Path
        ``labels_he.npz`` (sparse label matrix of the render).
    scalef : float
        Full-res pixel -> render pixel scale factor.
    crop_origin : tuple[float, float]
        Full-res (x, y) of the render's top-left corner.

    Returns
    -------
    DataFrame
        Index = nucleus label id; ``nucleus_x_fullres_px``, ``nucleus_y_fullres_px``.
    """
    lab = sp.load_npz(npz_path).tocoo()
    ids = lab.data.astype(np.int64)
    n = np.bincount(ids)
    sx = np.bincount(ids, weights=lab.col.astype(float))
    sy = np.bincount(ids, weights=lab.row.astype(float))
    keep = np.flatnonzero(n)
    return pd.DataFrame({"nucleus_x_fullres_px": sx[keep] / n[keep] / scalef + crop_origin[0],
                         "nucleus_y_fullres_px": sy[keep] / n[keep] / scalef + crop_origin[1]}, index=keep)


def build_cells(adata, io: SampleIO, scalef: float, crop_origin):
    """Sum bins into cells on raw and on destriped counts, and attach morphometrics.

    Parameters
    ----------
    adata : AnnData
        Bin object with ``labels_joint`` / ``labels_joint_source`` / ``labels_he`` and the
        raw counts in ``layers['counts_raw']``.
    io : SampleIO
        Paths of the sample (the H&E label file gives nucleus pixel areas and centroids).
    scalef, crop_origin
        Render scale factor and full-res crop origin (for the nucleus centroids).

    Returns
    -------
    AnnData
        One row per reconstructed object. ``X`` and ``layers['counts']`` are RAW INTEGER bin
        sums; ``layers['counts_destriped']`` is kept for reference only.
        ``obs['has_nucleus']`` separates nucleus-anchored cells from GEX-only salvaged objects.
    """
    keys = ["spatial", "spatial_cropped_150_buffer"]
    destriped_X = adata.X
    adata.X = adata.layers["counts_raw"]          # raw integer counts first (the primary matrix)
    cdata = b2c.bin_to_cell(adata, labels_key="labels_joint", spatial_keys=keys)
    cdata.layers["counts"] = sp.csr_matrix(cdata.X, dtype=np.float32)
    assert np.allclose(cdata.layers["counts"].data, np.round(cdata.layers["counts"].data)), \
        "raw-count aggregation produced non-integers"
    adata.X = destriped_X
    cdata_d = b2c.bin_to_cell(adata, labels_key="labels_joint", spatial_keys=keys)
    assert (cdata_d.obs_names == cdata.obs_names).all()
    cdata.layers["counts_destriped"] = sp.csr_matrix(cdata_d.X, dtype=np.float32)
    cdata.X = cdata.layers["counts"].copy()
    del cdata_d
    log(f"  {cdata.n_obs:,} cells x {cdata.n_vars:,} genes; median bins/cell {np.median(cdata.obs.bin_count):.0f}")

    src = cdata.obs["labels_joint_source"].astype(str)
    cdata.obs["has_nucleus"] = (src == "primary").to_numpy()
    n_prim = int(cdata.obs.has_nucleus.sum())
    log(f"  primary (H&E nucleus-anchored) {n_prim:,}; secondary (GEX-only) {cdata.n_obs - n_prim:,}")

    nuc_bins = adata.obs.loc[adata.obs.labels_he > 0].groupby("labels_he").size()
    oid = cdata.obs["object_id"].to_numpy(int)
    cdata.obs["n_bins_nucleus"] = pd.Series(nuc_bins).reindex(oid).fillna(0).to_numpy(int)
    cdata.obs.loc[~cdata.obs.has_nucleus, "n_bins_nucleus"] = 0
    na = nucleus_pixel_areas(io.npz_he, MPP).reindex(oid).to_numpy(float)
    na[~cdata.obs.has_nucleus.to_numpy()] = np.nan
    cdata.obs["nucleus_area_um2"] = na.astype(np.float32)
    cen = nucleus_centroids_fullres(io.npz_he, scalef, crop_origin).reindex(oid)
    for col in ("nucleus_x_fullres_px", "nucleus_y_fullres_px"):
        v = cen[col].to_numpy(float)
        v[~cdata.obs.has_nucleus.to_numpy()] = np.nan
        cdata.obs[col] = v
    log(f"  nucleus centroids: {int(np.isfinite(cdata.obs.nucleus_x_fullres_px).sum()):,} cells")

    morph = footprint_morphometrics(adata, "labels_joint").set_index("object_id").reindex(oid)
    for col in ("cell_area_um2", "cell_perimeter_um", "cell_circularity", "eccentricity",
                "footprint_bbox_rows", "footprint_bbox_cols"):
        cdata.obs[col] = morph[col].to_numpy()
    assert np.allclose(cdata.obs.cell_area_um2, cdata.obs.bin_count * BIN_SIDE_UM ** 2), \
        "footprint area disagrees with bin_count x 4 um^2"
    mpp_full = adata.uns["spatial"][list(adata.uns["spatial"])[0]]["scalefactors"]["microns_per_pixel"]
    cdata.obs["cell_area_px2"] = (cdata.obs.cell_area_um2 / mpp_full ** 2).astype(np.float32)
    with np.errstate(divide="ignore", invalid="ignore"):
        cdata.obs["nucleus_cell_area_ratio"] = (cdata.obs.nucleus_area_um2 / cdata.obs.cell_area_um2).astype(np.float32)
    tot_raw = np.asarray(cdata.layers["counts"].sum(axis=1)).ravel()
    tot_des = np.asarray(cdata.layers["counts_destriped"].sum(axis=1)).ravel()
    cdata.obs["total_counts_raw"] = tot_raw
    cdata.obs["total_counts_destriped"] = tot_des
    with np.errstate(divide="ignore", invalid="ignore"):
        cdata.obs["destripe_ratio"] = np.where(tot_raw > 0, tot_des / np.where(tot_raw > 0, tot_raw, 1), np.nan)
    cdata.obs["sample"] = io.sample
    cdata.obs["condition"] = io.condition
    cdata.obs["segmentation_method"] = "bin2cell_0.3.4_stardist_he_expanded"
    cdata.obs_names = [f"b2c_{io.sample}_{i}" for i in cdata.obs["object_id"]]
    return cdata


def registration_checks(adata, cdata, he_crop: np.ndarray, scalef: float, crop_origin, io: SampleIO) -> pd.DataFrame:
    """Verify the three coordinate systems agree and that cells land on stained tissue.

    Parameters
    ----------
    adata, cdata : AnnData
        Bin object and cell object.
    he_crop : ndarray
        The H&E render StarDist segmented (rows x cols x 3).
    scalef : float
        Full-res pixel -> render pixel scale factor.
    crop_origin : tuple[float, float]
        Full-res (x, y) of the render's top-left corner.
    io : SampleIO
        Paths of the sample.

    Returns
    -------
    DataFrame
        One row per check with ``status`` PASS / FAIL / INFO. Any FAIL blocks integration.
    """
    rows = []

    def chk(name, ok, value, expected):
        rows.append(dict(check=name, status="PASS" if ok else "FAIL", value=value, expected=expected))
        log(f"    [{'PASS' if ok else 'FAIL'}] {name:60s} {value}")

    lib = list(adata.uns["spatial"])[0]
    sf = adata.uns["spatial"][lib]["scalefactors"]
    hires = adata.uns["spatial"][lib]["images"]["hires"]
    full_h = hires.shape[0] / sf["tissue_hires_scalef"]; full_w = hires.shape[1] / sf["tissue_hires_scalef"]
    with Image.open(io.source_image) as im:
        W, H = im.size
    chk("full-res image size == hires / tissue_hires_scalef",
        abs(W - full_w) / W < 0.01 and abs(H - full_h) / H < 0.01,
        f"tif {W}x{H}; inferred {full_w:.0f}x{full_h:.0f}", "within 1%")
    xy = adata.obsm["spatial"]
    n_out = int(((xy[:, 0] >= W) | (xy[:, 1] >= H) | (xy < 0).any(axis=1)).sum())
    chk("bins inside the full-res image (capture area may overhang the image edge)",
        n_out / len(xy) <= MAX_FRAC_BINS_OUTSIDE_IMAGE,
        f"{n_out} of {len(xy):,} outside ({100 * n_out / len(xy):.3f}%); image {W}x{H}",
        f"<= {100 * MAX_FRAC_BINS_OUTSIDE_IMAGE:g}% outside")
    r = adata.obs["array_row"].to_numpy(); c = adata.obs["array_col"].to_numpy()
    one_row = r == pd.Series(r).mode()[0]
    xs = xy[one_row]; cs = c[one_row]; o = np.argsort(cs)
    pitch = np.abs(np.diff(xs[o, 0])[np.diff(cs[o]) == 1])
    expected_pitch = BIN_SIDE_UM / sf["microns_per_pixel"]
    chk("bin pitch along a row == 2 um / microns_per_pixel",
        abs(np.median(pitch) - expected_pitch) / expected_pitch < 0.02,
        f"{np.median(pitch):.3f} px (expected {expected_pitch:.3f})", "within 2%")
    xyc = adata.obsm["spatial_cropped_150_buffer"] * scalef
    inside = (xyc[:, 0] >= 0) & (xyc[:, 0] < he_crop.shape[1]) & (xyc[:, 1] >= 0) & (xyc[:, 1] < he_crop.shape[0])
    chk("bins inside the cropped render", float((~inside).mean()) <= MAX_FRAC_BINS_OUTSIDE_IMAGE,
        f"{int((~inside).sum())} outside; render {he_crop.shape[1]}x{he_crop.shape[0]} px",
        f"<= {100 * MAX_FRAC_BINS_OUTSIDE_IMAGE:g}% outside")
    chk("crop origin consistent (spatial - cropped == origin)",
        bool(np.allclose(adata.obsm["spatial"][:, 0] - adata.obsm["spatial_cropped_150_buffer"][:, 0], crop_origin[0]) and
             np.allclose(adata.obsm["spatial"][:, 1] - adata.obsm["spatial_cropped_150_buffer"][:, 1], crop_origin[1])),
        f"origin {crop_origin}", "constant offset")
    rng = np.random.default_rng(RANDOM_SEED)
    sub = rng.choice(cdata.n_obs, min(50000, cdata.n_obs), replace=False)
    cxy = cdata.obsm["spatial_cropped_150_buffer"][sub] * scalef
    px = he_crop[np.clip(cxy[:, 1].round().astype(int), 0, he_crop.shape[0] - 1),
                 np.clip(cxy[:, 0].round().astype(int), 0, he_crop.shape[1] - 1)].mean(axis=1)
    px_sw = he_crop[np.clip(cxy[:, 0].round().astype(int), 0, he_crop.shape[0] - 1),
                    np.clip(cxy[:, 1].round().astype(int), 0, he_crop.shape[1] - 1)].mean(axis=1)
    on = float((px < TISSUE_PIXEL_MAX_MEAN_RGB).mean()); on_sw = float((px_sw < TISSUE_PIXEL_MAX_MEAN_RGB).mean())
    chk("cell centroids land on stained tissue (mean RGB < 225)", on >= 0.90 and on > on_sw + 0.2,
        f"{100 * on:.2f}% (swapped-axes control {100 * on_sw:.2f}%)", ">= 90% and >> control")
    lab = sp.load_npz(io.npz_he).tocoo()
    take = rng.choice(lab.nnz, min(2_000_000, lab.nnz), replace=False)
    nuc_px = he_crop[lab.row[take], lab.col[take]].mean(axis=1)
    bg = he_crop[rng.integers(0, he_crop.shape[0], 2_000_000), rng.integers(0, he_crop.shape[1], 2_000_000)]
    bg_tissue = bg.mean(axis=1)[bg.mean(axis=1) < TISSUE_PIXEL_MAX_MEAN_RGB]
    chk("StarDist nucleus pixels darker than random tissue pixels",
        float(np.median(nuc_px)) < float(np.median(bg_tissue)) - 15,
        f"nucleus median RGB {np.median(nuc_px):.1f} vs tissue {np.median(bg_tissue):.1f}", "nucleus darker by > 15")
    chk("label mask shape == render shape", lab.shape == he_crop.shape[:2],
        f"{lab.shape} vs {he_crop.shape[:2]}", "equal")
    prim = cdata.obs["has_nucleus"].to_numpy()
    nxy = cdata.obs.loc[prim, ["nucleus_x_fullres_px", "nucleus_y_fullres_px"]].to_numpy(float)
    from scipy.spatial import cKDTree
    _, nearest = cKDTree(adata.obsm["spatial"]).query(nxy, k=1)
    own = adata.obs["labels_he"].to_numpy(int)[nearest] == cdata.obs.loc[prim, "object_id"].to_numpy(int)
    dist_um = np.hypot(*(nxy - cdata.obsm["spatial"][prim]).T) * sf["microns_per_pixel"]
    chk("bin nearest each nucleus centroid lies in that nucleus (centroid transform)",
        float(own.mean()) >= MIN_FRAC_NUCLEUS_CENTROID_IN_OWN_NUCLEUS,
        f"{100 * own.mean():.2f}%; nucleus-to-cell centroid median {np.median(dist_um):.2f} um, p99 {np.percentile(dist_um, 99):.2f} um",
        f">= {100 * MIN_FRAC_NUCLEUS_CENTROID_IN_OWN_NUCLEUS:g}%")
    ac = adata.uns["bin2cell"]["array_check"]
    rows.append(dict(check="array_row flipped vs image y (bin2cell check_array_coordinates)",
                     status="INFO", value=ac["row"]["flipped"], expected="recorded"))
    rows.append(dict(check="array_col flipped vs image x", status="INFO", value=ac["col"]["flipped"],
                     expected="recorded"))
    return pd.DataFrame(rows)


def qc_tables(adata, cdata, io: SampleIO) -> None:
    """Write the unfiltered QC tables of one sample into ``io.out / 'qc'``.

    Writes ``bin_assignment_summary.csv`` (where the bins and UMIs went),
    ``cell_metrics_unfiltered.csv`` (one row per object), ``cell_metrics_summary.csv``
    (distribution quantiles) and ``qc_gate_preview.csv`` -- a PREVIEW of the downstream QC
    gates for review only; notebook 02 applies the authoritative gates.
    """
    import scanpy as sc
    out = io.out
    cdata.var["mt"] = cdata.var_names.str.match(r"^mt-", case=False)
    sc.pp.calculate_qc_metrics(cdata, qc_vars=["mt"], percent_top=None, log1p=True, inplace=True)
    n = adata.n_obs
    umi = adata.obs.n_counts.to_numpy(float)
    o = adata.obs

    def row(quantity, mask):
        return dict(quantity=quantity, n=int(mask.sum()), frac_of_bins=float(mask.sum()) / n,
                    umi=umi[mask].sum(), frac_of_umi=umi[mask].sum() / umi.sum())

    def count_row(quantity, value):
        return dict(quantity=quantity, n=int(value), frac_of_bins=np.nan, umi=np.nan, frac_of_umi=np.nan)

    rows = [
        row("bins under tissue (2 um filtered matrix)", np.ones(n, bool)),
        row("bins inside a StarDist nucleus", (o.labels_he > 0).to_numpy()),
        row(f"bins in an expanded H&E cell (<= {MAX_BIN_DISTANCE} bins)", (o.labels_he_expanded > 0).to_numpy()),
        row("bins in a GEX-derived object", (o.labels_gex > 0).to_numpy()),
        row("bins assigned in labels_joint (H&E cells + salvaged GEX objects)", (o.labels_joint > 0).to_numpy()),
        row("  of which salvaged from GEX (secondary)", (o.labels_joint_source == "secondary").to_numpy()),
        row("bins UNASSIGNED", (o.labels_joint == 0).to_numpy()),
        count_row("cells (objects) total", cdata.n_obs),
        count_row("cells with a nucleus (primary)", cdata.obs.has_nucleus.sum()),
        count_row("StarDist nuclei detected (H&E)", len(np.unique(sp.load_npz(io.npz_he).data))),
        count_row("StarDist nuclei that received >= 1 bin", (o.labels_he[o.labels_he > 0]).nunique()),
    ]
    pd.DataFrame(rows).to_csv(out / "qc" / "bin_assignment_summary.csv", index=False)
    for r in rows:
        log(f"    {r['quantity']:70s} {r['n']:>12,}" + (f"  ({100 * r['frac_of_bins']:.1f}% bins, {100 * r['frac_of_umi']:.1f}% UMI)" if np.isfinite(r['frac_of_bins']) else ""))

    cols = ["object_id", "has_nucleus", "labels_joint_source", "bin_count", "n_bins_nucleus",
            "cell_area_um2", "cell_perimeter_um", "cell_circularity", "nucleus_area_um2",
            "nucleus_cell_area_ratio", "total_counts", "n_genes_by_counts", "pct_counts_mt",
            "total_counts_destriped", "destripe_ratio", "array_row", "array_col"]
    met = cdata.obs[cols].copy()
    met["x_fullres_px"] = cdata.obsm["spatial"][:, 0]; met["y_fullres_px"] = cdata.obsm["spatial"][:, 1]
    met.index.name = "cell_id"
    met.to_csv(out / "qc" / "cell_metrics_unfiltered.csv")

    summ = []
    for subset, msk in (("all", np.ones(cdata.n_obs, bool)), ("primary_nucleus_anchored", cdata.obs.has_nucleus.to_numpy())):
        for col in ("total_counts", "n_genes_by_counts", "bin_count", "cell_area_um2", "nucleus_area_um2",
                    "nucleus_cell_area_ratio", "pct_counts_mt", "cell_circularity", "destripe_ratio"):
            v = cdata.obs[col].to_numpy(float)[msk]; v = v[np.isfinite(v)]
            summ.append(dict(subset=subset, metric=col, n=len(v), mean=v.mean(), median=np.median(v),
                             p05=np.percentile(v, 5), p25=np.percentile(v, 25), p75=np.percentile(v, 75),
                             p95=np.percentile(v, 95), p99=np.percentile(v, 99), max=v.max()))
    pd.DataFrame(summ).round(4).to_csv(out / "qc" / "cell_metrics_summary.csv", index=False)

    keep = np.ones(cdata.n_obs, bool); led = []
    n0 = cdata.n_obs

    def gate(name, mask, detail):
        nonlocal keep
        rem = int((keep & ~mask).sum()); keep &= mask
        led.append(dict(criterion=name, n_removed=rem, pct_removed=round(100 * rem / n0, 3),
                        n_remaining=int(keep.sum()), detail=detail))

    gate("has_nucleus (primary, H&E-anchored object)", cdata.obs.has_nucleus.to_numpy(),
         "GEX-only salvaged objects have no nucleus; only nucleus-anchored cells enter notebook 02")
    gate(f"total_counts >= {PREVIEW_MIN_COUNTS}", cdata.obs.total_counts.to_numpy() >= PREVIEW_MIN_COUNTS, "QC_MIN_COUNTS (notebook 02)")
    gate(f"n_genes >= {PREVIEW_MIN_GENES}", cdata.obs.n_genes_by_counts.to_numpy() >= PREVIEW_MIN_GENES, "QC_MIN_GENES (notebook 02)")
    for m in ("log1p_total_counts", "log1p_n_genes_by_counts"):
        v = cdata.obs[m].to_numpy(float)[keep]; med = np.median(v); mad = np.median(np.abs(v - med))
        lo, hi = med - PREVIEW_MAD_NMADS * mad, med + PREVIEW_MAD_NMADS * mad
        vv = cdata.obs[m].to_numpy(float)
        gate(f"MAD({PREVIEW_MAD_NMADS}) {m}", (vv >= lo) & (vv <= hi), f"keep raw [{np.expm1(lo):.0f}, {np.expm1(hi):.0f}]")
    a = cdata.obs.cell_area_um2.to_numpy(float)
    gate(f"cell area in [{PREVIEW_MIN_AREA_UM2}, {PREVIEW_MAX_AREA_UM2}] um^2", (a >= PREVIEW_MIN_AREA_UM2) & (a <= PREVIEW_MAX_AREA_UM2), "MIN/MAX_CELL_AREA_UM2")
    led.append(dict(criterion="TOTAL (preview)", n_removed=n0 - int(keep.sum()), pct_removed=round(100 * (n0 - keep.sum()) / n0, 3),
                    n_remaining=int(keep.sum()), detail="PREVIEW ONLY -- authoritative ledger: notebook 02, tables/qc_filter_ledger.csv"))
    pd.DataFrame(led).to_csv(out / "qc" / "qc_gate_preview.csv", index=False)
    for r in led:
        log(f"    {r['criterion']:55s} -{r['n_removed']:>8,} -> {r['n_remaining']:>8,}")


def qc_figures(adata, cdata, he_crop: np.ndarray, scalef: float, io: SampleIO, crop_origin=(0.0, 0.0)) -> None:
    """Write the QC figures of one sample into ``io.out / 'qc'``.

    Figures: per-cell metric histograms (nucleus-anchored vs GEX-only), cells on the hires
    H&E coloured by UMI / area / provenance, ``N_ZOOM_REGIONS`` zoomed overlays (StarDist
    nucleus borders, cell footprint edges, centroids), a whole-section overview with the zoom
    boxes, and the area / counts / bins relationships.
    """
    fig_dir = io.out / "qc"
    prim = cdata.obs.has_nucleus.to_numpy()
    crop0 = np.asarray(crop_origin, dtype=float)

    fig, axes = plt.subplots(2, 4, figsize=(19, 8.5))
    panels = [("total_counts", "UMI per cell", True), ("n_genes_by_counts", "genes per cell", True),
              ("bin_count", "2 um bins per cell", False), ("cell_area_um2", "cell area (um^2)", False),
              ("nucleus_area_um2", "nucleus area (um^2)", False), ("nucleus_cell_area_ratio", "nucleus : cell area", False),
              ("pct_counts_mt", "% mitochondrial", False), ("destripe_ratio", "destriped / raw total", False)]
    for ax, (col, lab, logx) in zip(axes.ravel(), panels):
        for msk, name, c in ((prim, "nucleus-anchored", "#2a78d6"), (~prim, "GEX-only salvaged", "#b04a2a")):
            raw = cdata.obs[col].to_numpy(float)[msk]; raw = raw[np.isfinite(raw)]
            if len(raw) == 0:
                continue
            v = np.log10(raw + 1) if logx else raw
            hi = np.percentile(v, 99.5)
            ax.hist(v, bins=np.linspace(v.min(), hi, 60), histtype="stepfilled", alpha=0.45, color=c,
                    label=f"{name} (n={len(v):,}, median {np.median(raw):.1f})")
        ax.set_xlabel(("log10(1 + " + lab + ")") if logx else lab); ax.set_ylabel("cells"); ax.legend(fontsize=6.5)
    fig.suptitle(f"{io.sample} Bin2Cell cells -- unfiltered per-cell metrics", fontsize=12)
    fig.tight_layout(); fig.savefig(fig_dir / "cell_metrics_unfiltered_hist.png", dpi=170); plt.close(fig)

    lib = list(adata.uns["spatial"])[0]
    hs = adata.uns["spatial"][lib]["scalefactors"]["tissue_hires_scalef"]
    hires = np.asarray(Image.open(io.image_dir / "tissue_hires_image.png").convert("RGB"))
    xy = cdata.obsm["spatial"] * hs
    fig, axes = plt.subplots(1, 3, figsize=(21, 7.5))
    for ax, (col, lab, cmap) in zip(axes, (("total_counts", "UMI per cell", "viridis"), ("cell_area_um2", "cell area (um^2)", "cividis"),
                                           ("has_nucleus", "nucleus-anchored (1) vs GEX-only (0)", "coolwarm"))):
        v = cdata.obs[col].to_numpy(float)
        ax.imshow(hires); vmin, vmax = np.nanpercentile(v, [2, 98]) if col != "has_nucleus" else (0, 1)
        s_ = ax.scatter(xy[:, 0], xy[:, 1], c=v, s=0.15, cmap=cmap, vmin=vmin, vmax=vmax, linewidths=0, rasterized=True)
        fig.colorbar(s_, ax=ax, fraction=0.03, pad=0.01, label=lab); ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"{lab}  (n={cdata.n_obs:,} cells on hires H&E)", fontsize=9)
    fig.tight_layout(); fig.savefig(fig_dir / "cells_spatial_on_hires.png", dpi=200); plt.close(fig)

    img, r, c = label_image_2um(adata, "labels_joint")
    lab = img[r, c]

    def nb(dr, dc):
        rr, cc = r + dr, c + dc
        ok = (rr >= 0) & (rr < img.shape[0]) & (cc >= 0) & (cc < img.shape[1])
        v = np.full(len(r), -1); v[ok] = img[rr[ok], cc[ok]]
        return v

    boundary = (lab > 0) & ((nb(1, 0) != lab) | (nb(-1, 0) != lab) | (nb(0, 1) != lab) | (nb(0, -1) != lab))
    bxy = adata.obsm["spatial_cropped_150_buffer"] * scalef
    src = adata.obs.labels_joint_source.astype(str).to_numpy()
    cxy = cdata.obsm["spatial_cropped_150_buffer"] * scalef
    side = int(ZOOM_SIDE_UM / MPP)
    rng = np.random.default_rng(RANDOM_SEED)
    qs = [(0.15, 0.5), (0.5, 0.2), (0.5, 0.8), (0.85, 0.5)]
    centres = [(np.quantile(cxy[:, 0], qx), np.quantile(cxy[:, 1], qy)) for qx, qy in qs]
    pick = rng.choice(cdata.n_obs, N_ZOOM_REGIONS - len(qs), replace=False)
    centres += [tuple(cxy[i]) for i in pick]
    crops = []
    for k, (cx, cy) in enumerate(centres):
        x0 = int(max(0, cx - side // 2)); y0 = int(max(0, cy - side // 2))
        x1 = min(he_crop.shape[1], x0 + side); y1 = min(he_crop.shape[0], y0 + side)
        crops.append((x0, y0, x1, y1))
        nuc = b2c.view_labels(image_path=str(io.he_render), labels_npz_path=str(io.npz_he),
                              crop=(x0, y0, x1, y1), fill=False, border=True, border_color=[255, 255, 0])
        fig, axes = plt.subplots(1, 3, figsize=(21, 7.4))
        axes[0].imshow(he_crop[y0:y1, x0:x1]); axes[0].set_title(f"region {k + 1}: H&E ({ZOOM_SIDE_UM:.0f} um square)")
        axes[1].imshow(nuc); axes[1].set_title("StarDist nucleus borders (yellow)")
        axes[2].imshow(he_crop[y0:y1, x0:x1])
        m = (bxy[:, 0] >= x0) & (bxy[:, 0] < x1) & (bxy[:, 1] >= y0) & (bxy[:, 1] < y1)
        mb = m & boundary
        axes[2].scatter(bxy[mb & (src == "primary"), 0] - x0, bxy[mb & (src == "primary"), 1] - y0, s=1.6, c="#00d4ff", linewidths=0, label="cell footprint edge (nucleus-anchored)")
        axes[2].scatter(bxy[mb & (src == "secondary"), 0] - x0, bxy[mb & (src == "secondary"), 1] - y0, s=1.6, c="#ff7f0e", linewidths=0, label="footprint edge (GEX-only)")
        mn = m & (adata.obs.labels_he.to_numpy() > 0)
        axes[2].scatter(bxy[mn, 0] - x0, bxy[mn, 1] - y0, s=0.5, c="#ffff00", alpha=0.5, linewidths=0, label="bins inside nucleus")
        mc = (cxy[:, 0] >= x0) & (cxy[:, 0] < x1) & (cxy[:, 1] >= y0) & (cxy[:, 1] < y1)
        axes[2].scatter(cxy[mc, 0] - x0, cxy[mc, 1] - y0, s=6, c="#ff00ff", marker="+", linewidths=0.7, label="cell centroid")
        nuc_c = (cdata.obs[["nucleus_x_fullres_px", "nucleus_y_fullres_px"]].to_numpy(float) - crop0) * scalef
        mnc = np.isfinite(nuc_c[:, 0]) & (nuc_c[:, 0] >= x0) & (nuc_c[:, 0] < x1) & (nuc_c[:, 1] >= y0) & (nuc_c[:, 1] < y1)
        axes[2].scatter(nuc_c[mnc, 0] - x0, nuc_c[mnc, 1] - y0, s=5, c="#00ff00", marker="x", linewidths=0.6, label="nucleus centroid")
        axes[2].set_xlim(0, x1 - x0); axes[2].set_ylim(y1 - y0, 0)
        axes[2].set_title(f"bins -> cells: {int(mc.sum())} centroids in view"); axes[2].legend(fontsize=6.5, loc="lower right")
        for ax in axes:
            ax.set_xticks([]); ax.set_yticks([])
        um_px = 50 / MPP
        axes[0].plot([15, 15 + um_px], [side - 20, side - 20], color="k", lw=3); axes[0].text(15 + um_px / 2, side - 28, "50 um", ha="center", fontsize=8)
        fig.suptitle(f"{io.sample} -- Bin2Cell segmentation overlay, region {k + 1} "
                     f"(render {MPP} um/px; crop x[{x0},{x1}) y[{y0},{y1}) px)", fontsize=11)
        fig.tight_layout(); fig.savefig(fig_dir / f"overlay_region_{k + 1}.png", dpi=OVERLAY_DPI); plt.close(fig)
    pd.DataFrame(crops, columns=["x0", "y0", "x1", "y1"]).assign(units=f"px at {MPP} um/px in {io.he_render.name}").to_csv(fig_dir / "overlay_regions.csv", index=False)

    ds = 8
    small = he_crop[::ds, ::ds]
    fig, ax = plt.subplots(figsize=(12, 12 * small.shape[0] / small.shape[1]))
    ax.imshow(small)
    sub = rng.choice(cdata.n_obs, min(40000, cdata.n_obs), replace=False)
    ax.scatter(cxy[sub, 0] / ds, cxy[sub, 1] / ds, s=0.3, c="#00a0ff", alpha=0.5, linewidths=0)
    for k, (x0, y0, x1, y1) in enumerate(crops):
        ax.add_patch(plt.Rectangle((x0 / ds, y0 / ds), (x1 - x0) / ds, (y1 - y0) / ds, fill=False, color="red", lw=1.2))
        ax.text(x0 / ds, y0 / ds - 3, str(k + 1), color="red", fontsize=9)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_title(f"{io.sample}: 40k cell centroids on the segmented H&E render; red = zoom regions")
    fig.savefig(fig_dir / "overview_centroids_and_regions.png", dpi=150, bbox_inches="tight"); plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    ok = prim & np.isfinite(cdata.obs.nucleus_area_um2.to_numpy(float))
    axes[0].hexbin(cdata.obs.nucleus_area_um2.to_numpy(float)[ok], cdata.obs.cell_area_um2.to_numpy(float)[ok], gridsize=60, bins="log", cmap="Greys")
    axes[0].set_xlabel("nucleus area (um^2, StarDist mask)"); axes[0].set_ylabel("cell area (um^2, bins x 4)"); axes[0].set_title("expansion check: cell vs nucleus area")
    axes[1].hexbin(cdata.obs.cell_area_um2.to_numpy(float), np.log10(cdata.obs.total_counts.to_numpy(float) + 1), gridsize=60, bins="log", cmap="Greys")
    axes[1].set_xlabel("cell area (um^2)"); axes[1].set_ylabel("log10(UMI + 1)"); axes[1].set_title("counts vs area (technical coupling)")
    axes[2].hist(cdata.obs.bin_count.to_numpy()[prim], bins=np.arange(0, 120, 2), color="#2a78d6", alpha=0.7, label="nucleus-anchored")
    axes[2].hist(cdata.obs.bin_count.to_numpy()[~prim], bins=np.arange(0, 120, 2), color="#b04a2a", alpha=0.7, label="GEX-only")
    axes[2].set_xlabel("bins per cell"); axes[2].set_ylabel("cells"); axes[2].legend(); axes[2].set_title("bins per cell")
    fig.tight_layout(); fig.savefig(fig_dir / "area_counts_bins_relationships.png", dpi=170); plt.close(fig)


def run_bin2cell_sample(sample: str, out_base: Path, force: bool = False) -> int:
    """Full Bin2Cell reconstruction of one sample (7 stages, ~10 min on one GPU).

    Skipped entirely -- nothing is written -- when the final cell object already exists and
    ``force`` is False. The expensive intermediates (H&E render, both StarDist label files,
    the GEX image) are reused when present, so an interrupted run restarts where it stopped.

    Parameters
    ----------
    sample : str
        Key of ``SAMPLE_INPUTS``.
    out_base : Path
        Root output directory (``B2C_DIR``).
    force : bool
        Recompute every intermediate and overwrite the cell object.

    Returns
    -------
    int
        Number of failed registration checks (0 = safe to integrate).
    """
    io = sample_io(sample, out_base)
    if io.cells_h5ad.exists() and not force:
        log(f"  {sample}: {io.cells_h5ad.name} exists -- skipped (FORCE=True recomputes)")
        return 0
    (io.out / "qc").mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    banner(f"Bin2Cell: {sample}  ->  {io.out}")
    if not io.source_image.exists() or not (io.bin_dir / "filtered_feature_bc_matrix.h5").exists():
        raise FileNotFoundError(f"inputs missing: {io.source_image} / {io.bin_dir}")
    env = record_environment(io.out)
    log(f"  env: bin2cell {env['bin2cell']}, stardist {env['stardist']}, tensorflow {env['tensorflow']}, GPUs {env['gpus']}")
    params = dict(PARAMS, INPUTS={k: str(v) for k, v in SAMPLE_INPUTS[sample].items()})
    (io.out / "params.json").write_text(json.dumps(params, indent=2))

    banner("1/7  read 2 um bins")
    adata = read_bins(io)
    lib = list(adata.uns["spatial"])[0]

    banner(f"2/7  H&E render at {MPP} um/px (cropped to capture area + {CROP_BUFFER} px)")
    # scaled_he_image also checks the full-res image dimensions against the hires scale
    # factor (raises on mismatch) and stores the cropped coordinates in obsm.
    reuse_he = io.he_render.exists() and not force
    b2c.scaled_he_image(adata, mpp=MPP, crop=True, buffer=CROP_BUFFER, save_path=None if reuse_he else str(io.he_render))
    log(f"  {'reusing' if reuse_he else 'wrote'} {io.he_render.name}")
    scalef = float(adata.uns["spatial"][lib]["scalefactors"][f"tissue_{MPP}_mpp_{CROP_BUFFER}_buffer_scalef"])
    crop_origin = (float((adata.obsm["spatial"][:, 0] - adata.obsm["spatial_cropped_150_buffer"][:, 0])[0]),
                   float((adata.obsm["spatial"][:, 1] - adata.obsm["spatial_cropped_150_buffer"][:, 1])[0]))
    he_crop = b2c.load_image(str(io.he_render))
    log(f"  render {he_crop.shape[1]}x{he_crop.shape[0]} px; scale factor {scalef:.6f}; crop origin {crop_origin}")
    del adata.uns["spatial"][lib]["images"][f"{MPP}_mpp_{CROP_BUFFER}_buffer"]   # not carried in the h5ad

    banner("3/7  destripe")
    safe_destripe(adata)

    banner("4/7  StarDist nuclear segmentation of the H&E render")
    run_stardist(io.he_render, io.npz_he, force, stardist_model=STARDIST_HE_MODEL, prob_thresh=STARDIST_HE_PROB_THRESH)
    b2c.insert_labels(adata, str(io.npz_he), basis="spatial", spatial_key="spatial_cropped_150_buffer",
                      mpp=MPP, labels_key="labels_he")
    log(f"  bins inside a nucleus: {int((adata.obs.labels_he > 0).sum()):,} / {adata.n_obs:,}")
    b2c.expand_labels(adata, labels_key="labels_he", expanded_labels_key="labels_he_expanded",
                      algorithm="max_bin_distance", max_bin_distance=MAX_BIN_DISTANCE, k=EXPAND_K)
    log(f"  bins in an expanded cell (<= {MAX_BIN_DISTANCE} bins, expression-PCA tie break): "
        f"{int((adata.obs.labels_he_expanded > 0).sum()):,}")

    banner("5/7  GEX-based segmentation (salvage of nucleus-free objects)")
    b2c.check_array_coordinates(adata)            # needed by insert_labels(basis="array")
    log(f"  array_check: {adata.uns['bin2cell']['array_check']}")
    if not io.gex_render.exists() or force:
        b2c.grid_image(adata, "n_counts_adjusted", mpp=MPP, sigma=GEX_SIGMA, save_path=str(io.gex_render))
    run_stardist(io.gex_render, io.npz_gex, force, stardist_model=STARDIST_GEX_MODEL,
                 prob_thresh=STARDIST_GEX_PROB_THRESH, nms_thresh=STARDIST_GEX_NMS_THRESH)
    b2c.insert_labels(adata, str(io.npz_gex), basis="array", mpp=MPP, labels_key="labels_gex")
    b2c.salvage_secondary_labels(adata, primary_label="labels_he_expanded", secondary_label="labels_gex",
                                 labels_key="labels_joint")
    bins_tab = adata.obs[["in_tissue", "array_row", "array_col", "n_counts", "n_genes", "n_counts_adjusted",
                          "destripe_factor", "labels_he", "labels_he_expanded", "labels_gex", "labels_joint",
                          "labels_joint_source"]].copy()
    bins_tab["x_fullres_px"] = adata.obsm["spatial"][:, 0]; bins_tab["y_fullres_px"] = adata.obsm["spatial"][:, 1]
    bins_tab.to_parquet(io.out / "bins_labels.parquet")
    log(f"  wrote bins_labels.parquet ({len(bins_tab):,} rows)")

    banner("6/7  bins -> cells")
    cdata = build_cells(adata, io, scalef, crop_origin)
    cfg = SAMPLE_INPUTS[sample]
    cdata.uns["bin2cell"] = dict(
        params=params, environment=env, source_bins=str(io.bin_dir), source_image=str(io.source_image),
        he_render=str(io.he_render), he_render_scalef=scalef, he_render_mpp=MPP,
        crop_origin_fullres_px=list(crop_origin), labels_he_npz=str(io.npz_he), labels_gex_npz=str(io.npz_gex),
        coordinate_note="obsm['spatial'] = mean of member-bin centres in FULL-RESOLUTION H&E pixels "
                        "(x = column, y = row, origin top-left); obsm['spatial_cropped_150_buffer'] = "
                        "spatial - crop_origin, x scalef indexes the render",
        array_check=adata.uns["bin2cell"]["array_check"],
        counts_note="X and layers['counts'] are RAW INTEGER bin sums; layers['counts_destriped'] are "
                    "the b2c.destripe-adjusted sums (non-integer)")
    cdata.uns["spatial"][lib]["metadata"] = dict(
        sample=sample, chemistry_description=cfg["chemistry"], spaceranger_version=cfg["spaceranger"],
        assay=cfg["assay"], genome=cfg["genome"],
        source_image_path=str(io.image_dir / "tissue_hires_image.png"),
        fullres_image_path=str(io.source_image), segmentation="bin2cell 0.3.4")
    cdata.uns["spatial"][lib]["images"] = {}

    banner("7/7  registration checks, QC, overlays")
    checks = registration_checks(adata, cdata, he_crop, scalef, crop_origin, io)
    checks.to_csv(io.out / "qc" / "registration_checks.csv", index=False)
    n_fail = int((checks.status == "FAIL").sum())
    qc_tables(adata, cdata, io)
    try:
        qc_figures(adata, cdata, he_crop, scalef, io, crop_origin)
    except Exception as exc:  # noqa: BLE001
        log(f"  QC FIGURES FAILED: {exc!r}")
        n_fail += 1
    cdata.write_h5ad(io.cells_h5ad, compression="gzip")
    log(f"  wrote {io.cells_h5ad} ({cdata.n_obs:,} cells x {cdata.n_vars:,} genes)")
    log(f"  done in {(time.time() - t0) / 60:.1f} min; registration checks failed: {n_fail}")
    if n_fail:
        warn_banner(f"{sample}: REGISTRATION / QC FAILURE -- do not integrate until resolved")
    return n_fail


def bin2cell_status(out_base: Path, samples: list[str]) -> pd.DataFrame:
    """One row per sample: is the cell object there, how many cells, did the checks pass.

    Parameters
    ----------
    out_base : Path
        Root Bin2Cell output directory (``B2C_DIR``).
    samples : list[str]
        Samples to report.
    """
    rows = []
    for s in samples:
        io = sample_io(s, out_base)
        r = dict(sample=s, condition=io.condition, spaceranger=SAMPLE_INPUTS[s]["spaceranger"],
                 genome=SAMPLE_INPUTS[s]["genome"], cell_object=io.cells_h5ad.exists())
        p = io.out / "qc" / "bin_assignment_summary.csv"
        if p.exists():
            t = pd.read_csv(p).set_index("quantity")["n"]
            r["n_objects"] = int(t.get("cells (objects) total", np.nan))
            r["n_nucleus_anchored"] = int(t.get("cells with a nucleus (primary)", np.nan))
        p = io.out / "qc" / "registration_checks.csv"
        if p.exists():
            r["registration_checks_failed"] = int((pd.read_csv(p).status == "FAIL").sum())
        p = io.out / "qc" / "qc_gate_preview.csv"
        if p.exists():
            r["n_after_gate_preview"] = int(pd.read_csv(p).n_remaining.iloc[-1])
        rows.append(r)
    return pd.DataFrame(rows)


def bin2cell_table(out_base: Path, samples: list[str], name: str) -> pd.DataFrame:
    """Concatenate one ``qc/<name>.csv`` across samples (e.g. ``registration_checks``).

    Registration checks were renamed after ``Ctrl_1A2`` was first run (2026-09-11); old names
    are mapped through ``CHECK_NAME_ALIASES`` (the old checks were stricter, so an old PASS
    implies the new one) and the text as written is kept in ``check_as_written``.
    """
    frames = []
    for s in samples:
        p = sample_io(s, out_base).out / "qc" / f"{name}.csv"
        if p.exists():
            frames.append(pd.read_csv(p).assign(sample=s))
    df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    if "check" in df:
        df["check_as_written"] = df["check"]
        df["check"] = df["check"].replace(CHECK_NAME_ALIASES)
    return df


if __name__ == "__main__":
    failed = {}
    for s in SAMPLES:
        failed[s] = run_bin2cell_sample(s, B2C_DIR, force=FORCE)
    print("registration checks failed per sample:", failed, flush=True)
    sys.exit(1 if any(failed.values()) else 0)
