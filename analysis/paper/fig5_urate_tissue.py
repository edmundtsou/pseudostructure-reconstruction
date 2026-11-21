"""Prepare current COAST Fig. 5 units and tissue-level URAT1/NPT4 expression flags.

This is a descriptive spatial view of the four sections in the hvm COAST fit. Coordinates,
labels, and counts all come from that fit, so no separate segmentation or older tubule matrix
is joined here.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = next(p for p in Path(__file__).resolve().parents if (p / 'pseudospace').is_dir())
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
if str(Path(__file__).resolve().parent) not in sys.path:
    sys.path.insert(0, str(Path(__file__).resolve().parent))

import fig4_screen as fs  # noqa: E402
import fig5_screen as f5  # noqa: E402
from fig3_common import RUNS  # noqa: E402

TARGETS = ('Slc17a3', 'Slc22a12', 'Slc17a1')


def _checked_xy(qc_obs: pd.DataFrame, p1_obs: pd.DataFrame) -> pd.DataFrame:
    """Use QC centroids after checking their coordinates against the current pass-1 units."""
    cols = ['centroid_x', 'centroid_y']
    common = qc_obs.index.intersection(p1_obs.index)
    if len(common) == 0:
        raise ValueError('No shared unit names between current QC and pass-1 objects')
    a = qc_obs.loc[common, cols].to_numpy(float)
    b = p1_obs.loc[common, cols].to_numpy(float)
    finite = np.isfinite(a).all(axis=1) & np.isfinite(b).all(axis=1)
    if not finite.all() or np.any(np.abs(a[finite] - b[finite]) > 1.0):
        raise ValueError('QC/pass-1 unit centroids disagree by more than the 1-pixel tolerance')
    return qc_obs[cols]


def tissue_units(results_root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return current hvm PT units and section background centroids.

    PT rows are indexed by COAST observation names and include native-image centroids, the
    Fig. 5 mouse-reference quantile, Fig. 4 PT segment code, non-mito library, mito fraction,
    nearest Podocyte distance, and the three raw target-gene counts. The second table contains
    the current pseudospace-reference unit centroids used as each section's tissue background.
    """
    import anndata as ad
    from scipy import sparse
    from scipy.spatial import cKDTree

    root = Path(results_root)
    run = root / RUNS['hvm']
    qc = ad.read_h5ad(run / '02_qc/units_qc.h5ad')
    ref = pd.read_parquet(run / '04_pseudospace/pseudospace.parquet')
    ref = ref.set_index('obs_name') if 'obs_name' in ref.columns else ref
    if not ref.index.is_unique or not qc.obs_names.is_unique:
        raise ValueError('Current hvm observation names must be unique')

    p1 = ad.read_h5ad(run / '03_integration/pass1_labelled.h5ad', backed='r')
    xy = _checked_xy(qc.obs, p1.obs)
    p1_obs = p1.obs

    sections = f5.SECTIONS
    pt = ref[ref['sample'].astype(str).isin(sections)
             & ref['segment_class'].astype(str).isin(fs.PT)].copy()
    if not pt.index.isin(qc.obs_names).all():
        raise ValueError('Current pseudospace PT units are missing from current QC observations')
    pt['sample'] = pt['sample'].astype(str)
    sample = pt['sample'].to_numpy()
    seg = pt['segment_class'].astype(str).map(fs.PT)
    if seg.isna().any():
        raise ValueError('A PT label has no Fig. 4 PT segment code')

    counts = sparse.csc_matrix(qc[pt.index.to_numpy()].layers['counts'])
    mt = qc.var['mt'].to_numpy(bool)
    total = np.asarray(counts.sum(axis=1)).ravel().astype(float)
    library = np.asarray(counts[:, ~mt].sum(axis=1)).ravel().astype(float)
    mito = np.divide(np.asarray(counts[:, mt].sum(axis=1)).ravel(), total,
                     out=np.zeros_like(total), where=total > 0)
    gene_idx = qc.var_names.get_indexer(TARGETS)
    if (gene_idx < 0).any():
        missing = [TARGETS[i] for i in np.flatnonzero(gene_idx < 0)]
        raise ValueError(f'Target genes absent from current hvm count matrix: {missing}')
    target_counts = counts[:, gene_idx].toarray()

    positions = fs.quantile_axis(pt['pseudospace'].to_numpy(float), np.isin(sample, f5.MOUSE))
    pt_xy = xy.loc[pt.index, ['centroid_x', 'centroid_y']].to_numpy(float)
    glom_distance = np.full(len(pt), np.nan)
    for section in sections:
        glomeruli = p1_obs[(p1_obs['sample'].astype(str) == section)
                           & (p1_obs['segment_class'].astype(str) == 'Podocyte')]
        ix = np.flatnonzero(sample == section)
        if len(glomeruli) and len(ix):
            tree = cKDTree(glomeruli[['centroid_x', 'centroid_y']].to_numpy(float))
            glom_distance[ix] = tree.query(pt_xy[ix])[0]

    units = pd.DataFrame({
        'sample': sample,
        'x': pt_xy[:, 0],
        'y': pt_xy[:, 1],
        'position': positions,
        'segment': seg.to_numpy(dtype=int),
        'library': library,
        'mito': mito,
        'glom_distance': glom_distance,
        **{gene: target_counts[:, i] for i, gene in enumerate(TARGETS)},
    }, index=pt.index)
    units.index.name = pt.index.name

    background_index = ref.index.intersection(xy.index)
    background = pd.DataFrame({
        'sample': ref.loc[background_index, 'sample'].astype(str).to_numpy(),
        'x': xy.loc[background_index, 'centroid_x'].to_numpy(float),
        'y': xy.loc[background_index, 'centroid_y'].to_numpy(float),
    }, index=background_index)
    p1.file.close()
    return units, background


def depth_strata(library: pd.Series) -> pd.Series:
    """Return stable, near-equal-count library strata numbered 0–9, aligned to the input."""
    if library.index.has_duplicates:
        raise ValueError('Library index must be unique')
    values = library.to_numpy(float)
    if not np.isfinite(values).all():
        raise ValueError('Library must be finite')
    strata = np.empty(len(library), dtype=np.int8)
    order = np.argsort(values, kind='stable')
    for number, group in enumerate(np.array_split(order, 10)):
        strata[group] = number
    return pd.Series(strata, index=library.index, name='depth_stratum')


def high_expression(units: pd.DataFrame, fraction: float = 0.25) -> pd.DataFrame:
    """Flag within-decile high expressors for NPT4 and URAT1 in one section.

    Units are stably sorted by non-mito library and split into ten near-equal groups. Within
    each group, the top `fraction` by raw counts per non-mito library is flagged separately
    for each transcript, up to the number with positive counts. Stable sorting breaks ties.
    """
    if not 0 < fraction <= 1:
        raise ValueError('fraction must be in (0, 1]')
    required = {'library', 'Slc17a3', 'Slc22a12'}
    missing = required.difference(units.columns)
    if missing:
        raise KeyError(f'Missing required unit columns: {sorted(missing)}')
    if units.index.has_duplicates:
        raise ValueError('Unit index must be unique')
    if (units['library'] <= 0).any() or not np.isfinite(units['library']).all():
        raise ValueError('Non-mito library must be finite and positive')

    out = pd.DataFrame(False, index=units.index, columns=['Slc17a3', 'Slc22a12'], dtype=bool)
    strata = depth_strata(units['library'])
    groups = [np.flatnonzero(strata.to_numpy() == k) for k in range(10)]
    for gene in out.columns:
        rate = units[gene].to_numpy(float) / units['library'].to_numpy(float)
        flags = np.zeros(len(units), dtype=bool)
        for group in groups:
            if not len(group):
                continue
            order = group[np.argsort(-rate[group], kind='stable')]
            n_top = int(np.ceil(fraction * len(group)))
            chosen = order[:n_top]
            flags[chosen[rate[chosen] > 0]] = True
        out[gene] = flags
    return out
