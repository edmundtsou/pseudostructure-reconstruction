"""Fig. 4c: why tubules. COAST run on single cells of our own Visium HD sections, exactly as on tubules.

Cells are bin2cell reconstructions (StarDist nuclei on the H&E, expanded by <= 4 bins; analysis/bin2cell/bin2cell_run.py,
the same parameters for every section). Mouse: the two control sections after the cell QC of the bin2cell bundle
(nucleus-anchored, >= 200 UMI, >= 150 genes, 3 MADs of log1p counts and genes, <= 25 % mitochondrial, 10-600 um^2).
Human: both healthy-cortex slices of the one donor, through the same QC here, genes renamed to their shared 1:1 mouse
orthologs as in the hvm tubule fit. Each species is fitted on its own by coast_fit.two_pass, unchanged from the tubule
fits (gene filter, normalisation, Harmony, Leiden + marker labels, nephron units refitted, ElPiGraph curve oriented
from the early to the late anchors). Cells have no polygon, so the tubule doublet step has no cell counterpart.
The cell axis and the tubule axis are scored by the Fig. 3 measure (per-section Spearman rho of position against the
nephron rank of the label group), cells by their own label and by the label of the tubule their nucleus sits in.

  python analysis/paper/fig4_sc.py [--data-root <dir>] [--results-root <dir>] [--species mouse,human] [--labels-only]

Reads <data-root>/bin2cell/mouse_bundle/02_qc_integration_annotation_composition/checkpoints/ (mouse) and
<results-root>/bin2cell/<slice>/<slice>_bin2cell_cells_unfiltered.h5ad (human). Writes
<results-root>/paper_figures/fig4/cells/<species>/: cells.parquet (sample, labels, t, UMAP), nodes.npy (curve in UMAP),
summary.json (cell counts, rho of cells and tubules per section).
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = next(p for p in Path(__file__).resolve().parents if (p / 'pseudospace').is_dir())
sys.path.insert(0, str(Path(__file__).resolve().parent))
MOUSE = {'Ctrl_1A2': 'Ctrl1A2', 'Ctrl_1A4': 'Ctrl1A4', 'IR_2A2': 'IR2A2', 'IR_2A4': 'IR2A4'}
CTRL = ['Ctrl1A2', 'Ctrl1A4']
HUMAN = {'HUK1_COR_1': 'HUK1_COR1', 'HUK1_MED_1': 'HUK1_MED1'}
BUNDLE = 'bin2cell/mouse_bundle/02_qc_integration_annotation_composition/checkpoints'
# cell QC of the bin2cell bundle (run_config_02_qc.json)
MIN_COUNTS, MIN_GENES, N_MADS, MAX_PCT_MT, AREA_UM2 = 200, 150, 3.0, 25.0, (10.0, 600.0)
SEED = 0


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--data-root', default=os.environ.get('PSEUDOSPACE_DATA_ROOT', str(REPO / 'data')))
    p.add_argument('--results-root', default=os.environ.get('PSEUDOSPACE_RESULTS_ROOT', str(REPO / 'results')))
    p.add_argument('--runs', default='aki,hvm')
    p.add_argument('--out', default='paper_figures/fig4/cells', help='relative to the results root')
    p.add_argument('--labels-only', action='store_true', help='only add the tubule labels to an existing fit')
    return p.parse_args(argv)


def cell_qc(a):
    """Boolean mask of the bundle's cell gates on raw counts (a: one section, nucleus-anchored cells)."""
    from scipy import sparse
    X = sparse.csr_matrix(a.layers['counts'])
    total = np.asarray(X.sum(1)).ravel()
    genes = np.asarray((X > 0).sum(1)).ravel()
    mt = a.var_names.str.match(r'^mt-', case=False)
    pct_mt = 100 * np.asarray(X[:, mt].sum(1)).ravel() / np.maximum(total, 1)
    keep = (total >= MIN_COUNTS) & (genes >= MIN_GENES)
    for v in (np.log1p(total), np.log1p(genes)):   # MADs over all nucleus-anchored cells, as in the bundle
        med = np.median(v)
        mad = np.median(np.abs(v - med))
        keep &= (v >= med - N_MADS * mad) & (v <= med + N_MADS * mad)
    if mt.sum() >= 10:
        keep &= pct_mt <= MAX_PCT_MT
    area = a.obs['cell_area_um2'].to_numpy(float)
    return keep & (area >= AREA_UM2[0]) & (area <= AREA_UM2[1])


def load_mouse(data_root: Path, sections):
    import anndata as ad
    parts = {}
    for raw, s in MOUSE.items():
        if s not in sections:
            continue
        a = ad.read_h5ad(data_root / BUNDLE / f'{raw}_filtered.h5ad')
        a = ad.AnnData(a.layers['counts'].astype(np.float32), obs=a.obs[['cell_area_um2']].copy(), var=a.var[[]])
        a.obs['cell_uid'] = a.obs_names + '-' + raw
        parts[s] = a
    a = ad.concat(parts, label='sample', index_unique='_')
    meta = pd.read_parquet(data_root / BUNDLE / 'cell_metadata.parquet', columns=['cell_uid', 'cell_type'])
    a.obs['bundle_type'] = a.obs['cell_uid'].map(meta.set_index('cell_uid')['cell_type']).astype(str)
    return a


def load_human(results_root: Path):
    import anndata as ad
    shared = pd.read_csv(results_root / 'package_pseudospace/runs/hvm/00_ortholog_space/shared_orthologs.csv')
    parts, ledger = {}, {}
    for raw, s in HUMAN.items():
        c = ad.read_h5ad(results_root / 'bin2cell' / raw / f'{raw}_bin2cell_cells_unfiltered.h5ad')
        c.var_names_make_unique()
        c = c[c.obs['has_nucleus'].to_numpy(bool)].copy()
        keep = cell_qc(c)
        ledger[s] = {'nucleus_anchored': int(c.n_obs), 'after_qc': int(keep.sum())}
        c = c[keep]
        have = shared[shared.human_source.isin(c.var_names)]
        a = ad.AnnData(c[:, have.human_source.to_numpy()].layers['counts'].astype(np.float32),
                       obs=c.obs[['cell_area_um2']].copy(), var=pd.DataFrame(index=have.mouse_symbol.to_numpy()))
        parts[s] = a
    a = ad.concat(parts, label='sample', index_unique='_')
    a.obs['bundle_type'] = 'NA'
    return a, ledger


def load_run(run: str, data_root: Path, results_root: Path):
    """The cells of one tubule fit, in that fit's gene space and conditions: aki = the four mouse sections (Control,
    IR); hvm = the two control mice and the two human slices on the shared 1:1 orthologs (mouse, human)."""
    import anndata as ad
    groups = pd.read_csv(REPO / 'analysis/package_configs' / {'aki': 'aki', 'hvm': 'human_vs_mouse'}[run] / 'samples.csv',
                         comment='#').set_index('sample')['group']
    if run == 'aki':
        a, extra = load_mouse(data_root, list(groups.index)), {}
    else:
        m = load_mouse(data_root, CTRL)
        h, ledger = load_human(results_root)
        a, extra = ad.concat({'m': m, 'h': h}, join='inner', index_unique=None), {'qc': ledger}
    a.obs['condition'] = a.obs['sample'].astype(str).map(groups).astype('category')
    return a, extra


def fit_run(a, run: str, out: Path, extra: dict):
    """COAST on cells exactly as on tubules (coast_fit.two_pass with the Fig. 1 pre-Harmony UMAP: filter genes,
    normalise, Harmony on section with variable genes per condition, Leiden + marker labels, nephron cells refitted,
    ElPiGraph curve oriented early -> late). Saves pass-1 embeddings (before / after Harmony) and the pass-2 cells."""
    from scipy.spatial import cKDTree
    from coast_fit import panels, two_pass

    P = panels(run)
    a.layers['counts'] = a.X.copy()
    a.obs['sample'] = a.obs['sample'].astype('category')
    p1, p2, _ = two_pass(a, P, pre_umap=True)
    pre, post = p1.obsm['X_umap_pre_harmony'], p1.obsm['X_umap']
    pd.DataFrame({'sample': p1.obs['sample'].astype(str), 'condition': p1.obs['condition'].astype(str),
                  'label': p1.obs['segment_class'].astype(str), 'pre1': pre[:, 0], 'pre2': pre[:, 1],
                  'post1': post[:, 0], 'post2': post[:, 1]}, index=p1.obs_names.copy()).to_parquet(out / 'pass1.parquet')
    np.save(out / 'pass1_pca.npy', p1.obsm['X_pca'][:, :50].astype(np.float32))
    np.save(out / 'pass1_harmony.npy', p1.obsm['X_harmony'][:, :50].astype(np.float32))
    umap, nodes = p2.obsm['X_umap'], p2.uns['coast']['nodes']
    nearest = cKDTree(nodes).query(p2.obsm['X_harmony'])[1]
    np.save(out / 'nodes.npy', np.array([np.median(umap[nearest == k], 0) for k in range(len(nodes))
                                         if (nearest == k).any()]))
    cells = pd.DataFrame({'sample': p2.obs['sample'].astype(str), 'condition': p2.obs['condition'].astype(str),
                          'label': p2.obs['segment_class'].astype(str), 'bundle_type': p2.obs['bundle_type'].astype(str),
                          't': p2.obs['pseudostructure'].astype(float), 'umap1': umap[:, 0], 'umap2': umap[:, 1]},
                         index=p2.obs_names.copy())
    cells.to_parquet(out / 'cells.parquet')
    (out / 'summary.json').write_text(json.dumps({'n_cells_qc': a.obs['sample'].value_counts().to_dict(),
                                                  'n_cells_nephron': cells['sample'].value_counts().to_dict(), **extra},
                                                 indent=1))


def tubule_labels(cells: pd.DataFrame, run: str, data_root: Path, results_root: Path) -> pd.Series:
    """COAST segment label of the tubule each cell's nucleus sits in (NaN outside a COAST tubule). Mouse: the bundle's
    nucleus -> v4fix polygon assignment; human: nucleus centroid in a HUK1 v2 polygon (both in full-res H&E pixels).
    A unit <sample>_unit_N of the tubule fits is polygon N of its segmentation."""
    import anndata as ad
    units = pd.read_parquet(results_root / 'coast_runs' / run / '04_pseudospace/pseudospace.parquet')['segment_class']
    sample = cells['sample'].astype(str)
    fi = np.full(len(cells), np.nan)
    mouse = sample.isin(MOUSE.values()).to_numpy()
    if mouse.any():
        raw = {v: k for k, v in MOUSE.items()}
        uid = [n.removesuffix('_' + s) + '-' + raw[s] for n, s in zip(cells.index[mouse], sample[mouse])]
        c2t = pd.read_parquet(data_root / BUNDLE / 'cell_to_tubule.parquet', columns=['cell_uid', 'feature_index'])
        fi[mouse] = pd.Series(uid).map(c2t.set_index('cell_uid')['feature_index']).to_numpy(float)
    if (~mouse).any():
        from shapely import STRtree, points
        from shapely.geometry import shape
        for raw, s in HUMAN.items():
            k = (sample == s).to_numpy()
            if not k.any():
                continue
            obs = ad.read_h5ad(results_root / 'bin2cell' / raw / f'{raw}_bin2cell_cells_unfiltered.h5ad', backed='r').obs
            names = cells.index[k].str.removesuffix('_' + s)
            xy = obs.loc[names, ['nucleus_x_fullres_px', 'nucleus_y_fullres_px']].to_numpy(float)
            feats = json.loads((data_root / f'{s}_v2.geojson').read_text())['features']
            tree = STRtree([shape(f['geometry']) for f in feats])
            hit = tree.query(points(xy), predicate='within')   # (point, polygon) pairs
            first = pd.Series(hit[1]).groupby(hit[0]).min()     # a nucleus on two touching polygons takes the first
            v = np.full(k.sum(), np.nan)
            v[first.index.to_numpy()] = first.to_numpy()
            fi[k] = v
    name = [f'{s}_unit_{int(f)}' if np.isfinite(f) else '' for s, f in zip(sample, np.asarray(fi, float))]
    return pd.Series(pd.Series(name).map(units).to_numpy(), index=cells.index, name='tubule_label')


def label_summary(out: Path, run: str, data_root: Path, results_root: Path):
    """Order of the cell axis and of the tubule axis, by the Fig. 3 measure: per section, Spearman rho of position
    against the nephron rank of the label group (PTS1 ... CD). Cells are ranked by their own COAST label and by the
    label of the tubule their nucleus sits in; tubules by their own label."""
    from fig3_common import GROUP, RANK, load_positions, section_rho
    cells = pd.read_parquet(out / 'cells.parquet')
    cells['tubule_label'] = tubule_labels(cells, run, data_root, results_root)
    cells['group'] = cells.label.map(GROUP)
    cells['tubule_group'] = cells.tubule_label.map(GROUP)
    cells.to_parquet(out / 'cells.parquet')
    tub = load_positions(results_root, run)
    secs = sorted(cells['sample'].unique())
    tub = tub[tub['sample'].isin(secs)]
    rho = lambda coord, table: section_rho(coord, table, secs).set_index('sample')['rho'].round(3).to_dict()
    summ = json.loads((out / 'summary.json').read_text())
    summ.update({
        'rho_cells_own_label': rho(cells.t, cells.assign(rank=cells.group.map(RANK))),
        'rho_cells_tubule_label': rho(cells.t, cells.assign(rank=cells.tubule_group.map(RANK))),
        'rho_tubules': rho(tub.t, tub),
        'cells_in_labelled_tubule': cells.tubule_group.notna().groupby(cells['sample']).sum().astype(int).to_dict(),
        'cells_by_group': cells.group.value_counts().to_dict(),
        'median_t_cells_by_tubule_group': cells.groupby('tubule_group')['t'].median().round(3).to_dict(),
        'median_t_tubules_by_group': tub.groupby('group')['t'].median().round(3).to_dict()})
    if (cells.bundle_type != 'NA').any():
        summ['bundle_vs_coast'] = pd.crosstab(cells.label, cells.bundle_type).to_dict(orient='index')
    (out / 'summary.json').write_text(json.dumps(summ, indent=1))
    print(run, {k: v for k, v in summ.items() if k.startswith(('rho', 'median', 'cells_in'))}, flush=True)


def main(argv=None):
    args = parse_args(argv)
    root = Path(args.results_root) / args.out
    for run in args.runs.split(','):
        out = root / run
        out.mkdir(parents=True, exist_ok=True)
        if not args.labels_only:
            a, extra = load_run(run, Path(args.data_root), Path(args.results_root))
            fit_run(a, run, out, extra)
        label_summary(out, run, Path(args.data_root), Path(args.results_root))


if __name__ == '__main__':
    main()
