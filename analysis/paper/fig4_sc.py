"""Fig. 4c: COAST on single cells of our own Visium HD sections, the same method as on tubules with cells as units.

Cells are bin2cell reconstructions (StarDist nuclei on the H&E, expanded by <= 4 bins; analysis/bin2cell/bin2cell_run.py,
the same parameters for every section). Mouse: the two control sections after the cell QC of the bin2cell bundle
(nucleus-anchored, >= 200 UMI, >= 150 genes, 3 MADs of log1p counts and genes, <= 25 % mitochondrial, 10-600 um^2).
Human: both healthy-cortex slices of the one donor, through the same QC here, with genes renamed to their shared 1:1
mouse orthologs as in the hvm tubule fit. Each species is fitted on its own by coast_fit.two_pass, unchanged from the
tubule fits: pass 1 labels every Leiden cluster with the marker panels, pass 2 refits the nephron cells, and the
elastic curve is oriented from the early to the late anchor genes. PT cells (pass-1 labels PT-S1/S2/S3) are then
rank-scaled to [0, 1] along the curve, as the tubule PT axis of Fig. 4 is, and the Fig. 4b genes are fitted along it with
COAST's penalised NB GAM (coast.tl.fit_curves). For mouse, the bundle's own marker cell types (PT_S1/S2/S3 ...)
are carried along as an independent check of the COAST cell labels.

  python analysis/paper/fig4_sc.py [--data-root <dir>] [--results-root <dir>] [--species mouse,human]

Reads <data-root>/bin2cell/mouse_bundle/02_qc_integration_annotation_composition/checkpoints/ (mouse) and
<results-root>/bin2cell/<slice>/<slice>_bin2cell_cells_unfiltered.h5ad (human). Writes
<results-root>/paper_figures/fig4/cells/<species>/: cells.parquet (sample, labels, t, PT axis, UMAP), nodes.npy (curve
in UMAP), curves.csv (sample x gene x PT position), summary.json.
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
MOUSE = {'Ctrl_1A2': 'Ctrl1A2', 'Ctrl_1A4': 'Ctrl1A4'}
HUMAN = {'HUK1_COR_1': 'HUK1_COR1', 'HUK1_MED_1': 'HUK1_MED1'}
BUNDLE = 'bin2cell/mouse_bundle/02_qc_integration_annotation_composition/checkpoints'
PT = ['PT-S1', 'PT-S2', 'PT-S3']
GENES = ['Pigr', 'Slc22a7', 'Cyp4a12a']
# cell QC of the bin2cell bundle (run_config_02_qc.json)
MIN_COUNTS, MIN_GENES, N_MADS, MAX_PCT_MT, AREA_UM2 = 200, 150, 3.0, 25.0, (10.0, 600.0)
SEED = 0


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--data-root', default=os.environ.get('PSEUDOSPACE_DATA_ROOT', str(REPO / 'data')))
    p.add_argument('--results-root', default=os.environ.get('PSEUDOSPACE_RESULTS_ROOT', str(REPO / 'results')))
    p.add_argument('--species', default='mouse,human')
    p.add_argument('--jobs', type=int, default=8)
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


def load_mouse(data_root: Path):
    import anndata as ad
    parts = {}
    for raw, s in MOUSE.items():
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


def fit_species(a, run: str, out: Path, jobs: int, extra: dict):
    import coast
    import scanpy as sc
    from scipy.spatial import cKDTree
    from scipy.stats import spearmanr
    from coast_fit import _pass, exclude_regex, panels, two_pass

    P = panels(run)
    a.layers['counts'] = a.X.copy()
    a.obs['sample'] = a.obs['sample'].astype('category')
    a.obs['condition'] = pd.Categorical(['all'] * a.n_obs)
    p1, p2, _ = two_pass(a, P)
    # the whole-nephron curve, recorded only: on cells it does not order the nephron (2026-11-15, mouse rho 0.45-0.54,
    # PT running S3 -> S1), so the panel uses a PT curve, as the earlier single-nucleus panel did
    order = {l: i for i, l in enumerate(P['segment_order'])}
    t_neph = p2.obs['pseudostructure'].to_numpy(float)
    rank_neph = p2.obs['segment_class'].astype(str).map(order).to_numpy(float)
    s_neph = p2.obs['sample'].astype(str).to_numpy()
    rho_neph = {s: float(spearmanr(t_neph[(s_neph == s) & np.isfinite(rank_neph)],
                                   rank_neph[(s_neph == s) & np.isfinite(rank_neph)]).statistic) for s in np.unique(s_neph)}

    # PT pass: the pass-1 PT cells through the same COAST steps, curve oriented from the PT-S1 to the PT-S3 panel
    pt_cells = p1.obs_names[p1.obs['segment_class'].astype(str).isin(PT)]
    units = a[pt_cells].copy()
    units.X = units.layers['counts'].copy()
    p3 = _pass(units, P, 'condition', exclude_regex(), 'leiden_pt')
    p3.obs['segment_class'] = p1.obs.loc[p3.obs_names, 'segment_class'].astype(str).astype('category')
    coast.tl.pseudostructure(p3, early_genes=P['label_panels']['PT-S1'], late_genes=P['label_panels']['PT-S3'],
                             seed=SEED)
    t = p3.obs['pseudostructure'].to_numpy(float)
    label = p3.obs['segment_class'].astype(str).to_numpy()
    sample = p3.obs['sample'].astype(str).to_numpy()

    sc.pp.neighbors(p3, use_rep='X_harmony', n_neighbors=15, random_state=SEED)
    sc.tl.umap(p3, random_state=SEED)
    umap = p3.obsm['X_umap']
    nearest = cKDTree(p3.uns['coast']['nodes']).query(p3.obsm['X_harmony'])[1]
    np.save(out / 'nodes.npy', np.array([np.median(umap[nearest == k], 0) for k in range(len(p3.uns['coast']['nodes']))
                                         if (nearest == k).any()]))

    # PT axis rank-scaled to [0, 1], as the tubule PT axis of Fig. 4
    x = (pd.Series(t).rank(method='average').to_numpy() - 0.5) / len(t)
    cells = pd.DataFrame({'sample': sample, 'label': label, 'bundle_type': p3.obs['bundle_type'].astype(str).to_numpy(),
                          't': t, 'x_pt': x, 'umap1': umap[:, 0], 'umap2': umap[:, 1]}, index=p3.obs_names)
    cells.to_parquet(out / 'cells.parquet')

    src = a[p3.obs_names]
    genes = [g for g in GENES if g in src.var_names]
    mt = src.var_names.str.match(r'^mt-', case=False)
    b = src[:, genes].copy()
    b.obs['x'] = x
    b.obs['lib'] = np.asarray(src.layers['counts'][:, ~mt].sum(1)).ravel()
    curves = coast.tl.fit_curves(b, genes, x='x', sample_key='sample', library='lib', limits=(0, 1), n_points=201,
                                 n_jobs=jobs)
    curves.to_csv(out / 'curves.csv')

    rank = pd.Series(label).map({l: i for i, l in enumerate(PT)}).to_numpy(float)
    pt = np.ones(len(t), bool)
    summary = {
        'n_cells_pass1': int(p1.n_obs), 'n_cells_nephron': int(p2.n_obs), 'n_pt_cells': int(len(t)),
        'per_sample_label': pd.crosstab(cells['sample'], cells['label']).to_dict(orient='index'),
        'rho_pt_by_sample': {s: float(spearmanr(t[sample == s], rank[sample == s]).statistic)
                             for s in np.unique(sample)},
        'whole_nephron_curve_rho_by_sample': rho_neph,
        'median_x_pt': cells.groupby('label')['x_pt'].median().round(3).to_dict(), **extra}
    if (cells.bundle_type != 'NA').any():   # COAST cluster labels against the bundle's own cell types (PT cells)
        summary['bundle_vs_coast_pt'] = pd.crosstab(cells.label[pt], cells.bundle_type[pt]).to_dict(orient='index')
    (out / 'summary.json').write_text(json.dumps(summary, indent=1))
    print(json.dumps({k: v for k, v in summary.items() if k != 'per_sample_label'}, indent=1), flush=True)


def tubule_labels(cells: pd.DataFrame, species: str, data_root: Path, results_root: Path) -> pd.Series:
    """COAST segment label of the tubule each cell's nucleus sits in (NaN outside a COAST tubule). Mouse: the bundle's
    nucleus -> v4fix polygon assignment; human: nucleus centroid in a HUK1 v2 polygon (both in full-res H&E pixels).
    A unit <sample>_unit_N of the tubule fits is polygon N of its segmentation."""
    import anndata as ad
    run = 'aki' if species == 'mouse' else 'hvm'
    units = pd.read_parquet(results_root / 'coast_runs' / run / '04_pseudospace/pseudospace.parquet')['segment_class']
    sample = cells['sample'].astype(str)
    if species == 'mouse':
        raw = {v: k for k, v in MOUSE.items()}
        uid = [n.removesuffix('_' + s) + '-' + raw[s] for n, s in zip(cells.index, sample)]
        c2t = pd.read_parquet(data_root / BUNDLE / 'cell_to_tubule.parquet', columns=['cell_uid', 'feature_index'])
        fi = pd.Series(uid).map(c2t.set_index('cell_uid')['feature_index']).to_numpy()
    else:
        from shapely import STRtree, points
        from shapely.geometry import shape
        fi = np.full(len(cells), np.nan)
        for raw, s in HUMAN.items():
            k = (sample == s).to_numpy()
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


def label_summary(out: Path, species: str, data_root: Path, results_root: Path):
    """Add the tubule labels to cells.parquet and the PT order of the cell axis against them to summary.json."""
    from scipy.stats import spearmanr
    cells = pd.read_parquet(out / 'cells.parquet')
    cells['tubule_label'] = tubule_labels(cells, species, data_root, results_root)
    cells.to_parquet(out / 'cells.parquet')
    rank = cells.tubule_label.map({l: i for i, l in enumerate(PT)})
    ok = rank.notna()
    summ = json.loads((out / 'summary.json').read_text())
    summ['cells_in_pt_tubule'] = cells[ok].groupby('sample').size().to_dict()
    summ['cells_by_tubule_label'] = cells.tubule_label.fillna('none').value_counts().to_dict()
    summ['rho_tubule_pt_by_sample'] = {s: float(spearmanr(g.x_pt, rank[g.index]).statistic)
                                      for s, g in cells[ok].groupby('sample')}
    summ['median_x_pt_by_tubule_label'] = cells[ok].groupby('tubule_label')['x_pt'].median().round(3).to_dict()
    (out / 'summary.json').write_text(json.dumps(summ, indent=1))
    print(species, {k: summ[k] for k in ('cells_in_pt_tubule', 'rho_tubule_pt_by_sample', 'median_x_pt_by_tubule_label')},
          flush=True)


def main(argv=None):
    args = parse_args(argv)
    root = Path(args.results_root) / 'paper_figures/fig4/cells'
    for species in args.species.split(','):
        out = root / species
        out.mkdir(parents=True, exist_ok=True)
        if not args.labels_only:
            if species == 'mouse':
                fit_species(load_mouse(Path(args.data_root)), 'aki', out, args.jobs, {})
            else:
                a, ledger = load_human(Path(args.results_root))
                fit_species(a, 'hvm', out, args.jobs, {'qc': ledger})
        label_summary(out, species, Path(args.data_root), Path(args.results_root))


if __name__ == '__main__':
    main()
