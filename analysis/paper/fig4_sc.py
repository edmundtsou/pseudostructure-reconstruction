"""Fig. 4c: COAST on single-nucleus RNA-seq of mouse PT, without tissue positions.

Input: CELLxGENE Census PT nuclei (data/external/census_pt_cells/mouse_pt_cells.h5ad; S1/S2/S3 cell types from the
authors), restricted to the adult male mice (3-6 months), matching the two male control mice. COAST runs as on the
tubules: normalise, Harmony over donors, elastic curve oriented from the PT-S1 to the PT-S3 marker panel of
analysis/package_configs/aki/marker_panels.yaml. The Fig. 4b genes are fitted along it per donor with COAST's
penalised NB GAM (coast.tl.fit_curves). Cyp4a12a has no 1:1 human ortholog and is not in the Census panel.

  python analysis/paper/fig4_sc.py [--data-root <dir>] [--results-root <dir>]

Writes <results-root>/paper_figures/fig4/sc/: nuclei.parquet (donor, label, t, UMAP), nodes.npy (curve in UMAP),
curves.csv (donor x gene x position), summary.json (nuclei per donor and label, rho of t with the label order).
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

REPO = next(p for p in Path(__file__).resolve().parents if (p / 'pseudospace').is_dir())
STAGE, SEX = '3-6 month-old mature stage', 'male'
LABEL = {'epithelial cell of proximal tubule segment 1': 'PTS1', 'epithelial cell of proximal tubule segment 2': 'PTS2',
         'epithelial cell of proximal tubule segment 3': 'PTS3'}
GENES = ['Pigr', 'Slc22a7']
SEED = 0


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--data-root', default=os.environ.get('PSEUDOSPACE_DATA_ROOT', str(REPO / 'data')))
    p.add_argument('--results-root', default=os.environ.get('PSEUDOSPACE_RESULTS_ROOT', str(REPO / 'results')))
    p.add_argument('--jobs', type=int, default=8)
    return p.parse_args(argv)


def main(argv=None):
    import anndata as ad
    import coast
    import yaml
    from scipy.spatial import cKDTree
    from scipy.stats import spearmanr

    args = parse_args(argv)
    out = Path(args.results_root) / 'paper_figures/fig4/sc'
    out.mkdir(parents=True, exist_ok=True)
    P = yaml.safe_load((REPO / 'analysis/package_configs/aki/marker_panels.yaml').read_text())
    early, late = P['label_panels']['PT-S1'], P['label_panels']['PT-S3']

    a = ad.read_h5ad(Path(args.data_root) / 'external/census_pt_cells/mouse_pt_cells.h5ad')
    a = a[(a.obs.development_stage == STAGE) & (a.obs.sex == SEX)].copy()
    a.obs['label'] = a.obs.cell_type.map(LABEL).astype('category')
    a.obs['donor'] = a.obs.donor_id.astype(str).astype('category')
    a = coast.pp.filter_genes(a, keep=GENES + early + late)
    coast.pp.normalize(a)
    coast.pp.integrate(a, batch_key='donor', random_state=SEED)
    coast.tl.pseudostructure(a, early_genes=early, late_genes=late, seed=SEED)
    t = a.obs['pseudostructure'].to_numpy(float)

    # curve nodes placed in UMAP as coast.pl.curve does: median UMAP position of the nuclei nearest each node
    fitted, nodes = a.obsm[a.uns['coast']['use_rep']], a.uns['coast']['nodes']
    nearest = cKDTree(nodes).query(fitted)[1]
    umap = a.obsm['X_umap']
    np.save(out / 'nodes.npy', np.array([np.median(umap[nearest == k], 0) for k in range(len(nodes))
                                         if (nearest == k).any()]))
    pd.DataFrame({'donor': a.obs.donor.astype(str), 'label': a.obs.label.astype(str), 't': t,
                  'umap1': umap[:, 0], 'umap2': umap[:, 1]}, index=a.obs_names).to_parquet(out / 'nuclei.parquet')

    curves = coast.tl.fit_curves(a, GENES, sample_key='donor', limits=(0, 1), n_points=101, n_jobs=args.jobs)
    curves.to_csv(out / 'curves.csv')

    rank = a.obs.label.map({'PTS1': 0, 'PTS2': 1, 'PTS3': 2}).astype(float).to_numpy()
    rho = {d: float(spearmanr(t[k], rank[k]).statistic) for d in a.obs.donor.cat.categories
           if (k := (a.obs.donor == d).to_numpy()).any()}
    summary = {'n_nuclei': int(a.n_obs), 'n_genes': int(a.n_vars),
               'per_donor_label': pd.crosstab(a.obs.donor, a.obs.label).to_dict(orient='index'),
               'rho_by_donor': rho,
               'median_t': a.obs.groupby('label', observed=True)['pseudostructure'].median().round(3).to_dict()}
    (out / 'summary.json').write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))


if __name__ == '__main__':
    main()
