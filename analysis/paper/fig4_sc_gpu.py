"""fig4_sc.py on the GPU: the same COAST steps, with PCA, Harmony, the neighbour graph, UMAP and Leiden computed by
rapids-singlecell instead of scanpy / harmonypy. Everything else (gene filter, normalisation, variable genes per
condition, marker labels, nephron selection, the ElPiGraph curve, tubule labels) is COAST and fig4_sc.py unchanged.
Numerics differ slightly from the CPU path (GPU PCA and Harmony, RAPIDS UMAP and Leiden); runs in the rapids-sc env.

  python analysis/paper/fig4_sc_gpu.py [--runs aki,hvm] [--out paper_figures/fig4/cells_gpu] [fig4_sc.py flags]
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))


def patch():
    """Route the five heavy scanpy / harmonypy calls COAST makes to rapids-singlecell, same arguments."""
    import anndata as ad
    import cupy as cp
    import harmonypy
    import rapids_singlecell as rsc
    import scanpy as sc

    def pca(adata, n_comps=50, mask_var=None, random_state=0, **_):
        cols = adata.var[mask_var].to_numpy(bool) if isinstance(mask_var, str) else slice(None)
        tmp = ad.AnnData(adata.X[:, cols])
        rsc.get.anndata_to_GPU(tmp)
        rsc.pp.pca(tmp, n_comps=n_comps, random_state=random_state)
        adata.obsm['X_pca'] = np.asarray(tmp.obsm['X_pca'].get() if hasattr(tmp.obsm['X_pca'], 'get')
                                         else tmp.obsm['X_pca'])
        del tmp
        cp.get_default_memory_pool().free_all_blocks()

    class _Harmony:
        def __init__(self, z):
            self.Z_corr = z

    def run_harmony(data, meta, key, theta=2, lamb=1, max_iter_harmony=10, random_state=0, **_):
        tmp = ad.AnnData(obs=meta[[key]].copy(), obsm={'X_pca': np.asarray(data, dtype=np.float32)})
        rsc.pp.harmony_integrate(tmp, key=key, basis='X_pca', theta=theta, ridge_lambda=lamb,
                                 max_iter_harmony=max_iter_harmony, random_state=random_state)
        z = tmp.obsm['X_pca_harmony']
        return _Harmony(np.asarray(z.get() if hasattr(z, 'get') else z, dtype=float))

    def neighbors(adata, n_neighbors=15, use_rep=None, random_state=0, **_):
        rsc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep=use_rep, random_state=random_state)

    def umap(adata, min_dist=0.5, random_state=0, **_):
        rsc.tl.umap(adata, min_dist=min_dist, random_state=random_state)

    def leiden(adata, resolution=1.0, random_state=0, key_added='leiden', n_iterations=-1, **_):
        rsc.tl.leiden(adata, resolution=resolution, random_state=random_state, key_added=key_added,
                      n_iterations=n_iterations if n_iterations > 0 else 100)

    sc.pp.pca, sc.pp.neighbors, sc.tl.umap, sc.tl.leiden = pca, neighbors, umap, leiden
    harmonypy.run_harmony = run_harmony


if __name__ == '__main__':
    patch()
    import fig4_sc
    argv = sys.argv[1:]
    if not any(a.startswith('--out') for a in argv):
        argv += ['--out', 'paper_figures/fig4/cells_gpu']
    fig4_sc.main(argv)
