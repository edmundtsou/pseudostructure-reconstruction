"""Load the matched, four-section Fig. 5 changed-gene profiles."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from fig5_urate import SECTIONS


def load_global_shapes(fig5_root: Path) -> tuple[np.ndarray, pd.DataFrame, np.ndarray]:
    """Return shared x, changed-gene metadata, and matched scaled curves (section, gene, x).

    Profiles are log2(CP10k + 0.5), centered separately within each gene and section,
    then divided by the largest section SD for that gene. Thus all four profiles share
    one scale while retaining section-specific shape and relative amplitude.
    """
    fig5_root = Path(fig5_root)
    shape_dir = fig5_root / 'shape_clusters'
    screen_dir = fig5_root / 'screen'
    metadata = pd.read_csv(shape_dir / 'genes.csv').set_index('gene')
    metadata = metadata.loc[metadata['changed'].astype(bool)].copy()
    metadata = metadata.sort_values(['cluster', 'rms_condition'], ascending=[True, False], kind='mergesort')
    names = metadata.index.to_numpy(str)

    with np.load(shape_dir / 'shapes.npz', allow_pickle=True) as saved, \
            np.load(screen_dir / 'all_curves.npz', allow_pickle=True) as fitted:
        x = np.asarray(saved['x'], dtype=float)
        fit_grid = np.asarray(fitted['grid'], dtype=float)
        if (x.ndim != 1 or x.size < 2 or fit_grid.ndim != 1 or not np.isfinite(x).all()
                or not (np.diff(x) > 0).all()):
            raise ValueError('The shape grid must contain at least two increasing finite positions')
        grid_idx = pd.Index(fit_grid).get_indexer(x)
        if (grid_idx < 0).any() or len(np.unique(grid_idx)) != len(grid_idx):
            raise ValueError('Saved shape grid does not exactly match unique screen grid positions')
        shape_names = np.asarray(saved['genes']).astype(str)
        fit_names = np.asarray(fitted['genes']).astype(str)
        shape_idx = pd.Index(shape_names).get_indexer(names)
        fit_idx = pd.Index(fit_names).get_indexer(names)
        if (shape_idx < 0).any() or (fit_idx < 0).any():
            raise ValueError('Changed genes do not exactly match both shape and screen archives')

        raw = []
        for section in SECTIONS:
            values = np.asarray(fitted[f'curve|{section}'][fit_idx][:, grid_idx], dtype=float)
            if not np.isfinite(values).all():
                raise ValueError(f'Nonfinite fitted values for changed genes in {section}')
            raw.append(np.logaddexp2(values, -1.0))
        raw = np.stack(raw)
        centered = raw - raw.mean(axis=2, keepdims=True)
        section_sd = centered.std(axis=2, ddof=0)
        scale = section_sd.max(axis=0)
        if not np.isfinite(scale).all() or (scale <= 0).any():
            raise ValueError('Every changed gene must have positive finite support SD')
        curves = centered / scale[None, :, None]

        # Ensure the archives came from the same fit, in both gene and coordinate order.
        if not np.array_equal(shape_names[shape_idx], names) or not np.array_equal(fit_names[fit_idx], names):
            raise ValueError('Gene order mismatch while matching changed-gene profiles')
        for i, section in enumerate(SECTIONS):
            z = np.asarray(saved[f'z|{section}'][shape_idx], dtype=float)
            if z.shape[1] != len(x):
                raise ValueError('Saved shape grid length mismatch')
            reconstructed = centered[i] / section_sd[i, :, None]
            np.testing.assert_allclose(reconstructed, z, rtol=1e-6, atol=1e-7,
                                       err_msg=f'{section} screen/shape values disagree')

    if not np.isfinite(curves).all():
        raise ValueError('Scaled profiles contain nonfinite values')
    return x, metadata, curves
