"""Describe Fig. 5 urate transporter profiles in the study's Visium HD sections.

The inputs and interpretation are documented in docs/paper/fig5-significance-outline.md.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
GENES = ['Slc17a3', 'Slc22a12', 'Slc17a1']
SECTIONS = ['Ctrl1A2', 'Ctrl1A4', 'HUK1_COR1', 'HUK1_MED1']


def section_profiles(screen: Path) -> pd.DataFrame:
    """Positions and correlations for the exporter and reabsorber on the single COAST axis."""
    genes = pd.read_csv(screen / 'genes.csv', index_col='gene')
    axes = pd.read_csv(screen / 'axes.csv').squeeze(axis=0)
    rows = []
    # The existing screen writes gene symbols as a NumPy object array.
    with np.load(screen / 'all_curves.npz', allow_pickle=True) as curves:
        grid = curves['grid']
        idx = pd.Index(curves['genes']).get_indexer(GENES)
        if (idx < 0).any():
            raise ValueError('A target gene is missing from the Fig. 5 curve archive')
        for section in SECTIONS:
            fitted = curves[f'curve|{section}'][idx]
            bins = curves[f'bins|{section}'][idx[:2]]
            supported = ((grid >= axes[f'x_lo_{section}']) & (grid <= axes[f'x_hi_{section}'])
                         & np.isfinite(fitted[:2]).all(axis=0))
            paired_curves = fitted[:2, supported]
            raw = genes.loc[GENES, f'R_{section}']
            rows.append({
                'section': section,
                **raw.to_dict(),
                'NPT4 peak': grid[supported][paired_curves[0].argmax()],
                'URAT1 peak': grid[supported][paired_curves[1].argmax()],
                'curve r': np.corrcoef(paired_curves)[0, 1],
                'bin r': np.corrcoef(bins[:, np.isfinite(bins).all(axis=0)])[0, 1],
            })
    return pd.DataFrame(rows).set_index('section')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results-root', type=Path,
                        default=os.environ.get('PSEUDOSPACE_RESULTS_ROOT', REPO / 'results'))
    args = parser.parse_args()

    print('COAST raw endpoint contrasts, supported maxima and profile correlations')
    print(section_profiles(args.results_root / 'paper_figures/fig5/screen').round(3).to_string())


if __name__ == '__main__':
    main()
