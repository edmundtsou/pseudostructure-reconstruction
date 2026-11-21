"""Describe Fig. 5 urate transporter profiles and check their directions in public data.

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


def lake_endpoint_contrasts(cache: Path, min_nuclei: int) -> pd.DataFrame:
    """Paired S3 minus S1 log2(CPM + 1) in the normal cortical donor groups from notebook 38."""
    with np.load(cache) as lake:
        meta = pd.DataFrame(
            [group.split('|') for group in lake['groups']],
            columns=['donor', 'sex', 'segment', 'region'],
        )
        keep = meta.region.eq('cortex') & (lake['counts'] >= min_nuclei)
        rates = pd.DataFrame(lake['sums'], columns=lake['genes'])[GENES]
        rates = np.log2(rates.div(lake['totals'], axis=0) * 1e6 + 1)
        by_segment = rates.loc[keep].set_index([meta.donor[keep], meta.segment[keep]]).unstack('segment')
    return (by_segment.xs('S3', axis=1, level='segment')
            - by_segment.xs('S1', axis=1, level='segment')).dropna()


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
    parser.add_argument('--data-root', type=Path, default=os.environ.get('PSEUDOSPACE_DATA_ROOT', REPO / 'data'))
    parser.add_argument('--results-root', type=Path,
                        default=os.environ.get('PSEUDOSPACE_RESULTS_ROOT', REPO / 'results'))
    parser.add_argument('--lake-pseudobulk', type=Path, required=True, help='notebook 38 donor pseudobulk NPZ')
    args = parser.parse_args()

    print('COAST raw endpoint contrasts, supported maxima and profile correlations')
    print(section_profiles(args.results_root / 'paper_figures/fig5/screen').round(3).to_string())
    for minimum in (50, 25):
        contrast = lake_endpoint_contrasts(args.lake_pseudobulk, minimum)
        print(f'\nLake normal cortical donors with at least {minimum} nuclei in each endpoint (n = {len(contrast)})')
        print(pd.DataFrame({'median S3−S1': contrast.median().round(3), 'negative': (contrast < 0).sum(),
                            'positive': (contrast > 0).sum()}).to_string())

    mouse = pd.read_csv(args.data_root / 'external/segment_microdissection_sex/GSE212213_TPM_Male_Female.txt.gz',
                        sep='\t', index_col='Gene_symbol')
    print('\nChen microdissected mouse PT, mean TPM by sex and segment')
    print(mouse.loc[GENES, ['PTS1', 'PTS2', 'PTS3', 'fPTS1', 'fPTS2', 'fPTS3']].to_string())


if __name__ == '__main__':
    main()
