"""Render the Fig. 5 urate profiles and their external endpoint checks."""
from __future__ import annotations

import argparse
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.transforms import ScaledTranslation

from fig5_urate import GENES, REPO, SECTIONS, lake_endpoint_contrasts

COLORS = ('#17668C', '#AD481B')
STYLES = ('-', '--')
MARKERS = ('o', 'D')
SEGMENTS = ['PTS1', 'PTS2', 'PTS3']
WIDTH_MM = 183


def render(screen: Path, lake_cache: Path, data_root: Path):
    """Plot section-resolved COAST profiles, mouse segments and paired human donors."""
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 8,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.linewidth': 0.7,
        'pdf.fonttype': 42,
        'svg.fonttype': 'none',
        'legend.frameon': False,
    })
    fig, axs = plt.subplots(2, 2, figsize=(WIDTH_MM / 25.4, 129 / 25.4),
                            gridspec_kw={'height_ratios': [1.22, 1]})
    fig.subplots_adjust(left=0.11, right=0.977, top=0.85, bottom=0.12, wspace=0.33, hspace=0.69)
    fig.text(0.11, 0.974, 'Urate transporter expression along the proximal tubule',
             weight='bold', size=12, va='top')
    fig.legend([Line2D([], [], color=c, lw=2) for c in COLORS],
               ['NPT4 · urate export', 'URAT1 · urate reuptake'],
               loc='upper left', bbox_to_anchor=(0.106, 0.932), ncol=2, frameon=False)

    axes = pd.read_csv(screen / 'axes.csv').iloc[0]
    with np.load(screen / 'all_curves.npz', allow_pickle=True) as saved:
        grid = saved['grid']
        rows = pd.Index(saved['genes']).get_indexer(GENES[:2])
        if (rows < 0).any():
            raise ValueError('A urate transporter is missing from the Fig. 5 curve archive')
        for ax, sections in zip(axs[0], (SECTIONS[:2], SECTIONS[2:])):
            for section, linestyle, marker in zip(sections, STYLES, MARKERS):
                fits = saved[f'curve|{section}'][rows]
                bins = saved[f'bins|{section}'][rows]
                bin_x = saved[f'bin_x|{section}'][rows]
                supported = ((grid >= axes[f'x_lo_{section}']) & (grid <= axes[f'x_hi_{section}'])
                             & np.isfinite(fits).all(axis=0))
                for curve, obs, xpos, color in zip(fits, bins, bin_x, COLORS):
                    peak = curve[supported].max()
                    ax.plot(grid[supported], curve[supported] - peak, color=color, ls=linestyle, lw=1.55)
                    ax.scatter(xpos, obs - peak, color=color, marker=marker, s=8, alpha=0.53,
                               linewidths=0, zorder=3)
            ax.axhline(0, color='0.83', lw=0.6, zorder=0)
            ax.set(xlim=(0, 1), ylim=(-2.7, 0.22), yticks=[-2, -1, 0],
                   xticks=[0, 0.25, 0.5, 0.75, 1], xticklabels=['0', '0.25', '0.5', '0.75', '1'],
                   xlabel='PT position · mouse-reference quantile')
    axs[0, 0].set_ylabel('Expression relative to peak (log₂)')

    mouse = pd.read_csv(data_root / 'external/segment_microdissection_sex/GSE212213_TPM_Male_Female.txt.gz',
                        sep='\t', index_col='Gene_symbol').loc[GENES[:2]]
    ax = axs[1, 0]
    for prefix, style, marker in zip(('', 'f'), STYLES, MARKERS):
        tpm = mouse[[prefix + name for name in SEGMENTS]].to_numpy(float)
        if not np.all(tpm > 0):
            raise ValueError('Mouse segment TPMs must be positive for the log fold change')
        for rel, color in zip(np.log2(tpm / tpm[:, :1]), COLORS):
            ax.plot(range(3), rel, color=color, ls=style, marker=marker, lw=1.4, ms=3.1)
    ax.set(xlim=(-0.2, 2.2), ylim=(-0.28, 4.13), xticks=[0, 1, 2], xticklabels=['S1', 'S2', 'S3'],
           ylabel='Expression relative to S1 (log₂)')
    ax.legend([Line2D([], [], color='0.35', ls=style, marker=marker, ms=3) for style, marker in zip(STYLES, MARKERS)],
              ['Male', 'Female'], loc='upper left', bbox_to_anchor=(0.02, 1.015),
              fontsize=7, borderaxespad=0, labelspacing=0.25)

    ax = axs[1, 1]
    lake = lake_endpoint_contrasts(lake_cache, min_nuclei=50)[GENES[:2]].to_numpy()
    ax.axhline(0, color='0.73', lw=0.7, zorder=0)
    for x, y in lake:
        ax.plot([0, 1], [x, y], color='0.78', lw=0.62, zorder=1)
    for j, (color, marker) in enumerate(zip(COLORS, ('o', 's'))):
        ax.scatter(np.full(len(lake), j), lake[:, j], s=17, color=color, marker=marker, alpha=0.68,
                   linewidths=0, zorder=2)
    ax.set(xlim=(-0.32, 1.32), xticks=[0, 1], xticklabels=['NPT4', 'URAT1'],
           ylabel='S3 − S1 expression (log₂)')

    titles = ['Mouse COAST · two kidneys', 'Human COAST · two sections, one donor',
              'Mouse microdissection', f'Human kidney atlas · {len(lake)} donors']
    for letter, ax, title in zip('abcd', axs.flat, titles):
        ax.set_title(title, loc='left', fontsize=8, pad=7, weight='bold')
        offset = ScaledTranslation(-16 / 72, 7 / 72, fig.dpi_scale_trans)
        ax.text(0, 1, letter, transform=ax.transAxes + offset, weight='bold', size=9, va='bottom')
    fig.text(0.11, 0.021, 'COAST lines show fits for the two sections per species. '
             'Circles and diamonds mark their observed positional bins (20 per section).',
             color='0.28', fontsize=7, va='bottom')
    return fig


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path, default=os.environ.get('PSEUDOSPACE_DATA_ROOT', REPO / 'data'))
    parser.add_argument('--results-root', type=Path,
                        default=os.environ.get('PSEUDOSPACE_RESULTS_ROOT', REPO / 'results'))
    parser.add_argument('--lake-pseudobulk', type=Path, required=True)
    parser.add_argument('--audit-layout', action='store_true',
                        help='run the layout audit with the nature-figure scripts on PYTHONPATH')
    args = parser.parse_args()
    root = args.results_root / 'paper_figures/fig5'
    export = root / 'urate'
    export.mkdir(parents=True, exist_ok=True)
    fig = render(root / 'screen', args.lake_pseudobulk, args.data_root)
    if args.audit_layout:
        from audit_panel_alignment import require_matplotlib_panel_alignment

        require_matplotlib_panel_alignment(
            fig, json_out=export / 'fig5_urate.alignment.json', require_panel_labels=True, strict=True,
        )
    fig.savefig(export / 'fig5_urate.pdf')
    fig.savefig(export / 'fig5_urate.svg')
    fig.savefig(export / 'fig5_urate.png', dpi=300)
    fig.savefig(export / 'fig5_urate.tiff', dpi=600, pil_kwargs={'compression': 'tiff_lzw'})


if __name__ == '__main__':
    main()
