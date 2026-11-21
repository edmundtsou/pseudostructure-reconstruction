"""Render the Fig. 5 urate evidence from the study's four Visium HD sections."""
from __future__ import annotations

import argparse
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.transforms import ScaledTranslation

from fig5_urate import GENES, REPO, SECTIONS
from fig5_urate_tissue import HIGH_FRACTIONS, overlap_table, raw_endpoint_contrasts, tissue_units

GENE_COLORS = ('#17668C', '#AD481B')
MAP_COLORS = (*GENE_COLORS, '#714E90')
SECTION_COLORS = ('#304555', '#304555', '#64763E', '#64763E')
SECTION_NAMES = ('Mouse 1', 'Mouse 2', 'Human 1', 'Human 2')
STYLES = ('-', '--', '-', '--')
MARKERS = ('o', 'D', 'o', 'D')
TRANSPORTERS = {'Slc17a3': 'NPT4', 'Slc17a1': 'NPT1', 'Slc22a12': 'URAT1',
                'Slc22a6': 'OAT1', 'Slc22a7': 'OAT2', 'Slc22a8': 'OAT3'}
WIDTH_MM = 183


def ortholog_panel(ax, genes):
    """Global context for the focused transporter comparison."""
    eligible = genes.loc[genes.assessable]
    x = eligible[[f'E_{s}' for s in SECTIONS[:2]]].mean(axis=1)
    y = eligible[[f'E_{s}' for s in SECTIONS[2:]]].mean(axis=1)
    ax.axhline(0, color='0.83', lw=0.6, zorder=0)
    ax.axvline(0, color='0.83', lw=0.6, zorder=0)
    ax.scatter(x, y, color='#A6AFB6', alpha=0.38, s=3, lw=0)
    ax.scatter(x[eligible.catR], y[eligible.catR], color='#704F8A', s=7, lw=0)
    for gene, color, offset, label in zip(GENES[:2], GENE_COLORS, [(9, -12), (8, 12)], ('NPT4', 'URAT1')):
        xy = (x[gene], y[gene])
        ax.scatter(*xy, s=22, color=color, linewidth=0.7, edgecolor='white', zorder=5)
        ax.annotate(label, xy, xytext=offset, textcoords='offset points', color=color, fontsize=6.5,
                    arrowprops={'arrowstyle': '-', 'color': color, 'lw': 0.5})
    ax.set(xlim=(-9, 9), ylim=(-9, 9), xticks=[-8, -4, 0, 4, 8], yticks=[-8, -4, 0, 4, 8],
           xlabel='Mouse late − early (log₂)', ylabel='Human late − early (log₂)')
    ax.legend([Line2D([], [], marker='o', ls='', color='#704F8A', ms=3)],
              [f'{int(eligible.catR.sum())} reversals'], loc='upper left', fontsize=6.5, handletextpad=0.15)


def profile_panels(axs, screen, axes):
    """Fitted and observed transcript profiles on their supported COAST intervals."""
    with np.load(screen / 'all_curves.npz', allow_pickle=True) as saved:
        grid = saved['grid']
        rows = pd.Index(saved['genes']).get_indexer(GENES[:2])
        if (rows < 0).any():
            raise ValueError('A urate transporter is missing from the Fig. 5 curve archive')
        for ax, indices in zip(axs, (range(2), range(2, 4))):
            for i in indices:
                section = SECTIONS[i]
                fits = saved[f'curve|{section}'][rows]
                bins = saved[f'bins|{section}'][rows]
                bin_x = saved[f'bin_x|{section}'][rows]
                supported = ((grid >= axes[f'x_lo_{section}']) & (grid <= axes[f'x_hi_{section}'])
                             & np.isfinite(fits).all(axis=0))
                for curve, obs, xpos, color in zip(fits, bins, bin_x, GENE_COLORS):
                    peak = curve[supported].max()
                    ax.plot(grid[supported], curve[supported] - peak, color=color, ls=STYLES[i], lw=1.25)
                    ax.scatter(xpos, obs - peak, color=color, marker=MARKERS[i], s=6, alpha=0.53, lw=0, zorder=3)
            ax.axhline(0, color='0.83', lw=0.6, zorder=0)
            ax.set(xlim=(0, 1), ylim=(-2.65, 0.18), yticks=[-2, -1, 0],
                   xticks=[0, 0.5, 1], xticklabels=['0', '0.5', '1'], xlabel='COAST position')
    axs[0].set_ylabel('Relative to fitted peak (log₂)')


def tissue_panel(ax, units, background, displayed):
    """Overlay high-expression calls on native tissue centroids, omitting COAST coordinates."""
    dims = ax.get_position()
    aspect = (dims.width * ax.figure.get_figwidth()) / (dims.height * ax.figure.get_figheight())
    width = min(0.216, 0.82 / aspect)
    stretch = np.array([width, -width * aspect])
    for i, section in enumerate(SECTIONS):
        bg = background.loc[background['sample'].eq(section), ['x', 'y']].to_numpy()
        pt = units.loc[displayed[section].index, ['x', 'y']].to_numpy()
        lo, hi = bg.min(0), bg.max(0)
        center, scale = (lo + hi) / 2, (hi - lo).max()
        center_ax = np.array([(i + 0.5) / 4, 0.55])
        bg = (bg - center) / scale * stretch + center_ax
        pt = (pt - center) / scale * stretch + center_ax
        ax.scatter(bg[:, 0], bg[:, 1], s=0.23, color='#DDDEE0', lw=0, rasterized=True)
        high = displayed[section]
        classes = (high[GENES[0]] & ~high[GENES[1]], high[GENES[1]] & ~high[GENES[0]], high.both)
        for mask, color in zip(classes, MAP_COLORS):
            ax.scatter(pt[mask, 0], pt[mask, 1], s=1.5, color=color, lw=0, rasterized=True)
        ax.text(center_ax[0], 0.07, SECTION_NAMES[i], fontsize=7, ha='center', va='top')
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.set_axis_off()
    ax.legend([Line2D([], [], marker='o', ls='', color=c, ms=4) for c in MAP_COLORS],
              ['NPT4 high only', 'URAT1 high only', 'Both high'], loc='upper right', fontsize=7,
              ncol=3, bbox_to_anchor=(1.01, 1.13), handletextpad=0.15, columnspacing=0.9, borderaxespad=0)


def overlap_panel(ax, overlap):
    ax.axhline(1, color='0.65', lw=0.65, zorder=0)
    for i, section in enumerate(SECTIONS):
        s = overlap.loc[overlap.section.eq(section), 'observed_expected']
        ax.plot(range(len(HIGH_FRACTIONS)), s, color=SECTION_COLORS[i], ls=STYLES[i], marker=MARKERS[i], ms=3, lw=1)
    ax.set(xlim=(-0.2, 2.2), xticks=range(3), xticklabels=['20%', '25%', '33%'],
           xlabel='High-expression threshold', ylabel='Joint high · observed / expected')


def transporter_panel(ax, genes):
    """One raw log contrast per gene and section; values use the same end cores as the screen."""
    cores = genes.loc[list(TRANSPORTERS), [f'R_{s}' for s in SECTIONS]]
    values = cores.to_numpy()
    clevel = np.ceil(np.abs(values).max())
    ax.imshow(values, cmap='RdBu_r', vmin=-clevel, vmax=clevel, aspect='auto', interpolation='nearest')
    for i, j in np.ndindex(values.shape):
        value = values[i, j]
        ax.text(j, i, f'{value:+.1f}', ha='center', va='center', fontsize=6,
                color='white' if abs(value) > clevel * 0.60 else '#222222')
    ax.axvline(1.5, color='white', lw=1.1)
    ax.set(xticks=range(4), xticklabels=['M1', 'M2', 'H1', 'H2'],
           yticks=range(len(cores)), yticklabels=list(TRANSPORTERS.values()), xlabel='Raw late − early (log₂)')
    ax.tick_params(length=0, pad=2)
    ax.spines[['left', 'bottom']].set_visible(False)


def robustness_table(units, genes, heldout, axes):
    """Contrast of raw endpoint contrasts, URAT1 minus NPT4, under four definitions."""
    raw = raw_endpoint_contrasts(units, axes, GENES[:2])
    rows = []
    for section in SECTIONS:
        values = {
            'Cores': genes.loc[GENES[:2], f'R_{section}'],
            'All labels': raw[f'all_labels_{section}'],
            'Held out': heldout.loc[GENES[:2], f'R_{section}'],
            'QC subset': raw[f'qc_subset_{section}'],
        }
        rows.append({'section': section, **{name: val.iloc[1] - val.iloc[0] for name, val in values.items()}})
    return pd.DataFrame(rows).set_index('section')


def robustness_panel(ax, robustness):
    ax.axhline(0, color='0.7', lw=0.65, zorder=0)
    for i, section in enumerate(SECTIONS):
        ax.plot(np.arange(4) + (i % 2 - 0.5) * 0.10, robustness.loc[section], color=SECTION_COLORS[i],
                ls=STYLES[i], marker=MARKERS[i], ms=3, lw=1)
    ax.set(xticks=range(4), xticklabels=['Cores', 'All\nlabels', 'Held\nout', 'QC\nsubset'],
           ylabel='URAT1 − NPT4 contrast (log₂)')


def render(results_root: Path):
    plt.rcParams.update({
        'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'], 'font.size': 7,
        'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': 0.6, 'legend.frameon': False,
        'pdf.fonttype': 42, 'svg.fonttype': 'none',
    })
    root = results_root / 'paper_figures/fig5'
    screen = root / 'screen'
    genes = pd.read_csv(screen / 'genes.csv', index_col='gene')
    axes = pd.read_csv(screen / 'axes.csv').iloc[0]
    heldout = pd.read_csv(root / 'screen_genehalf_2026-11-18/genes.csv', index_col='gene')
    units, background = tissue_units(results_root)
    overlap, high = overlap_table(units)
    robustness = robustness_table(units, genes, heldout, axes)

    fig = plt.figure(figsize=(WIDTH_MM / 25.4, 178 / 25.4))
    gs = fig.add_gridspec(3, 3, height_ratios=[1, 1.27, 1],
                          left=0.084, right=0.985, top=0.86, bottom=0.15, wspace=0.43, hspace=0.65)
    a, b, c = [fig.add_subplot(gs[0, i]) for i in range(3)]
    d = fig.add_subplot(gs[1, :])
    e, f, g = [fig.add_subplot(gs[2, i]) for i in range(3)]
    fig.text(0.084, 0.975, 'Urate transporter expression along the proximal tubule',
             weight='bold', fontsize=12, va='top')
    fig.legend([Line2D([], [], color=color, lw=2) for color in GENE_COLORS],
               ['NPT4 · urate export', 'URAT1 · urate reuptake'],
               loc='upper left', bbox_to_anchor=(0.08, 0.930), ncol=2, fontsize=8)
    fig.legend([Line2D([], [], color=color, ls=style, marker=marker, ms=3)
                for color, style, marker in zip(SECTION_COLORS, STYLES, MARKERS)],
               SECTION_NAMES, loc='upper right', bbox_to_anchor=(0.985, 0.935), ncol=2, fontsize=7)
    ortholog_panel(a, genes)
    profile_panels((b, c), screen, axes)
    tissue_panel(d, units, background, high)
    overlap_panel(e, overlap)
    transporter_panel(f, genes)
    robustness_panel(g, robustness)
    titles = ['Ortholog endpoint contrasts', 'Mouse expression profiles', 'Human expression profiles',
              'High-expressing tubule profiles in tissue', 'Tubule-level overlap',
              'Urate transport components', 'Endpoint sensitivity']
    for letter, ax, title in zip('abcdefg', [a, b, c, d, e, f, g], titles):
        ax.set_title(title, loc='left', fontsize=7.5, pad=7, weight='bold')
        offset = ScaledTranslation(-14 / 72, 7 / 72, fig.dpi_scale_trans)
        ax.text(0, 1, letter, transform=ax.transAxes + offset, weight='bold', size=8.5, va='bottom')
    fig.text(0.084, 0.026, 'Mouse 1–2 are separate kidneys. Human 1–2 are two cortical sections from one kidney.\n'
             'Points in b and c are observed counts per library in 20 position bins per section.', fontsize=7,
             color='0.28', va='bottom')
    return fig, overlap, robustness


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results-root', type=Path,
                        default=os.environ.get('PSEUDOSPACE_RESULTS_ROOT', REPO / 'results'))
    parser.add_argument('--audit-layout', action='store_true',
                        help='run the layout audit with the nature-figure scripts on PYTHONPATH')
    args = parser.parse_args()
    export = args.results_root / 'paper_figures/fig5/urate'
    export.mkdir(parents=True, exist_ok=True)
    fig, overlap, robustness = render(args.results_root)
    if args.audit_layout:
        from audit_panel_alignment import require_matplotlib_panel_alignment

        require_matplotlib_panel_alignment(
            fig, json_out=export / 'fig5_urate.alignment.json', require_panel_labels=True, strict=True,
        )
    overlap.to_csv(export / 'tubule_overlap.csv', index=False)
    robustness.to_csv(export / 'endpoint_sensitivity.csv')
    fig.savefig(export / 'fig5_urate.pdf')
    fig.savefig(export / 'fig5_urate.svg')
    fig.savefig(export / 'fig5_urate.png', dpi=300)
    fig.savefig(export / 'fig5_urate.tiff', dpi=600, pil_kwargs={'compression': 'tiff_lzw'})


if __name__ == '__main__':
    main()
