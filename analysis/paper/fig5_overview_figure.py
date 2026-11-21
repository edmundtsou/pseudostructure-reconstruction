"""Render an overview of cross-species PT shapes and two transport examples for Fig. 5."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.transforms import ScaledTranslation

from fig5_global_shapes import SECTIONS, load_global_shapes
from fig5_urate import REPO

SPECIES = ('Mouse', 'Human')
SECTION_NAMES = ('Kidney 1', 'Kidney 2', 'Cortex 1', 'Cortex 2')
SPECIES_COLORS = ('#31597C', '#687D3D')
GENE_COLORS = ('#17688F', '#AE491F')
STYLES = ('-', '--', '-', '--')
MARKERS = ('o', 'D', 'o', 'D')
CLASS_NAMES = {
    1: 'Early → late', 2: 'Late → early',
    3: 'Middle → late', 4: 'Early → middle',
}
PAIR_PANELS = (
    ('Urate export and reuptake', ('Slc17a3', 'Slc22a12'), ('NPT4', 'URAT1')),
    ('Organic-anion uptake and export', ('Slc22a8', 'Abcc2'), ('OAT3', 'MRP2')),
)
WIDTH_MM = 183


def panel_label(ax, letter, title, *, pad=7):
    ax.set_title(title, loc='left', pad=pad, fontsize=8, weight='bold')
    offset = ScaledTranslation(-15 / 72, pad / 72, ax.figure.dpi_scale_trans)
    ax.text(0, 1, letter, transform=ax.transAxes + offset, weight='bold', fontsize=9, va='bottom')


def heatmaps(fig, axs, x, genes, curves):
    """Show the same ordered genes in all four sections on the shared fitting interval."""
    cmap = LinearSegmentedColormap.from_list('centered_expression', ['#234A7B', '#FAF8F3', '#AD401F'])
    counts = genes.groupby('cluster').size()
    boundaries = np.r_[0, counts.cumsum().to_numpy()]
    half_step = (x[1] - x[0]) / 2
    for i, ax in enumerate(axs):
        heat = ax.imshow(curves[i], origin='upper', aspect='auto', interpolation='nearest',
                         extent=(x[0] - half_step, x[-1] + half_step, len(genes), 0),
                         cmap=cmap, vmin=-2, vmax=2, rasterized=True)
        for bound in boundaries[1:-1]:
            ax.axhline(bound, color='white', linewidth=0.8)
        ax.set(xlim=(x[0] - half_step, x[-1] + half_step), xticks=[0.10, 0.50, 0.90], yticks=[])
        ax.set_xticklabels(['0.1', '0.5', '0.9'])
        ax.tick_params(axis='x', length=2, pad=2)
        ax.spines[['left', 'bottom']].set_visible(False)
        dim = ax.get_position()
        fig.text(dim.x0 + dim.width / 2, dim.y1 + 0.007, SECTION_NAMES[i], ha='center', fontsize=7)
    for species, left, right, color in zip(SPECIES, axs[::2], axs[1::2], SPECIES_COLORS):
        l, r = left.get_position(), right.get_position()
        label = species if species == 'Mouse' else 'Human · one donor'
        fig.text((l.x0 + r.x1) / 2, l.y1 + 0.024, label, ha='center', color=color, weight='bold', fontsize=8)
    left = axs[0]
    for i, (cluster, count) in enumerate(counts.items()):
        mid = (boundaries[i] + boundaries[i + 1]) / 2
        left.annotate(f'{CLASS_NAMES[cluster]}\n{count} genes', xy=(0, mid),
                      xycoords=('axes fraction', 'data'), xytext=(-6, 0),
                      textcoords='offset points', va='center', ha='right', fontsize=7)
    top = left.get_position()
    fig.text(0.051, top.y1 + 0.027, 'a', weight='bold', fontsize=9)
    fig.text(0.051, top.y1 + 0.011, 'Matched genes', fontsize=7, weight='bold')
    end = axs[-1].get_position()
    color_ax = fig.add_axes([end.x1 + 0.013, end.y0 + end.height * .235, 0.011, end.height * .53])
    cb = fig.colorbar(heat, cax=color_ax, ticks=[-2, 0, 2], extend='both')
    cb.outline.set_visible(False)
    cb.ax.tick_params(length=2, pad=1.5, labelsize=6.5)
    cb.set_label('Centered expression', labelpad=1, fontsize=7)
    return color_ax


def class_profiles(axs, x, genes, curves):
    for i, (cluster, label) in enumerate(CLASS_NAMES.items()):
        ax = axs[i]
        mean = curves[:, genes.cluster.eq(cluster).to_numpy(), :].mean(1)
        ax.axhline(0, lw=0.55, color='0.85')
        for j in range(len(SECTIONS)):
            ax.plot(x, mean[j], color=SPECIES_COLORS[j // 2], ls=STYLES[j], lw=1.1)
        ax.set(xlim=(x[0], x[-1]), ylim=(-1.5, 1.15), xticks=[0.10, 0.50, 0.90], yticks=[-1, 0, 1])
        ax.set_xticklabels(['0.1', '0.5', '0.9'])
        ax.set_title(label, loc='left', pad=5, weight='bold')
        if i != 0:
            ax.tick_params(axis='y', labelleft=False)
    axs[0].set_ylabel('Mean profile')
    panel_label(axs[0], 'b', CLASS_NAMES[1], pad=5)


def paired_profiles(fig, axs, x, screen):
    """Colored fits and raw position bins on a within-gene scale shared across sections."""
    with np.load(screen / 'all_curves.npz', allow_pickle=True) as data:
        grid_idx = pd.Index(data['grid']).get_indexer(x)
        for i, (title, genes, labels) in enumerate(PAIR_PANELS):
            left, right = axs[2 * i:2 * i + 2]
            rows = pd.Index(data['genes']).get_indexer(genes)
            if (rows < 0).any() or (grid_idx < 0).any():
                raise ValueError('Missing fitted pair or supported coordinate')
            fits = np.stack([data[f'curve|{section}'][rows][:, grid_idx] for section in SECTIONS])
            fitted_peaks = fits.max(axis=-1)
            scale = np.ptp(fits, axis=-1).max(axis=0)
            if not np.isfinite(fits).all() or (scale <= 0).any():
                raise ValueError('Transporter profiles need finite, nonflat fits')
            for j, ax in enumerate((left, right)):
                for k in range(2 * j, 2 * j + 2):
                    section = SECTIONS[k]
                    values = data[f'bins|{section}'][rows]
                    positions = data[f'bin_x|{section}'][rows]
                    for g, color in enumerate(GENE_COLORS):
                        ax.plot(x, (fits[k, g] - fitted_peaks[k, g]) / scale[g],
                                color=color, ls=STYLES[k], lw=1.05)
                        valid = np.isfinite(values[g]) & (positions[g] >= x[0]) & (positions[g] <= x[-1])
                        ax.scatter(positions[g, valid], (values[g, valid] - fitted_peaks[k, g]) / scale[g],
                                   color=color, s=5, lw=0, alpha=.62, marker=MARKERS[k], zorder=3)
                ax.axhline(0, lw=0.55, color='0.85')
                ax.set(xlim=(x[0], x[-1]), ylim=(-1.22, .22), xticks=[0.10, .50, .90], yticks=[-1, -.5, 0])
                ax.set_xticklabels(['0.1', '0.5', '0.9'])
                ax.set_title(SPECIES[j], loc='left', color=SPECIES_COLORS[j], pad=4, fontsize=7)
                if ax is not axs[0]:
                    ax.tick_params(axis='y', labelleft=False)
            l, r = left.get_position(), right.get_position()
            fig.text(l.x0 - 0.028, l.y1 + .037, 'cd'[i], weight='bold', fontsize=9)
            fig.text(l.x0, l.y1 + .037, title, weight='bold', fontsize=8)
            fig.legend([Line2D([], [], color=color, lw=1.2) for color in GENE_COLORS], labels,
                       loc='lower right', bbox_to_anchor=(r.x1, l.y1 + .015), borderaxespad=0,
                       ncol=2, handlelength=1.25, handletextpad=.3, columnspacing=.7, fontsize=7)
    axs[0].set_ylabel('Peak-relative profile')


def observed_overlap(ax, pairs):
    """Plot the same-tubule expression check at a representative threshold."""
    rows = pairs[pairs.fraction.eq(.25)].set_index(['gene_a', 'gene_b', 'section'])
    ax.axvline(1, color='0.6', lw=0.7)
    for i, (_, genes, _) in enumerate(PAIR_PANELS):
        ax.axhline(i, color='0.88', lw=.5)
        for j, section in enumerate(SECTIONS):
            ratio = rows.loc[(*genes, section), 'matched_mito_observed_expected']
            ax.scatter(ratio, i + (j - 1.5) * .15, s=18, marker=MARKERS[j], color=SPECIES_COLORS[j // 2], lw=0)
    ax.set(xlim=(.62, 1.62), ylim=(1.5, -.5), yticks=[0, 1], yticklabels=['NPT4–URAT1', 'OAT3–MRP2'],
           xlabel='Joint-high tubules, observed / expected')
    ax.tick_params(axis='y', length=0)
    panel_label(ax, 'e', 'Raw counts matched by depth and mitochondrial fraction')


def render(fig5_root):
    plt.rcParams.update({
        'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'], 'font.size': 7.5,
        'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': .6,
        'pdf.fonttype': 42, 'svg.fonttype': 'none', 'legend.frameon': False,
        'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
    })
    x, genes, curves = load_global_shapes(fig5_root)
    summary = json.loads((fig5_root / 'shape_clusters/summary.json').read_text())
    if len(genes) != summary['genes_changed']['condition']:
        raise ValueError('Summary and plotted changed-gene count differ')
    pairs = pd.read_csv(fig5_root / 'story_comparison/pair_profiles.csv')
    fig = plt.figure(figsize=(WIDTH_MM / 25.4, 203 / 25.4))
    grid = fig.add_gridspec(4, 4, height_ratios=[3.25, 1, 1.34, .84],
                            left=.173, right=.929, top=.865, bottom=.109, wspace=.30, hspace=.66)
    top = [fig.add_subplot(grid[0, i]) for i in range(4)]
    classes = [fig.add_subplot(grid[1, i]) for i in range(4)]
    examples = [fig.add_subplot(grid[2, i]) for i in range(4)]
    bottom = fig.add_subplot(grid[3, :])
    fig.text(.051, .974, 'Proximal tubule expression patterns in mouse and human',
             weight='bold', fontsize=11.5, va='top')
    fig.text(.051, .935, f"{len(genes)} of {summary['n_eligible']} eligible orthologs changed profile shape",
             fontsize=8.5, va='top')
    both = int((genes.graded_A & genes.graded_B).sum())
    fig.text(.929, .935, f'{both} varied ≥1.4-fold in both species',
             fontsize=7.5, ha='right', va='top', color='.29')
    color_ax = heatmaps(fig, top, x, genes, curves)
    class_profiles(classes, x, genes, curves)
    paired_profiles(fig, examples, x, fig5_root / 'screen')
    observed_overlap(bottom, pairs)
    fig.text(.173, .057, 'Heatmap rows share a gene-specific scale across all four sections.\n'
             'COAST position increases from early to late. Dots in c and d show observed expression bins.\n'
             'Circles and solid lines mark section 1. Diamonds and dashed lines mark section 2.',
             fontsize=6.8, va='top', color='.28')
    return fig, color_ax, genes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results-root', type=Path,
                        default=os.environ.get('PSEUDOSPACE_RESULTS_ROOT', REPO / 'results'))
    parser.add_argument('--audit-layout', action='store_true',
                        help='run the layout audit with the nature-figure scripts on PYTHONPATH')
    args = parser.parse_args()
    root = args.results_root / 'paper_figures/fig5'
    out = root / 'overview'
    out.mkdir(parents=True, exist_ok=True)
    fig, color_ax, genes = render(root)
    if args.audit_layout:
        from audit_panel_alignment import require_matplotlib_panel_alignment

        require_matplotlib_panel_alignment(
            fig, exclude_axes=[color_ax], json_out=out / 'fig5_overview.alignment.json', strict=True,
        )
    genes.to_csv(out / 'heatmap_gene_order.csv')
    fig.savefig(out / 'fig5_overview.pdf')
    fig.savefig(out / 'fig5_overview.svg')
    fig.savefig(out / 'fig5_overview.png', dpi=300)
    fig.savefig(out / 'fig5_overview.tiff', dpi=600, pil_kwargs={'compression': 'tiff_lzw'})


if __name__ == '__main__':
    main()
