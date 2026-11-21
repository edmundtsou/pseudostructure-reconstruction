"""Draft Fig. 5 from the full ortholog universe, twofold positional filters and raw count examples."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, to_rgb
from matplotlib.lines import Line2D
from sklearn.cluster import KMeans

from fig5_curve_audit import WINDOW, classify, curve_arrays
from fig5_overview_figure import SPECIES_COLORS
from fig5_urate import REPO, SECTIONS

PROGRAMS = {
    'SLC transport': 'R-HSA-425407',
    'Amino-acid transport': 'R-HSA-352230',
    'Protein digestion': 'KEGG_2019_Mouse::Protein digestion and absorption',
    'Amino-acid metabolism': 'R-HSA-71291',
    'Carbohydrate metabolism': 'R-HSA-71387',
    'Fatty-acid metabolism': 'R-HSA-8978868',
    'Lysosome': 'KEGG_2019_Mouse::Lysosome',
    'Peroxisome': 'KEGG_2019_Mouse::Peroxisome',
}
CLUSTER_COLORS = ('#A14C31', '#2F7875', '#6F6CA0', '#7F6840')
CLUSTER_NAMES = ('Early +\nLate −', 'Early −\nLate +', 'Early dip', 'Middle dip')
STATUS_COLORS = ('#537741', '#868789', '#AE4D39')
HEAT = LinearSegmentedColormap.from_list('centered_expression', ['#234A7B', '#F9F7F1', '#AD401F'])
WIDTH_MM = 183
HEIGHT_MM = 195


def load_profiles(root):
    meta = pd.read_csv(root / 'screen/genes.csv').set_index('gene')
    genes = meta.index[meta.screened]
    x = np.arange(WINDOW[0], WINDOW[1] + .001, .005).round(3)
    fit, observed = curve_arrays(root / 'screen', genes, x)
    held_fit, held_observed = curve_arrays(root / 'screen_genehalf_2026-11-18', genes, x)
    if not np.isfinite(fit).all() or len(genes) != meta.screened.sum():
        raise ValueError('The matched heatmap needs all screened gene profiles on common support')
    quality = classify(fit, observed, held_fit, held_observed, genes, amplitude=1.)
    ranges = np.ptp(fit, axis=2).reshape(2, 2, len(genes)).min(1)
    centered = fit - fit.mean(axis=-1, keepdims=True)
    means = centered.reshape(2, 2, len(genes), len(x)).mean(1)
    table = quality.assign(mouse_min_range=ranges[0], human_min_range=ranges[1])
    return x, fit, means, table


def cluster_profiles(fit, means, table):
    """Cluster whole human-minus-mouse difference curves on a shared within-gene scale."""
    retained = table.mouse_both_grade | table.human_both_grade
    scale = fit[:, retained].std(axis=-1).max(0)
    contrast = (means[1, retained] - means[0, retained]) / scale[:, None]
    model = KMeans(n_clusters=4, n_init=50, random_state=2026).fit(contrast)
    centers = model.cluster_centers_
    early, late = centers[:, :25].mean(1), centers[:, -25:].mean(1)
    early_gain = np.argmax(early - late)
    late_gain = np.argmax(late - early)
    middle_dip = np.argmin(centers[:, 45:95].mean(1) - (early + late) / 2)
    order = [early_gain, late_gain]
    order.extend(c for c in range(len(centers)) if c not in (early_gain, late_gain, middle_dip))
    order.append(middle_dip)
    if len(set(order)) != len(centers):
        raise ValueError('Centroids do not match the draft contrast labels')
    label = np.argsort(order)[model.labels_] + 1
    distance = np.square(contrast - centers[model.labels_]).mean(1)
    ranked = np.lexsort((distance, label))
    return table.loc[retained].assign(cluster=label).iloc[ranked]


def program_counts(root, retained):
    annotations = json.loads((root.parent / 'gene_sets.json').read_text())
    selections = {'All retained': retained.index}
    for name, key in PROGRAMS.items():
        source, = [members for label, members in annotations.items() if label.endswith(key)]
        selections[name] = retained.index.intersection(source)
    rows = []
    for name, genes in selections.items():
        members = retained.loc[genes]
        both = members.reliable_both.sum()
        conserved = members.conserved.sum()
        divergent = members.divergent.sum()
        rows.append({
            'program': name, 'retained': len(members),
            **{f'cluster_{i}': int(members.cluster.eq(i).sum()) for i in range(1, 5)},
            'both_supported': both, 'conserved': conserved, 'intermediate': both - conserved - divergent,
            'divergent': divergent,
        })
    return pd.DataFrame(rows).set_index('program')


def heading(ax, letter, title, y=1.18):
    ax.text(0, y, letter, transform=ax.transAxes, weight='bold', fontsize=8.8, va='bottom')
    ax.annotate(title, xy=(0, y), xycoords='axes fraction', xytext=(12, 0), textcoords='offset points',
                weight='bold', fontsize=7.7, va='bottom')


def paired_heatmap(ax, x, means, clusters=None):
    """The same gene order and absolute centered log2 scale within each pair."""
    extent = len(means[0])
    ax.set_facecolor('none')
    ax.set(xlim=(0, 2.1 if clusters is None else 2.53), ylim=(extent, 0), yticks=[])
    half_pixel = (x[1] - x[0]) / 2
    edges = (x[0] - half_pixel, x[-1] + half_pixel)
    tick_x, tick_labels = [], []
    for i, label in enumerate(('Mouse', 'Human')):
        offset = 1.1 * i
        heat = ax.imshow(means[i], extent=(offset, offset + 1, extent, 0), aspect='auto',
                         interpolation='none', origin='upper', cmap=HEAT, vmin=-1.5, vmax=1.5)
        tick_x.extend(offset + (np.array([.2, .5, .8]) - edges[0]) / (edges[1] - edges[0]))
        tick_labels.extend(['.2', '.5', '.8'])
        ax.annotate(label, (offset + .5, 1), xycoords=('data', 'axes fraction'), xytext=(0, 4),
                    textcoords='offset points', ha='center', color=SPECIES_COLORS[i], weight='bold')
    if clusters is not None:
        counts = clusters.value_counts().sort_index()
        bounds = np.r_[0, counts.cumsum()]
        for i, count in enumerate(counts):
            ax.axhspan(bounds[i], bounds[i + 1], xmin=2.15 / 2.53, xmax=2.21 / 2.53,
                       color=CLUSTER_COLORS[i], lw=0)
            ax.text(2.27, (bounds[i] + bounds[i + 1]) / 2, f'{i + 1} · {count}', va='center', fontsize=6)
        for bound in bounds[1:-1]:
            ax.axhline(bound, xmax=2.1 / 2.53, color='white', lw=.65)
    ax.set_xticks(tick_x, labels=tick_labels)
    ax.tick_params(axis='x', length=0, pad=2)
    ax.spines[['bottom', 'left']].set_visible(False)
    return heat


def range_filter(ax, table):
    m, h = table.mouse_both_grade, table.human_both_grade
    retained = m | h
    xx, yy = table.mouse_min_range, table.human_min_range
    ax.set(xlim=(0, 6), ylim=(0, 5), xlabel='Mouse range (log₂)', ylabel='Human range (log₂)')
    ax.set_xscale('asinh', linear_width=.22)
    ax.set_yscale('asinh', linear_width=.22)
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_ticks([0, .5, 1, 2, 4], labels=['0', '0.5', '1', '2', '4'])
        axis.set_minor_locator(plt.NullLocator())
    ax.axvline(1, color='.44', lw=.65, ls='--')
    ax.axhline(1, color='.44', lw=.65, ls='--')
    ax.scatter(xx[~retained], yy[~retained], s=.8, color='#717580', alpha=.14, edgecolor='none', rasterized=True)
    ax.scatter(xx[retained], yy[retained], s=2.4, color='#464A65', alpha=.55, edgecolor='none', rasterized=True)
    for number, name, x, y, align in [
        ((~m & h).sum(), 'Human only', .04, .95, 'left'),
        ((m & h).sum(), 'Both', .96, .95, 'right'),
        ((m & ~h).sum(), 'Mouse only', .96, .38, 'right'),
    ]:
        ax.text(x, y, f'{number}\n{name}', transform=ax.transAxes, ha=align, va='top', fontsize=6.2,
                bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': .94, 'pad': 1})


def counts_grid(ax, values, column_labels, colors, row_labels=()):
    data = values.to_numpy()
    totals = data.sum(axis=1, keepdims=True)
    fraction = np.divide(data, totals, out=np.zeros(data.shape), where=totals > 0)
    palette = np.asarray([to_rgb(color) for color in colors])
    tint = .65 * np.sqrt(fraction)[:, :, None]
    ax.imshow(1 - tint * (1 - palette[None, :, :]), origin='upper', aspect='auto', interpolation='nearest')
    for i, counts in enumerate(data):
        for j, count in enumerate(counts):
            ax.text(j, i, str(count), ha='center', va='center', color='#202126' if count else '.58', fontsize=6.4)
    ax.set_xticks(np.arange(len(column_labels)), labels=column_labels)
    ax.set_yticks(np.arange(len(row_labels)), labels=row_labels)
    ax.xaxis.tick_top()
    ax.tick_params(axis='both', length=0, pad=2, labelsize=6.4)
    for spine in ax.spines.values():
        spine.set_visible(False)
    for i in range(1, len(data)):
        ax.axhline(i - .5, color='.87', lw=.45)
    for j in range(1, len(column_labels)):
        ax.axvline(j - .5, color='white', lw=.8)


def examples(root, axs, x, fit, order):
    with np.load(root / 'curve_audit/reaggregated_counts.npz') as raw:
        indices = pd.Index(raw['genes']).get_indexer(order)
        if (indices < 0).any() or raw['sections'].tolist() != SECTIONS:
            raise ValueError('The count audit does not match Fig. 5 genes or sections')
        xx = raw['centers']
        yy = raw['profiles'][:, indices]
    for j, (ax, name, role) in enumerate(zip(axs, order, (
        'Conserved transport', 'Conserved peptidase', 'Earlier in human', 'Later in human',
    ))):
        offset = fit[:, j].mean(-1)
        for i in range(4):
            color = SPECIES_COLORS[i // 2]
            ax.plot(x, fit[i, j] - offset[i], color=color, lw=.95, ls=('-', '--')[i % 2])
            ax.scatter(xx, yy[i, j] - offset[i], marker=('o', 'D')[i % 2], color=color,
                       s=4, alpha=.70, lw=0, zorder=3)
        ax.axhline(0, lw=.55, color='.82')
        ax.set(xlim=(.14, .93), xticks=[.2, .5, .8])
        ax.set_xticklabels(['0.2', '0.5', '0.8'])
        heading(ax, 'fghi'[j], name.upper(), y=1.14)
        ax.set_title(role, loc='left', fontsize=6.5, pad=4.5)
    axs[0].set_ylabel('Centered log₂ rate')


def render(root):
    x, fit, means, table = load_profiles(root)
    retained = cluster_profiles(fit, means, table)
    program = program_counts(root, retained)
    order = np.argsort(-means[0] @ (x - x.mean()), kind='stable')
    keep = table.index.get_indexer(retained.index)
    names = ['Slc6a19', 'Dpep1', 'Mme', 'Pck1']
    if not pd.Index(names).isin(retained.index).all():
        raise ValueError('Draft examples must pass the requested twofold filter')
    plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 7, 'axes.labelsize': 6.6,
                         'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': .55,
                         'xtick.major.size': 2, 'ytick.major.size': 2, 'pdf.fonttype': 42, 'svg.fonttype': 'none'})
    fig = plt.figure(figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4))
    grid = fig.add_gridspec(3, 24, height_ratios=(1.68, 1.12, 1), top=.855, bottom=.123,
                            left=.08, right=.976, hspace=.67, wspace=.55)
    a, b, c = [fig.add_subplot(spec) for spec in (grid[0, :6], grid[0, 8:13], grid[0, 15:])]
    d, e = [fig.add_subplot(spec) for spec in (grid[1, 5:14], grid[1, 15:])]
    bottom = [fig.add_subplot(grid[2, i:i + 5]) for i in (0, 6, 12, 18)]
    for ax, letter in zip([a, b, c, d, e, *bottom], 'abcdefghi'):
        ax.set_gid(letter)
    fig.text(.054, .969, 'Proximal tubule expression in mouse and human', weight='bold', fontsize=11, va='top')
    fig.text(.054, .928, 'Two mouse kidneys and two cortical sections from one human donor', fontsize=7.3, color='.26')
    heading(a, 'a', f'{len(table):,} screened orthologs', y=1.115)
    paired_heatmap(a, x, means[:, order])
    heading(b, 'b', 'Twofold positional filter', y=1.115)
    range_filter(b, table)
    b.set_title(f'{len(table) - len(retained):,} below cutoff', fontsize=6.7, loc='left', pad=5)
    heading(c, 'c', f'{len(retained)} retained in four clusters', y=1.115)
    heat = paired_heatmap(c, x, means[:, keep], retained.cluster)
    base = a.get_position()
    color_ax = fig.add_axes([base.x0 + .052, base.y0 - .038, .106, .008])
    cb = fig.colorbar(heat, cax=color_ax, orientation='horizontal', ticks=[-1.5, 0, 1.5], extend='both')
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=5.7, length=1.4, pad=1.2)
    cb.set_ticklabels(['−1.5', '0', '+1.5'])
    cb.set_label('Centered log₂ expression', fontsize=5.9, labelpad=1)

    heading(d, 'd', 'Program counts by cluster', y=1.26)
    counts_grid(d, program.iloc[:, 1:5], [f'{i + 1}  {name}' for i, name in enumerate(CLUSTER_NAMES)],
                CLUSTER_COLORS, program.index)
    heading(e, 'e', f'Shape agreement in both groups ({retained.reliable_both.sum()})', y=1.26)
    counts_grid(e, program[['conserved', 'intermediate', 'divergent']],
                ['Conserved', 'Intermediate', 'Divergent'], STATUS_COLORS)

    examples(root, bottom, x, fit[:, table.index.get_indexer(names)], names)
    pairs = pd.read_csv(root / 'curve_audit/observed_program_pairs.csv')
    pair = pairs.loc[pairs.gene_a.eq('Mme') & pairs.gene_b.eq('Dpep1'), SECTIONS].to_numpy().ravel()
    if pair.size != 4 or not np.isfinite(pair).all():
        raise ValueError('Missing MME–DPEP1 observed-bin correlations')
    r = [f'{v:+.2f}'.replace('-', '−') for v in pair]
    fig.legend([Line2D([], [], color=color, lw=1.2) for color in SPECIES_COLORS], ['Mouse', 'Human'],
               bbox_to_anchor=(.94, .092), loc='upper right', ncol=2, borderaxespad=0,
               handlelength=1.4, frameon=False, fontsize=6.6)
    fig.text(.08, .083, 'COAST position · Points show observed position-bin rates.', fontsize=6.5)
    fig.text(.08, .062, f'MME–DPEP1 across bins · Mouse r {r[0]}, {r[1]} · Human r {r[2]}, {r[3]}', fontsize=6.5)
    fig.text(.08, .031, 'The cutoff requires a twofold range in both sections of either group.\n'
             'Cluster labels describe the human-minus-mouse contrast.', color='.28', fontsize=6.2)
    table = table.assign(cluster=retained.cluster,
                         all_order=pd.Series(np.arange(1, len(order) + 1), index=table.index[order]),
                         cluster_order=pd.Series(np.arange(1, len(keep) + 1), index=retained.index))
    return fig, color_ax, table, program


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results-root', type=Path,
                        default=os.environ.get('PSEUDOSPACE_RESULTS_ROOT', REPO / 'results'))
    parser.add_argument('--audit-layout', action='store_true')
    args = parser.parse_args()
    root = args.results_root / 'paper_figures/fig5'
    out = root / 'three_row'
    out.mkdir(parents=True, exist_ok=True)
    fig, color_ax, table, program = render(root)
    if args.audit_layout:
        from audit_panel_alignment import require_matplotlib_panel_alignment

        require_matplotlib_panel_alignment(fig, exclude_axes=[color_ax], strict=True,
                                          json_out=out / 'fig5_three_row.alignment.json')
    fig.savefig(out / 'fig5_three_row.png', dpi=300)
    fig.savefig(out / 'fig5_three_row.svg', dpi=600)
    fig.savefig(out / 'fig5_three_row.pdf', dpi=600)
    fig.savefig(out / 'fig5_three_row.tiff', dpi=600, pil_kwargs={'compression': 'tiff_lzw'})
    table.to_csv(out / 'ortholog_profiles.csv')
    program.to_csv(out / 'program_counts.csv')


if __name__ == '__main__':
    main()
