"""Render the Fig. 5 biological discovery draft from local COAST results."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from sklearn.cluster import KMeans

from fig5_curve_audit import WINDOW, curve_arrays, shape
from fig5_three_row_figure import HEAT, CLUSTER_COLORS
from fig5_urate import REPO, SECTIONS

WIDTH_MM = 183
HEIGHT_MM = 201
SPECIES = ('Mouse', 'Human')
COLORS = ('#275D83', '#8B3E6C')
SECTION_STYLES = ('-', '--', '-', '--')
S2_CLASSES = {
    'S1 high, off at S2': ('Off', '#47798B'),
    'S2 peak': ('Peak', '#B16428'),
    'S2 trough': ('Trough', '#999999'),
    'on at S2': ('On', '#65527E'),
}
SWITCH_COLORS = {
    'early in mouse, delayed in human (model)': ('#B75520', 'o', 'Delayed early'),
    'late in both (model)': ('#3D6489', 'o', 'Late in both'),
    'early in both (conserved early)': ('#597E3C', 'D', 'Early in both'),
}
GAIN_CLASSES = {
    4: ('Early drop', '#736791'),
    1: ('Early high', '#476581'),
    2: ('Early–mid', '#6E8B7B'),
    3: ('Plateau', '#A87945'),
    0: ('Late rise', '#A64C3E'),
}


def axes_mm(fig, x, y, w, h):
    return fig.add_axes([x / WIDTH_MM, y / HEIGHT_MM, w / WIDTH_MM, h / HEIGHT_MM])


def heading(ax, letter, title, lift=12):
    ax.annotate(letter, (0, 1), xycoords='axes fraction', xytext=(0, lift),
                textcoords='offset points', fontweight='bold', fontsize=9, va='bottom')
    ax.annotate(title, (0, 1), xycoords='axes fraction', xytext=(13, lift + .3),
                textcoords='offset points', fontweight='bold', fontsize=7.5, va='bottom')


def paired_map(ax, values, x, band=None, guides=(), band_width=.36):
    """Draw identically ordered log-rate deviations for both species."""
    n = values.shape[1]
    gap = .15
    edge = (x[1] - x[0]) / 2
    left, right = x[0] - edge, x[-1] + edge
    x_local = (np.array([.25, .5, .75]) - left) / (right - left)
    for i, (name, color) in enumerate(zip(SPECIES, COLORS)):
        offset = i * (1 + gap)
        ax.imshow(values[i], extent=(offset, offset + 1, n, 0), cmap=HEAT,
                  vmin=-1.5, vmax=1.5, aspect='auto', interpolation='nearest', rasterized=True)
        ax.annotate(name, (offset + .5, 1), xycoords=('data', 'axes fraction'),
                    xytext=(0, 2.5), textcoords='offset points',
                    ha='center', color=color, fontsize=6.5, fontweight='bold')
        for guide in guides:
            ax.axvline(offset + (guide - left) / (right - left), color='.28',
                       linewidth=.45, dashes=(2, 2), alpha=.8)
    stop = 2 + gap
    if band:
        start = 0
        for label, count, color in band:
            ax.fill_betweenx([start, start + count], stop + .035, stop + .072, color=color, lw=0)
            ax.text(stop + .095, start + count / 2, label, fontsize=6, ha='left', va='center')
            start += count
            if start < n:
                ax.hlines(start, 0, stop, color='white', lw=.55)
        if start != n:
            raise ValueError('Heatmap row groups must cover every displayed gene')
    ax.set(xlim=(0, stop + (band_width if band else 0)), ylim=(n, 0), yticks=[],
           xticks=np.r_[x_local, x_local + 1 + gap], xticklabels=['.25', '.50', '.75'] * 2)
    ax.tick_params(axis='x', length=0, pad=2.3, labelsize=6.1)
    for spine in ax.spines.values():
        spine.set_visible(False)


def switch_panel(ax, switches):
    for label, group in switches.groupby('class'):
        color, marker, _ = SWITCH_COLORS.get(label, ('.67', 'o', 'Other'))
        ax.scatter(group.mouse_switch, group.human_switch_plot,
                   s=19 if marker == 'D' else 11, marker=marker, color=color,
                   edgecolors='white', linewidths=.35, alpha=.95, zorder=3)
    limits = (.07, .91)
    ax.plot(limits, limits, color='.66', linestyle='--', lw=.65)
    ax.set(xlim=limits, ylim=limits, xticks=[.2, .4, .6, .8], yticks=[.2, .4, .6, .8],
           xlabel='Mouse switch position', ylabel='Human switch position')
    handles = [Line2D([], [], marker=m, color=c, ls='', markersize=3.8, label=label)
               for c, m, label in SWITCH_COLORS.values()]
    handles.append(Line2D([], [], marker='o', color='.67', ls='', markersize=3.8, label='Other'))
    ax.legend(handles=handles, loc='upper center', bbox_to_anchor=(.50, -.155), ncol=4,
              fontsize=5.5, handletextpad=.25, handlelength=.7, columnspacing=.7)


def score_panel(ax, tubules):
    paths = []
    for j, section in enumerate(SECTIONS):
        frame = tubules.loc[tubules['sample'].eq(section), ['E', 'L']].copy()
        frame['progress_bin'] = pd.qcut(frame.mean(axis=1).rank(method='first'), 20, labels=False)
        avg = frame.groupby('progress_bin')[['E', 'L']].mean()
        paths.append(avg.assign(section=section))
        ax.plot(avg.E, avg.L, color=COLORS[j // 2], ls=SECTION_STYLES[j],
                linewidth=1.28, marker='o', markevery=[0, -1], markersize=2.25, zorder=3)
    limits = (-2.15, 2.75)
    ax.plot(limits, limits, color='.68', lw=.65, ls=':')
    ax.set(xlim=limits, ylim=limits, xticks=[-2, -1, 0, 1, 2], yticks=[-2, -1, 0, 1, 2],
           xlabel='Early-program score', ylabel='Late-program score')
    return pd.concat(paths)


def gene_order(gene_meta, profiles, x, labels):
    """Group mouse-defined S2 profiles by class and event position in mouse."""
    ordered, bands = [], []
    for name, (label, color) in labels.items():
        group = gene_meta.index[gene_meta['class'].eq(name)]
        z = profiles.loc[group].to_numpy()
        if name == 'S2 peak':
            position = x[z.argmax(axis=-1)]
        elif name == 'S2 trough':
            position = x[z.argmin(axis=-1)]
        else:
            halfway = (z.min(-1) + z.max(-1)) / 2
            crossed = z >= halfway[:, None] if name == 'on at S2' else z <= halfway[:, None]
            position = x[crossed.argmax(axis=-1)]
        ordered.extend(group[np.argsort(position, kind='stable')])
        bands.append((f'{label} {len(group)}', len(group), color))
    if len(ordered) != len(gene_meta):
        raise ValueError('Unexpected S2 gene classes')
    return ordered, bands


def lysosome_panel(ax, lysosome, early):
    other = lysosome.loc[~lysosome.index.isin(early)]
    selected = lysosome.loc[early]
    xlim, ylim = (-4.2, 1.75), (-1.4, 2.82)
    ax.plot(ylim, ylim, color='.67', ls='--', lw=.65)
    for frame, color, size in ((other, '.67', 10), (selected, '#A45824', 17)):
        ax.scatter(frame.observed_mouse, frame.observed_human, color=color, s=size,
                   edgecolors='white', linewidths=.35, zorder=3)
    for gene, offset in {'Lgmn': (-14, -7), 'Ctsb': (4, 4), 'Ctsh': (-12, 7)}.items():
        entry = lysosome.loc[gene]
        ax.annotate(gene.upper(), (entry.observed_mouse, entry.observed_human), xytext=offset,
                    textcoords='offset points', fontsize=5.7, color='#7D431C', fontstyle='italic')
    n = int((lysosome.observed_human > lysosome.observed_mouse).sum())
    ax.text(.06, .94, f'{n} of {len(lysosome)} more early-weighted', transform=ax.transAxes,
            va='top', fontsize=6.2, color='.25')
    ax.set(xlim=xlim, ylim=ylim, xticks=[-4, -2, 0, 1], yticks=[-1, 0, 1, 2],
           xlabel='Mouse early − late', ylabel='Human early − late')


def human_gains(genes, fit, x):
    """Recover the five exploratory human profile groups on the fixed study screen."""
    mouse_range, human_range = np.ptp(fit, axis=-1)
    retained = ((human_range >= .5) & (human_range >= 1.65 * np.maximum(.1, mouse_range))
                & (human_range - mouse_range >= .25))
    human = fit[1, retained]
    modes = KMeans(n_clusters=5, n_init=40, random_state=21).fit_predict(shape(human)[:, ::4])
    return pd.DataFrame({'mode': modes, 'human_peak': x[human.argmax(-1)],
                         'mouse_range': mouse_range[retained], 'human_range': human_range[retained]},
                        index=genes[retained])


def lysosome_contrasts(genes, memberships, raw, raw_x):
    keep = genes.isin(memberships['KEGG_2019_Mouse::Lysosome'])
    early_late = raw[:, keep][:, :, raw_x < .34].mean(-1) - raw[:, keep][:, :, raw_x > .73].mean(-1)
    return pd.DataFrame(early_late.T, columns=['observed_mouse', 'observed_human'], index=genes[keep])


def metabolic_panel(ax, x, raw_x, pooled_fit, pooled_obs, genes):
    keys = ['Gpt', 'Gpt2']
    gene_index = genes.get_indexer(keys)
    if (gene_index < 0).any():
        raise ValueError('Missing GPT or GPT2 in the ortholog screen')
    for i, color in enumerate(COLORS):
        for j, ls in enumerate(('-', '--')):
            ax.plot(x, pooled_fit[i, gene_index[j]], color=color, ls=ls, lw=1.25)
            ax.scatter(raw_x, pooled_obs[i, gene_index[j]], color=color,
                       marker='o' if j == 0 else 's', s=5, lw=0, alpha=.55)
    ax.axhline(0, c='.70', lw=.55, zorder=0)
    legend = [Line2D([], [], ls=ls, c='.23', lw=1, label=key.upper())
              for ls, key in zip(('-', '--'), keys)]
    ax.legend(handles=legend, loc='upper left', fontsize=5.7, handlelength=1.55, labelspacing=.3)
    ax.set(xlim=WINDOW, xticks=[.2, .5, .8], yticks=[-1, 0, 1],
           xlabel='COAST position', ylabel='Centered log₂ rate')


def render(root, draft, s2):
    meta = pd.read_csv(root / 'three_row/ortholog_profiles.csv', index_col='gene')
    genes = meta.index
    x = np.arange(WINDOW[0], WINDOW[1] + .001, .005).round(3)
    fit, _ = curve_arrays(root / 'screen', genes, x)
    with np.load(root / 'curve_audit/reaggregated_counts.npz') as source:
        rows = pd.Index(source['genes']).get_indexer(genes)
        if (rows < 0).any():
            raise ValueError('Observed bins must cover every screened ortholog')
        raw_x = source['centers']
        inside = (raw_x >= x[0]) & (raw_x <= x[-1])
        raw = source['profiles'][:, rows][:, :, inside]
        raw_x = raw_x[inside]
    fit = fit.reshape(2, 2, len(genes), -1).mean(1)
    raw = raw.reshape(2, 2, len(genes), -1).mean(1)
    fit_offset = fit.mean(-1, keepdims=True)
    centered_fit, centered_raw = fit - fit_offset, raw - fit_offset

    switches = pd.read_csv(s2 / 'switches.csv')
    tubules = pd.read_csv(s2 / 'tubules.csv')
    s2_genes = pd.read_csv(s2 / 's2_genes.csv', index_col='gene')
    s2_order, s2_bands = gene_order(s2_genes, pd.DataFrame(fit[0], index=genes), x, S2_CLASSES)
    gain = human_gains(genes, fit, x)
    memberships = json.loads((root.parent / 'gene_sets.json').read_text())
    lysosome = lysosome_contrasts(genes, memberships, raw, raw_x)
    early_lysosome = gain.index[gain['mode'].isin([1, 4]) & gain.index.isin(lysosome.index)]
    retained = meta.loc[meta['cluster_order'].notna()].sort_values('cluster_order')
    order = meta.sort_values('all_order').index
    cutoff = retained.groupby('cluster', sort=True).size()
    cutoff_bands = [(f'C{i + 1}  {n}', n, CLUSTER_COLORS[i]) for i, n in enumerate(cutoff)]
    gain_order, gain_bands = [], []
    for mode, (label, color) in GAIN_CLASSES.items():
        g = gain.loc[gain['mode'].eq(mode)].sort_values('human_peak')
        gain_order.extend(g.index)
        gain_bands.append((f'{label} {len(g)}', len(g), color))

    gain.loc[gain_order].to_csv(draft / 'panel_f_human_gain_genes.csv')
    lysosome.assign(highlighted=lysosome.index.isin(early_lysosome)).to_csv(draft / 'panel_g_lysosome.csv')
    s2_genes.loc[s2_order].to_csv(draft / 'panel_e_s2_gene_order.csv')

    fig = plt.figure(figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4), facecolor='white')
    axes = [axes_mm(fig, l, 151, 77, 34) for l in (11, 98)]
    axes += [axes_mm(fig, l, y, 48, 48) for y in (89, 22) for l in (11, 69, 127)]
    for ax, letter in zip(axes, 'abcdefgh'):
        ax.set_gid(letter)
    a, b, c, d, e, f, g, h = axes
    titles = ('All screened orthologs', 'Genes with twofold positional variation',
              'Gene switch positions', 'Early and late programs', 'Mouse S2 genes',
              'Stronger gradients in human', 'Early lysosomal expression', 'GPT and GPT2 profiles')
    for ax, letter, title in zip(axes, 'abcdefgh', titles):
        heading(ax, letter, title)

    def heatmap(ax, labels, **kwargs):
        j = genes.get_indexer(labels)
        if (j < 0).any():
            raise ValueError('Heatmap rows are missing from the matched ortholog set')
        paired_map(ax, centered_fit[:, j], x, **kwargs)

    heatmap(a, order)
    a.set_ylabel(f'{len(order):,} orthologs', labelpad=3)
    heatmap(b, retained.index, band=cutoff_bands, band_width=.37)
    b.set_ylabel(f'{len(retained)} orthologs', labelpad=3)
    switch_panel(c, switches)
    score_panel(d, tubules).to_csv(draft / 'panel_d_binned_scores.csv')
    heatmap(e, s2_order, band=s2_bands, guides=(.36, .65), band_width=.54)
    e.set_ylabel(f'{len(s2_order)} orthologs', labelpad=3)
    heatmap(f, gain_order, band=gain_bands, band_width=1.04)
    f.set_ylabel(f'{len(gain_order)} orthologs', labelpad=3)
    lysosome_panel(g, lysosome, early_lysosome)
    metabolic_panel(h, x, raw_x, centered_fit, centered_raw, genes)

    color_axis = axes_mm(fig, 30, 8.0, 23, 2.0)
    scale = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(-1.5, 1.5), cmap=HEAT)
    cb = fig.colorbar(scale, cax=color_axis, orientation='horizontal', ticks=[-1.5, 0, 1.5])
    cb.ax.tick_params(labelsize=5.5, length=1.6, pad=1.2)
    cb.outline.set_linewidth(.45)
    fig.text(11 / WIDTH_MM, 9 / HEIGHT_MM, 'Heatmaps', va='center', fontsize=6.2, weight='bold')
    fig.text(56 / WIDTH_MM, 9 / HEIGHT_MM, 'Centered log₂ rate', va='center', fontsize=6.2)
    fig.text(11 / WIDTH_MM, 16.3 / HEIGHT_MM, 'COAST position', va='center', fontsize=6.4)
    fig.legend([Line2D([], [], c=color, lw=1.4) for color in COLORS],
               ['Mouse · 2 kidneys', 'Human · 2 sections, 1 donor'],
               loc='center right', bbox_to_anchor=(175 / WIDTH_MM, 9 / HEIGHT_MM),
               ncol=2, borderaxespad=0, handlelength=1.6, fontsize=6.2)
    return fig, axes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results-root', type=Path,
                        default=os.environ.get('PSEUDOSPACE_RESULTS_ROOT', REPO / 'results'))
    parser.add_argument('--s2-input-dir', type=Path)
    parser.add_argument('--audit-layout', action='store_true')
    args = parser.parse_args()
    root = args.results_root / 'paper_figures/fig5'
    draft = root / 'story_draft'
    draft.mkdir(parents=True, exist_ok=True)
    mpl.rcParams.update({
        'font.family': 'sans-serif', 'font.sans-serif': ['DejaVu Sans'],
        'font.size': 6.5, 'pdf.fonttype': 42, 'svg.fonttype': 'none',
        'axes.linewidth': .56, 'axes.spines.top': False, 'axes.spines.right': False,
        'axes.labelpad': 2.4, 'xtick.major.size': 2, 'ytick.major.size': 2,
        'xtick.major.pad': 2, 'ytick.major.pad': 2, 'xtick.major.width': .55,
        'ytick.major.width': .55, 'legend.frameon': False,
    })
    fig, axes = render(root, draft, args.s2_input_dir or draft / 's2_inputs')
    if args.audit_layout:
        from audit_panel_alignment import require_matplotlib_panel_alignment

        require_matplotlib_panel_alignment(
            fig, axes=axes, strict=True,
            row_groups=[list('ab'), list('cde'), list('fgh')],
            column_groups=[['c', 'f'], ['d', 'g'], ['e', 'h']],
            json_out=draft / 'fig5_story_draft.alignment.json',
        )
    fig.savefig(draft / 'fig5_story_draft.pdf', dpi=600)
    fig.savefig(draft / 'fig5_story_draft.png', dpi=300)
    fig.savefig(draft / 'fig5_story_draft.svg', dpi=600)
    fig.savefig(draft / 'fig5_story_draft.tiff', dpi=600, pil_kwargs={'compression': 'tiff_lzw'})
    plt.close(fig)


if __name__ == '__main__':
    main()
