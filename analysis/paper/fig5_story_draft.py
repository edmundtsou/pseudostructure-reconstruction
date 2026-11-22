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
from fig5_three_row_figure import HEAT
from fig5_urate import REPO, SECTIONS

WIDTH_MM = 183
HEIGHT_MM = 182
SPECIES = ('Mouse', 'Human')
COLORS = ('#275D83', '#8B3E6C')
SECTION_STYLES = ('-', '--', '-', '--')
S2_CLASSES = {
    'S1 high, off at S2': ('Off', '#47798B'),
    'S2 peak': ('Peak', '#B16428'),
    'S2 trough': ('Trough', '#999999'),
    'on at S2': ('On', '#65527E'),
}
SWITCH_DISPLACEMENT = .10
SWITCH_COLORS = ('#BFC3C8', '#AF602A')
S2_BAND = (.36, .65)
LATE_GENES = ('Gpt', 'Dcxr', 'Ggact', 'Pter', 'Abhd14a')


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


def screen_panel(ax, meta):
    mouse, human = meta.mouse_min_range.ge(1), meta.human_min_range.ge(1)
    names = ('Mouse only', 'Both', 'Human only')
    counts = [int(test.sum()) for test in (mouse & ~human, mouse & human, ~mouse & human)]
    y = np.array([2., 1., 0.])
    ax.barh(y, counts, height=.45, color=[COLORS[0], '.48', COLORS[1]], linewidth=0)
    for yy, label, n in zip(y, names, counts):
        ax.text(-7, yy, label, ha='right', va='center', fontsize=6.2)
        ax.text(n + 5, yy, str(n), ha='left', va='center', fontsize=6.3, fontweight='bold')
    ax.text(-88, 3.25, f'{len(meta):,} orthologs screened', fontsize=7, fontweight='bold')
    ax.text(-88, 2.81, f'{sum(counts)} with ≥2-fold positional range', fontsize=6.1, color='.32')
    ax.set(xlim=(-88, 229), ylim=(-.45, 3.7))
    ax.set_axis_off()


def switch_panel(ax, switches):
    bins = np.linspace(0, 1, 21)
    distance = (switches.human_switch_plot - switches.mouse_switch).abs().round(3)
    moved = distance.ge(SWITCH_DISPLACEMENT)
    groups = (~moved, moved)
    # Bars point away from one baseline so both species can share one positional scale.
    for sign, col in ((1, 'mouse_switch'), (-1, 'human_switch_plot')):
        baseline = np.zeros(len(bins) - 1)
        for group, color in zip(groups, SWITCH_COLORS):
            counts, _ = np.histogram(switches.loc[group, col], bins=bins)
            ax.bar(bins[:-1], sign * counts, width=np.diff(bins)*.91, align='edge',
                   bottom=sign * baseline, color=color, lw=0)
            baseline += counts
    ax.axhline(0, color='.27', lw=.6)
    ax.set(xlim=(0, 1), ylim=(-26, 26), xticks=[0, .25, .5, .75, 1],
           yticks=[-24, -12, 0, 12, 24], yticklabels=['24', '12', '0', '12', '24'],
           ylabel='Genes per bin', xlabel='Switch position')
    ax.text(.04, .84, 'Mouse', transform=ax.transAxes, c=COLORS[0], fontsize=6.5, fontweight='bold')
    ax.text(.04, .09, 'Human', transform=ax.transAxes, c=COLORS[1], fontsize=6.5, fontweight='bold')
    cutoff = f'{SWITCH_DISPLACEMENT:.2f}'
    labels = (f'Similar (<{cutoff}), {(~moved).sum()} genes', f'Moved (≥{cutoff}), {moved.sum()} genes')
    handles = [Line2D([], [], marker='s', color=color, ls='', markersize=3.7, label=label)
               for color, label in zip(SWITCH_COLORS, labels)]
    ax.legend(handles=handles, loc='upper center', bbox_to_anchor=(.50, -.195), ncol=2,
              fontsize=5.7, handletextpad=.4, handlelength=.8, columnspacing=.9)


def shade_s2(ax):
    ax.axvspan(*S2_BAND, color='#F0EBDF', lw=0, zorder=0)


def stage_panel(ax, tubules):
    rows = []
    shade_s2(ax)
    for j, section in enumerate(SECTIONS):
        frame = tubules.loc[tubules['sample'].eq(section), ['E', 'L', 'xq']].copy()
        frame['bin'] = pd.qcut(frame.xq, 20, labels=False)
        frame['early_only'] = frame.E.gt(0) & frame.L.le(0)
        summary = frame.groupby('bin').agg(position=('xq', 'mean'), fraction=('early_only', 'mean'),
                                           tubules=('early_only', 'size'))
        rows.append(summary.assign(section=section))
        ax.plot(summary.position, 100*summary.fraction, ls=SECTION_STYLES[j],
                color=COLORS[j // 2], lw=1.45)
    ax.text(np.mean(S2_BAND), 104, 'Mouse S2', ha='center', va='top', color='#6B5638', fontsize=6)
    ax.set(xlim=(0, 1), ylim=(0, 106), xticks=[0, .25, .5, .75, 1], yticks=[0, 25, 50, 75, 100],
           xlabel='COAST position', ylabel='Tubules (%)')
    return pd.concat(rows)


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


def lysosome_panel(ax, lysosome):
    delta = (lysosome.observed_human - lysosome.observed_mouse).sort_values()
    ranks = np.arange(1, len(delta) + 1)
    ax.axhline(0, color='.45', lw=.55, zorder=0)
    ax.bar(ranks, delta.values, color=np.where(delta.gt(0), COLORS[1], '.62'), lw=0)
    ax.text(.12, .84, f'{delta.gt(0).sum()} of {len(delta)}\nmore early-weighted',
            transform=ax.transAxes, va='top', fontsize=6.2, color=COLORS[1])
    ax.set(xlim=(0, len(ranks) + 1), ylim=(-1, 4.5), xticks=[1, 25, 50], yticks=[-1, 0, 2, 4],
           xlabel='Lysosomal genes, ranked', ylabel='Change in early bias (log₂)')


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


def profile_panel(ax, x, raw_x, fit, raw, genes, faint=False, s2=False):
    if s2:
        shade_s2(ax)
    ax.axhline(0, color='.68', lw=.55, zorder=1)
    for i, color in enumerate(COLORS):
        if faint:
            ax.plot(x, fit[i, genes].T, color=color, lw=.50, alpha=.20)
        ax.plot(x, fit[i, genes].mean(0), color=color, lw=1.45)
        ax.plot(raw_x, raw[i, genes].mean(0), ls='', marker='o', color=color,
                markersize=1.6, markeredgewidth=0, alpha=.68)
    ax.set(xlim=WINDOW, xticks=[.25, .5, .75], xticklabels=['.25', '.50', '.75'])


def metabolic_panel(ax, x, raw_x, pooled_fit, pooled_obs, genes):
    keys = ['Gpt', 'Gpt2']
    gene_index = genes.get_indexer(keys)
    if (gene_index < 0).any():
        raise ValueError('Missing GPT or GPT2 in the ortholog screen')
    for i, color in enumerate(COLORS):
        for j, ls in enumerate(('-', '--')):
            ax.plot(x, pooled_fit[i, gene_index[j]], color=color, ls=ls, lw=1.23)
            ax.plot(raw_x, pooled_obs[i, gene_index[j]], color=color, ls='',
                    marker='o' if j == 0 else 's', markersize=1.6, markeredgewidth=0, alpha=.60)
    ax.axhline(0, c='.68', lw=.55, zorder=0)
    legend = [Line2D([], [], ls=ls, c='.23', lw=.95, label=key.upper())
              for ls, key in zip(('-', '--'), keys)]
    ax.legend(handles=legend, loc='upper left', fontsize=5.4, handlelength=1.35, labelspacing=.24)
    ax.set(xlim=WINDOW, xticks=[.25, .5, .75], xticklabels=['.25', '.50', '.75'],
           yticks=[-1, 0, 1])


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
    gain.to_csv(draft / 'panel_human_gain_genes.csv')
    lysosome.assign(highlighted=lysosome.index.isin(early_lysosome)).to_csv(draft / 'panel_f_lysosome.csv')
    s2_genes.loc[s2_order].to_csv(draft / 'panel_d_s2_gene_order.csv')

    fig = plt.figure(figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4), facecolor='white')
    a = axes_mm(fig, 13, 130, 41, 35)
    b = axes_mm(fig, 64, 130, 54, 35)
    c = axes_mm(fig, 131, 130, 43, 35)
    d = axes_mm(fig, 13, 66, 74, 42)
    modules = [axes_mm(fig, l, 66, 20, 42) for l in (104, 129, 154)]
    f, g, h, i = bottom = [axes_mm(fig, l, 17, 32, 33) for l in (13, 56, 99, 142)]
    axes = [a, b, c, d, *modules, *bottom]
    for ax, letter in zip(axes, ('a', 'b', 'c', 'd', 'e1', 'e2', 'e3', 'f', 'g', 'h', 'i')):
        ax.set_gid(letter)
    for ax, letter, title in zip((a, b, c, d, modules[0], f, g, h, i), 'abcdefghi',
                                 ('Positional variation', 'Switches along the proximal tubule',
                                  'Early program only', 'Mouse S2 genes', 'Mouse-defined S2 profiles',
                                  'Lysosomal shift', 'Early lysosomal genes', 'Late metabolic genes',
                                  'GPT and GPT2')):
        heading(ax, letter, title, lift=13 if ax in (d, modules[0]) else 8)

    screen_panel(a, meta)
    switch_panel(b, switches)
    stage_panel(c, tubules).to_csv(draft / 'panel_c_early_only.csv')
    order = genes.get_indexer(s2_order)
    if (order < 0).any():
        raise ValueError('S2 heatmap rows are missing from the matched ortholog set')
    paired_map(d, centered_fit[:, order], x, band=s2_bands, guides=S2_BAND, band_width=.44)
    d.set_ylabel(f'{len(s2_order)} orthologs', labelpad=3)
    d.set_xlabel('COAST position', labelpad=2)
    for ax, cls in zip(modules, ('S1 high, off at S2', 'S2 peak', 'on at S2')):
        group = s2_genes.index[s2_genes['class'].eq(cls)]
        name, color = S2_CLASSES[cls]
        profile_panel(ax, x, raw_x, centered_fit, centered_raw, genes.get_indexer(group), s2=True)
        ax.set_title(f'{name} · {len(group)}', fontsize=6.3, color=color, pad=2)
        ax.set(ylim=(-1.3, 1.45), yticks=[-1, 0, 1])
    modules[0].set_ylabel('Centered log₂ rate')
    for ax in modules[1:]:
        ax.tick_params(labelleft=False)
    modules[1].set_xlabel('COAST position', labelpad=2)

    lysosome_panel(f, lysosome)
    profile_panel(g, x, raw_x, centered_fit, centered_raw, genes.get_indexer(early_lysosome))
    profile_panel(h, x, raw_x, centered_fit, centered_raw, genes.get_indexer(LATE_GENES), faint=True)
    metabolic_panel(i, x, raw_x, centered_fit, centered_raw, genes)
    for ax, count in ((g, len(early_lysosome)), (h, len(LATE_GENES))):
        ax.text(.08, .93, f'{count} genes', transform=ax.transAxes, va='top', color='.29', fontsize=6)
    for ax in (g, h, i):
        ax.set(xlabel='COAST position', ylim=(-1.52, 1.88), yticks=[-1, 0, 1],
               ylabel='Centered log₂ rate' if ax is g else '')
    for ax in (h, i):
        ax.tick_params(labelleft=False)

    color_axis = axes_mm(fig, 28, 6.2, 22, 1.8)
    scale = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(-1.5, 1.5), cmap=HEAT)
    cb = fig.colorbar(scale, cax=color_axis, orientation='horizontal', ticks=[-1.5, 0, 1.5])
    cb.ax.tick_params(labelsize=5.5, length=1.5, pad=1.0)
    cb.outline.set_linewidth(.45)
    fig.text(13 / WIDTH_MM, 7.2 / HEIGHT_MM, 'Heatmap', va='center', fontsize=6, weight='bold')
    fig.text(52 / WIDTH_MM, 7.2 / HEIGHT_MM, 'Centered log₂ rate', va='center', fontsize=6)
    fig.legend([Line2D([], [], c=color, lw=1.4) for color in COLORS],
               ['Mouse · 2 kidneys', 'Human · 2 sections, 1 donor'],
               loc='center right', bbox_to_anchor=(174 / WIDTH_MM, 7.2 / HEIGHT_MM),
               ncol=2, borderaxespad=0, handlelength=1.5, fontsize=6.2)
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
            row_groups=[['a', 'b', 'c'], ['d', 'e1', 'e2', 'e3'], list('fghi')],
            exemptions=[
                {'panels': ['a', 'b', 'c'], 'checks': ['panel-width', 'horizontal-gutter'],
                 'reason': 'The shared-axis switch histogram receives extra width.'},
                {'panels': ['d', 'e1', 'e2', 'e3'], 'checks': ['panel-width', 'horizontal-gutter'],
                 'reason': 'One paired heatmap beside an equally spaced triptych.'},
            ],
            json_out=draft / 'fig5_story_draft.alignment.json',
        )
    fig.savefig(draft / 'fig5_story_draft.pdf', dpi=600)
    fig.savefig(draft / 'fig5_story_draft.png', dpi=300)
    fig.savefig(draft / 'fig5_story_draft.svg', dpi=600)
    fig.savefig(draft / 'fig5_story_draft.tiff', dpi=600, pil_kwargs={'compression': 'tiff_lzw'})
    plt.close(fig)


if __name__ == '__main__':
    main()
