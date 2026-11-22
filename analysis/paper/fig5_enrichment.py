"""Compare cluster ORA with preranked GSEA of Fig. 5 positional profiles."""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.stats import false_discovery_control

from fig5_curve_audit import WINDOW, cross, curve_arrays, shape
from fig5_urate import REPO, SECTIONS

sys.path.insert(0, str(REPO))
from pseudospace.modules import enrich_modules  # noqa: E402
from pseudospace.pathway_pipelines import rank_first_gsea  # noqa: E402

MAX_SIZE = 300
WIDTH_MM = 183
HEIGHT_MM = 76
RANGE_PROGRAMS = {
    'Peroxisome': 'KEGG_2019_Mouse::Peroxisome',
    'Peroxisomal protein import': 'R-HSA-9033241',
    'Fatty-acid metabolism': 'R-HSA-8978868',
    'SLC transport': 'R-HSA-425407',
    'Amino-acid transport': 'R-HSA-352230',
    'Protein digestion': 'KEGG_2019_Mouse::Protein digestion and absorption',
    'Lysosome': 'KEGG_2019_Mouse::Lysosome',
}


def local_sets(sets, background, minimum):
    background = set(background)
    mapped = {key: sorted(set(members) & background) for key, members in sets.items()}
    return {key: members for key, members in mapped.items()
            if minimum <= len(members) <= MAX_SIZE and len(members) < len(background)}


def run_ora(table, fit, sets):
    """Test eligible cluster backgrounds separately from shape or large amplitude losses."""
    cluster = table['cluster'].dropna()
    both = table.loc[table.reliable_both]
    ranges = np.ptp(fit, axis=-1)
    loss = ranges[2:].max(0) <= .5 * ranges[:2].min(0)
    gain = ranges[:2].max(0) <= .5 * ranges[2:].min(0)
    screens = [
        ('clusters', cluster.index,
         {f'Cluster {i}': cluster.index[cluster.eq(i)] for i in sorted(cluster.astype(int).unique())}, 5),
        ('shape_class', both.index, {'Conserved': both.index[both.conserved],
                                    'Divergent': both.index[both.divergent]}, 3),
        ('twofold_range_loss', cluster.index,
         {'Shallower in human': table.index[loss & table.mouse_both_grade].intersection(cluster.index),
          'Shallower in mouse': table.index[gain & table.human_both_grade].intersection(cluster.index)}, 5),
    ]
    results = []
    for family, universe, foregrounds, minimum in screens:
        test = enrich_modules(foregrounds, local_sets(sets, universe, minimum), background=universe)
        test.insert(0, 'family', family)
        test.insert(1, 'background_size', len(universe))
        results.append(test)
    return pd.concat(results, ignore_index=True)


def range_adjusted_loss(mouse, human):
    """Compare human amplitude with the typical amplitude at a similar mouse range."""
    bands = np.array_split(np.argsort(mouse, kind='stable'), max(1, min(32, len(mouse) // 20)))
    reference = np.array([np.median(mouse[band]) for band in bands])
    quartiles = np.array([np.quantile(human[band], [.25, .5, .75]) for band in bands])
    expected = np.interp(mouse, reference, quartiles[:, 1])
    spread = np.interp(mouse, reference, quartiles[:, 2] - quartiles[:, 0])
    if (spread <= 0).any() or (np.diff(reference) <= 0).any():
        raise ValueError('The mouse amplitude strata need increasing centers and a positive human IQR')
    return (expected - human) / spread


def rank_scores(fit, x):
    amplitude = np.ptp(fit, axis=-1)
    mouse, human = amplitude[:2].mean(0), amplitude[2:].mean(0)
    centered = fit - fit.mean(axis=-1, keepdims=True)
    shared_scale = centered.std(axis=-1).max(0)
    w = (x.mean() - x) / np.abs(x - x.mean()).sum()
    delta = centered[2:].mean(0) - centered[:2].mean(0)
    unit = shape(fit)
    return {
        'human_range_loss': mouse - human,
        'human_early_bias': (delta @ w) / shared_scale,
        'shape_conservation': cross(unit[:2], unit[2:]).min(0),
        'mouse_range_adjusted_loss': range_adjusted_loss(mouse, human),
        'unit_variance_early_bias': (unit[2:].mean(0) - unit[:2].mean(0)) @ w,
    }


def run_gsea(table, shaped, inputs, sets, permutations, threads):
    retained = table.cluster.notna()
    screens = [
        ('all_screened', 'human_range_loss', np.ones(len(table), bool), 10),
        ('all_screened', 'mouse_range_adjusted_loss', np.ones(len(table), bool), 10),
        ('twofold_either', 'human_range_loss', retained, 5),
        ('twofold_either', 'human_early_bias', retained, 5),
        ('twofold_both', 'shape_conservation', table.reliable_both, 3),
        ('graded_both', 'unit_variance_early_bias', shaped, 5),
    ]
    tests = []
    ranks = []
    for source, grid, profiles in inputs:
        scores = rank_scores(profiles, grid)
        for universe, statistic, keep, minimum in screens:
            score = pd.Series(scores[statistic], index=table.index).loc[keep]
            test = rank_first_gsea(score, local_sets(sets, score.index, minimum), metric='rank_normal',
                                   min_size=minimum, max_size=MAX_SIZE, permutation_num=permutations,
                                   seed=2026, threads=threads)
            test = test.drop(columns=['p_one_sided', 'call'])
            test['leading_edge'] = test.leading_edge.map('; '.join)
            test = test.assign(source=source, universe=universe, statistic=statistic)
            tests.append(test)
            ranks.append(score.rename('score').reset_index().assign(source=source, universe=universe,
                                                                   statistic=statistic))
            print(source, universe, statistic, len(test), 'sets; within-rank q ≤ 0.05:', int(test.q.le(.05).sum()),
                  flush=True)
    tests = pd.concat(tests, ignore_index=True)
    # A zero nominal p means no sampled statistic exceeded the observed value.
    # The floor also makes the joint BH calculation finite at the permutation resolution.
    tests['joint_bh_q'] = tests.groupby('source').p.transform(
        lambda values: false_discovery_control(values.clip(lower=1 / (permutations + 1))))
    return tests, pd.concat(ranks, ignore_index=True)


def plot_amplitude_results(gsea, output, *, audit_layout=False):
    """Show results for a named set of programs across three profile sources."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 7, 'pdf.fonttype': 42,
                         'svg.fonttype': 'none', 'axes.spines.right': False,
                         'axes.spines.top': False, 'axes.linewidth': .55})
    fig, axes = plt.subplots(1, 2, figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4), sharex=True, sharey=True,
                             gridspec_kw={'left': .27, 'right': .975, 'top': .725, 'bottom': .25, 'wspace': .24})
    names = list(RANGE_PROGRAMS)
    markers = {'fit': ('o', -.22, 'Fitted profiles'),
               'gene_held_out': ('D', 0, 'Gene-held-out fits'),
               'observed': ('^', .22, 'Observed bins')}
    ranks = ('human_range_loss', 'mouse_range_adjusted_loss')
    captions = ('Absolute range loss', 'Range loss at similar mouse strength')
    all_genes = gsea.loc[gsea.universe.eq('all_screened')]
    color = '#31597C'
    for ax, rank, caption in zip(axes, ranks, captions):
        ax.set_title(caption, loc='left', weight='bold', fontsize=7.8, pad=6)
        ax.axvline(0, color='.55', lw=.65, ls='--')
        for y, key in enumerate(RANGE_PROGRAMS.values()):
            if y:
                ax.axhline(y - .5, color='.89', lw=.5)
            rows = all_genes.loc[all_genes.statistic.eq(rank) & all_genes.pathway_id.str.endswith(key)]
            if rows.source.duplicated().any() or set(rows.source) != set(markers):
                raise ValueError(f'Expected one pathway result per profile source for {key}')
            for row in rows.itertuples(index=False):
                marker, offset, _ = markers[row.source]
                ax.scatter(row.NES, y + offset, marker=marker, s=13, lw=.65,
                           facecolors=color if row.joint_bh_q <= .05 else 'white', edgecolors=color, zorder=3)
        ax.set(xlim=(-3.1, 3.1), xticks=[-3, -2, -1, 0, 1, 2, 3], ylim=(len(names) - .5, -.5))
        ax.spines['left'].set_visible(False)
        ax.tick_params(axis='y', length=0, pad=4)
    axes[0].set_yticks(np.arange(len(names)), names)
    fig.text(.055, .954, 'Which programs retain positional variation?', weight='bold', fontsize=10.7, va='top')
    fig.text(.055, .825, 'Ranked comparison across all 3,912 screened orthologs', fontsize=7.2, color='.28')
    fig.text(.621, .149, 'Normalized enrichment score', fontsize=7, ha='center')
    fig.text(.621, .103, '← More retained                          Greater loss →', fontsize=6.4, ha='center')
    fig.legend(handles=[Line2D([], [], color=color, marker=mark, lw=0, markersize=4, label=label)
                        for mark, _, label in markers.values()], loc='lower left',
               bbox_to_anchor=(.048, .035), frameon=False, ncol=3, fontsize=6.6,
               columnspacing=.9, handletextpad=.25)
    fig.text(.967, .062, 'Filled markers · joint q ≤ 0.05', fontsize=6.6, ha='right')
    if audit_layout:
        from audit_panel_alignment import require_matplotlib_panel_alignment

        require_matplotlib_panel_alignment(fig, strict=True, json_out=output / 'range_loss_comparison.alignment.json')
    fig.savefig(output / 'range_loss_comparison.pdf')
    fig.savefig(output / 'range_loss_comparison.svg')
    fig.savefig(output / 'range_loss_comparison.png', dpi=300)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results-root', type=Path,
                        default=os.environ.get('PSEUDOSPACE_RESULTS_ROOT', REPO / 'results'))
    parser.add_argument('--permutations', type=int, default=4000)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--audit-layout', action='store_true')
    args = parser.parse_args()
    if args.permutations < 1 or args.threads < 1:
        parser.error('Use positive permutation and thread counts')
    root = args.results_root / 'paper_figures/fig5'
    output = root / 'enrichment'
    output.mkdir(parents=True, exist_ok=True)
    sets = json.loads((root.parent / 'gene_sets.json').read_text())
    table = pd.read_csv(root / 'three_row/ortholog_profiles.csv', index_col='gene')
    prior = pd.read_csv(root / 'curve_audit/genes.csv', index_col='gene').loc[table.index]
    x = np.arange(WINDOW[0], WINDOW[1] + .001, .005).round(3)
    fit, _ = curve_arrays(root / 'screen', table.index, x)
    held, _ = curve_arrays(root / 'screen_genehalf_2026-11-18', table.index, x)
    with np.load(root / 'curve_audit/reaggregated_counts.npz') as raw:
        rows = pd.Index(raw['genes']).get_indexer(table.index)
        if (rows < 0).any() or raw['sections'].tolist() != SECTIONS:
            raise ValueError('The count audit does not match Fig. 5 genes or sections')
        xx = raw['centers']
        keep = (xx >= x[0]) & (xx <= x[-1])
        observed = raw['profiles'][:, rows][:, :, keep].astype(float)
    observed = gaussian_filter1d(observed, sigma=.75, axis=-1)
    run_ora(table, fit, sets).to_csv(output / 'ora.csv', index=False)
    gsea, ranks = run_gsea(table, prior.reliable_both,
                           [('fit', x, fit), ('gene_held_out', x, held), ('observed', xx[keep], observed)],
                           sets, args.permutations, args.threads)
    gsea.to_csv(output / 'gsea.csv', index=False)
    ranks.to_csv(output / 'gene_scores.csv', index=False)
    plot_amplitude_results(gsea, output, audit_layout=args.audit_layout)


if __name__ == '__main__':
    main()
