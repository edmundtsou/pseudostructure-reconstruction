"""Compare specified proximal-tubule biological stories using the study's four Visium HD sections."""
from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd

from fig5_urate import REPO, SECTIONS
from fig5_urate_tissue import overlap_table, raw_endpoint_contrasts, tissue_units

PROGRAMS = {
    'Urate': ('Slc17a3', 'Slc22a12', 'Slc17a1', 'Slc22a6', 'Slc22a7', 'Slc22a8', 'Slc22a11', 'Slc2a9'),
    'Glucuronate–xylulose': ('Akr1a1', 'Cryl1', '4931406C07Rik', 'Dcxr', 'Sord', 'Xylb'),
    'Kynurenine–NAD': ('Kmo', 'Kynu', 'Haao', 'Acmsd', 'Qprt', 'Aldh8a1'),
    'Apical peptide processing': ('Mme', 'Dpp4', 'Anpep', 'Cubn', 'Lrp2'),
    'Prostaglandin clearance': ('Slco2a1', 'Hpgd', 'Ptgr1', 'Ptgr2'),
    'Glucose': ('Slc5a2', 'Slc5a1', 'Pdzk1ip1', 'Slc2a2'),
    'Phosphate': ('Slc34a1', 'Slc34a3', 'Slc20a2', 'Xpr1'),
    'Neutral amino acids': ('Slc6a19', 'Slc6a18', 'Slc36a2'),
    'Lactate': ('Slc5a12', 'Slc5a8', 'Slc16a1'),
    'Organic anions': ('Slc22a6', 'Slc22a8', 'Abcc2', 'Abcc4', 'Slc22a2'),
}

# Pairs connect opposing transport roles, consecutive enzymes, or a substrate and branch step.
PAIRS = {
    'Urate': [('Slc17a3', 'Slc22a12')],
    'Glucuronate–xylulose': [('Dcxr', 'Sord'), ('Sord', 'Xylb'), ('Dcxr', 'Xylb')],
    'Kynurenine–NAD': [('Haao', 'Acmsd'), ('Acmsd', 'Qprt')],
    'Apical peptide processing': [('Mme', 'Dpp4'), ('Mme', 'Anpep')],
    'Prostaglandin clearance': [('Slco2a1', 'Hpgd'), ('Hpgd', 'Ptgr1')],
    'Glucose': [('Slc5a2', 'Slc5a1'), ('Slc5a2', 'Pdzk1ip1'), ('Slc5a2', 'Slc2a2')],
    'Phosphate': [('Slc34a1', 'Slc34a3'), ('Slc34a1', 'Xpr1')],
    'Neutral amino acids': [('Slc6a19', 'Slc6a18')],
    'Lactate': [('Slc5a12', 'Slc5a8')],
    'Organic anions': [('Slc22a6', 'Abcc2'), ('Slc22a8', 'Abcc2')],
}

PLOT_PAIRS = {
    ('Slc17a3', 'Slc22a12'): 'NPT4–URAT1',
    ('Slc22a8', 'Abcc2'): 'OAT3–MRP2',
    ('Mme', 'Anpep'): 'MME–ANPEP',
    ('Mme', 'Dpp4'): 'MME–DPP4',
    ('Haao', 'Acmsd'): 'HAAO–ACMSD',
    ('Dcxr', 'Sord'): 'DCXR–SORD',
    ('Slc5a2', 'Slc5a1'): 'SGLT2–SGLT1',
}
WIDTH_MM = 183


def comparison(results_root):
    root = results_root / 'paper_figures/fig5'
    screen = root / 'screen'
    genes = pd.read_csv(screen / 'genes.csv', index_col='gene')
    heldout = pd.read_csv(root / 'screen_genehalf_2026-11-18/genes.csv', index_col='gene')
    axes = pd.read_csv(screen / 'axes.csv').iloc[0]
    membership = pd.DataFrame((story, gene) for story, members in PROGRAMS.items() for gene in members)
    membership.columns = ['story', 'gene']
    available = sorted(set(membership.gene) & set(genes.index))
    units, _ = tissue_units(results_root, targets=available)
    contrasts = raw_endpoint_contrasts(units, axes, available)
    endpoint_cols = [f'R_{s}' for s in SECTIONS]
    np.testing.assert_allclose(
        contrasts[[f'cores_{s}' for s in SECTIONS]], genes.loc[available, endpoint_cols], atol=1e-5,
    )
    heldout_stats = heldout[endpoint_cols + ['assessable', 'catR']].add_prefix('heldout_')
    primary = genes[[f'detected_{s}' for s in SECTIONS] + ['assessable', 'catR']]
    gene_table = membership.join(pd.concat([primary, heldout_stats, contrasts], axis=1), on='gene')

    overlap_rows = []
    with np.load(screen / 'all_curves.npz', allow_pickle=True) as saved:
        curve_index = pd.Index(saved['genes'])
        for story, pairs in PAIRS.items():
            for pair in pairs:
                if any(gene not in available for gene in pair):
                    continue
                overlaps, _ = overlap_table(units, genes=pair)
                matched, _ = overlap_table(units, genes=pair, match_mito=True)
                keys = ['section', 'fraction']
                overlaps = overlaps.merge(
                    matched[keys + ['observed_expected']].rename(
                        columns={'observed_expected': 'matched_mito_observed_expected'},
                    ), on=keys, validate='one_to_one',
                )
                rows = curve_index.get_indexer(pair)
                for section in SECTIONS:
                    observed = saved[f'bins|{section}'][rows]
                    supported = np.isfinite(observed).all(axis=0)
                    values = observed[:, supported]
                    r = np.corrcoef(values)[0, 1] if values.shape[1] >= 3 and np.all(values.std(axis=1) > 0) else np.nan
                    mask = overlaps.section.eq(section)
                    overlaps.loc[mask, 'bin_pearson'] = r
                    overlaps.loc[mask, 'finite_bins'] = supported.sum()
                overlap_rows.append(overlaps.assign(story=story, gene_a=pair[0], gene_b=pair[1]))
    return gene_table, pd.concat(overlap_rows, ignore_index=True)


def render(pairs):
    """Show raw position-bin correlation and depth/mitochondrial-matched tubule overlap."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.transforms import ScaledTranslation

    from fig5_urate_figure import MARKERS, SECTION_COLORS, SECTION_NAMES

    plt.rcParams.update({
        'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'], 'font.size': 8,
        'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': 0.6,
        'legend.frameon': False, 'pdf.fonttype': 42, 'svg.fonttype': 'none',
    })
    fig, axs = plt.subplots(1, 2, figsize=(WIDTH_MM / 25.4, 110 / 25.4))
    fig.subplots_adjust(left=0.19, right=0.98, top=0.78, bottom=0.22, wspace=0.24)
    fig.text(0.19, 0.96, 'Cross-species comparison of candidate gene pairs', fontsize=12, weight='bold', va='top')
    fig.legend([Line2D([], [], color=color, ls='', marker=marker, ms=4)
                for color, marker in zip(SECTION_COLORS, MARKERS)],
               SECTION_NAMES, ncol=4, fontsize=7, loc='upper left', bbox_to_anchor=(0.185, 0.90),
               columnspacing=1.1, handletextpad=0.2)
    points = pairs[pairs.fraction.eq(0.25)].set_index(['gene_a', 'gene_b', 'section'])
    panels = [('bin_pearson', 'Observed position bins', 'Pearson correlation', 0, (-1.06, 1.06)),
              ('matched_mito_observed_expected', 'Joint-high tubule profiles', 'Observed / expected', 1, (0, 1.78))]
    y = np.arange(len(PLOT_PAIRS))
    for letter, ax, (column, title, xlabel, zero, limits) in zip('ab', axs, panels):
        ax.axvline(zero, color='0.6', lw=0.65, zorder=0)
        for row in y:
            ax.axhline(row, color='0.91', lw=0.45, zorder=0)
        for i, section in enumerate(SECTIONS):
            values = [points.loc[(a, b, section), column] for a, b in PLOT_PAIRS]
            ax.scatter(values, y + (i - 1.5) * 0.15, color=SECTION_COLORS[i], marker=MARKERS[i], s=15, lw=0)
        ax.set(xlim=limits, ylim=(len(y) - 0.45, -0.55), yticks=y, xlabel=xlabel)
        ax.tick_params(axis='y', length=0)
        ax.set_title(title, loc='left', fontsize=8.5, weight='bold', pad=8)
        offset = ScaledTranslation(-14 / 72, 8 / 72, fig.dpi_scale_trans)
        ax.text(0, 1, letter, transform=ax.transAxes + offset, weight='bold', fontsize=9, va='bottom')
    axs[0].set_yticklabels(list(PLOT_PAIRS.values()))
    axs[1].set_yticklabels([])
    fig.text(0.19, 0.035, 'Panel a uses bins with expression detected for both transcripts.\n'
             'Panel b uses the top 25% within depth and mitochondrial-fraction strata.', fontsize=7, va='bottom')
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results-root', type=Path,
                        default=os.environ.get('PSEUDOSPACE_RESULTS_ROOT', REPO / 'results'))
    parser.add_argument('--audit-layout', action='store_true',
                        help='run the layout audit with the nature-figure scripts on PYTHONPATH')
    args = parser.parse_args()
    export = args.results_root / 'paper_figures/fig5/story_comparison'
    export.mkdir(parents=True, exist_ok=True)
    genes, pairs = comparison(args.results_root)
    genes.to_csv(export / 'gene_contrasts.csv', index=False)
    pairs.to_csv(export / 'pair_profiles.csv', index=False)
    fig = render(pairs)
    if args.audit_layout:
        from audit_panel_alignment import require_matplotlib_panel_alignment

        require_matplotlib_panel_alignment(
            fig, json_out=export / 'story_comparison.alignment.json', require_panel_labels=True, strict=True,
        )
    fig.savefig(export / 'story_comparison.pdf')
    fig.savefig(export / 'story_comparison.svg')
    fig.savefig(export / 'story_comparison.png', dpi=300)
    fig.savefig(export / 'story_comparison.tiff', dpi=600, pil_kwargs={'compression': 'tiff_lzw'})


if __name__ == '__main__':
    main()
