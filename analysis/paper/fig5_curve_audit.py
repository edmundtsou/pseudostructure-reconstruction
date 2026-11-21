"""Compare full proximal tubule profiles across the four Fig. 5 sections."""
from __future__ import annotations

import argparse
import json
import os
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from scipy.ndimage import gaussian_filter1d

from fig5_urate import REPO, SECTIONS

WINDOW = (.19, .875)
COLORS = ('#31597C', '#31597C', '#687D3D', '#687D3D')
STYLES = ('-', '--', '-', '--')
PROGRAMS = {
    'Early amino-acid transport': ('Slc6a19', 'Slc7a9', 'Slc3a2', 'Slc7a8', 'Slc7a7'),
    'Additional amino-acid transport': ('Slc36a1', 'Slc16a10', 'Slc3a1', 'Slc6a18', 'Slc6a6'),
    'Peptide processing': ('Mme', 'Dpp4', 'Dpep1', 'Enpep', 'Anpep'),
    'Carbohydrate metabolism': ('Pck1', 'Dcxr', 'Miox', 'Fbp1', 'Aldob'),
    'Peroxisome-associated genes': ('Cat', 'Acox2', 'Hao2', 'Amacr', 'Idh1'),
}


def curve_arrays(path, genes, x):
    with np.load(path / 'all_curves.npz', allow_pickle=True) as data:
        rows = pd.Index(data['genes']).get_indexer(genes)
        cols = pd.Index(data['grid'].round(3)).get_indexer(x)
        if (rows < 0).any() or (cols < 0).any():
            raise ValueError(f'Missing study genes or grid positions in {path}')
        fits = np.asarray([data[f'curve|{s}'][rows][:, cols] for s in SECTIONS], dtype=float)
        counts = np.asarray([data[f'bins|{s}'][rows] for s in SECTIONS], dtype=float)
        positions = np.asarray([data[f'bin_x|{s}'][rows] for s in SECTIONS], dtype=float)
    fits, counts = np.logaddexp2(fits, -1), np.logaddexp2(counts, -1)
    bins = np.asarray([[np.interp(x, p, c, left=np.nan, right=np.nan) for p, c in zip(ps, cs)]
                       for ps, cs in zip(positions, counts)])
    return fits, gaussian_filter1d(bins, sigma=5, axis=-1)


def shape(a):
    centered = a - a.mean(axis=-1, keepdims=True)
    sd = centered.std(axis=-1, keepdims=True)
    return np.divide(centered, sd, out=np.full_like(centered, np.nan), where=sd > 0)


def cross(a, b):
    return np.array([(u * v).mean(axis=-1) for u in a for v in b])


def within(a):
    return np.array([(a[0] * a[1]).mean(axis=-1), (a[2] * a[3]).mean(axis=-1)])


def classify(fit, raw, held_fit, held_raw, index, amplitude=.5):
    """Separate reliable shape differences from weak or inconsistent profiles."""
    amp = np.ptp(fit, axis=-1)
    strong = amp >= amplitude
    zf, zr, zh, zhr = map(shape, (fit, raw, held_fit, held_raw))
    wf, wr = within(zf).min(axis=0), within(zr).min(axis=0)
    eligible = strong.all(axis=0) & (wf >= .65) & (wr >= .55)
    cf, cr = cross(zf[:2], zf[2:]), cross(zr[:2], zr[2:])
    ch, held_raw_r = cross(zh[:2], zh[2:]), cross(zhr[:2], zhr[2:])
    same = eligible & (cf.min(0) >= .70) & (cr.min(0) >= .50)
    different = eligible & (cf.max(0) <= .30) & (cr.max(0) <= .45)
    conserved = same & (ch.min(0) >= .65) & (held_raw_r.min(0) >= .45)
    divergent = different & (ch.max(0) <= .40) & (held_raw_r.max(0) <= .50)
    return pd.DataFrame({
        'mouse_amplitude': amp[:2].mean(0), 'human_amplitude': amp[2:].mean(0),
        'mouse_both_grade': strong[:2].all(0), 'human_both_grade': strong[2:].all(0),
        'within_fit_min_r': wf, 'within_observed_min_r': wr, 'reliable_both': eligible,
        'cross_fit_min_r': cf.min(0), 'cross_fit_max_r': cf.max(0),
        'cross_observed_min_r': cr.min(0), 'cross_observed_max_r': cr.max(0),
        'cross_held_fit_min_r': ch.min(0), 'cross_held_fit_max_r': ch.max(0),
        'cross_held_observed_min_r': held_raw_r.min(0), 'cross_held_observed_max_r': held_raw_r.max(0),
        'primary_conserved': same, 'primary_divergent': different,
        'conserved': conserved, 'divergent': divergent,
    }, index=index)


def fixed_width_counts(results_root, genes):
    """Re-aggregate study counts on the same COAST axis, also within the lower human mito tertiles."""
    import anndata as ad
    from scipy import sparse

    from fig3_common import RUNS, load_positions

    run = results_root / RUNS['hvm']
    ref = load_positions(results_root, 'hvm')
    pt = ref.loc[ref['sample'].isin(SECTIONS) & ref['group'].isin(('PTS1', 'PTS2', 'PTS3'))]
    qc = ad.read_h5ad(run / '02_qc/units_qc.h5ad')
    counts = sparse.csr_matrix(qc[pt.index.to_numpy()].layers['counts'])
    gene_idx = qc.var_names.get_indexer(genes)
    if (gene_idx < 0).any():
        raise ValueError('Profile genes missing from COAST counts')
    mito_genes = qc.var['mt'].to_numpy(bool)
    lib = np.asarray(counts[:, ~mito_genes].sum(1)).ravel()
    mito_counts = np.asarray(counts[:, mito_genes].sum(1)).ravel()
    mito = np.divide(mito_counts, lib + mito_counts, out=np.zeros(len(lib)), where=lib + mito_counts > 0)
    y = counts[:, gene_idx]
    samples = pt['sample'].to_numpy()
    # Use the single-axis screen's mouse-reference quantiles for these study units.
    mouse = pt['sample'].isin(SECTIONS[:2])
    mouse_pos = pt.loc[mouse, 't'].sort_values()
    ranks = (mouse_pos.rank(method='average') - .5) / len(mouse_pos)
    x = np.interp(pt['t'], mouse_pos, ranks)
    ranges = np.array([np.quantile(x[samples == section], [.05, .95]) for section in SECTIONS])
    left, right = ranges[:, 0].max(), ranges[:, 1].min()
    if left >= right:
        raise ValueError('Four sections have no overlapping COAST support')
    edges = np.linspace(left, right, 21)
    centers = (edges[:-1] + edges[1:]) / 2
    profiles, n_units = [], []
    strata = [samples == section for section in SECTIONS]
    for human in strata[2:]:
        strata.append(human & (mito <= np.quantile(mito[human], 1 / 3)))
    for section in strata:
        sizes, rate = [], []
        for a, b in zip(edges[:-1], edges[1:]):
            mask = section & (x >= a) & ((x <= b) if b == edges[-1] else (x < b))
            sizes.append(mask.sum())
            total = lib[mask].sum()
            observed = np.asarray(y[mask].sum(0)).ravel()
            rate.append(np.log2(1e4 * observed / total + .5) if total else np.full(len(genes), np.nan))
        profiles.append(np.asarray(rate).T)
        n_units.append(sizes)
    return centers, np.asarray(profiles), np.asarray(n_units)


def annotate(table, memberships):
    fields = ['reliable_both', 'conserved', 'divergent']
    rows = []
    for name, members in memberships.items():
        subset = table.loc[table.index.intersection(members)]
        rows.append({'set': name, 'screened': len(subset), **subset[fields].sum().to_dict(),
                     **{f'{key}_genes': ','.join(subset.index[subset[key]]) for key in fields[1:]}})
    return pd.DataFrame(rows).sort_values(['reliable_both', 'set'], ascending=[False, True])


def observed_program_pairs(genes, full, low):
    rows = []
    for name, members in PROGRAMS.items():
        for a, b in combinations(members, 2):
            i, j = genes.get_loc(a), genes.get_loc(b)
            full_r = (full[:4, i] * full[:4, j]).mean(-1)
            low_r = (low[4:, i] * low[4:, j]).mean(-1)
            rows.append({'program': name, 'gene_a': a, 'gene_b': b,
                         **dict(zip(SECTIONS, full_r)),
                         **dict(zip((f'{s}_lowmito' for s in SECTIONS[2:]), low_r))})
    return pd.DataFrame(rows)


def overview(out, table, x, fit):
    """Display every well-supported gene, ordered by the shape of the full curves."""
    from sklearn.cluster import KMeans

    z = shape(fit)
    means = np.stack([z[:2].mean(0), z[2:].mean(0)])
    eligible = np.flatnonzero(table.reliable_both)
    joined = means[:, eligible, ::5].transpose(1, 0, 2).reshape(len(eligible), -1)
    km = KMeans(n_clusters=10, n_init=30, random_state=2026).fit(joined)
    group_mean = km.cluster_centers_.reshape(10, 2, -1)
    ranks = np.lexsort((group_mean[:, 1].argmax(1), group_mean[:, 0].argmax(1)))
    labels = np.argsort(ranks)[km.labels_]
    distance = np.square(joined - km.cluster_centers_[km.labels_]).mean(1)
    order = np.lexsort((distance, labels))
    members = table.iloc[eligible[order]].assign(shape_group=labels[order] + 1)
    members.to_csv(out / 'shape_group_order.csv')
    boundaries = np.flatnonzero(np.diff(labels[order])) + 1
    fig, axs = plt.subplots(1, 3, figsize=(7.6, 5.5), gridspec_kw={'width_ratios': [1, 1, .11]},
                            layout='constrained')
    fig.suptitle(f'Full-curve patterns in {len(eligible)} orthologs with reproducible position profiles')
    for i, ax in enumerate(axs[:2]):
        heat = ax.imshow(means[i, eligible[order]], aspect='auto', interpolation='nearest', origin='upper',
                         extent=(x[0], x[-1], len(order), 0), cmap='RdBu_r', vmin=-2, vmax=2)
        ax.set(title=('Mouse', 'Human')[i], xlabel='COAST position', yticks=[])
    flags = members[['conserved', 'divergent']].to_numpy(int) @ np.array([1, -1])
    status_colors = ('#8E1764', '#E3E3E3', '#3E682E')
    axs[2].imshow(flags[:, None], aspect='auto', interpolation='nearest', origin='upper',
                  extent=(0, 1, len(order), 0), vmin=-1, vmax=1, cmap=ListedColormap(status_colors))
    axs[2].set(title='Class', xticks=[], yticks=[])
    for ax in axs:
        for boundary in boundaries:
            ax.axhline(boundary, color='white', linewidth=.8)
    fig.colorbar(heat, ax=axs[:2], orientation='horizontal', shrink=.7, label='Standardized profile shape')
    fig.legend(handles=[Patch(color=color, label=label) for color, label in
                        zip(status_colors, ('Divergent', 'Intermediate', 'Conserved'))],
               loc='outside lower center', ncol=3, frameon=False)
    fig.savefig(out / 'supported_shapes.png', dpi=220)
    plt.close(fig)
    return members


def gene_card(ax, gene, row, x, fit, centers, raw, counts):
    offset = fit[:, row].mean(axis=-1)
    for i in range(4):
        ax.plot(x, fit[i, row] - offset[i], STYLES[i], color=COLORS[i], lw=1.1)
        ax.scatter(centers, raw[i, row] - offset[i], color=COLORS[i], s=5, alpha=.5,
                   marker='o' if i % 2 == 0 else 'D', edgecolor='none')
    for i in range(2):
        supported = counts[i + 4] >= 20
        low = raw[i + 4, row]
        full = raw[i + 2, row]
        aligned = low - low[supported].mean() + (full[supported] - offset[i + 2]).mean()
        ax.scatter(centers[supported], aligned[supported], facecolors='none', edgecolors=COLORS[i + 2],
                   s=9, linewidths=.5, marker='o' if i == 0 else 'D')
    ax.axhline(0, color='.85', lw=.55)
    ax.set(title=gene.upper(), xlim=(.14, .93), xticks=[.2, .5, .8])


def program_profiles(out, table, x, fit, centers, raw, counts, atlas_order):
    plt.rcParams.update({'font.size': 7.5, 'axes.spines.top': False, 'axes.spines.right': False,
                         'pdf.fonttype': 42})
    fig, axs = plt.subplots(len(PROGRAMS), 5, figsize=(10, 8.4), squeeze=False, layout='constrained')
    fig.suptitle('Conserved and rearranged proximal tubule expression profiles', fontsize=12)
    for (program, genes), axes in zip(PROGRAMS.items(), axs):
        axes[0].set_ylabel(program + '\nCentered log₂ rate')
        for gene, ax in zip(genes, axes):
            gene_card(ax, gene, table.index.get_loc(gene), x, fit, centers, raw, counts)
    fig.supxlabel('COAST position  |  Blue, mouse  ·  Green, human\n'
                  'Lines, fitted profiles  ·  Filled markers, counts  ·  Open markers, lower-mito human counts')
    for suffix in ('pdf', 'png'):
        fig.savefig(out / f'program_profiles.{suffix}', dpi=240)
    plt.close(fig)

    with PdfPages(out / 'supported_gene_atlas.pdf') as pdf:
        for start in range(0, len(atlas_order), 24):
            page = atlas_order[start:start + 24]
            fig, axs = plt.subplots(6, 4, figsize=(9.5, 11), layout='constrained')
            fig.suptitle(f'Supported orthologs {start + 1}–{start + len(page)}')
            for ax in axs.flat:
                ax.set_visible(False)
            for gene, ax in zip(page, axs.flat):
                ax.set_visible(True)
                gene_card(ax, gene, table.index.get_loc(gene), x, fit, centers, raw, counts)
            pdf.savefig(fig)
            plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results-root', type=Path,
                        default=os.environ.get('PSEUDOSPACE_RESULTS_ROOT', REPO / 'results'))
    args = parser.parse_args()
    root = args.results_root / 'paper_figures/fig5'
    out = root / 'curve_audit'
    out.mkdir(parents=True, exist_ok=True)
    metadata = pd.read_csv(root / 'screen/genes.csv').set_index('gene')
    genes = metadata.index[metadata.screened]
    x = np.arange(WINDOW[0], WINDOW[1] + .001, .005).round(3)
    arrays = (*curve_arrays(root / 'screen', genes, x),
              *curve_arrays(root / 'screen_genehalf_2026-11-18', genes, x))
    table = classify(*arrays, genes)
    for section, fits in zip(SECTIONS, arrays[0]):
        table[f'peak_{section}'] = np.where(np.isfinite(fits).all(-1), x[fits.argmax(-1)], np.nan)
    centers, raw, counts = fixed_width_counts(args.results_root, genes)
    full = (centers >= WINDOW[0]) & (centers <= WINDOW[1])
    low = full & (counts[4:] >= 20).all(0)
    if low.sum() < 3:
        raise ValueError('Fewer than three common position bins have adequate human low-mito coverage')
    raw_z = shape(raw[:, :, full])
    low_z = shape(raw[:, :, low])
    fixed = cross(raw_z[:2], raw_z[2:4])
    low_cross = cross(low_z[:2], low_z[4:])
    table = table.assign(
        rebin_cross_min_r=fixed.min(0), rebin_cross_max_r=fixed.max(0),
        lowmito_cross_min_r=low_cross.min(0), lowmito_cross_max_r=low_cross.max(0),
        lowmito_vs_full_min_r=(low_z[2:4] * low_z[4:]).mean(-1).min(0),
    )
    table.to_csv(out / 'genes.csv')
    np.savez_compressed(out / 'reaggregated_counts.npz', genes=genes.to_numpy(str), sections=SECTIONS,
                        centers=centers, profiles=raw[:4], bin_units=counts[:4],
                        lowmito_profiles=raw[4:], lowmito_bin_units=counts[4:])
    memberships = json.loads((root.parent / 'gene_sets.json').read_text())
    annotate(table, memberships).to_csv(out / 'annotated_programs.csv', index=False)
    observed_program_pairs(genes, raw_z, low_z).to_csv(out / 'observed_program_pairs.csv', index=False)
    supported = overview(out, table, x, arrays[0])
    program_profiles(out, table, x, arrays[0], centers, raw, counts, supported.index)
    print(f'{len(table)} screened, {len(supported)} reliable in all sections, '
          f'{table.conserved.sum()} conserved and {table.divergent.sum()} divergent.')
    print('Minimum unit count per fixed bin and section:', counts.min(1).tolist())
    print('Common low-mito checkpoints:', centers[low].round(3).tolist())
    for amplitude in (.35, .5, .75):
        check = classify(*arrays, genes, amplitude)
        mouse, human = check.mouse_both_grade, check.human_both_grade
        print(amplitude, (mouse & ~human).sum(), (human & ~mouse).sum(),
              check[['reliable_both', 'conserved', 'divergent']].sum().tolist())


if __name__ == '__main__':
    main()
