"""Fig. 4 panels drawn from the screen (fig4_screen.py), sized for the figure layout (make_fig4.py).

 gam  how one curve is fitted, refitted from the raw counts: tubules, spline basis, NB fit, LR test vs decoys
 ex   bulk, clusters and pseudostructure on one plot for three example genes
 lines  every gene the clusters miss, by pattern (median curve against the cluster staircase; no longer in Fig. 4)
 sc   COAST on mouse PT single-nucleus RNA-seq (fig4_sc.py): nuclei and curve in UMAP by the authors' S1/S2/S3
      labels, label densities along the axis, and the Fig. 4b genes along it next to the tubule curves
 a  one gene three ways: bulk, S1/S2/S3 clusters, pseudostructure (the top category 1 gene)
 b  every gene that passes either category: curves next to the same genes painted with their cluster means, and bulk
 c  the top three category 1 genes (peak between segments)
 d  the top three category 2 genes (change within a segment)

Reads <results-root>/paper_figures/fig4/screen/ and writes PNGs to <results-root>/paper_figures/fig4/panels/.
Protocol: docs/paper/fig4-continuous-vs-discrete.md.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = next(p for p in Path(__file__).resolve().parents if (p / 'pseudospace').is_dir())
sys.path.insert(0, str(Path(__file__).resolve().parent))
from fig3_panels import COLOURS, INK, MUTED, SECTION_COLOURS, SECTION_NAMES, style  # noqa: E402

MICE = ['Ctrl1A2', 'Ctrl1A4']
SEGS = ['PTS1', 'PTS2', 'PTS3']
CLUSTER, BULK = '#555555', '#AAAAAA'
LABEL_LOG2 = 0.5
LINES = None                # straight lines through the clusters, set by load()
# Fig. 4b examples, chosen by hand after fig4_linear_gain.py from its top 10 (COAST fit: ranks 1, 5, 9): one per shape
# the straight lines miss, the cleanest in both mice. Peak at S2/S3; switch within S2; dip just past S2/S3.
OVERLAY = ['Pigr', 'Slc22a7', 'Cyp4a12a']
# panel boxes of make_fig4.py in SVG units (1 unit = 0.01 inch here)
SIZE = {'gam': (1160, 250), 'ex': (540, 760), 'sc': (580, 760), 'lines': (580, 760), 'a': (540, 292), 'b': (580, 967), 'c': (540, 282), 'd': (540, 282)}


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--results-root', default=os.environ.get('PSEUDOSPACE_RESULTS_ROOT', str(REPO / 'results')))
    p.add_argument('--out', default='paper_figures/fig4', help='relative to the results root')
    return p.parse_args(argv)


def load(root: Path):
    genes = pd.read_csv(root / 'screen/genes.csv').set_index('gene')
    z = np.load(root / 'screen/curves.npz')
    names = list(z['genes'])
    curves = {m: {part: dict(zip(names, z[f'{part}|{m}'])) for part in ('curve', 'clusters', 'bulk', 'bins', 'bin_x')}
              for m in MICE}
    # Post-hoc rule (after the screen, before the figure was final; recorded in the protocol doc): a call is shown
    # only if it holds on raw counts within single labels, by >= LABEL_LOG2 in the same direction in both mice.
    # Category 2 examples are ranked by that within-label change: the curve's change can be inflated where a
    # spline dives into near-zero counts (Dpys in late S3).
    chk = pd.read_csv(root / 'screen/label_check.csv').set_index('gene')
    a, b = MICE
    within = chk[[f'within_label_change_{a}', f'within_label_change_{b}']]
    same = (np.sign(within[f'within_label_change_{a}']) == np.sign(chk[f'curve_change_{a}'])) & \
           (np.sign(within[f'within_label_change_{b}']) == np.sign(chk[f'curve_change_{b}']))
    ok2 = same & (within.abs().min(1) >= LABEL_LOG2)
    peak = chk[[f'rise_in_S2_{a}', f'rise_in_S2_{b}', f'fall_in_S3_{a}', f'fall_in_S3_{b}']].min(1)
    genes['label_ok'] = False
    genes.loc[chk.index, 'label_ok'] = np.where(chk.category == 1, peak >= LABEL_LOG2, ok2)
    genes['label_score'] = np.nan
    genes.loc[chk.index, 'label_score'] = np.where(chk.category == 1, peak, within.abs().min(1))
    genes['cat1'] &= genes.label_ok
    genes['cat2'] &= genes.label_ok
    global LINES
    LINES = pd.read_csv(root / 'screen/linear_gain.csv').set_index('gene')    # fig4_linear_gain.py
    return genes, z['grid'], curves


def segment_band(ax, b12, b23, y, h):
    from matplotlib.patches import Rectangle

    for (lo, hi), seg in zip(((0, b12), (b12, b23), (b23, 1)), SEGS):
        ax.add_patch(Rectangle((lo, y), hi - lo, h, color=COLOURS[seg], lw=0, clip_on=False))


def gene_axes(ax, gene, row, grid, curves, show=('bins', 'bulk', 'clusters', 'curve'), ylim=None):
    """One gene in both mice: bin means, bulk, cluster steps and the curve."""
    b12, b23 = row.b12, row.b23
    edges = ((0, b12), (b12, b23), (b23, 1))
    for m in MICE:
        c = SECTION_COLOURS[m]
        if 'bins' in show:
            ax.scatter(curves[m]['bin_x'][gene], curves[m]['bins'][gene], s=5, color=c, lw=0, alpha=0.8, zorder=3)
        if 'bulk' in show:
            ax.plot([0, 1], [curves[m]['bulk'][gene]] * 2, color=BULK, lw=1.0, ls=(0, (2, 1.5)), zorder=1)
        if 'clusters' in show:
            for (lo, hi), v in zip(edges, curves[m]['clusters'][gene]):
                ax.plot([lo, hi], [v, v], color=CLUSTER, lw=1.3, zorder=2, solid_capstyle='butt')
        if 'lines' in show:
            lg = LINES.loc[gene]
            yl = lg[f'lines_b0_{m}'] + lg[f'lines_b1_{m}'] * grid + lg[f'lines_b2_{m}'] * np.maximum(grid - lg.c2, 0)
            ax.plot(grid, (yl + np.log(1e4)) / np.log(2), color=c, lw=0.9, ls=(0, (1.5, 1.2)), zorder=3)
        if 'curve' in show:
            # stop at the last bin: past it the spline tail and the mislabelled end units take over (protocol doc)
            k = grid <= np.nanmax(curves[m]['bin_x'][gene])
            ax.plot(grid[k], curves[m]['curve'][gene][k], color=c, lw=1.4, zorder=4)
    lo, hi = ylim or ax.get_ylim()
    ax.set_ylim(lo, hi)
    pad = (hi - lo) * 0.06
    segment_band(ax, b12, b23, lo - 2.2 * pad, pad)
    for b in (b12, b23):
        ax.axvline(b, color='#CCCCCC', lw=0.5, zorder=0)
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['', ''])
    ax.tick_params(axis='x', length=0)


def ylim_of(gene, curves, pad=0.5):
    """The data's range (bin means and cluster means); a spline that dives where counts run out is clipped."""
    v = np.concatenate([np.r_[curves[m]['bins'][gene], curves[m]['clusters'][gene]] for m in MICE])
    v = v[np.isfinite(v)]
    return v.min() - pad, v.max() + pad


def examples(genes):
    c1 = genes[genes.cat1].sort_values('cat1_score', ascending=False).index.tolist()
    c2 = genes[genes.cat2 & ~genes.cat1].sort_values('label_score', ascending=False).index.tolist()
    return c1, c2


# ---------------------------------------------------------------------------------------------------- a
def panel_a(plt, genes, grid, curves, out):
    c1, c2 = examples(genes)
    gene = (c1 or c2)[0]
    fig = plt.figure(figsize=(SIZE['a'][0] / 100, SIZE['a'][1] / 100))
    lim = ylim_of(gene, curves)
    titles = ('bulk', 'S1 / S2 / S3 clusters', 'pseudostructure')
    for i, show in enumerate((('bins', 'bulk'), ('bins', 'clusters'), ('bins', 'curve'))):
        ax = fig.add_axes([0.08 + i * 0.31, 0.16, 0.27, 0.70])
        gene_axes(ax, gene, genes.loc[gene], grid, curves, show, lim)
        ax.set_title(titles[i], fontsize=7.5, color=INK, pad=3)
        if i:
            ax.set_yticklabels([])
        else:
            ax.set_ylabel('log2 counts per 10k')
            ax.text(0.02, 0.97, gene, transform=ax.transAxes, fontsize=8, style='italic', color=INK, va='top')
        if i == 2:
            for j, m in enumerate(MICE):
                ax.text(0.03, 0.96 - j * 0.08, m, transform=ax.transAxes, fontsize=6.5, color=SECTION_COLOURS[m],
                        va='top')
    fig.savefig(out / 'f4_a.png', facecolor='white')
    plt.close(fig)
    return gene


# ---------------------------------------------------------------------------------------------------- b
def panel_b(plt, genes, grid, curves, out, mark):
    from matplotlib.colors import LinearSegmentedColormap

    c1 = genes[genes.cat1].copy()
    c2 = genes[genes.cat2 & ~genes.cat1].copy()
    c1['key'] = c1['peak_at_Ctrl1A2'] * 10 + c1[['peak_x_Ctrl1A2', 'peak_x_Ctrl1A4']].mean(1)
    sign = np.sign(c2['change_Ctrl1A2'])
    c2['key'] = c2['change_segment_Ctrl1A2'] * 10 + (sign < 0) * 5 + c2['cat2_score'].rank(ascending=False) / 1e4
    rows = c1.sort_values('key').index.tolist() + c2.sort_values('key').index.tolist()
    n1 = len(c1)
    curve = np.stack([np.mean([curves[m]['curve'][g] for m in MICE], 0) for g in rows])
    steps = np.zeros_like(curve)
    bulk = np.zeros(len(rows))
    for i, g in enumerate(rows):
        b12, b23 = genes.loc[g, ['b12', 'b23']]
        seg = np.digitize(grid, [b12, b23])
        steps[i] = np.mean([curves[m]['clusters'][g] for m in MICE], 0)[seg]
        bulk[i] = np.mean([curves[m]['bulk'][g] for m in MICE])
    lo, hi = curve.min(1, keepdims=True), curve.max(1, keepdims=True)
    scale = lambda v: (v - lo) / (hi - lo)  # noqa: E731
    cmap = LinearSegmentedColormap.from_list('seq', ['#FFFFFF', '#C9D3EE', '#5B7DC8', '#1D2F6F'])
    fig = plt.figure(figsize=(SIZE['b'][0] / 100, SIZE['b'][1] / 100))
    top, bottom = 0.955, 0.04
    widths = [(0.12, 0.38), (0.54, 0.38), (0.95, 0.035)]
    titles = ('pseudostructure', 'S1 / S2 / S3 clusters', 'bulk')
    for (x0, w), mat, title in zip(widths, (scale(curve), scale(steps), scale(bulk[:, None])), titles):
        ax = fig.add_axes([x0, bottom, w, top - bottom])
        ax.set_title(title, fontsize=7.5, color=INK, pad=3)
        ax.imshow(np.clip(mat, 0, 1), aspect='auto', interpolation='nearest', cmap=cmap, vmin=0, vmax=1,
                  extent=(0, 1, len(rows), 0))
        ax.axhline(n1, color=INK, lw=0.8)
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        if w > 0.1:
            band = fig.add_axes([x0, bottom - 0.014, w, 0.009])
            b12, b23 = genes.loc[rows, 'b12'].mean(), genes.loc[rows, 'b23'].mean()
            segment_band(band, b12, b23, 0, 1)
            band.set_xlim(0, 1)
            band.set_ylim(0, 1)
            band.axis('off')
    lab = fig.add_axes([0.0, bottom, 0.115, top - bottom])
    lab.set_xlim(0, 1)
    lab.set_ylim(len(rows), 0)
    lab.axis('off')
    for g in mark:
        if g in rows:
            lab.text(0.95, rows.index(g) + 0.5, g, ha='right', va='center', fontsize=6, style='italic', color=INK)
    fig.savefig(out / 'f4_b.png', facecolor='white')
    plt.close(fig)
    return {'n_rows': len(rows), 'n_cat1': n1}


# ---------------------------------------------------------------------------------------------------- b (lines)
GROUPS = [('peak at S2/S3', lambda g: g.cat1),
          ('rises within S2', lambda g: ~g.cat1 & g.cat2 & (g.change_segment_Ctrl1A2 == 1) & (g.change_Ctrl1A2 > 0)),
          ('falls within S2', lambda g: ~g.cat1 & g.cat2 & (g.change_segment_Ctrl1A2 == 1) & (g.change_Ctrl1A2 < 0)),
          ('falls within S3', lambda g: ~g.cat1 & g.cat2 & (g.change_segment_Ctrl1A2 == 2) & (g.change_Ctrl1A2 < 0))]
MISSED = '#F4B6B0'


def centred(gene, genes, grid, curves):
    """Curve and cluster staircase of one gene, both mice averaged, as log2 fold from its bulk PT mean."""
    b12, b23 = genes.loc[gene, ['b12', 'b23']]
    seg = np.digitize(grid, [b12, b23])
    bulk = np.mean([curves[m]['bulk'][gene] for m in MICE])
    curve = np.mean([curves[m]['curve'][gene] for m in MICE], 0) - bulk
    steps = np.mean([curves[m]['clusters'][gene] for m in MICE], 0)[seg] - bulk
    # the spline can dive where counts run out; keep it inside the range of this gene's own bin means
    bins = np.concatenate([curves[m]['bins'][gene] for m in MICE]) - bulk
    bins = bins[np.isfinite(bins)]
    return np.clip(curve, bins.min() - 0.25, bins.max() + 0.25), steps


def panel_b_lines(plt, genes, grid, curves, out):
    fig = plt.figure(figsize=(SIZE['lines'][0] / 100, SIZE['lines'][1] / 100))
    rows = [(name, genes.index[rule(genes)].tolist()) for name, rule in GROUPS]
    rows = [(n, g) for n, g in rows if g]
    h, gap, top = 0.78 / len(rows), 0.20 / len(rows), 0.96
    counts = {}
    for i, (name, members) in enumerate(rows):
        ax = fig.add_axes([0.14, top - (i + 1) * h - i * gap + gap * 0.35, 0.82, h])
        pairs = [centred(g, genes, grid, curves) for g in members]
        cur, stp = np.stack([c for c, _ in pairs]), np.stack([s for _, s in pairs])
        for c in cur:
            ax.plot(grid, c, color=SECTION_COLOURS['Ctrl1A4'], lw=0.5, alpha=0.35, zorder=1)
        mc, ms = np.median(cur, 0), np.median(stp, 0)
        ax.fill_between(grid, mc, ms, color=MISSED, lw=0, zorder=2)
        ax.plot(grid, ms, color=CLUSTER, lw=2.0, zorder=3, solid_capstyle='butt')
        ax.plot(grid, mc, color=SECTION_COLOURS['Ctrl1A2'], lw=2.2, zorder=4)
        ax.axhline(0, color=BULK, lw=0.8, ls=(0, (2, 1.5)), zorder=0)
        both = np.r_[mc, ms]
        b12, b23 = genes.loc[members, 'b12'].mean(), genes.loc[members, 'b23'].mean()
        gene_band(ax, b12, b23, min(-2, np.floor(both.min() - 1)), max(2, np.ceil(both.max() + 1)))
        ax.set_title(f'{name}  ·  {len(members)} genes', fontsize=7.5, color=INK, pad=3, loc='left')
        ax.set_ylabel('log2 fold from PT mean')
        counts[name] = len(members)
    fig.savefig(out / 'f4_b.png', facecolor='white')
    plt.close(fig)
    return counts


def gene_band(ax, b12, b23, lo, hi):
    ax.set_ylim(lo, hi)
    pad = (hi - lo) * 0.05
    segment_band(ax, b12, b23, lo - 2.4 * pad, pad)
    for b in (b12, b23):
        ax.axvline(b, color='#CCCCCC', lw=0.5, zorder=0)
    ax.set_xlim(0, 1)
    ax.set_xticks([])


# ---------------------------------------------------------------------------------------------------- c, d
def panel_examples(plt, key, gene_list, genes, grid, curves, out):
    fig = plt.figure(figsize=(SIZE[key][0] / 100, SIZE[key][1] / 100))
    for i, g in enumerate(gene_list[:3]):
        ax = fig.add_axes([0.08 + i * 0.315, 0.13, 0.25, 0.74])
        gene_axes(ax, g, genes.loc[g], grid, curves, ('bins', 'clusters', 'curve'), ylim_of(g, curves))
        ax.set_title(g, fontsize=8, style='italic', color=INK, pad=3)
        if i == 0:
            ax.set_ylabel('log2 counts per 10k')
    fig.savefig(out / f'f4_{key}.png', facecolor='white')
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------- GAM anatomy
def refit(results_root, gene, mouse, genes):
    """Refit one gene exactly as fig4_screen.py does (COAST's penalised NB GAM on both mice, shared knots, smoothing
    and dispersion); returns the pieces the figure shows for one mouse. Uses COAST's internal basis and fitter
    (coast.tl._cr_basis, coast.tl._nb_gam) because the figure draws the coefficients, which fit_curves does not
    return."""
    import anndata as ad
    from scipy import sparse
    from coast.tl import _cr_basis, _nb_gam
    from fig4_screen import GRID, N_BINS, N_SPLITS, PT, log2_cp10k, pattern_axis
    from fig3_common import RUNS

    root = Path(results_root)
    qc = ad.read_h5ad(root / RUNS['aki'] / '02_qc/units_qc.h5ad', backed='r')
    ref = pd.read_parquet(root / RUNS['aki'] / '04_pseudospace/pseudospace.parquet')
    ref = ref.set_index('obs_name') if 'obs_name' in ref.columns else ref
    pt = ref[ref['sample'].astype(str).isin(MICE) & ref.segment_class.astype(str).isin(PT)]
    units = pt.index.to_numpy()
    rows = qc.obs_names.get_indexer(units)
    counts = sparse.csr_matrix(qc.layers['counts'][np.sort(rows)])[np.argsort(np.argsort(rows))]
    mt = qc.var['mt'].to_numpy(bool)
    lib = np.asarray(counts[:, ~mt].sum(1)).ravel()
    y_all = np.asarray(counts[:, qc.var_names.get_loc(gene)].todense()).ravel()
    axes = {(k, h): pd.read_parquet(root / 'paper_figures/fig3/robustness/aki/reruns' / f'gene_half_s{k}_{h}.parquet')['t']
            .reindex(units).to_numpy(float) for k in range(N_SPLITS) for h in 'AB'}
    x, ok = pattern_axis(axes, genes.loc[gene, 'pattern'])
    sample = pt['sample'].astype(str).to_numpy()
    xs, ys, ls, ss = x[ok], y_all[ok], lib[ok], sample[ok]
    knots = np.quantile(xs, np.linspace(0, 1, 10))          # coast.tl.fit_curves defaults: k = 10 at quantiles
    basis, S = _cr_basis(xs, knots)
    blocks = [(np.flatnonzero(ss == m), basis[ss == m], np.log(ls[ss == m])) for m in sorted(set(ss))]
    beta = _nb_gam(ys, blocks, S)[sorted(set(ss)).index(mouse)]
    on_grid = _cr_basis(GRID, knots)[0]
    k = ss == mouse
    label = pt.segment_class.astype(str).map(PT).to_numpy()[ok][k]
    xm, ym, lm = xs[k], ys[k], ls[k]
    inside = (GRID >= xm.min()) & (GRID <= xm.max())
    bins = np.array_split(np.argsort(xm), N_BINS)
    return {'x': xm, 'y': ym, 'lib': lm, 'label': label, 'grid': GRID, 'basis': on_grid, 'beta': beta,
            'knots': knots, 'curve': np.where(inside, log2_cp10k(on_grid @ beta), np.nan),
            'bin_x': np.array([xm[b].mean() for b in bins]),
            'bins': np.array([np.log2(1e4 * ym[b].sum() / lm[b].sum()) for b in bins])}


def log2_scale(b):
    """A knot coefficient (log rate per unit of library) in log2 counts per 10k."""
    return (b + np.log(1e4)) / np.log(2)


def panel_gam(plt, args, root, genes, grid, curves, gene, out, mouse='Ctrl1A2'):
    """How one curve is fitted, on real data: tubules, spline basis, NB fit, test against the clusters."""
    r = refit(args.results_root, gene, mouse, genes)
    stored = curves[mouse]['curve'][gene]
    drift = float(np.nanmax(np.abs(stored - r['curve'])))
    assert drift < 0.05, f'refit of {gene} differs from the screen by {drift:.3f} log2'
    row = genes.loc[gene]
    seg_col = [COLOURS[s] for s in SEGS]
    fig = plt.figure(figsize=(SIZE['gam'][0] / 100, SIZE['gam'][1] / 100))
    w, x0s = 0.19, (0.055, 0.30, 0.545, 0.80)
    axs = [fig.add_axes([x0, 0.20, w, 0.66]) for x0 in x0s]
    titles = ('tubules', 'spline basis', 'penalised NB GAM fit', 'curve vs clusters')
    # 1 tubules: counts per 10k of each PT tubule, coloured by its S1/S2/S3 label
    ax = axs[0]
    cp = np.log2(1 + 1e4 * r['y'] / r['lib'])
    jitter = np.random.default_rng(0).uniform(-0.08, 0.08, len(cp)) * (r['y'] == 0)
    ax.scatter(r['x'], cp + jitter, s=1.2, lw=0, c=[seg_col[l] for l in r['label']], alpha=0.6, rasterized=True)
    ax.set_ylabel('log2(1 + counts per 10k)')
    ax.text(0.03, 0.97, gene, transform=ax.transAxes, fontsize=8, style='italic', color=INK, va='top')
    # 2 basis: the cubic regression spline's basis functions, each scaled by its coefficient (the curve's log rate
    # at that knot), on the log2 scale; they sum to the fit
    ax = axs[1]
    cmap = plt.get_cmap('viridis')
    nb = r['basis'].shape[1]
    for j in range(nb):
        ax.plot(r['grid'], log2_scale(r['beta'][j]) * r['basis'][:, j], color=cmap(j / nb), lw=1.1)
    ax.axhline(0, color=BULK, lw=0.6, ls=(0, (2, 1.5)))
    for kn in r['knots']:
        ax.axvline(kn, color='#DDDDDD', lw=0.5, zorder=0)
    ax.set_ylabel('βₖ Bₖ(t)  (log2 per 10k)')
    # 3 fit: the sum, with bin means of the same tubules
    ax = axs[2]
    ax.scatter(r['bin_x'], r['bins'], s=6, color=SECTION_COLOURS[mouse], lw=0, zorder=3)
    ax.plot(r['grid'], r['curve'], color=SECTION_COLOURS[mouse], lw=1.6, zorder=4)
    ax.set_ylabel('log2 counts per 10k')
    ax.text(0.03, 0.97, 'log μ = log lib + Σ βₖ Bₖ(t)\npenalty λ ∫ f″(t)² dt', transform=ax.transAxes,
            fontsize=6.5, color=INK, va='top', linespacing=1.3)
    # 4 test: likelihood ratio of clusters + spline over clusters, every gene against its label-shuffled decoy
    ax = axs[3]
    lr, dlr = genes[f'lr_{mouse}'].dropna(), genes[f'decoy_lr_{mouse}'].dropna()
    edges = np.geomspace(0.1, max(lr.max(), dlr.max()) * 1.1, 45)
    ax.hist(dlr.clip(lower=0.1), edges, color='#BBBBBB', lw=0, label='label-shuffled decoy')
    ax.hist(lr.clip(lower=0.1), edges, histtype='step', color=SECTION_COLOURS[mouse], lw=1.0, label='genes')
    ax.set_xscale('log')
    import json
    thr = json.loads((root / 'screen/summary.json').read_text())['thresholds'][mouse]
    ax.axvline(thr, color=INK, lw=0.8, ls=(0, (3, 2)))
    ax.text(thr, 0.97, ' 5% FDR', transform=ax.get_xaxis_transform(), fontsize=6, color=INK, va='top')
    g_lr = row[f'lr_{mouse}']
    ax.annotate(gene, (g_lr, 0), xytext=(g_lr, ax.get_ylim()[1] * 0.35), fontsize=6.5, style='italic',
                ha='center', color=INK, arrowprops=dict(arrowstyle='-|>', lw=0.6, color=INK))
    ax.set_xlabel('likelihood ratio, curve over clusters', labelpad=1)
    ax.set_ylabel('genes')
    ax.legend(loc='upper right', bbox_to_anchor=(1, 0.80), fontsize=5.5, frameon=False, handlelength=1.2)
    for i, ax in enumerate(axs):
        ax.set_title(titles[i], fontsize=7.5, color=INK, pad=3)
        if i < 3:
            lo, hi_ = ax.get_ylim()
            pad = (hi_ - lo) * 0.05
            segment_band(ax, row.b12, row.b23, lo - 2.4 * pad, pad)
            ax.set_ylim(lo, hi_)
            ax.set_xlim(0, 1)
            ax.set_xticks([0, 1])
            ax.set_xticklabels(['', ''])
            ax.tick_params(axis='x', length=0)
            ax.set_xlabel('pseudostructure', labelpad=11)
    fig.savefig(out / 'f4_gam.png', facecolor='white')
    plt.close(fig)
    return {'gene': gene, 'mouse': mouse, 'refit_max_drift_log2': drift, 'lr': float(row[f'lr_{mouse}']),
            'threshold': thr}


def panel_overlay(plt, gene_list, genes, grid, curves, out):
    """Each example gene, both mice averaged: one bar for bulk, one bar per cluster, the curve on top."""
    from matplotlib.lines import Line2D

    fig = plt.figure(figsize=(SIZE['ex'][0] / 100, SIZE['ex'][1] / 100))
    n = len(gene_list)
    h = 0.80 / n
    curve_c = SECTION_COLOURS['Ctrl1A2']
    for i, g in enumerate(gene_list):
        ax = fig.add_axes([0.13, 0.95 - (i + 1) * h - i * (0.17 / n) + 0.02, 0.84, h])
        row = genes.loc[g]
        mean = lambda part: np.mean([curves[m][part][g] for m in MICE], 0)
        bulk, clusters, curve = mean('bulk'), mean('clusters'), mean('curve')
        ax.scatter(mean('bin_x'), mean('bins'), s=6, color=MUTED, lw=0, alpha=0.7, zorder=3)
        ax.plot([0, 1], [bulk, bulk], color=BULK, lw=6, alpha=0.45, solid_capstyle='butt', zorder=1)
        for (lo, hi), v in zip(((0, row.b12), (row.b12, row.b23), (row.b23, 1)), clusters):
            ax.plot([lo, hi], [v, v], color=CLUSTER, lw=3.5, solid_capstyle='butt', zorder=2)
        # stop at the last bin: past it the spline tail and the mislabelled end units take over (protocol doc)
        k = grid <= min(np.nanmax(curves[m]['bin_x'][g]) for m in MICE)
        ax.plot(grid[k], curve[k], color=curve_c, lw=2, zorder=4)
        v = np.r_[mean('bins'), clusters, curve[k]]
        lo, hi = np.nanmin(v) - 0.4, np.nanmax(v) + 0.4
        ax.set_ylim(lo, hi)
        pad = (hi - lo) * 0.06
        segment_band(ax, row.b12, row.b23, lo - 2.2 * pad, pad)
        for x in (row.b12, row.b23):
            ax.axvline(x, color='#CCCCCC', lw=0.5, zorder=0)
        ax.set_xlim(0, 1)
        ax.set_xticks([])
        ax.set_title(g, fontsize=8, style='italic', color=INK, pad=3, loc='left')
        ax.set_ylabel('log2 counts per 10k')
        if i == 1:
            keys = [Line2D([], [], color=BULK, lw=6, alpha=0.45), Line2D([], [], color=CLUSTER, lw=3.5),
                    Line2D([], [], color=curve_c, lw=2)]
            ax.legend(keys, ['bulk', 'S1 / S2 / S3 clusters', 'pseudostructure'], loc='lower right', fontsize=6.5,
                      frameon=False, handlelength=2.2, labelspacing=0.4)
    fig.savefig(out / 'f4_ex.png', facecolor='white')
    plt.close(fig)


def panel_sc(plt, args, genes, curves, out):
    """Cells against tubules (fig4_sc.py): COAST run the same way on bin2cell cells as on tubules. Per species, the
    cells in UMAP with the curve, coloured by the segment of the tubule each nucleus sits in; where each segment falls
    along the cell axis and along the tubule axis. The per-section rho of both axes is in summary.json, not drawn."""
    import json
    from scipy.stats import gaussian_kde
    from fig3_common import GROUP_ORDER, load_positions

    cell_dir = Path(args.results_root) / args.out / 'cells'
    fig = plt.figure(figsize=(SIZE['sc'][0] / 100, SIZE['sc'][1] / 100))
    summ = {}
    for species, name, run, y0 in [('mouse', 'Mouse control', 'aki', 0.53), ('human', 'Human', 'hvm', 0.04)]:
        d = cell_dir / species
        cells = pd.read_parquet(d / 'cells.parquet')
        nodes = np.load(d / 'nodes.npy')
        summ[species] = json.loads((d / 'summary.json').read_text())
        tub = load_positions(args.results_root, run)
        tub = tub[tub['sample'].isin(cells['sample'].unique())]
        groups = [g for g in GROUP_ORDER if (cells.tubule_group == g).sum() >= 50 and (tub.group == g).sum() >= 50]
        # UMAP coloured by the cells' own COAST cluster labels (Leiden + marker panels, as for tubules), each label
        # written at the median of its cells; the densities on the right are by the tubule each nucleus sits in
        ax = fig.add_axes([0.0, y0 + 0.04, 0.36, 0.36])
        grp = cells.group.fillna('').to_numpy()
        order = np.random.default_rng(0).permutation(len(cells))
        ax.scatter(cells.umap1.to_numpy()[order], cells.umap2.to_numpy()[order], s=0.15, lw=0, rasterized=True,
                   c=[COLOURS.get(g, '#E3E3E3') for g in grp[order]])
        ax.plot(*nodes.T, '-o', color='black', ms=1.4, lw=0.7, mfc='white', mew=0.4, zorder=3)
        for lab, part in cells.groupby('label'):
            ax.text(part.umap1.median(), part.umap2.median(), lab, fontsize=6, weight='bold', ha='center',
                    va='center', color=INK, zorder=4,
                    bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none', alpha=0.75))
        ax.set_axis_off()
        ax.set_title(f"{name}: {len(cells):,} cells, by cell cluster", fontsize=7, color=INK, pad=1)
        grid = np.linspace(0, 1, 300)
        for j, (title, t, g) in enumerate([('cells, by tubule', cells.t, cells.tubule_group), ('tubules', tub.t, tub.group)]):
            ax = fig.add_axes([0.42 + j * 0.30, y0 + 0.03, 0.26, 0.38])
            for i, grp in enumerate(groups):
                base = len(groups) - 1 - i
                yv = gaussian_kde(t[g == grp].to_numpy(float), bw_method=0.15)(grid)
                ax.fill_between(grid, base, base + 0.9 * yv / yv.max(), color=COLOURS[grp], alpha=0.5, lw=0)
                ax.plot(grid, base + 0.9 * yv / yv.max(), color=COLOURS[grp], lw=0.7)
                if j == 0:
                    ax.text(-0.04, base + 0.15, grp, fontsize=5.5, color=COLOURS[grp], ha='right', weight='bold')
            ax.set_xlim(0, 1)
            ax.set_ylim(0, len(groups))
            ax.set_yticks([])
            ax.spines['left'].set_visible(False)
            ax.set_xticks([0, 1])
            ax.set_xlabel('pseudostructure', labelpad=0)
            ax.set_title(title, fontsize=7, color=INK, pad=1)
    fig.savefig(out / 'f4_sc.png', facecolor='white')
    plt.close(fig)
    return summ


def main(argv=None):
    args = parse_args(argv)
    root = Path(args.results_root) / args.out
    out = root / 'panels'
    out.mkdir(parents=True, exist_ok=True)
    plt = style()
    genes, grid, curves = load(root)
    c1, c2 = examples(genes)
    a = panel_a(plt, genes, grid, curves, out)
    # presentation choice made after the screen: panel a shows the top gene, so c and d skip it (no duplicate)
    c1, c2 = [g for g in c1 if g != a][:3], [g for g in c2 if g != a][:3]
    b = panel_b_lines(plt, genes, grid, curves, out)
    gam = panel_gam(plt, args, root, genes, grid, curves, a, out)
    top10 = LINES.index[:10]
    assert all(g in top10 and g in curves[MICE[0]]['curve'] for g in OVERLAY), 'Fig. 4b examples left the top 10'
    ex = OVERLAY
    panel_overlay(plt, ex, genes, grid, curves, out)
    print('gam:', gam, '| overlay:', ex, '| sc rho:', {k: (v['rho_tubules'], v['rho_cells_tubule_label']) for k, v in panel_sc(plt, args, genes, curves, out).items()})
    if c1:
        panel_examples(plt, 'c', c1, genes, grid, curves, out)
    if c2:
        panel_examples(plt, 'd', c2, genes, grid, curves, out)
    pd.Series({'panel_a': a, 'cat1_examples': ','.join(c1), 'cat2_examples': ','.join(c2), **b}).to_csv(
        out / 'panels.csv', header=False)
    print('a:', a, '| c:', c1, '| d:', c2, '| b rows:', b)


if __name__ == '__main__':
    main()
