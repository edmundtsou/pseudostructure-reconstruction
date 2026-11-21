"""Prepare current COAST Fig. 5 units and tissue-level URAT1/NPT4 expression flags.

This is a descriptive spatial view of the four sections in the hvm COAST fit. Coordinates,
labels, and counts all come from that fit, so no separate segmentation or older tubule matrix
is joined here.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = next(p for p in Path(__file__).resolve().parents if (p / 'pseudospace').is_dir())
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
if str(Path(__file__).resolve().parent) not in sys.path:
    sys.path.insert(0, str(Path(__file__).resolve().parent))

import fig4_screen as fs  # noqa: E402
import fig5_screen as f5  # noqa: E402
from fig3_common import RUNS  # noqa: E402

TARGETS = ('Slc17a3', 'Slc22a12', 'Slc17a1')
HIGH_FRACTIONS = (0.20, 0.25, 1 / 3)


def _checked_xy(qc_obs: pd.DataFrame, p1_obs: pd.DataFrame) -> pd.DataFrame:
    """Use QC centroids after checking their coordinates against the current pass-1 units."""
    cols = ['centroid_x', 'centroid_y']
    common = qc_obs.index.intersection(p1_obs.index)
    if len(common) == 0:
        raise ValueError('No shared unit names between current QC and pass-1 objects')
    a = qc_obs.loc[common, cols].to_numpy(float)
    b = p1_obs.loc[common, cols].to_numpy(float)
    finite = np.isfinite(a).all(axis=1) & np.isfinite(b).all(axis=1)
    if not finite.all() or np.any(np.abs(a[finite] - b[finite]) > 1.0):
        raise ValueError('QC/pass-1 unit centroids disagree by more than the 1-pixel tolerance')
    return qc_obs[cols]


def tissue_units(results_root: Path, targets=TARGETS) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return current hvm PT units and section background centroids.

    PT rows are indexed by COAST observation names and include native-image centroids, the
    Fig. 5 mouse-reference quantile, Fig. 4 PT segment code, non-mito library, mito fraction,
    nearest Podocyte distance, and raw target-gene counts. The second table contains
    the current pseudospace-reference unit centroids used as each section's tissue background.
    """
    import anndata as ad
    from scipy import sparse
    from scipy.spatial import cKDTree

    root = Path(results_root)
    run = root / RUNS['hvm']
    qc = ad.read_h5ad(run / '02_qc/units_qc.h5ad')
    ref = pd.read_parquet(run / '04_pseudospace/pseudospace.parquet')
    ref = ref.set_index('obs_name') if 'obs_name' in ref.columns else ref
    if not ref.index.is_unique or not qc.obs_names.is_unique:
        raise ValueError('Current hvm observation names must be unique')

    p1 = ad.read_h5ad(run / '03_integration/pass1_labelled.h5ad', backed='r')
    xy = _checked_xy(qc.obs, p1.obs)
    p1_obs = p1.obs

    sections = f5.SECTIONS
    pt = ref[ref['sample'].astype(str).isin(sections)
             & ref['segment_class'].astype(str).isin(fs.PT)].copy()
    if not pt.index.isin(qc.obs_names).all():
        raise ValueError('Current pseudospace PT units are missing from current QC observations')
    pt['sample'] = pt['sample'].astype(str)
    sample = pt['sample'].to_numpy()
    seg = pt['segment_class'].astype(str).map(fs.PT)
    if seg.isna().any():
        raise ValueError('A PT label has no Fig. 4 PT segment code')

    counts = sparse.csc_matrix(qc[pt.index.to_numpy()].layers['counts'])
    mt = qc.var['mt'].to_numpy(bool)
    total = np.asarray(counts.sum(axis=1)).ravel().astype(float)
    library = np.asarray(counts[:, ~mt].sum(axis=1)).ravel().astype(float)
    mito = np.divide(np.asarray(counts[:, mt].sum(axis=1)).ravel(), total,
                     out=np.zeros_like(total), where=total > 0)
    gene_idx = qc.var_names.get_indexer(targets)
    if (gene_idx < 0).any():
        missing = [targets[i] for i in np.flatnonzero(gene_idx < 0)]
        raise ValueError(f'Target genes absent from current hvm count matrix: {missing}')
    target_counts = counts[:, gene_idx].toarray()

    coordinate = pt['pseudospace'].to_numpy(float)
    if hasattr(fs, 'quantile_axis'):
        positions = fs.quantile_axis(coordinate, np.isin(sample, f5.MOUSE))
    else:
        # Earlier Fig. 5 code defines the same reference mapping inside its screen.
        positions = f5.mouse_quantile(coordinate, sample, np.ones(len(pt), dtype=bool))
    pt_xy = xy.loc[pt.index, ['centroid_x', 'centroid_y']].to_numpy(float)
    glom_distance = np.full(len(pt), np.nan)
    for section in sections:
        glomeruli = p1_obs[(p1_obs['sample'].astype(str) == section)
                           & (p1_obs['segment_class'].astype(str) == 'Podocyte')]
        ix = np.flatnonzero(sample == section)
        if len(glomeruli) and len(ix):
            tree = cKDTree(glomeruli[['centroid_x', 'centroid_y']].to_numpy(float))
            glom_distance[ix] = tree.query(pt_xy[ix])[0]

    units = pd.DataFrame({
        'sample': sample,
        'x': pt_xy[:, 0],
        'y': pt_xy[:, 1],
        'position': positions,
        'segment': seg.to_numpy(dtype=int),
        'library': library,
        'mito': mito,
        'glom_distance': glom_distance,
        **{gene: target_counts[:, i] for i, gene in enumerate(targets)},
    }, index=pt.index)
    units.index.name = pt.index.name

    background_index = ref.index.intersection(xy.index)
    background = pd.DataFrame({
        'sample': ref.loc[background_index, 'sample'].astype(str).to_numpy(),
        'x': xy.loc[background_index, 'centroid_x'].to_numpy(float),
        'y': xy.loc[background_index, 'centroid_y'].to_numpy(float),
    }, index=background_index)
    p1.file.close()
    return units, background


def depth_strata(library: pd.Series, mito: pd.Series | None = None) -> pd.Series:
    """Return ten depth strata, optionally split into mitochondrial-fraction tertiles."""
    if library.index.has_duplicates:
        raise ValueError('Library index must be unique')
    values = library.to_numpy(float)
    if not np.isfinite(values).all():
        raise ValueError('Library must be finite')
    groups = np.array_split(np.argsort(values, kind='stable'), 10)
    if mito is not None:
        if not mito.index.equals(library.index) or not np.isfinite(mito).all():
            raise ValueError('Mitochondrial fractions must be finite and aligned with library sizes')
        mito_values = mito.to_numpy(float)
        groups = [sub for group in groups for sub in
                  np.array_split(group[np.argsort(mito_values[group], kind='stable')], 3)]
    strata = np.empty(len(library), dtype=np.int8)
    for number, group in enumerate(groups):
        strata[group] = number
    return pd.Series(strata, index=library.index, name='depth_stratum')


def high_expression(units: pd.DataFrame, fraction: float = 0.25, *, genes=TARGETS[:2],
                    strata: pd.Series | None = None) -> pd.DataFrame:
    """Flag within-decile high expressors in one section.

    Units are stably sorted by non-mito library and split into ten near-equal groups. Within
    each group, the top `fraction` by raw counts per non-mito library is flagged separately
    for each transcript, up to the number with positive counts. Stable sorting breaks ties.
    """
    if not 0 < fraction <= 1:
        raise ValueError('fraction must be in (0, 1]')
    required = {'library', *genes}
    missing = required.difference(units.columns)
    if missing:
        raise KeyError(f'Missing required unit columns: {sorted(missing)}')
    if units.index.has_duplicates:
        raise ValueError('Unit index must be unique')
    if (units['library'] <= 0).any() or not np.isfinite(units['library']).all():
        raise ValueError('Non-mito library must be finite and positive')

    out = pd.DataFrame(False, index=units.index, columns=list(genes), dtype=bool)
    if strata is None:
        strata = depth_strata(units['library'])
    if not strata.index.equals(units.index) or strata.isna().any():
        raise ValueError('Strata must cover all units in the same order')
    groups = [np.flatnonzero(strata.to_numpy() == k) for k in strata.unique()]
    for gene in out.columns:
        rate = units[gene].to_numpy(float) / units['library'].to_numpy(float)
        flags = np.zeros(len(units), dtype=bool)
        for group in groups:
            if not len(group):
                continue
            order = group[np.argsort(-rate[group], kind='stable')]
            n_top = int(np.ceil(fraction * len(group)))
            chosen = order[:n_top]
            flags[chosen[rate[chosen] > 0]] = True
        out[gene] = flags
    return out


def overlap_table(units, genes=TARGETS[:2], *, match_mito=False):
    """Compare joint-high counts with their expected values under within-stratum independence."""
    rows, displayed = [], {}
    left, right = genes
    for section in f5.SECTIONS:
        block = units[units['sample'].eq(section)]
        strata = depth_strata(block.library, block.mito if match_mito else None)
        for fraction in HIGH_FRACTIONS:
            high = high_expression(block, fraction, genes=genes, strata=strata)
            high['both'] = high[left] & high[right]
            grouped = high.groupby(strata)
            marginal = grouped[[left, right]].mean()
            expected = (marginal.prod(axis=1) * grouped.size()).sum() / len(high)
            observed = high.both.mean()
            rows.append({'section': section, 'fraction': fraction, 'n': len(high),
                         'joint_high': int(high.both.sum()), 'observed': observed, 'expected': expected,
                         'observed_expected': observed / expected if expected > 0 else np.nan})
            if fraction == 0.25:
                displayed[section] = high
    return pd.DataFrame(rows), displayed


def raw_endpoint_contrasts(units, axes, genes):
    """Return gene-wise log₂ late-minus-early rates across the three Fig. 5 count contrasts."""
    genes = list(genes)
    contrasts = {}
    for section in f5.SECTIONS:
        block = units[units['sample'].eq(section)]
        early, late = block.segment.eq(0), block.segment.eq(2)
        early_core = early & block.position.le(axes['b12'] - fs.MARGIN)
        late_core = late & block.position.ge(axes['b23'] + fs.MARGIN)
        if section in f5.MOUSE:
            cutoff = block.loc[early_core, 'glom_distance'].quantile(0.90)
            restricted = early_core, late_core & block.glom_distance.le(cutoff)
        else:
            low_mito = block.mito.le(block.mito.quantile(1 / 3))
            restricted = early_core & low_mito, late_core & low_mito
        for name, masks in (('all_labels', (early, late)), ('cores', (early_core, late_core)),
                            ('qc_subset', restricted)):
            rates = []
            for mask in masks:
                library = block.loc[mask, 'library'].sum()
                if library <= 0:
                    raise ValueError(f'{section}: {name} has an endpoint without expression coverage')
                rates.append((block.loc[mask, genes].sum() + 0.5) / library)
            contrasts[f'{name}_{section}'] = np.log2(rates[1] / rates[0])
    return pd.DataFrame(contrasts)
