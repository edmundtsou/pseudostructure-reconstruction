"""Synthetic checks for the Fig. 5 screen (analysis/paper/fig5_screen.py)."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'analysis/paper'))
pytest.importorskip('pandas')


def test_mouse_quantile_maps_human_units_through_the_mice():
    from fig5_screen import mouse_quantile

    mean = np.array([0.1, 0.2, 0.3, 0.4, 0.15, 0.9])
    sample = np.array(['Ctrl1A2', 'Ctrl1A4', 'Ctrl1A2', 'Ctrl1A4', 'HUK1_COR1', 'HUK1_MED1'])
    x = mouse_quantile(mean, sample, np.ones(6, bool))
    assert np.allclose(x[:4], [0.125, 0.375, 0.625, 0.875])
    assert x[4] == pytest.approx(0.25) and x[5] == pytest.approx(0.875)


def test_end_contrast_reads_only_inside_the_section_range():
    from fig5_screen import GRID, end_contrast

    e = end_contrast(2.0 * GRID, 0.3, 0.7, (0.0, 1.0))
    s1, s3 = GRID[GRID <= 0.24].mean(), GRID[GRID >= 0.76].mean()
    assert e == pytest.approx(2.0 * (s3 - s1))
    assert np.isnan(end_contrast(2.0 * GRID, 0.3, 0.7, (0.0, 0.7)))      # no S3 core inside the range


def test_reversal_needs_both_species_and_is_blind_under_relabeling():
    from fig5_screen import SECTIONS, SPLITS, category_r

    # gene 0: mouse rises, human falls; gene 1: one human section disagrees; gene 2: the sections disagree within
    # each species, in the pattern of relabeling 1 (rows: Ctrl1A2, Ctrl1A4, COR1, MED1; columns: genes)
    E = dict(zip(SECTIONS, np.array([[1.0, 1.0, 1.0], [1.0, 1.0, -1.0], [-1.0, 1.0, 1.0], [-1.0, -1.0, -1.0]])))
    R = {m: v.copy() for m, v in E.items()}
    ok = np.ones(3, bool)
    calls = category_r(E, R, ok, *SPLITS['species'])
    assert calls.tolist() == [True, False, False]
    # gene 2 reverses between the relabeling-1 arms {COR1, Ctrl1A2} and {MED1, Ctrl1A4}
    assert category_r(E, R, ok, *SPLITS['relabeling 1']).tolist() == [False, False, True]
    R['Ctrl1A2'] = -R['Ctrl1A2']                                            # raw counts disagree with the curve
    assert not category_r(E, R, ok, *SPLITS['species'])[0]


def test_steeper_in_human_needs_every_pair():
    from fig5_screen import SPLITS, category_h

    E = {'Ctrl1A2': np.array([0.3, 0.3]), 'Ctrl1A4': np.array([0.4, 1.0]),
         'HUK1_COR1': np.array([1.5, 1.5]), 'HUK1_MED1': np.array([1.2, 1.2])}
    R = {m: v.copy() for m, v in E.items()}
    assert category_h(E, R, np.ones(2, bool), *SPLITS['species']).tolist() == [True, False]   # 1.2 - 1.0 < 0.5


def test_flattening_residual_removes_a_uniform_attenuation():
    from fig5_screen import flattening_residuals

    rng = np.random.default_rng(0)
    m = rng.normal(0, 1, 500)
    E = {'Ctrl1A2': m, 'Ctrl1A4': m + rng.normal(0, 0.01, 500), 'HUK1_COR1': 0.4 * m, 'HUK1_MED1': 0.4 * m}
    E['HUK1_COR1'][0] = E['HUK1_MED1'][0] = -0.4 * m[0]                    # one reversed gene
    d, pairs, lines = flattening_residuals(E, np.ones(500, bool), 'species')
    assert all(abs(s - 0.4) < 0.02 for _, s in lines)
    assert abs(d[0]) > 0.5 * abs(m[0]) and np.median(np.abs(d[1:])) < 0.02


def test_urate_high_expression_is_depth_stratified_and_requires_detection():
    import pandas as pd

    from fig5_urate_tissue import high_expression

    units = pd.DataFrame({
        'library': np.repeat(np.arange(1, 11) * 100, 4),
        'Slc17a3': np.tile([0, 2, 0, 0], 10),
        'Slc22a12': np.tile([1, 2, 3, 4], 10),
    })
    flags = high_expression(units, fraction=0.5)
    assert flags.index.equals(units.index)
    assert flags.Slc17a3.tolist() == [False, True, False, False] * 10
    assert flags.Slc22a12.tolist() == [False, False, True, True] * 10


def test_tissue_matching_splits_depth_strata_by_mitochondrial_fraction():
    import pandas as pd

    from fig5_urate_tissue import depth_strata, high_expression

    units = pd.DataFrame({
        'library': np.repeat(np.arange(1, 11) * 100, 6),
        'mito': np.tile([0.2, 0.1, 0.3, 0.1, 0.2, 0.3], 10),
        'a': np.tile([1, 2, 3, 4, 5, 6], 10),
        'b': np.tile([6, 5, 4, 3, 2, 1], 10),
    })
    strata = depth_strata(units.library, units.mito)
    assert strata.iloc[:6].tolist() == [1, 0, 2, 0, 1, 2]
    high = high_expression(units, fraction=0.5, genes=('a', 'b'), strata=strata)
    assert high.sum().tolist() == [30, 30]
    assert not high.all(axis=1).any()


def test_global_heatmap_keeps_gene_order_and_between_section_amplitude(tmp_path):
    import pandas as pd

    from fig5_global_shapes import SECTIONS, load_global_shapes

    shape_dir = tmp_path / 'shape_clusters'
    screen_dir = tmp_path / 'screen'
    shape_dir.mkdir()
    screen_dir.mkdir()
    pd.DataFrame({'gene': ['a', 'b'], 'changed': [True, True], 'cluster': [2, 1],
                  'rms_condition': [1.4, 1.1]}).to_csv(shape_dir / 'genes.csv', index=False)
    x = np.array([0.1, 0.5, 0.9])
    pattern = np.array([[-1., 0., 1.], [1., 0., -1.]])
    shapes, fits = {}, {}
    for i, section in enumerate(SECTIONS):
        amplitude = 1 if i < 2 else 0.25
        signal = 3 + pattern * amplitude
        shapes[f'z|{section}'] = pattern / pattern.std(axis=1, keepdims=True)
        fits[f'curve|{section}'] = np.log2(np.exp2(signal) - 0.5)
    np.savez(shape_dir / 'shapes.npz', x=x, genes=['a', 'b'], **shapes)
    np.savez(screen_dir / 'all_curves.npz', grid=x, genes=['a', 'b'], **fits)

    grid, genes, curves = load_global_shapes(tmp_path)
    np.testing.assert_array_equal(grid, x)
    assert genes.index.tolist() == ['b', 'a']
    np.testing.assert_allclose(curves[0], shapes[f'z|{SECTIONS[0]}'][::-1])
    np.testing.assert_allclose(curves[2], 0.25 * curves[0], atol=1e-12)


def test_full_curve_audit_separates_reversal_from_flat_signal():
    import pandas as pd

    from fig5_curve_audit import classify

    x = np.linspace(-1., 1., 15)
    fitted = np.array([
        [x, x, x],
        [x, x, x],
        [x, -x, .1 * x],
        [x, -x, .1 * x],
    ])
    result = classify(fitted, fitted, fitted, fitted, pd.Index(['same', 'reversed', 'weak']))
    assert result['conserved'].tolist() == [True, False, False]
    assert result['divergent'].tolist() == [False, True, False]
    assert result['reliable_both'].tolist() == [True, True, False]


def test_profile_center_separates_relocation_from_attenuation():
    from fig5_enrichment import position_centers

    x = np.linspace(0, 1, 21)
    profiles = np.vstack([x, 3 + .01 * x, 1 - x, np.ones_like(x)])
    np.testing.assert_allclose(position_centers(profiles, x), [2 / 3, 2 / 3, 1 / 3, np.nan], atol=.001)
