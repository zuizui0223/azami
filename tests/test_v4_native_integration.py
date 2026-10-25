import numpy as np
import pytest
from analysis.v4.native_integration import summaries, covariance_levels, association, analyse


def test_within_and_between_covariance_match_explicit_equal_species_calculation():
    values = np.array([[0., 1.], [2., 3.], [10., -1.], [12., -3.], [14., -5.]])
    groups = np.array([1, 1, 2, 2, 2])
    means, cov, _ = summaries(values, groups)
    within, between = covariance_levels(means, cov, [1, 1])
    expected = (np.cov(values[:2], rowvar=False, bias=True)+
                np.cov(values[2:], rowvar=False, bias=True))/2
    assert within == pytest.approx(expected)
    assert between == pytest.approx(np.cov(means, rowvar=False, bias=True))


def test_scalar_rv_retains_signed_correlation_separately():
    rv, signed = association(np.array([[1., -.6], [-.6, 1.]]), [0], [1])
    assert rv == pytest.approx(.36)
    assert signed == pytest.approx(-.6)


def test_hue_rv_invariant_to_rotating_circular_origin():
    rng = np.random.default_rng(58)
    values = rng.normal(size=(100, 3))
    angle = .9
    rotation = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    altered = values.copy()
    altered[:, :2] = values[:, :2]@rotation
    first = association(np.cov(values, rowvar=False), [0, 1], [2])[0]
    second = association(np.cov(altered, rowvar=False), [0, 1], [2])[0]
    assert first == pytest.approx(second)


def test_identical_raw_and_adjusted_have_exactly_zero_paired_contrasts():
    rng = np.random.default_rng(3)
    values = rng.normal(size=(150, 2))
    result = analyse(values, values.copy(), np.repeat(np.arange(30), 5),
                     ["orientation_image_vertical_angle", "corolla_lab_lightness"],
                     [("angle", "lightness")], rng, replicates=20, equal_n=True)
    for level in ("within", "between"):
        contrast = result["pairs"][0]["estimates"]["raw_minus_adjusted/"+level]
        assert contrast["estimate"] == 0.
        assert contrast["interval"]["lower"] == contrast["interval"]["upper"] == 0.
