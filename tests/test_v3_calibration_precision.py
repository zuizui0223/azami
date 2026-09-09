import itertools

import numpy as np
import pytest

from analysis.v3.calibration_precision import (
    BET_FRACTIONS, _log_lower_capital, amendment, bounded_mean_sequence, precision_decision,
)
from analysis.v3.multicoordinate_calibration import wilson_interval


def test_budget_counts_metrics_as_well_as_both_grids_and_scenarios():
    spec = amendment()
    assert spec["confidence_sequences"] == 4 * 2 * 2
    assert spec["alpha_per_sequence"] * 16 == .05


def test_fixed_bet_factor_has_expectation_one_for_any_bounded_distribution():
    x = np.array([0., .2, .8, 1.])
    probabilities = np.array([.15, .2, .35, .3])
    mean = x @ probabilities
    for f in BET_FRACTIONS:
        plus = 1 - f + f * x / mean
        minus = 1 - f + f * (1 - x) / (1 - mean)
        assert min(plus.min(), minus.min()) >= 0
        assert plus @ probabilities == pytest.approx(1.)
        assert minus @ probabilities == pytest.approx(1.)


def test_exact_bernoulli_enumeration_bounds_anytime_exclusion_probability():
    alpha = .2
    for mean in (.1, .5, .9):
        crossing_probability = 0.
        for sequence in itertools.product([0., 1.], repeat=8):
            data = np.asarray(sequence)
            crosses = any(max(_log_lower_capital(data[:n], mean),
                              _log_lower_capital(1-data[:n], 1-mean)) >= np.log(2/alpha)
                          for n in range(1, 9))
            if crosses:
                crossing_probability += mean ** data.sum() * (1-mean) ** (8-data.sum())
        assert crossing_probability <= alpha + 1e-12


def test_duplicating_correlated_coordinates_does_not_create_more_outers():
    indicators = np.ones((100, 351))
    indicators[:10] = 0
    means = indicators.mean(axis=1)
    duplicated = np.tile(indicators, (1, 10)).mean(axis=1)
    original = bounded_mean_sequence(means, alpha=.003125)
    assert original == bounded_mean_sequence(duplicated, alpha=.003125)
    assert original["independent_outer_datasets"] == 100
    # The old pooled Wilson denominator falsely suggests much more precision.
    assert wilson_interval(int(indicators.sum()), int(indicators.size))[2] < .015
    assert original["half_width"] > .015


def test_boundary_samples_have_nonzero_uncertainty_and_reflection_symmetry():
    zero = bounded_mean_sequence(np.zeros(400), alpha=.003125)
    one = bounded_mean_sequence(np.ones(400), alpha=.003125)
    assert zero["low"] == 0 and zero["high"] > 0
    assert one["high"] == 1 and one["low"] < 1
    assert one["low"] == pytest.approx(1-zero["high"])
    varying = np.array([.93, .81, .95, .91, 1.])
    result = bounded_mean_sequence(varying, alpha=.003125)
    assert result["low"] <= varying.mean() <= result["high"]


@pytest.mark.parametrize("bad", [[], [[.4]], [.1, np.nan], [-.1], [1.1]])
def test_invalid_independent_samples_fail_closed(bad):
    with pytest.raises(ValueError):
        bounded_mean_sequence(bad, alpha=.003125)


def test_minimum_maximum_and_coverage_failure_never_rescued():
    assert not precision_decision([0]*25, [1.]*25)["precision_satisfied"]
    assert precision_decision([0]*400, [1.]*400)["precision_satisfied"]
    assert not precision_decision([0]*400, [.8]*400)["admission_satisfied"]
    with pytest.raises(ValueError, match="maximum"):
        precision_decision([0]*401, [1.]*401)
    with pytest.raises(ValueError, match="Aligned"):
        precision_decision([.5], [.95])
