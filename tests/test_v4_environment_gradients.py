import numpy as np
import pandas as pd
import pytest

from analysis.v4.environment_gradients import exposure_support, gradient_points
from analysis.v4.environment_units import interval_per_increment


def test_between_support_does_not_count_abundant_species_repeatedly():
    rows = pd.DataFrame({"accepted_key": [1]*100+[2]*2,
                         "rain__between": [10]*100+[100]*2,
                         "rain__within": [0]*100+[-2,2]})
    support = exposure_support(rows, ["rain"])
    assert support["rain__between"]["median"] == 55
    assert support["rain__between"]["n"] == 2
    assert support["rain__within"]["n"] == 102


def test_change_interval_reverses_below_reference_and_zero_is_exact():
    curve = gradient_points(
        {"raw_units_per_environment_sd": 2, "raw_ci95": [1,3]},
        {"q05": 0, "median": 10, "q95": 20}, 10, count=3)
    np.testing.assert_allclose(curve["change"], [-2,0,2])
    np.testing.assert_allclose(curve["lower"], [-3,0,1])
    np.testing.assert_allclose(curve["upper"], [-1,0,3])
    with pytest.raises(ValueError):
        gradient_points({}, {}, 0)


def test_physical_increment_undoes_only_internal_predictor_standardization():
    effect = {"raw_units_per_environment_sd": 2, "raw_ci95": [1,3]}
    beta, interval = interval_per_increment(effect, sd=500, increment=100)
    assert beta == .4
    np.testing.assert_allclose(interval,[.2,.6])
