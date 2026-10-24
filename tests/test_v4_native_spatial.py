import numpy as np
import pandas as pd
import pytest
from analysis.v4.native_spatial import cells, balanced_indices, conditional_residuals


def test_grid_wraps_dateline_and_bounds_poles():
    ids = cells([0, 0, -90, 90], [-180, 180, 0, 0])
    assert ids[0] == ids[1]
    assert np.all((ids >= 0) & (ids < 72*40))
    with pytest.raises(ValueError):
        cells([91], [0])


def test_balanced_selection_is_identity_stable_and_retains_species_cells():
    frame = pd.DataFrame({"obs_id": [8, 2, 1, 5, 3],
                          "accepted_key": [1, 1, 1, 2, 2], "_cell": [1, 1, 2, 1, 1]})
    first = frame.loc[balanced_indices(frame)]
    shuffled = frame.iloc[::-1]
    second = shuffled.loc[balanced_indices(shuffled)]
    assert sorted(first.obs_id) == sorted(second.obs_id)
    assert not first.duplicated(["accepted_key", "_cell"]).any()
    assert set(first.accepted_key) == {1, 2}


def test_residuals_subtract_species_blup_not_whole_species_mean():
    sample = pd.DataFrame({"accepted_key": [1, 1, 2, 2], "trait": [1., 2., 3., 4.]})
    sample.attrs["endpoint"] = "trait"
    component = {"reference_prediction": 0., "variance_ratio": 1., "effects": {}}
    result = conditional_residuals(sample, component, [], {})
    assert result == pytest.approx([0., 1., 3.-7/3, 4.-7/3])
