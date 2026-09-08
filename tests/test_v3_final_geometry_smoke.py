import numpy as np

from analysis.v3.multicoordinate_calibration import generate
from analysis.v3.run_final_geometry_smoke import definition


def test_final_geometry_smoke_contract_is_bounded_and_non_authorizing():
    spec, prelim = definition()
    assert spec["outer_replicates_per_scenario"] == 1
    assert spec["bootstrap_replicates"] == 199
    assert spec["grid_degrees"] == [2, 5]
    assert spec["module_dimensions"] == {"orientation": 1, "visible_colour": 8, "gross_shape": 4}
    assert spec["predictors"] == 9
    assert spec["nuisance_columns"] == 13
    assert spec["family_slots"] == 36
    assert list(prelim["process_indices"]) == ["wetting_moisture", "radiation", "heat_drying", "mechanical"]
    assert spec["empirical_trait_environment_values_read"] == 0
    assert spec["ecological_models_executed"] == 0
    assert not spec["ecological_fitting_authorized"]


def test_smoke_generator_exercises_final_multicoordinate_shapes_without_empirical_values():
    data = generate("crossed_spatial_null", 2026090891)
    assert data["predictors"].shape[1] == 9
    assert data["nuisance"].shape[1] == 13
    assert data["responses"]["orientation"].shape[1] == 1
    assert data["responses"]["visible_colour"].shape[1] == 8
    assert data["responses"]["gross_shape"].shape[1] == 4
    assert np.isfinite(data["predictors"]).all()
    assert np.isfinite(data["nuisance"]).all()
    assert len(np.unique(data["components"])) < len(data["components"])
