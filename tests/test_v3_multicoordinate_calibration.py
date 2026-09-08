import numpy as np

from analysis.v3.dependence_resampling import source_partition
from analysis.v3.module_ecology import fit_matched_module
from analysis.v3.multicoordinate_calibration import (
    MODULES,
    definition,
    generate,
    precision_decision,
    wilson_interval,
)


def test_final_calibration_contract_matches_frozen_9_by_13_by_1_8_4_geometry():
    contract, env = definition()
    assert contract["bootstrap_replicates"] == 999
    assert contract["predictors"] == 9 == len(env["variables"])
    assert contract["nuisance_columns"] == 13
    assert contract["module_dimensions"] == {"orientation": 1, "visible_colour": 8, "gross_shape": 4}
    assert contract["family_slots"] == 36
    assert not contract["ecological_fitting_authorized"]


def test_generated_final_geometry_is_deterministic_and_exercises_frozen_nuisance():
    a = generate("scale_difference", 202609081)
    b = generate("scale_difference", 202609081)
    assert a["predictors"].shape[1] == 9
    assert a["nuisance"].shape[1] == 13
    assert len(a["predictor_names"]) == 9
    assert len(a["nuisance_names"]) == 13
    assert np.array_equal(a["predictors"], b["predictors"])
    assert np.array_equal(a["nuisance"], b["nuisance"])
    assert np.array_equal(a["components"], b["components"])
    for module, d in (("orientation", 1), ("visible_colour", 8), ("gross_shape", 4)):
        assert a["responses"][module].shape == (len(a["taxa"]), d)
        assert a["truths"][module].shape == (3, d, 9)
        assert np.array_equal(a["responses"][module], b["responses"][module])
        assert np.array_equal(a["truths"][module], b["truths"][module])
        # Scale-difference truth is radiation-only for the direct difference.
        direct = a["truths"][module][2]
        assert np.all(direct[:, :4] == 0)
        assert np.any(direct[:, 4] != 0)
        assert np.all(direct[:, 5:] == 0)
    assert len(np.unique(a["components"])) < len(a["components"])


def test_all_three_multicoordinate_point_estimators_are_identifiable_on_generated_case():
    data = generate("crossed_spatial_null", 202609082)
    for module in MODULES:
        fit = fit_matched_module(
            data["responses"][module], data["predictors"], data["taxa"], data["nuisance"]
        )
        assert fit["coefficients"].shape == data["truths"][module].shape
        assert fit["observations"] == len(data["taxa"])
        assert len(fit["taxa"]) == 64
        assert np.isfinite(fit["coefficients"]).all()
    source = source_partition(
        np.arange(len(data["taxa"])), data["taxa"], data["components"],
        data["latitude"], data["longitude"], grid_degrees=2,
    )
    assert len(source.positions) == len(data["taxa"])


def test_wilson_precision_rule_cannot_stop_before_minimum_and_fails_bad_calibration():
    low, high, half = wilson_interval(5, 100)
    assert 0 <= low <= 0.05 <= high <= 1
    assert half > 0
    early = precision_decision(
        outer_with_false_family=0, outer_total=50,
        covered_coefficients=10000, coefficient_total=10000,
    )
    assert early["next_batch_required"]
    assert not early["precision_satisfied"]
    bad = precision_decision(
        outer_with_false_family=80, outer_total=400,
        covered_coefficients=8000, coefficient_total=10000,
    )
    assert bad["stop_at_max_without_admission"]
    assert not bad["admission_satisfied"]
