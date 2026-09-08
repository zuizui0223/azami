import json

import numpy as np

from analysis.v3.environment_model import test_family as _test_family
from analysis.v3.module_ecology import SCALES
from analysis.v3.simulate_full_family_tail_pilot import (
    SPEC,
    _slot_rows,
    _specs,
    generated_modules,
)


def test_tail_pilot_contract_resolves_first_holm_threshold_without_authorizing_ecology():
    spec, prelim = _specs()
    assert spec["family_slots"] == 36 == len(_test_family())
    assert spec["bootstrap_replicates"] == 999
    assert spec["minimum_plus_one_probability"] == 0.001
    assert spec["minimum_plus_one_probability"] <= 0.05 / 36
    assert spec["processes"] == list(prelim["process_indices"])
    assert not spec["ecological_fitting_authorized"]
    assert spec["empirical_trait_environment_values_read"] == 0


def test_generated_three_module_geometry_and_truth_are_deterministic():
    a = generated_modules("scale_difference", 123456)
    b = generated_modules("scale_difference", 123456)
    responses_a, truths_a, x_a, g_a, nuisance_a, components_a, lat_a, lon_a = a
    responses_b, truths_b, x_b, g_b, nuisance_b, components_b, lat_b, lon_b = b
    assert tuple(responses_a) == ("orientation", "visible_colour", "gross_shape")
    for module in responses_a:
        assert responses_a[module].shape == (len(x_a), 1)
        assert truths_a[module].shape == (3, 1, 9)
        assert np.array_equal(responses_a[module], responses_b[module])
        assert np.array_equal(truths_a[module], truths_b[module])
    assert np.array_equal(x_a, x_b)
    assert np.array_equal(g_a, g_b)
    assert np.array_equal(nuisance_a, nuisance_b)
    assert np.array_equal(components_a, components_b)
    assert np.array_equal(lat_a, lat_b)
    assert np.array_equal(lon_a, lon_b)
    expected = 0.6 * truths_a["orientation"] - 0.4 * truths_a["visible_colour"]
    assert np.allclose(truths_a["gross_shape"], expected)
    assert not np.allclose(
        responses_a["gross_shape"],
        0.6 * responses_a["orientation"] - 0.4 * responses_a["visible_colour"],
    )


def test_full_36_slot_assembly_reserves_every_module_process_scale_key():
    spec = json.loads(SPEC.read_text())
    _, prelim = _specs()
    blocks = prelim["process_indices"]
    _, truths, *_ = generated_modules("iid_null", 9988)
    summaries = {}
    for module in spec["modules"]:
        candidate_tests = []
        for scale in SCALES:
            for process, indices in blocks.items():
                candidate_tests.append({
                    "scale": scale,
                    "process": process,
                    "coefficients_tested": len(indices),
                    "status": "candidate_tail_estimated_not_calibrated",
                    "rank": 1,
                    "candidate_probability": 0.01,
                })
        summaries[module] = {"candidate_tests": candidate_tests}
    rows = _slot_rows(summaries, truths, blocks)
    assert len(rows) == 36
    keys = {(row["module"], row["process"], row["question"]) for row in rows}
    assert keys == set(_test_family())
    assert all(row["holm_probability"] is not None for row in rows)
