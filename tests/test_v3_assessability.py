import json

import numpy as np
import pandas as pd
import pytest

from analysis.v3.assessability import audit


def source():
    return pd.DataFrame({"obs_id": ["a", "b", "c", "d"], "accepted_key": ["A", "A", "A", "B"],
                         "native_range_status": ["native"] * 4, "pr_month": [1, 2, 3, np.nan]})


def states():
    return pd.DataFrame({"obs_id": ["a", "b", "c", "d"], "endpoint": ["angle"] * 4,
                         "state": ["usable", "usable", "qc_unusable", "download_failed"]})


def test_source_taxon_weights_are_not_recomputed_among_survivors():
    report = audit(source(), states(), {"orientation": ["angle"]})
    retained = report["joint_modules"]["orientation"]
    assert retained["fraction_of_source"] == .5
    assert retained["source_taxon_standardized_fraction"] == pytest.approx(1 / 3)
    assert report["endpoints"]["angle"]["milestones"]["attempted"]["observations"] == 4
    assert report["endpoints"]["angle"]["milestones"]["downloaded"]["observations"] == 3
    assert report["ecological_fitting_authorized"] is False
    json.dumps(report, allow_nan=False)


def test_missing_ledger_rows_are_unreported_not_failures():
    report = audit(source(), states().iloc[:2], {"orientation": ["angle"]})
    assert report["status"] == "COVERAGE_AUDIT_WITH_UNREPORTED_STATES"
    assert report["unreported_slots"] == 2
    assert report["endpoints"]["angle"]["terminal_states"] == {"unreported": 2, "usable": 2}
    assert report["expected_observation_endpoint_slots"] == 4


def test_joint_module_has_own_denominator_not_global_complete_case_filter():
    extra = pd.DataFrame({"obs_id": ["b"], "endpoint": ["chroma"], "state": ["usable"]})
    report = audit(source(), pd.concat([states(), extra]),
                   {"orientation": ["angle"], "colour": ["chroma"], "joint": ["angle", "chroma"]})
    assert report["joint_modules"]["orientation"]["observations"] == 2
    assert report["joint_modules"]["colour"]["observations"] == 1
    assert report["joint_modules"]["joint"]["observations"] == 1


def test_missing_exposure_and_no_exposure_are_not_passed_bias_tests():
    report = audit(source(), states(), {"orientation": ["angle"]})
    env = report["environmental_support"]
    assert env["vpd_month"]["status"] == "not_evaluable_exposure_absent"
    assert env["pr_month"]["strata"][-1]["stratum"] == "missing_exposure"
    assert env["pr_month"]["strata"][-1]["source_observations"] == 1


def test_constant_and_all_missing_exposure_are_explicit():
    x = source().assign(pr_month=7, vpd_month=np.nan)
    report = audit(x, states(), {"orientation": ["angle"]})
    assert report["environmental_support"]["pr_month"]["source_bin_edges"] == [7]
    assert report["environmental_support"]["vpd_month"]["status"] == "not_evaluable_no_finite_source_exposure"
    json.dumps(report, allow_nan=False)


@pytest.mark.parametrize("bad", ["non_native", "duplicate_source", "duplicate_state", "outside_source", "unknown_state", "unknown_endpoint"])
def test_invalid_baselines_or_state_ledger_fail_closed(bad):
    x, s = source(), states()
    if bad == "non_native":
        x.loc[0, "native_range_status"] = "introduced"
    elif bad == "duplicate_source":
        x.loc[0, "obs_id"] = "b"
    elif bad == "duplicate_state":
        s = pd.concat([s, s.iloc[:1]])
    elif bad == "outside_source":
        s.loc[0, "obs_id"] = "outside"
    elif bad == "unknown_state":
        s.loc[0, "state"] = "passed"
    else:
        s.loc[0, "endpoint"] = "unknown"
    with pytest.raises(ValueError):
        audit(x, s, {"orientation": ["angle"]})


def test_audit_does_not_emit_private_observation_identity():
    report = audit(source(), states(), {"orientation": ["angle"]})
    assert "obs_id" not in json.dumps(report)


def test_empty_disposition_ledger_preserves_full_source():
    report = audit(source(), states().iloc[:0], {"orientation": ["angle"]})
    assert report["unreported_slots"] == 4
    assert report["source_observations"] == 4
