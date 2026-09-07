"""Synthetic unit tests only; not empirical detector/trait validation."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import statsmodels.api as sm

SPEC = importlib.util.spec_from_file_location(
    "audit", Path(__file__).resolve().parents[1] / "analysis/v3/audit_revision_consistency.py")
audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(audit)


@pytest.mark.parametrize("p", [1, 2, 4, 9])
def test_hc3_matches_statsmodels(p):
    rng = np.random.default_rng(20260907 + p)
    x = rng.normal(size=(80, p))
    y = .7 * x[:, 0] + rng.normal(size=80) * (1 + np.abs(x[:, 0]))
    got = audit.ols_hc3(y, x)
    reference = sm.OLS(audit.zscore(y), sm.add_constant(audit.zscore(x))).fit(cov_type="HC3", use_t=False)
    assert got["beta"] == pytest.approx(reference.params[1], abs=1e-12)
    assert got["hc3_se"] == pytest.approx(reference.bse[1], abs=1e-12)
    assert [got["ci_low"], got["ci_high"]] == pytest.approx(reference.conf_int()[1], abs=1e-12)


def test_univariate_standardized_beta_is_correlation():
    x = np.arange(30.)
    y = x * x - 100
    assert audit.ols_hc3(y, x)["beta"] == pytest.approx(np.corrcoef(x, y)[0, 1])


def test_positive_rescaling_preserves_standardized_coefficient():
    x = np.arange(30.)
    y = np.sin(x) + .2 * x
    a = audit.ols_hc3(y, x)
    b = audit.ols_hc3(4 * y + 11, .1 * x - 273.15)
    assert a["beta"] == pytest.approx(b["beta"])
    assert a["hc3_se"] == pytest.approx(b["hc3_se"])


@pytest.mark.parametrize("x", [np.ones(10), np.r_[np.arange(9), np.nan]])
def test_bad_predictors_stop(x):
    with pytest.raises(ValueError):
        audit.ols_hc3(np.arange(10), x)


def test_singular_design_stops():
    x = np.arange(10.)
    with pytest.raises(ValueError, match="rank-deficient"):
        audit.ols_hc3(x*x, np.column_stack([x, x]))


@pytest.mark.parametrize("invalid", [None, "unknown", "", 2])
def test_unknown_boolean_stops(invalid):
    with pytest.raises(ValueError):
        audit.strict_bool(pd.Series([invalid]))


def test_false_string_is_false():
    assert not audit.strict_bool(pd.Series(["False"])).iloc[0]


def fixtures():
    scheduled = pd.DataFrame({"annotation_unit_id": ["h1", "h2"], "obs_id": [1, 2],
                              "photo_id": [11, 22], "taxon_name": ["A", "B"],
                              "orientation_angle_degrees": [40., 100.]})
    rows = []
    for _, head in scheduled.iterrows():
        for condition in sorted(audit.CONDITIONS):
            rows.append({"annotation_unit_id": head.annotation_unit_id, "obs_id": head.obs_id,
                         "photo_id": head.photo_id, "condition": condition,
                         "angle_deg": head.orientation_angle_degrees + (0 if condition == "baseline" else 3),
                         "usable": True, "head_clipped": False, "source_sha256": "a" * 64})
    return pd.DataFrame(rows), scheduled


def test_valid_records_and_nested_thresholds():
    rows, scheduled = fixtures()
    qc = audit.validate_linked_records(rows, scheduled)
    assert qc.stable_5.all() and qc.stable_10.all() and qc.stable_20.all()
    assert (qc.max_change == 3).all()


@pytest.mark.parametrize("change", ["missing_condition", "duplicate_condition", "missing_head", "bad_photo", "bad_hash", "different_hash", "bad_angle", "missing_column"])
def test_bad_records_fail_closed(change):
    rows, scheduled = fixtures()
    if change == "missing_condition": rows = rows.iloc[1:]
    elif change == "duplicate_condition": rows = pd.concat([rows, rows.iloc[:1]])
    elif change == "missing_head": rows = rows[rows.annotation_unit_id.ne("h1")]
    elif change == "bad_photo": rows.loc[0, "photo_id"] = 999
    elif change == "bad_hash": rows.loc[0, "source_sha256"] = "missing"
    elif change == "different_hash": rows.loc[0, "source_sha256"] = "b" * 64
    elif change == "bad_angle": rows.loc[0, "angle_deg"] = 181
    elif change == "missing_column": rows = rows.drop(columns="usable")
    with pytest.raises(ValueError):
        audit.validate_linked_records(rows, scheduled)


@pytest.mark.parametrize("change", ["unusable", "clipped", "baseline_mismatch", "large_shift", "failure"])
def test_rejected_head_is_not_stable(change):
    rows, scheduled = fixtures()
    if change == "unusable": rows.loc[0, "usable"] = False
    elif change == "clipped": rows.loc[0, "head_clipped"] = True
    elif change == "baseline_mismatch":
        rows.loc[(rows.annotation_unit_id == "h1") & (rows.condition == "baseline"), "angle_deg"] = 41
    elif change == "large_shift":
        rows.loc[(rows.annotation_unit_id == "h1") & (rows.condition == "x_plus_5pct"), "angle_deg"] = 61
    elif change == "failure":
        rows.loc[rows.annotation_unit_id == "h1", ["angle_deg", "source_sha256"]] = [np.nan, ""]
        rows.loc[rows.annotation_unit_id == "h1", "usable"] = False
    qc = audit.validate_linked_records(rows, scheduled).set_index("annotation_unit_id")
    assert not qc.loc["h1", "stable_20"]
    assert qc.loc["h2", "stable_5"]


def test_selection_is_invariant_to_head_order():
    n = 70
    h = pd.DataFrame({"annotation_unit_id": [str(i) for i in range(n)], "obs_id": range(n),
                       "photo_id": range(n), "taxon_name": ["A"] * n,
                       "orientation_status": ["usable"] * n, "orientation_angle_degrees": np.arange(n)})
    a, ah = audit.fixed_selection(h, ["A"])
    b, bh = audit.fixed_selection(h.sample(frac=1, random_state=1), ["A"])
    pd.testing.assert_frame_equal(a, b)
    pd.testing.assert_frame_equal(ah, bh)
    assert len(a) == 50
