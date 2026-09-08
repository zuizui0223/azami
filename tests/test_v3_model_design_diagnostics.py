import json

import numpy as np
import pandas as pd
import pytest

from analysis.v3.model_design_diagnostics import diagnose


def frame():
    rows = []
    for taxon in range(6):
        for j, offset in enumerate([-2., -1., 1., 2.]):
            rows.append({"obs_id": f"{taxon}-{j}", "accepted_key": str(taxon), "weight": 1.,
                         "pr_month": 10 * taxon + offset, "rsds_month": 10 * taxon - offset})
    return pd.DataFrame(rows)


def check(data, nuisance=None):
    return diagnose(data, ["pr_month", "rsds_month"], nuisance or [], weight_column="weight", cohort_id="synthetic_only")


def test_pooled_raw_correlation_cannot_substitute_for_either_scale():
    result = check(frame())
    assert result["raw"]["full_predictor_rank"] is True
    assert result["raw"]["correlation"][0][1] > .9
    assert result["within"]["correlation"][0][1] == pytest.approx(-1)
    assert result["among"]["correlation"][0][1] == pytest.approx(1)
    assert result["within"]["full_predictor_rank"] is False
    assert result["among"]["full_predictor_rank"] is False
    assert result["taxa_without_full_local_predictor_rank"] == 6
    assert result["ecological_models_executed"] == 0
    json.dumps(result, allow_nan=False)


def test_nuisance_alias_is_not_rescaled_roundoff_information():
    x = frame()
    x["timing"] = x["pr_month"]
    result = check(x, ["timing"])
    assert result["within_nuisance_residualized"]["matrix_rank"] == 0
    assert result["among_nuisance_residualized"]["matrix_rank"] == 0
    assert set(result["within_nuisance_residualized"]["constant_variables"]) == {"pr_month", "rsds_month"}


@pytest.mark.parametrize("bad", ["nan", "infinite", "zero_weight", "duplicate_id"])
def test_missingness_is_explicit_not_silently_dropped(bad):
    x = frame()
    if bad == "nan":
        x.loc[0, "pr_month"] = np.nan
    elif bad == "infinite":
        x.loc[0, "pr_month"] = np.inf
    elif bad == "zero_weight":
        x.loc[0, "weight"] = 0
    else:
        x.loc[0, "obs_id"] = x.loc[1, "obs_id"]
    with pytest.raises(ValueError):
        check(x)


def test_weighted_within_centering_uses_declared_weights():
    x = frame()
    x["weight"] = np.tile([1., 2., 3., 4.], 6)
    result = check(x)
    assert result["within"]["correlation"][0][1] == pytest.approx(-1)
    assert result["weight_column"] == "weight"


def test_unrequested_response_column_does_not_enter_diagnostics():
    x = frame()
    x["response"] = "do not convert me to a number"
    result = check(x)
    assert result["trait_values_read"] == 0
    assert "response" not in result["raw"]["variables"]


def test_unequal_taxon_sizes_cannot_turn_centering_dust_into_information():
    x = pd.DataFrame({"obs_id": [str(i) for i in range(70)],
                      "accepted_key": ["A"] * 30 + ["B"] * 40,
                      "weight": 1., "pr_month": .1})
    result = diagnose(x, ["pr_month"], [], weight_column="weight", cohort_id="constant")
    assert result["raw"]["matrix_rank"] == 0
    assert result["within"]["matrix_rank"] == 0
    assert result["within_nuisance_residualized"]["matrix_rank"] == 0


def test_empty_cohort_is_rejected_explicitly():
    with pytest.raises(ValueError, match="Empty realized model cohort"):
        check(frame().iloc[:0])


def test_among_alias_tolerance_is_invariant_to_predictor_origin():
    x = frame()
    original = diagnose(x, ["pr_month"], [], weight_column="weight", cohort_id="before")
    x["pr_month"] += 1e12
    shifted = diagnose(x, ["pr_month"], [], weight_column="weight", cohort_id="after")
    assert original["among_nuisance_residualized"]["matrix_rank"] == 1
    assert shifted["among_nuisance_residualized"]["matrix_rank"] == 1
