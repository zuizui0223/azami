import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from analysis.v3 import environment_matrix as env


ROOT = Path(__file__).resolve().parents[1]


def test_contract_has_month_aligned_exposures_and_no_climate_core():
    contract = env.load_contract(ROOT / "analysis" / "v3" / "environment_exposure_contract.json")
    ids = [row["id"] for row in contract["monthly_candidates"]]
    assert {"pr_month", "tas_month", "tasmax_month", "rsds_month", "vpd_month", "sfcWind_month"}.issubset(ids)
    assert contract["selection_is_phenotype_blind"] is True
    assert "core" not in contract
    assert contract["redundancy_rule"]["BIO18_vs_GSP"].startswith("compare both against observation-month precipitation")


def test_month_url_is_deterministic_and_month_checked():
    contract = env.load_contract()
    assert env.month_url(contract, "pr", 7).endswith("/pr/1981-2010/CHELSA_pr_07_1981-2010_V.2.1.tif")
    with pytest.raises(ValueError):
        env.month_url(contract, "pr", 13)


def test_primary_observation_validation_requires_native_status():
    frame = pd.DataFrame({"obs_id":["a"], "latitude":[35.0], "longitude":[135.0], "observation_month":[7]})
    with pytest.raises(ValueError):
        env.validate_observation_frame(frame, require_native=True)
    frame["native_range_status"] = "native"
    out = env.validate_observation_frame(frame, require_native=True)
    assert list(out["obs_id"]) == ["a"]


def test_duplicate_observations_are_rejected():
    frame = pd.DataFrame({
        "obs_id":["a","a"], "latitude":[35,36], "longitude":[135,136],
        "observation_month":[7,8], "native_range_status":["native","native"]
    })
    with pytest.raises(ValueError):
        env.validate_observation_frame(frame, require_native=True)


def test_environment_diagnostics_detect_redundancy_without_traits():
    rng = np.random.default_rng(7)
    n = 500
    rain = rng.normal(size=n)
    matrix = pd.DataFrame({
        "pr_month": rain,
        "BIO18": rain * 0.98 + rng.normal(scale=0.05, size=n),
        "GSP": rain * 0.97 + rng.normal(scale=0.07, size=n),
        "rsds_month": rng.normal(size=n),
        "vpd_month": rng.normal(size=n),
    })
    variables = list(matrix.columns)
    report, tables = env.diagnose_environment(matrix, variables, threshold=0.80)
    assert report["trait_columns_read"] == 0
    assert report["representation_frozen"] is False
    clusters = [set(group) for group in report["redundancy_components"]]
    assert any({"pr_month", "BIO18", "GSP"}.issubset(group) for group in clusters)
    assert tables["coverage"]["coverage"].eq(1.0).all()
    assert set(tables["vif"]["variable"]) == set(variables)
    assert report["matrix"]["matrix_rank"] == len(variables)


def test_missingness_is_reported_not_imputed():
    matrix = pd.DataFrame({"pr_month":[1.0, np.nan, 3.0], "vpd_month":[10.0, 11.0, 12.0]})
    coverage = env.coverage_table(matrix, ["pr_month", "vpd_month"]).set_index("variable")
    assert coverage.loc["pr_month", "n_finite"] == 2
    assert coverage.loc["vpd_month", "n_finite"] == 3
