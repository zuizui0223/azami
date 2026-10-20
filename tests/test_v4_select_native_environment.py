import numpy as np
import pandas as pd
from analysis.v4.select_native_environment import CANDIDATES, PROTECTED, select, make_levels


def test_redundant_vpd_removed_without_access_to_traits():
    rng = np.random.default_rng(42)
    values = pd.DataFrame(rng.normal(size=(2000, 6)), columns=CANDIDATES)
    values["chelsa_vpd_mean"] = values["chelsa_bio01"]+.01*rng.normal(size=2000)
    result = select({level: values for level in ("raw", "between", "within")})
    assert result["status"] == "selected"
    assert result["history"][0]["removed"] == "chelsa_vpd_mean"
    assert set(PROTECTED) <= set(result["selected"])


def test_protected_conflict_stops_instead_of_dropping_hypothesis():
    rng = np.random.default_rng(42)
    values = pd.DataFrame(rng.normal(size=(2000, 6)), columns=CANDIDATES)
    values["chelsa_rsds_mean"] = values["chelsa_bio12"]
    result = select({"raw": values})
    assert result["status"] == "unresolved"
    assert set(PROTECTED) <= set(result["selected"])


def test_means_precede_complete_environment_filter():
    frame = pd.DataFrame({v: [100., 200., 600., 400.] for v in CANDIDATES})
    frame["accepted_key"] = [1, 1, 1, 2]
    frame.loc[2, "chelsa_vpd_mean"] = np.nan
    levels = make_levels(frame)
    assert len(levels["raw"]) == 3
    assert levels["between"]["chelsa_bio12"].tolist() == [300., 400.]
    assert levels["within"]["chelsa_bio12"].tolist() == [-200., -100., 0.]
