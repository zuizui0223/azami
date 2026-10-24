"""Tests of source semantics and sign-free nuclear multiclimate statistic."""
import json
from pathlib import Path

import numpy as np
import pytest

from analysis.v4.t2_s1_s2_multiniche_phylo_bridge import (
    AXES, SOURCE_S1S2, SOURCE_S2, shrink_inverse,
    score_vectors, panel_inputs
)


def test_five_axes_match_existing_s1_without_preferred_direction():
    assert AXES == (
        "chelsa_bio12", "chelsa_bio15", "chelsa_bio01",
        "chelsa_bio04", "chelsa_rsds_mean"
    )
    s = json.loads(SOURCE_S1S2.read_text())
    t = json.loads(SOURCE_S2.read_text())
    assert s["S1"]["n_observations"] == 170066
    assert t["method"]["fixed_predictors"] == [
        "chelsa_bio01", "chelsa_bio04", "chelsa_bio12",
        "chelsa_bio15", "chelsa_rsds_mean"
    ]
    assert s["S2"]["source_species_with_sufficient_photo_observer_cell_support"] == 2


def test_shrink_inverse_symmetric_for_correlated_climate_and_sign_free_score():
    c = np.array([[1,0.999],[0.999,1]],float)
    inv = shrink_inverse(c)
    x = np.array([1.0,-0.3])
    center=np.zeros(2)
    assert np.isfinite(inv).all()
    assert score_vectors(x,center,inv)[0] == pytest.approx(
        score_vectors(-x,center,inv)[0])
    assert score_vectors(center,center,inv)[0] == pytest.approx(0)


def test_null_covariance_zero_axis_refuses_to_report_significance():
    with pytest.raises(ValueError,match="Degenerate"):
        shrink_inverse(np.diag([1.,0.]))
