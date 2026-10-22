"""Source-contract and synthetic numerical tests for whole-S1 WN/BM/OU1.

Synthetic matrices are unit tests ONLY; empirical results require frozen
170,066-photo authenticated source + 52 archived placement scenarios.
"""
import numpy as np
import pytest

from analysis.v4.t2_all244_ou_bm_orientation import (
    _aicc, _gls, evaluate_models, regression_design, ALPHA_GRID,
    COLLIDED_TIPS, MATRIX_CASES,
)


def synthetic(n=24, seed=20261022):
    rng=np.random.default_rng(seed)
    # Positive-definite shared-path tree-like matrix with multiple
    # nested groups, root shared stem and independent terminal lengths.
    cluster=np.repeat(np.arange(6),n//6)
    C=0.4*np.ones((n,n))+0.8*(cluster[:,None]==cluster[None,:])+0.6*np.eye(n)
    X=np.column_stack([np.ones(n),rng.standard_normal((n,5))])
    y=2*X[:,3]+rng.standard_normal(n)*0.4
    return y,X,C


def test_three_model_AICc_on_one_identical_species_cohort():
    y,X,C=synthetic()
    out=evaluate_models(y,X,C,float(y.std(ddof=1)),0.05)
    assert out["n_species"]==24
    assert out["best_model_AICc"] in {"WN","BM","OU1"}
    assert np.isfinite([out["WN_AICc"],out["BM_AICc"],out["OU1_AICc"]]).all()
    assert abs(out["OU1_minus_BM_AICc"]-(
        out["OU1_AICc"]-out["BM_AICc"]))<1e-9
    assert out["OU1_alpha_times_median_root_to_tip_branch_length"]>=0
    assert out["selected_model_BIO12_deg_per_sample_SD"] is not None
    assert out["OU1_is_single_optimum_no_optimum_shift_test"] is True


def test_gls_identity_agrees_with_OLS_for_fixed_design():
    y,X,C=synthetic()
    fit=_gls(y,X,np.eye(len(y)))
    assert fit is not None
    expected=np.linalg.lstsq(X,y,rcond=None)[0]
    np.testing.assert_allclose(fit["beta"],expected,rtol=1e-8)
    assert np.isfinite(_aicc(fit["loglik"],len(y),X.shape[1]+1))


def test_missing_climate_or_nonpositive_bm_covariance_fails_closed():
    y,X,C=synthetic()
    broken=C.copy()
    broken[0,0]=-1
    with pytest.raises(ValueError,match="positive rooted"):
        evaluate_models(y,X,broken,float(y.std()),0.05)
    badX=np.ones((len(y),3))
    with pytest.raises(ValueError,match="Bad original five-climate design"):
        regression_design([{
           "mean_S1_image_axis_deg":float(i),
           "S1_five_climate_means":{e:1.0 for e in
               ("chelsa_bio01","chelsa_bio04","chelsa_bio12",
                "chelsa_bio15","chelsa_rsds_mean")}
       } for i in range(len(y))],True)


def test_source_guard_does_not_recode_photo_angle_as_botanical_downward():
    assert COLLIDED_TIPS=={
        "Cirsium_parryi","Cirsium_setidens","Cirsium_yezoense"}
    assert len(MATRIX_CASES)==6
    assert set(p[0] for p in MATRIX_CASES)=={
        "source194","clean191","n10_143"}
    assert 0==ALPHA_GRID[0]
    assert ALPHA_GRID[-1]==250
