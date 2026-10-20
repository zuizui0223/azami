"""Synthetic mathematical regression gates for the exploratory v4 script."""
import numpy as np
import pandas as pd
from analysis.v4.functional_pilot import zarr, fit_r2, rv, shapley, fit_univariate_within, GROUPS

def test_multivariate_r2_nonnegative_and_perfect_fit():
    rng=np.random.default_rng(7)
    x=rng.normal(size=(50,2))
    assert 0.999999 < fit_r2(x,x) <= 1.
    assert fit_r2(np.empty((50,0)),x) == 0.

def test_shapley_contributions_sum_to_incremental_r2():
    rng=np.random.default_rng(4)
    columns=[x for cols in GROUPS.values() for x in cols]
    frame=pd.DataFrame(rng.normal(size=(70,len(columns))),columns=columns)
    frame["latitude"]=rng.normal(size=70)
    frame["longitude"]=rng.normal(size=70)
    frame["response"]=3*frame["chelsa_rsds_mean"]+rng.normal(size=70)
    for spatial in (False,True):
        result=shapley(frame,["response"],include_spatial=spatial)
        assert np.isclose(sum(result["block_shapley_r2"].values()),
                          result["full_r2"]-result["baseline_r2"],atol=1e-10)
        assert all(v >= -1e-9 for v in result["block_shapley_r2"].values())
        assert result["block_shapley_r2"]["radiation"] > 0

def test_rv_is_symmetric_bounded_and_does_not_use_signed_correlation():
    rng=np.random.default_rng(18)
    a=rng.normal(size=(60,2))
    b=np.column_stack([a[:,0],rng.normal(size=60)])
    assert np.isclose(rv(a,b),rv(b,a),atol=1e-12)
    assert 0 <= rv(a,b) <= 1.000000001

def test_equal_taxon_weighting_differs_from_observation_weighting():
    rng=np.random.default_rng(23)
    taxa=np.repeat(["A","B"],[200,10])
    x=np.r_[rng.normal(0,1,200),rng.normal(0,1,10)]
    y=np.r_[x[:200]+rng.normal(0,.2,200),-2*x[200:]+rng.normal(0,.2,10)]
    frame=pd.DataFrame({"taxon_name":taxa,"obs_id":np.arange(210).astype(str),
                        "chelsa_rsds_mean":x,"value":y})
    all_fit=fit_univariate_within(frame)
    without=fit_univariate_within(frame,["A"])
    assert all_fit["obs"]==210
    assert without["obs"]==10
    assert all_fit["beta_observation_weighted"] > 0
    assert all_fit["beta_equal_taxon_weight"] < all_fit["beta_observation_weighted"]
    assert without["beta_observation_weighted"] < 0
