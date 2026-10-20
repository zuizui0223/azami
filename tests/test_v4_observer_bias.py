"""Synthetic mathematical gates for observer and taxon source-bias sensitivity.
No private observation IDs, coordinates or user data are in the repository.
"""
import numpy as np
import pandas as pd
import pytest
from analysis.v4.observer_bias import cells, weights, fit, TRAITS, ENV

def fixture():
    rng=np.random.default_rng(123)
    out=[]
    for t in range(24):
        for u in range(5):
            for repeat in range(5 if t==0 and u==0 else 1):
                n=len(out);x=rng.normal()+t*.15
                row={"obs_id":str(n),"taxon":str(t),"user_id":f"fixture-user-{t}-{u}",
                     "cell":0 if u<3 else 1,"latitude":30+t,"longitude":-120+t*8,
                     "corolla_lab_chroma_median":x+.3*t+rng.normal(scale=.1),
                     "orientation_angle_degrees_median":x*4+t*.3}
                row.update({name:x+.05*i*rng.normal() for i,name in enumerate(ENV)})
                out.append(row)
    return pd.DataFrame(out)

def test_taxon_cell_observer_weights_neutralize_repeated_uploads():
    f=fixture()
    w=weights(f,"equal_taxon_cell_observer")
    assert np.isclose(w.sum(),1)
    bytax=f.assign(w=w).groupby("taxon").w.sum()
    assert np.allclose(bytax,bytax.iloc[0])
    repeated=f.taxon.eq('0') & f.user_id.eq('fixture-user-0-0')
    single=f.taxon.eq('0') & f.user_id.eq('fixture-user-0-1')
    assert np.isclose(w[repeated].sum(),w[single].sum())
    assert not np.isclose(weights(f,"observation")[repeated].sum(),
                          weights(f,"observation")[single].sum())

def test_missing_observer_fail_closed():
    f=fixture();f.loc[0,"user_id"]=None
    with pytest.raises(ValueError,match="Missing original user"):
        weights(f,"equal_taxon_cell_observer")

def test_regr_same_taxon_under_three_scales():
    f=fixture()
    for scale in ("pooled","within","among"):
        result=fit(f,"equal_taxon_cell_observer",TRAITS["colour"][0],ENV[0],scale)
        assert result["n_taxa"]==24 and result["n_observations"]==len(f)
        assert np.isfinite(result["beta"])
        assert np.isfinite(result["ci95"]).all()
    for scheme in ("observation","equal_taxon","equal_taxon_cell","equal_taxon_cell_observer"):
        result=fit(f,scheme,TRAITS["presentation"][0],ENV[2],"among")
        assert np.isfinite(result["beta"])

def test_dateline_regions_are_distinct_and_deterministic():
    k=cells(np.array([34,34]),np.array([179.9,-179.9]))
    assert k[0]!=k[1]
    assert (k==cells(np.array([34,34]),np.array([179.9,-179.9]))).all()


def test_observer_balanced_covariance_bootstrap_is_deterministic():
    from analysis.v4.observer_covariance import integration
    f=fixture()
    from analysis.v4.observer_bias import ALL_TRAITS
    f["corolla_lab_lightness_median"]=f["corolla_lab_chroma_median"]*0.3+60
    f["corolla_hue_sin_median"]=np.sin(f["corolla_lab_chroma_median"])
    f["corolla_hue_cos_median"]=np.cos(f["corolla_lab_chroma_median"])
    a=integration(f,reps=16,seed=219)
    b=integration(f,reps=16,seed=219)
    assert a==b and a["n_taxa"]==24
    assert -1<a["matched_taxon_bootstrap_95pct"][0]<1


def test_pooled_covariance_equals_exact_within_plus_between_parts():
    from analysis.v4.observer_effort_covariance import covariance_parts
    f=fixture()
    from analysis.v4.observer_effort_covariance import SCHEMES
    for scheme in SCHEMES:
        row=covariance_parts(f,scheme,"corolla_lab_chroma_median","chelsa_bio01")
        assert np.isclose(row["pooled_standardized_coefficient"],
                          row["between_contribution"]+row["within_contribution"],atol=1e-10)
        assert row["identity_residual"] < 1e-10
        if scheme=="observation":
            assert row["top_two_taxon_weight_mass"] > 2/len(f.taxon.unique())
