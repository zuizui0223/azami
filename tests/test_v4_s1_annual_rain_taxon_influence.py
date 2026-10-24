"""Original S1 annual rainfall continuous-angle robustness, no biological U/D binning."""
import numpy as np
import pytest

from analysis.v4 import s1_annual_rain_taxon_influence as rain
from analysis.v4.joint_s1_t2_phylogeny import EXPOSURES


def synthetic(n=82,seed=20261021):
    rng=np.random.default_rng(seed)
    v=rng.standard_normal((n,5))*np.array([4,150,260,13,4])+np.array([12,900,1150,48,10])
    y=60+0.50*(v[:,2]/100)-.03*(v[:,3]/10)+0.5*v[:,0]+rng.normal(0,.5,n)
    return [
       {"tip":f"Cirsium_test_{i:03}","accepted_key":f"k{i}",
        "n_eligible_images":7+((i*5)%29),
        "source_concept_collision":False,
        "image_mean_deg":float(y[i]),
        "image_median_deg":float(y[i]+rng.normal(0,.4)),
        "S1_five_climate_means":dict(zip(EXPOSURES,v[i].tolist()))}
       for i in range(n)]


def test_five_original_CHELSA_covariates_and_per_100mm_physical_slope():
    rows=synthetic()
    names,X,y,raw=rain.prepare(rows,"image_mean_deg")
    assert len(names)==82
    assert X.shape==(82,6)
    assert tuple(EXPOSURES)==tuple(rain.EXPOSURES)
    out=rain.ols_hc3(X,y)
    assert out["rank"]==6 and out["df"]==76
    estimate=out["beta"][3]
    assert estimate==pytest.approx(.50,abs=.07), "BIO12 effect expressed per 100 mm"
    assert out["CI_HC3"][3,0]>0
    assert np.isfinite(out["leverage"]).all()


def test_robust_single_taxon_leave_out_recomputes_original_five_predictors():
    rows=synthetic()
    report=rain.audit_cohort(rows,"image_mean_deg",reps=110)
    bio=report["rain_exposure_robustness"]["chelsa_bio12"]
    assert report["n_original_S1_taxon_keys"]==82
    assert report["taxon_equal_weighting"] is True
    assert bio["increment_original_units"]==100.
    assert bio["estimate_image_degrees_per_increment"]>0
    assert bio["species_leave_one_out_n"]==82
    assert bio["species_leave_one_out_effect_all_positive"] is True
    assert bio["species_pairs_bootstrap_completed"]>=100
    assert bio["species_pairs_bootstrap_percentile95"][0]>0
    assert bio["largest_single_taxon_effect_on_slope"]["taxon"].startswith("Cirsium_test")


def test_shift_angle_zero_and_temperature_scale_invariance():
    rows=synthetic()
    _,X,y,_=rain.prepare(rows,"image_mean_deg")
    b=rain.ols_hc3(X,y)
    for row in rows:
        row["image_mean_deg"]+=20
    _,X2,y2,_=rain.prepare(rows,"image_mean_deg")
    d=rain.ols_hc3(X2,y2)
    assert b["beta"][3]==pytest.approx(d["beta"][3],abs=1e-9)
    assert d["beta"][0]==pytest.approx(b["beta"][0]+20,abs=1e-8)
    for row in rows:
        row["S1_five_climate_means"]["chelsa_bio12"]+=250
    _,X3,y3,_=rain.prepare(rows,"image_mean_deg")
    assert rain.ols_hc3(X3,y3)["beta"][3]==pytest.approx(d["beta"][3],abs=1e-9)


def test_fails_with_nonindependent_or_missing_five_climate_inputs():
    rows=synthetic()
    rows[1]["tip"]=rows[0]["tip"]
    with pytest.raises(ValueError,match="Unverified/source-collided"):
        rain.prepare(rows,"image_mean_deg")
    rows=synthetic()
    rows[0]["S1_five_climate_means"]["chelsa_bio15"]=float("nan")
    with pytest.raises(ValueError,match="nonfinite"):
        rain.prepare(rows,"image_mean_deg")
    rows=synthetic()
    for row in rows:
        row["S1_five_climate_means"]["chelsa_bio04"]=100.
    with pytest.raises(ValueError,match="rank deficient"):
        rain.prepare(rows,"image_mean_deg")
