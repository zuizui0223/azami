import numpy as np
from analysis.v5.phenotype_water_occupancy import wquant,signed_profiles,climate_band_profile,RAIN,NAMES,K

def test_weighted_climate_quantile():
    np.testing.assert_allclose(wquant([0,10,20],[1,7,2],[.1,.5,.9]),[0,10,20])

def test_reversing_pairs_preserves_aligned_direction():
    # Same lower-to-higher climate contrast, opposite novelty direction.
    a=np.zeros((2,K));b=np.zeros((2,K))
    a[0,0]=1;b[0,100]=1
    a[1,100]=1;b[1,0]=1
    left=np.zeros((2,5));right=np.zeros((2,5))
    left[:,RAIN]=[100,200];right[:,RAIN]=[200,100]
    rows=signed_profiles(dict(hist_left=a,hist_right=b,left=left,right=right,blocks=np.array([1,2])),RAIN)
    aligned=[r for r in rows if r["direction"]=="aligned_lower_to_higher" and r["trait"]==NAMES[0]]
    np.testing.assert_allclose([r["toward_higher_climate"] for r in aligned],[-1,0,0,0,1])
    contrib=[r["contribution_to_total"] for r in rows if r["direction"]!="aligned_lower_to_higher" and r["trait"]==NAMES[0] and r["native_quintile"]==5 and "contribution_to_total" in r]
    assert abs(sum(contrib))<1e-12

def test_band_composition_normalizes_and_keeps_sparse():
    h=np.ones((4,3,5))/5
    x=np.zeros((4,5));x[:,RAIN]=[10,20,30,40]
    rows=climate_band_profile(dict(h=h,x=x,g=np.array([1,1,2,2])),RAIN,[15,25,35,45])
    for band in range(1,5):
        rr=[r for r in rows if r["climate_band"]==band and r["trait"]==NAMES[0]]
        assert all(r["status"]=="sparse" for r in rr)
        np.testing.assert_allclose(sum(r["frequency"] for r in rr),1)
