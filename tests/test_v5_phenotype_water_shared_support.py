import numpy as np
from analysis.v5.phenotype_water_shared_support import (
    shared_water_mask,subset_pairs,loo_diagnostic,RAIN,SEASON,ENDPOINTS,
)
from analysis.v5.phenotype_water_interaction import difference_features

def pair_data(a,b):
    a=np.asarray(a,float);b=np.asarray(b,float)
    n=len(a)
    return dict(left=a,right=b,blocks=np.arange(n),tail_left=np.zeros((n,2)),
                tail_right=np.zeros((n,2)),native_climate=a.copy(),
                design=difference_features(a,b,np.zeros(5),np.ones(5)))

def test_intersection_and_immutable_native_reference():
    rng=np.random.default_rng(18)
    a=rng.uniform(-1,1,(100,5))
    b=rng.uniform(-1,1,(100,5))
    ref=pair_data(a,b)
    ta=a[:2].copy();tb=b[:2].copy()
    ta[1,RAIN]=10
    target=pair_data(ta,tb)
    mask=shared_water_mask(target,ref)
    np.testing.assert_array_equal(mask,[True,False])
    selected=subset_pairs(target,mask)
    assert len(selected["left"])==1
    assert len(selected["native_climate"])==2
    assert len(target["left"])==2

def test_fail_closed_when_common_support_is_small():
    rng=np.random.default_rng(18)
    a=rng.uniform(-1,1,(7,5))
    b=rng.uniform(-1,1,(7,5))
    ref=pair_data(a,b)
    rows=loo_diagnostic([ref,ref],ENDPOINTS[0])
    assert len(rows)==9
    assert all(r["status"]=="insufficient_support" for r in rows)
