import numpy as np
from analysis.v5.phenotype_availability_enrichment import partition,group_profile

def test_partition_and_separate_outside_core():
    an=np.array([1,1,0,0,1],bool)
    ai=np.array([1,1,1,1,0],bool)
    on=np.array([1,0,0,1,0],bool)
    groups=partition(an,ai,on)
    np.testing.assert_array_equal(np.array(list(groups.values())),np.eye(5,dtype=bool))

def test_equal_cell_weight_and_native_tails():
    v=np.array([[0.,0.,0.],[0.,0.,0.],[10.,10.,10.]])
    ix=np.array([0,0,1])
    cuts=np.tile([2.,4.,6.,8.],(3,1))
    rows=group_profile(v,ix,2,np.arange(2),cuts,np.zeros(3),np.full(3,10.))
    for r in rows:
        assert r["cells"]==2 and r["records"]==3
        assert r["mean"]==5.
        assert r["existing_low"]==r["existing_high"]==.5
        assert r["outside"]==0
    assert group_profile(v,ix,2,np.array([],int),cuts,np.zeros(3),np.full(3,10.))==[]
