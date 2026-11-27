import numpy as np
from analysis.v5.phenotype_environment_explanation import (
    match_one,rematch,partition_density,local_weights,CALIPER
)
def test_matching_max_cardinality_not_greedy():
    # a0 can match b0 or b1, a1 can match b0 only.
    a=np.array([[0.],[.4]])
    b=np.array([[.1],[-.4]])
    l,r=match_one(a,b,(0,))
    assert len(l)==2
    assert set(zip(l,r))=={(0,1),(1,0)}
def test_caliper_does_not_force_match():
    l,r=match_one(np.array([[0.],[2.]]),np.array([[.1],[4.]]),(0,))
    np.testing.assert_array_equal(l,[0]);np.testing.assert_array_equal(r,[0])
def test_region_constraint():
    z=np.array([[0.],[2.],[2.],[0.]])
    a,b,g=rematch(np.array([0,1]),np.array([2,3]),np.array([1,2]),z,(0,))
    assert len(a)==len(b)==len(g)==0
def test_density_partition_including_no_overlap():
    a=np.array([[1.,2.],[3.,4.],[5.,6.]])
    b=np.array([[3.,2.],[4.,6.],[9.,8.]])
    for l,r in ((np.array([0,2]),np.array([0,1])),(np.array([],int),np.array([],int)),(np.arange(3),np.arange(3))):
        f,s,c,n=partition_density(a,b,l,r)
        np.testing.assert_allclose(f,c+n)
def test_common_region_weights():
    a=np.array([0.,.1,0.,.3,3.])
    b=np.array([.2,1.,.1,.4,0.])
    g=np.array([0,0,1,1,2])
    wa,wb,ea,eb,er,n=local_weights(a,b,g,0.)
    assert n==2
    np.testing.assert_allclose(wa.sum(),1.)
    np.testing.assert_allclose(wb.sum(),1.)
    for r in (0,1,2):
        np.testing.assert_allclose(wa[g==r].sum(),wb[g==r].sum())
    assert wa[4]==wb[4]==0.
