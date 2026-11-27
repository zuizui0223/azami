import numpy as np
from analysis.v5.phenotype_native_periphery import comparable_transitions,transition_features,subset_pairs

def test_transition_matching_uses_both_endpoints():
    # Identical direction is not sufficient when absolute climates differ.
    native=dict(x=np.array([[0.,0.],[.2,.2],[3.,3.],[3.2,3.2]]))
    intro=dict(x=np.array([[3.1,3.1],[3.3,3.3],[7.,7.],[7.2,7.2]]))
    p=(np.array([0,2]),np.array([1,3]),np.array([1,2]))
    a,b=comparable_transitions(native,p,intro,p,np.zeros(2),np.ones(2))
    np.testing.assert_array_equal(a,[1]);np.testing.assert_array_equal(b,[0])
    assert len(subset_pairs(p,a)[0])==1

def test_empty_native_pairs():
    d=dict(x=np.zeros((2,5)))
    a=(np.empty(0,int),)*3
    b=(np.array([0]),np.array([1]),np.array([1]))
    i,j=comparable_transitions(d,a,d,b,np.zeros(5),np.ones(5))
    assert len(i)==len(j)==0
