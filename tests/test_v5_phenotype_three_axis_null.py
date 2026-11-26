import numpy as np
from itertools import product
from analysis.v5.phenotype_three_axis_null import partition_null,summarize_null

def test_positive_mass_and_label_null():
    b=np.array([[-1.,1.,0.],[-1.,0.,1.]])
    signs=np.array(list(product([1.,-1.],repeat=2)))
    a=partition_null(b,np.array([True,True,False]),signs,2)
    np.testing.assert_allclose(a,[[.5,.5],[.5,0.],[0.,.5],[1.,0.]])
    r={x["metric"]:x for x in summarize_null(a)}
    assert r["total_redistribution"]["observed"]==1.
    assert r["total_redistribution"]["p"]==.5
    assert r["existing_share"]["observed"]==.5
    assert r["existing_share"]["p"]==.75

def test_undefined_ratio():
    r={x["metric"]:x for x in summarize_null(np.zeros((2,2)))}
    assert r["existing_share"]["status"]=="undefined"
    assert r["existing_share"]["p"]==1.
    assert r["total_redistribution"]["p"]==1.
