import numpy as np
from analysis.v5.phenotype_endpoint_scan import empirical_rank, cell_average_matrix, permute, circular_summary

def test_native_ranks_ties():
    r=empirical_rank(np.array([-1,0,1,2,3]),np.array([0.,1.,1.,2.]),np.ones(4))
    np.testing.assert_allclose(r,[0,.125,.5,.875,1])

def test_cell_equal():
    x=cell_average_matrix([1.,3.,5.],np.array([0,0,2]),3)
    np.testing.assert_allclose(x[:,0],[2,np.nan,5],equal_nan=True)

def test_exact_block_exchange():
    score,p,exact,n=permute(np.ones((3,2)),np.ones(3))
    assert exact and n==8 and p==.25 and score==1

def test_hue_seam():
    a=np.radians([359,1])
    m=np.array([np.sin(a).mean(),np.cos(a).mean()])
    result=circular_summary(m)
    assert min(result["mean_hue_degrees"],360-result["mean_hue_degrees"])<1e-6
