import numpy as np
from scipy.optimize._numdiff import approx_derivative
from analysis.v5.phenotype_curve_regional_structure import (
    contribution,regional_weights,weighted_fit,weighted_objective,tail_diagnostic,
    X,K,fit,sign_counts,
)

def test_score_contribution():
    c=contribution([1.,-1.,2.,-2.])
    assert c["cancellation"]==1 and c["net"]==0
    assert c["largest_absolute_share"]==1/3

def test_equal_region_weights():
    w=regional_weights(np.array([0,0,0,1]))
    np.testing.assert_allclose(w,[2/3,2/3,2/3,2.])

def test_weighted_gradient_and_repetition():
    rng=np.random.default_rng(30)
    a=rng.dirichlet(np.ones(K),size=12)
    b=rng.dirichlet(np.ones(K),size=12)
    w=np.repeat([1.,2.,3.],4)
    beta=rng.normal(0,.1,6)
    f,g=weighted_objective(beta,a,b,w)
    numeric=approx_derivative(lambda t: weighted_objective(t,a,b,w)[0],beta).ravel()
    np.testing.assert_allclose(g,numeric,atol=1e-6)
    u=fit(a,b,X,np.arange(6))[0]
    e=weighted_fit(a,b,w,u)
    q=fit(np.repeat(a,w.astype(int),axis=0),np.repeat(b,w.astype(int),axis=0),X,np.arange(6))[0]
    np.testing.assert_allclose(e,q,atol=1e-5,rtol=0)

def test_paired_and_regional_tail_gains_distinguished():
    # One pair gains only low tail; one only high tail. Both in same region.
    d=np.array([[.1,0,-.1,0,0],[0,0,-.1,0,.1]])
    p=sign_counts(d,np.ones(2),"cell_pair","equal_pair","chroma")
    r=sign_counts(d.mean(0,keepdims=True),np.ones(1),"region","equal_region","chroma")
    get=lambda q: next(x["weight_fraction"] for x in q if x["lower_sign"]==x["upper_sign"]==1)
    assert get(p)==0 and get(r)==1
