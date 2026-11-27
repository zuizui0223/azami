import numpy as np
from scipy.optimize._numdiff import approx_derivative
from analysis.v5.phenotype_water_curve_reconstruction import (
    design,profile,objective,fit_water,allocate,reconstruct,X,K,
    CHROMA,ANGLE,MODELS,
)

def test_signed_orientation_and_balance():
    rng=np.random.default_rng(71)
    a=rng.dirichlet(np.ones(K),size=8)
    b=rng.dirichlet(np.ones(K),size=8)
    x=rng.normal(size=(8,2));beta=np.array([.2,-.1,.04,.15])
    eta=design(x,True)@beta
    l,p=profile(a,b,eta)
    lr,pr=profile(b,a,design(-x,True)@beta)
    np.testing.assert_allclose(p,1-pr,atol=1e-10)
    np.testing.assert_allclose(l,lr,atol=1e-10)
    np.testing.assert_allclose(((a+b)*p).sum(1),1.,atol=1e-10)

def test_gradient():
    rng=np.random.default_rng(71)
    a=rng.dirichlet(np.ones(K),size=15)
    b=rng.dirichlet(np.ones(K),size=15)
    x=rng.normal(size=(15,2));z=design(x,True)
    w=rng.uniform(.5,1.5,size=15);beta=np.array([.1,.2,-.1,.05])
    f,g=objective(beta,a,b,z,w,4.)
    q=approx_derivative(lambda t:objective(t,a,b,z,w,4.)[0],beta).ravel()
    np.testing.assert_allclose(g,q,atol=1e-6)

def test_recover_water_tilt_and_shape():
    rng=np.random.default_rng(27)
    n=180
    x=rng.uniform(-2,2,size=(n,2))
    m=2*rng.dirichlet(np.ones(K),size=n)
    truth=np.array([.2,-.17,-.12,.16])
    _,p=profile(m/2,m/2,design(x,True)@truth)
    a,b=m*(1-p),m*p
    beta=fit_water(a,b,x,np.ones(n),True,4.)
    np.testing.assert_allclose(beta,truth,atol=.015,rtol=0)
    l,pa,pb=allocate(dict(a=a,b=b,x=x),"species_shape",beta)
    np.testing.assert_allclose(pa+pb,m,atol=1e-12)
    assert np.isfinite(l).all()

def test_exact_curve_reconstruction_when_allocations_equal_observed():
    rng=np.random.default_rng(55)
    n=100;m=2*rng.dirichlet(np.ones(K),size=n)
    true=np.array([.2,-.1,.03,-.15,.12,.01])
    _,p=profile(m/2,m/2,np.broadcast_to(X@true,m.shape))
    a,b=m*(1-p),m*p
    rows,beta=reconstruct(a,b,a,b,true)
    np.testing.assert_allclose(beta,true,atol=1e-7,rtol=0)
    np.testing.assert_allclose([r["observed_relative_enrichment"] for r in rows],
                               [r["reconstructed_relative_enrichment"] for r in rows],atol=1e-6)
