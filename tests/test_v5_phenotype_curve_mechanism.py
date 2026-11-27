import numpy as np
from scipy.optimize._numdiff import approx_derivative
from analysis.v5.phenotype_curve_mechanism import (
    eta_profile,objective,components,fit_model,X,FOCAL,K,variance_diagnostic,
)

def test_variance_partition():
    # Two regions, two cell means per region, no within-cell variance.
    # Four equal cells have values -2,-1,1,2 along all three axes.
    v=np.repeat(np.array([-2.,-1.,1.,2.])[:,None],3,axis=1)
    m=np.column_stack([v,v*v,v*v])
    from analysis.v5.phenotype_curve_mechanism import region_moments
    regions,n=region_moments(m,np.repeat([0,1],2))
    r=components(regions,n,np.ones((1,2)))[0]
    np.testing.assert_allclose(r,np.repeat(np.array([2.5,2.25,.25,0.])[:,None],3,axis=1))

def test_profile_gradient_and_original_model():
    rng=np.random.default_rng(16)
    a=rng.dirichlet(np.ones(K),size=12)
    b=rng.dirichlet(np.ones(K),size=12)
    z=rng.normal(size=(12,5))
    beta=rng.normal(0,.1,6+2*len(FOCAL))
    axes=(2,3)
    loss,grad=objective(beta,a,b,z,axes)
    numerical=approx_derivative(lambda t: objective(t,a,b,z,axes)[0],beta,method="3-point").ravel()
    np.testing.assert_allclose(grad,numerical,atol=1e-6)
    fit=fit_model(a,b,z,())
    from analysis.v5.phenotype_within_species_enrichment import profile_state
    ls,_,p=eta_profile(a,b,fit,z,())
    old,_,_=profile_state(a,b,X,fit)
    np.testing.assert_allclose(ls.sum(),old,atol=1e-10)
    np.testing.assert_allclose(((a+b)*p).sum(1),1.,atol=1e-10)

def test_climate_fit_recovers_direction_in_synthetic_pair_profiles():
    rng=np.random.default_rng(21)
    n=100
    z=np.zeros((n,5));z[:,2]=rng.uniform(-2,2,n)
    pooled=2*rng.dirichlet(np.ones(K),size=n)
    truth=np.zeros(6+len(FOCAL));truth[6]=.4
    # Balance pooled allocation using same profiled intercept.
    _,_,p=eta_profile(pooled/2,pooled/2,truth,z,(2,))
    a=pooled*(1-p);b=pooled*p
    est=fit_model(a,b,z,(2,))
    assert .20<est[6]<.41
