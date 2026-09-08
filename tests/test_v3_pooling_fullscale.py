"""Regression at the intended numerical scale, using generated data only."""
import numpy as np

from analysis.v3.hierarchical_ecology import spherical_basis
from analysis.v3.joint_partial_pooling import fit_joint_partial_pooling


def test_full_native_sized_synthetic_fit_meets_original_stationarity_threshold():
    rng = np.random.default_rng(2026090871)
    n, taxa, p = 319244, 354, 9
    weights = rng.lognormal(0,1.1,taxa-2)
    counts = np.r_[1,3,1+np.floor((n-taxa-2)*weights/weights.sum()).astype(int)]
    counts[-1] += n-int(counts.sum())
    groups = np.repeat(np.arange(taxa),counts)
    factor = rng.normal(size=n)
    x = np.sqrt(.65)*factor[:,None]+np.sqrt(.35)*rng.normal(size=(n,p))
    x[groups==1,1] = x[groups==1,0]
    lat, lon = rng.uniform(-70,70,n), rng.uniform(-170,170,n)
    nuisance = np.column_stack([.4*x[:,0]+rng.normal(size=n),rng.normal(size=n),spherical_basis(lat,lon)])
    mu = np.r_[0.,np.full(p-1,.3)]
    slopes = mu+rng.normal(0,.45,size=(taxa,p))
    y = rng.normal(size=taxa)[groups]+np.sum(x*slopes[groups],axis=1)+nuisance@np.linspace(.1,.4,nuisance.shape[1])+rng.normal(size=n)
    result = fit_joint_partial_pooling(y,x,groups,nuisance)
    assert len(result.taxa) == taxa
    assert result.prediction_error_covariance.shape == (taxa*p,taxa*p)
    assert result.support.n_observations.sum() == n
    assert all(a['projected_gradient_max'] <= 1e-4 for a in result.optimizer_attempts if a['accepted'])
    assert any(a['accepted'] for a in result.optimizer_attempts)
