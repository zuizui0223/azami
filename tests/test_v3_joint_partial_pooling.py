import numpy as np
import pytest
from scipy.linalg import block_diag, helmert

from analysis.v3.joint_partial_pooling import JointSlopeLikelihood, fit_joint_partial_pooling


def fixture(seed=981):
    rng = np.random.default_rng(seed)
    groups = np.repeat(['a','b','c','d'],[36,19,3,1])
    x = rng.normal(size=(len(groups),2))
    x[groups=='c',1] = x[groups=='c',0]  # local alias is not a deletion rule
    nuisance = np.column_stack([.5*x[:,0]+rng.normal(size=len(groups)),rng.normal(size=len(groups))])
    slopes = np.array([[1.4,-.3],[.7,.5],[2.,.2],[0.,1.]])
    y = np.sum(x*slopes[np.searchsorted(['a','b','c','d'],groups)],axis=1)+nuisance@np.array([.5,-.2])+rng.normal(0,.7,len(groups))
    return y,x,groups,nuisance


def dense_reference(y,x,groups,nuisance,likelihood,lam):
    """Independent full observation-contrast covariance; tiny test cases only."""
    contrast = block_diag(*[helmert(sum(groups==g),full=False) if sum(groups==g)>1 else np.zeros((0,1)) for g in likelihood.taxa])
    interaction = np.column_stack([x*(groups==g)[:,None] for g in likelihood.taxa])
    design = contrast@np.column_stack([x,nuisance@likelihood.nuisance_map])
    random = contrast@interaction
    yy = contrast@y
    prior = np.diag(np.tile(lam,len(likelihood.taxa)))
    v = np.eye(len(yy))+random@prior@random.T
    vi = np.linalg.inv(v)
    info = design.T@vi@design
    alpha = np.linalg.solve(info,design.T@vi@yy)
    residual = yy-design@alpha
    rss = residual@vi@residual
    sigma2 = rss/likelihood.df
    criterion = likelihood.df*np.log(sigma2)+np.linalg.slogdet(v)[1]+np.linalg.slogdet(info)[1]
    modes = prior@random.T@vi@residual
    selector = np.tile(np.column_stack([np.eye(x.shape[1]),np.zeros((x.shape[1],likelihood.nuisance_rank))]),(len(likelihood.taxa),1))
    mapping = selector-prior@random.T@vi@design
    covariance = sigma2*(prior-prior@random.T@vi@random@prior+mapping@np.linalg.inv(info)@mapping.T)
    return alpha, sigma2, criterion, modes.reshape(len(likelihood.taxa),-1), covariance


@pytest.mark.parametrize('lam',[[0.,0.],[.2,1.4],[0.,2.],[3.,0.]])
def test_block_likelihood_modes_and_full_prediction_covariance_match_dense(lam):
    y,x,g,b = fixture()
    likelihood = JointSlopeLikelihood(y,x,g,b)
    evaluation = likelihood.evaluate(lam)
    result = likelihood.result(evaluation)
    alpha,sigma2,criterion,modes,covariance = dense_reference(y,x,g,b,likelihood,np.array(lam))
    np.testing.assert_allclose(result.hypermean,alpha[:2],atol=1e-10)
    np.testing.assert_allclose(result.taxon_predictions,alpha[:2]+modes,atol=1e-10)
    np.testing.assert_allclose(result.prediction_error_covariance,covariance,atol=1e-10)
    assert result.residual_variance == pytest.approx(sigma2,rel=1e-10)
    assert result.reml_criterion == pytest.approx(criterion,rel=1e-10)
    assert np.linalg.eigvalsh(result.prediction_error_covariance).min() > -1e-10
    assert result.support.retained_in_joint_model.all()
    assert result.support.set_index('taxon').loc['c','within_predictor_rank']==1
    assert result.support.set_index('taxon').loc['d','no_local_slope_information']


def test_analytic_variance_gradient_matches_central_difference():
    likelihood = JointSlopeLikelihood(*fixture())
    lam = np.array([.4,1.2])
    exact = likelihood.evaluate(lam)['gradient']
    finite = []
    for j in range(2):
        step = np.eye(2)[j]*1e-5
        finite.append((likelihood.evaluate(lam+step)['criterion']-likelihood.evaluate(lam-step)['criterion'])/2e-5)
    np.testing.assert_allclose(exact,finite,rtol=1e-6,atol=1e-6)


def test_zero_heterogeneity_retains_hypermean_uncertainty_and_cross_taxon_dependence():
    likelihood = JointSlopeLikelihood(*fixture())
    result = likelihood.result(likelihood.evaluate([0,0]))
    for i in range(4):
        np.testing.assert_allclose(result.taxon_predictions[i],result.hypermean)
        for j in range(4):
            np.testing.assert_allclose(result.prediction_error_covariance[2*i:2*i+2,2*j:2*j+2],result.fixed_covariance[:2,:2])
    assert result.fixed_covariance[0,0] > 0


def test_no_information_taxon_receives_prediction_not_an_invented_local_estimate():
    likelihood = JointSlopeLikelihood(*fixture())
    result = likelihood.result(likelihood.evaluate([.3,.8]))
    np.testing.assert_allclose(result.taxon_predictions[-1],result.hypermean)
    expected = result.fixed_covariance[:2,:2]+np.diag(result.tau2)
    np.testing.assert_allclose(result.prediction_error_covariance[-2:,-2:],expected)


def test_exact_redundant_nuisance_is_compressed_not_counted_twice():
    y,x,g,b = fixture()
    first = JointSlopeLikelihood(y,x,g,b)
    second = JointSlopeLikelihood(y,x,g,np.column_stack([b,2*b[:,0],np.ones(len(y))]))
    a = first.result(first.evaluate([.5,.2]))
    c = second.result(second.evaluate([.5,.2]))
    assert first.nuisance_rank==second.nuisance_rank==2
    assert a.residual_df==c.residual_df
    np.testing.assert_allclose(a.hypermean,c.hypermean,atol=1e-10)
    np.testing.assert_allclose(a.prediction_error_covariance,c.prediction_error_covariance,atol=1e-10)


def test_global_alias_cannot_be_cured_by_shrinkage_or_silent_exposure_removal():
    y,x,g,b = fixture()
    with pytest.raises(ValueError,match='globally aliased'):
        JointSlopeLikelihood(y,x,g,x[:,:1])
    with pytest.raises(ValueError,match='globally aliased'):
        JointSlopeLikelihood(y,np.column_stack([x[:,0],x[:,0]]),g,b)


def test_missing_rows_not_dropped_and_variance_parameters_not_repaired():
    y,x,g,b = fixture()
    bad = y.copy(); bad[0] = np.nan
    with pytest.raises(ValueError,match='no silent row deletion'):
        JointSlopeLikelihood(bad,x,g,b)
    likelihood = JointSlopeLikelihood(y,x,g,b)
    for lam in ([1,-1],[1,np.nan],[1]):
        with pytest.raises(ValueError,match='nonnegative'):
            likelihood.evaluate(lam)


def test_joint_fit_reports_all_optimizer_attempts_and_sparse_taxa():
    result = fit_joint_partial_pooling(*fixture())
    assert len(result.optimizer_attempts)==4 and any(a['accepted'] for a in result.optimizer_attempts)
    assert len(result.support)==4 and result.taxon_predictions.shape==(4,2)
    assert not any(a['upper_bound_hit'] for a in result.optimizer_attempts if a['accepted'])
    assert result.uncertainty_scope=='conditional_working_model_not_spatially_calibrated'


def test_zero_variance_boundary_can_be_selected_exactly():
    x = np.tile([-1.,0.,1.],20)[:,None]
    g = np.repeat(np.arange(20),3)
    noise = np.tile([1.,-2.,1.],20)
    y = 1.3*x[:,0]+noise
    result = fit_joint_partial_pooling(y,x,g,np.zeros((len(y),0)))
    np.testing.assert_allclose(result.variance_ratios,0,atol=1e-8)
    np.testing.assert_allclose(result.hypermean,[1.3],atol=1e-9)


def test_common_unit_change_does_not_change_fitted_values_or_information():
    y,x,g,b = fixture()
    a = JointSlopeLikelihood(y,x,g,b)
    units = np.array([10.,.1])
    c = JointSlopeLikelihood(y,x*units,g,b)
    first = a.result(a.evaluate([.4,.8]))
    second = c.result(c.evaluate(np.array([.4,.8])/units**2))
    np.testing.assert_allclose(first.taxon_predictions,second.taxon_predictions*units,atol=1e-10)
    np.testing.assert_allclose(first.residuals,second.residuals,atol=1e-10)
