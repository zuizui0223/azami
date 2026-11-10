"""S2 X-only selection and interaction covariance invariants."""
import numpy as np
from analysis.v4.s2_xonly_native_intro_interactions import _selected_sets,_cluster_multiy

def test_selected_sets_require_blind_X():
    obj={"schema":"azami_s2_outcome_blind_common_climate_selection_v1",
         "original_source_observation_n":297622,"trait_values_not_loaded":True,
         "past_results_P_values_and_RF_rankings_not_loaded":True,
         "results":[{"species":"Cirsium arvense",
           "selected_shared_covariate_set_X_only":{"selected_variables":["chelsa_bio12"]}},
           {"species":"Cirsium vulgare",
           "selected_shared_covariate_set_X_only":{"selected_variables":["chelsa_bio01"]}}]}
    assert len(_selected_sets(obj))==2

def test_cluster_multiy_gives_unbiased_OLS():
    rng=np.random.default_rng(3)
    n=40*12;groups=np.repeat(np.arange(40),12)
    X=np.column_stack([np.ones(n),rng.normal(size=n),rng.normal(size=n)])
    Y=np.column_stack([.3+X[:,1]*.7-X[:,2]*.2+rng.normal(0,.1,n)])
    beta,scores,corr,condition,ncluster=_cluster_multiy(X,Y,np.ones(n),groups)
    assert abs(beta[1,0]-.7)<.04
    assert scores.shape==(40,3,1) and ncluster==40
