"""An actual EAzami exact-state-map CTMC, not a substitution with photo PIC/PGLS."""
import numpy as np
from analysis.v4 import eazami_ud_real_four_environment_transition_regime as m

def test_same_EAzami_frozen_genealogy_and_source_coverage():
    assert m.EA_CORE_GIT_BLOB=="09731b563907b3ec2778ea04086d729771376f7e"
    assert m.EA_TREES_BLOB=="d7f8c0b96ee9c257c815db06a5a0aa9fdae4841e"
    assert m.EA_STATES_BLOB=="f655544e77f1659323c1e915567dacbb51b01cc3"
    assert m.EA_JP_CSV_SHA=="d115635f57e11d8a9f5c417bf28992a3d2a00f2a4d02ebb710c96d892ad0a9b2"
    assert m.EA_TW_CSV_SHA=="63b43b670f55a16b5557e414e37297a9ee18942f06714094fc1b5788ec07c70c"
    assert m.AXES==("chelsa_bio18","chelsa_bio05","chelsa_rsds_mean","elevation_m")
    assert m.PANELS["n5_primary"]==("5",12,7,5,792)
    assert m.PANELS["n10_strict_sensitivity"]==("10",9,5,4,126)
    assert m.PANELS["n3_resolution_sensitivity"]==("3",13,7,6,1716)

def test_finite_conditional_fraction_never_called_biological_replicate_probability():
    out=m.rank(np.array([1.5,2.5,3.5]),2.5,True)
    assert out["extreme_state_maps"]==2
    assert out["all_fixed_count_state_maps"]==3
    assert out["conditional_exact_state_assignment_fraction"]==2/3
    assert out["not_a_biological_replicate_or_frequentist_P_value"]

def test_synthetic_true_eazami_signed_branch_direction(monkeypatch):
    class MockCore:
        def fit_symmetric_q(self,tree,states):return .12
        def ctmc_likelihood_and_edge_joint(self,tree,states,q,need_joint):
            return .25,{1:np.array([[.25,.4],[.1,.25]])}
    edges={key:{1:float(i+1)} for i,key in enumerate(m.AXES)}
    vector,changes,q=m.exact_transition_vector(MockCore(),None,{},edges)
    assert np.allclose(vector,np.array([.6,1.2,1.8,2.4]))
    assert changes==.5
    assert q==.12
