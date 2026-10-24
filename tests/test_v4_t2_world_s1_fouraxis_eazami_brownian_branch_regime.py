"""World original S1 0-180 photo angle, not made-up EAzami discrete U/D."""
from io import StringIO
import numpy as np
import pytest
from Bio import Phylo
from analysis.v4 import t2_world_s1_fouraxis_eazami_brownian_branch_regime as w

def test_exact_world_four_axis_source_and_noncategorical_transition_contract():
    assert w.FOUR==("chelsa_bio18","chelsa_bio05","chelsa_rsds_mean","elevation_m")
    assert w.ORIGINAL_SHA=="47b32972b20564ade163907ef6fdf3eb2a69f8b42c8d87b8cf4ea75ca1682dc1"
    assert w.FOUR_SHA=="823ecb8442eb03598253cb33f3e16bac86538c6c2652f48ca2dc0886d0e19626"
    assert w.ITS_SHA=="8fd10d299b4d520ce08f2552993b02d3cc1b10350f0fb68ea24874dbb9ccfe92"

def synthetic_tree():
    return Phylo.read(StringIO("(((a:1,b:1):2,(c:1,d:1):2):2,((e:1,f:1):2,(g:1,h:1):2):2);"),"newick")

def test_EAzami_branch_weighted_squared_change_reconstruction_is_linear():
    tree=synthetic_tree()
    labels=list("abcdefgh")
    D,inner=w.brownian_branch_operator(tree,labels)
    assert D.shape==(14,8)
    assert inner.sum()==6
    assert np.max(np.abs(D@np.ones(8)))<1e-8
    assert np.max(np.abs(D@np.arange(8,dtype=float)))>0
    four=np.column_stack([np.arange(8)/5,-np.arange(8)/9,
                          np.sin(np.arange(8)),np.cos(np.arange(8))])
    angle=np.array([10,11,20,21,75,76,90,91.])
    perm=np.vstack([angle,np.roll(angle,2),angle[::-1]])
    scores=w.branch_score(D,four,perm,np.ones(len(D),dtype=bool))
    assert scores.shape==(3,4)
    assert scores[0,0]>0 and scores[0,1]<0

def test_four_axis_maxT_and_one_sided_BIO18_999_null():
    sim=np.random.default_rng(20261024).normal(0,.5,(999,4))
    stats=w.rank(np.array([.4,.1,-.2,.6]),sim)
    assert stats["finite_tip_permutation_not_count_preserving_U_D_CTMP"]
    assert set(stats["four_axis_two_sided_fraction"])==set(w.FOUR)
    assert set(stats["four_axis_maxT_fraction"])==set(w.FOUR)
    assert 0<stats["BIO18_positive_direction_one_sided_fraction"]<=1
    assert stats["four_axis_maxT_fraction"]["chelsa_bio18"]>=stats["four_axis_two_sided_fraction"]["chelsa_bio18"]

def test_finite_8_tip_panel_yields_not_ancient_direction_event_count():
    tree=synthetic_tree()
    names=list("abcdefgh")
    source={name:{"mean_image_angle_deg":str(10+i*10),
        "n_photo_complete_four_axis":"12",
        "chelsa_bio18":str(100+i*15+(.1 if i%2 else .2)),
        "chelsa_bio05":str(10+i*i/5),
        "chelsa_rsds_mean":str(12+np.cos(i)),
        "elevation_m":str(100+3*i*i)} for i,name in enumerate(names)}
    r=w.evaluate_panel(source,names,[("constructed_synthetic_test_tree",tree)],199,20261024)
    assert r["n_original_source_names"]==8
    assert r["n_original_fouraxis_head_photos"]==96
    assert r["n_continuous_photo_tip_permutations"]==199
    assert r["not_Mk_transition_model_or_botanical_U_D"]
    assert not any("u_to_d_events" in key for key in r)

def test_fails_unsupported_too_few_tip_permutations():
    with pytest.raises(ValueError,match="Insufficient"):
        w.evaluate_panel({},[],[],19,20261024)
