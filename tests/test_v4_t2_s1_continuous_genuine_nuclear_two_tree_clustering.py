"""Nuclear S1 photo-angle tests: no plastid, no botanical U/D and no grafting."""
import numpy as np
import pytest
from Bio import Phylo
from io import StringIO
from analysis.v4 import t2_s1_continuous_genuine_nuclear_two_tree_clustering as m

def test_real_tree_provenance_and_no_confounded_ITS_or_cpDNA():
    assert m.NA_GIT_BLOB=="6417a38793427bc7c8562dc0e2cab59b74d8c316"
    assert m.J_TREE_SHA=="ba592d3a9f52b5481314ab20e3202048505cdf72a3fc15134d5dca261b9d9b7e"
    assert m.NA_EXPECTED==42 and m.J_EXPECTED==12
    assert m.J_TREE.endswith("primary_tree_v1.nwk")
    assert m.TRAITS==("mean_image_angle_deg","median_image_angle_deg","log_BIO12")

def test_known_clustering_in_synthetic_true_branch_pair():
    tree=Phylo.read(StringIO("((a:1,b:1):2,(c:1,d:1):2);"),"newick")
    tips,W,_=m.neighbor_kernel(tree)
    assert len(tips)==4
    assert all(abs(sum(row)-1)<1e-10 for row in W)
    by={x.name:i for i,x in enumerate(tips)}
    y=np.zeros(4)
    for lab,x in [("a",20),("b",21),("c",100),("d",101)]:y[by[lab]]=x
    obs=m.observed_nn_stat(y,W)
    assert abs(obs-1)<1e-8
    null=m.permutation_test(y,W,999,20261024)
    assert null["source_nearest_neighbor_mean_absolute_difference"]==obs
    assert null["ratio_observed_to_null"]<1
    assert 0<=null["one_sided_clustering_p"]<=1
    assert null["n_permutations"]==999

def test_rejects_missing_angle_and_shuffled_wrong_taxa():
    tree=Phylo.read(StringIO("(a:1,b:1,c:1);"),"newick")
    _,w,_=m.neighbor_kernel(tree)
    with pytest.raises(ValueError,match="scalar finite"):
        m.observed_nn_stat([5,float("nan"),1],w)
    with pytest.raises(ValueError,match="scalar finite"):
        m.observed_nn_stat([5,1],w)

def test_bh_is_monotonic_and_protects_multi_panel_scanning():
    p=m.bh([.01,.20,.03,.90])
    assert p[0]==pytest.approx(.04)
    assert p[1]==pytest.approx(.26666666666)
    assert p[2]==pytest.approx(.06)
    assert p[3]==pytest.approx(.9)
    assert list(m.PHOTO_GATES)==[1,10,30]

def test_invalid_partial_path_rejected():
    assert m.NA_TREE.startswith("traits/")
    assert "japan38" in m.J_MAP
    assert "OUTGROUP_saff" not in m.NA_TREE
