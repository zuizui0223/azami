"""S1 continuous image-angle ancestral reconstruction without binary posture coding."""
from io import StringIO
import numpy as np
import pytest
from Bio import Phylo
from analysis.v4 import s1_continuous_phylogenetic_history as s1h
from analysis.v4.native_phylogeny import tree_covariance


def tree_fixture():
    return Phylo.read(StringIO(
      "(((a:0.5,b:0.5):0.2,(c:0.4,d:0.4):0.3):0.5,"
      "((e:0.4,f:0.4):0.3,((g:0.2,h:0.2):0.2,i:0.4):0.3):0.5);"),
      "newick")


def test_native_S1_tree_covariance_matches_Brownian_shared_path():
    tree=tree_fixture()
    tips=[t.name for t in tree.get_terminals()]
    raw,ordered,path,parent=s1h.shared_path_matrix(tree,tips)
    expected=tree_covariance(tree,tips)
    assert len(ordered)>=len(tips)
    assert np.allclose(raw/np.mean(np.diag(raw)),expected,atol=1e-12)
    assert np.all(np.linalg.eigvalsh(raw)>0)
    assert all(np.isfinite(x).all() for x in path.values())
    assert tree.root not in parent


def test_continuous_ancestry_is_degree_scale_equivariant_not_discrete():
    tree=tree_fixture()
    tips=[t.name for t in tree.get_terminals()]
    y=np.array([45.,60.,70.,82.,100.,90.,120.,115.,135.])
    fit=s1h.continuous_BM(tree,tips,y,nugget_fraction=.05)
    assert fit["n_tips"]==9
    assert fit["n_edges"]==len(list(tree.find_clades()))-1
    assert fit["n_internal_edges"]<fit["n_edges"]
    assert fit["Brownian_rate_image_degree_sq_per_branch_unit"]>0
    assert fit["all_edges"]["5"]["conditional_expected_edges_increase_by_at_least_deg"]>=0
    assert fit["all_edges"]["10"]["conditional_expected_edges_increase_by_at_least_deg"]<=fit["all_edges"]["5"]["conditional_expected_edges_increase_by_at_least_deg"]
    assert fit["internal_edges_only"]["10"]["conditional_expected_edges_increase_by_at_least_deg"]<=fit["all_edges"]["10"]["conditional_expected_edges_increase_by_at_least_deg"]
    same=s1h.continuous_BM(tree,tips,y+10.,nugget_fraction=.05)
    assert fit["fitted_root_image_axis_deg"]+10==pytest.approx(same["fitted_root_image_axis_deg"],abs=1e-6)
    assert fit["Brownian_rate_image_degree_sq_per_branch_unit"]==pytest.approx(same["Brownian_rate_image_degree_sq_per_branch_unit"])
    with pytest.raises(ValueError,match="Unsupported"):
        s1h.continuous_BM(tree,tips,y,nugget_fraction=-.1)


def test_GLs_five_original_exposures_direct_effects_present():
    tree=tree_fixture(); tips=[t.name for t in tree.get_terminals()]
    n=len(tips); rng=np.random.default_rng(1)
    climate=rng.normal(size=(n,5))*np.array([10,1200,1000,12,5])
    y=65 + 0.02*climate[:,2] + rng.normal(scale=1,size=n)
    report=s1h.equal_species_five_climate_GLS(tree,tips,y,climate)
    assert set(report["coefficients_image_axis_deg_per_increment"])==set(s1h.EXPOSURES)
    assert report["n_tips"]==9
    assert report["df"]==3
    assert report["not_observation_weighted_original_S1"] is True
    assert report["not_a_causes_or_origins_test"] is True


def test_no_branch_length_fabrication_without_audit():
    tree=tree_fixture()
    names=[t.name for t in tree.get_terminals()]
    tree.get_terminals()[0].branch_length=0
    reg,meta=s1h.audited_pruned_tree(tree,names)
    assert meta["zero_length_edges_regularised"]==1
    assert 0<reg.get_terminals()[0].branch_length<.0001
    tree.get_terminals()[0].branch_length=-.5
    with pytest.raises(ValueError,match="Missing/negative"):
        s1h.audited_pruned_tree(tree,names)


def test_do_not_call_ancestral_degree_shifts_nodding_origins():
    code=s1h.continuous_BM.__doc__
    assert "NOT independently originated nodding" in code or "not independently originated nodding" in code.lower()
    assert s1h.SHIFT_DEGREES==(5.,10.,20.)
    assert s1h.ENDPOINT=="orientation_image_vertical_angle"
