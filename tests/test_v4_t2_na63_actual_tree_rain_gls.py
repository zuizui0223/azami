"""No-photo actual nuclear phylogeny uncertainty and PGLS contracts."""
from io import StringIO

from Bio import Phylo
import numpy as np
import pytest

from analysis.v4 import t2_na63_actual_tree_rain_gls as m


TREE = (
    "(((A:1,B:1):1,(C:1,D:1):1):1,"
    "((E:1,F:1):1,((G:1,H:1):0.5,(I:1,J:1):0.5):0.5):1);"
)


def small_tree():
    return Phylo.read(StringIO(TREE),"newick")


def test_source_brownian_covariance_carries_shared_ancestry():
    tree=small_tree()
    tips=[z.name for z in tree.get_terminals()]
    C=m.brownian_shared_ancestry(tree,tips)
    pos={x:i for i,x in enumerate(tips)}
    assert C.shape==(10,10)
    assert np.allclose(C,np.transpose(C))
    assert np.allclose(np.diag(C),3)
    assert C[pos["A"],pos["B"]]>C[pos["A"],pos["C"]]
    assert C[pos["A"],pos["C"]]>C[pos["A"],pos["H"]]
    assert np.linalg.eigvalsh(C).min()>0


def test_gls_lambda_zero_is_same_original_tip_ols_not_photo_angle():
    tree=small_tree()
    tips=[x.name for x in tree.get_terminals()]
    C=m.brownian_shared_ancestry(tree,tips)
    nodding=np.array([0,1,0,1,1,0,0,0,1,0.])
    temperature=np.array([7,9,3,8,5,12,6,4,10,11.])
    elevation=np.array([400,200,800,300,500,100,900,600,550,700.])
    X=np.column_stack((np.ones(10),nodding,
            (temperature-temperature.mean())/temperature.std(),
            (elevation-elevation.mean())/elevation.std()))
    rain=np.array([700,1200,950,1250,1700,550,1000,850,1850,760.])
    actual=m.gls(rain,X,C,0)
    ols=np.linalg.lstsq(X,rain,rcond=None)[0]
    assert np.allclose(actual["beta_nodding_BIO12_mm"],ols[1])
    assert 0<=actual["beta_nodding_one_sided_positive_t_p"]<=1
    assert actual["degrees_of_freedom"]==6
    full=m.gls(rain,X,C,1)
    assert full["beta_nodding_ci95_mm"][0]<full["beta_nodding_ci95_mm"][1]
    assert full["n_original_nuclear_tips"]==10


def test_estimated_lambda_profile_stays_inside_phylogenetic_covariance_domain():
    t=small_tree()
    C=m.brownian_shared_ancestry(t,[x.name for x in t.get_terminals()])
    trait=np.array([0,1,0,1,1,0,0,0,1,0.])
    values=np.array([700,1200,950,1250,1700,550,1000,850,1850,760.])
    X=np.column_stack((np.ones(10),trait))
    reports=m.fit_models(values,X,C)
    assert set(reports)==set(m.ORDERED_MODELS)
    assert reports["no_phylogeny"]["pagel_lambda"]==0
    assert reports["Brownian_complete_phylogeny"]["pagel_lambda"]==1
    assert 0<=reports["estimated_lambda"]["pagel_lambda"]<=1
    assert reports["estimated_lambda"]["loglikelihood"]>=(
        max(reports["no_phylogeny"]["loglikelihood"],
            reports["Brownian_complete_phylogeny"]["loglikelihood"])-1e-6)


def test_true_tree_Mk_ancestral_posterior_has_sensitivity_to_root_assumption():
    t=small_tree()
    D={"A":0,"B":1,"C":0,"D":1,"E":1,"F":0,"G":0,"H":0,"I":1,"J":0}
    results={z:m.internal_marginals(t,D,root_prior=z)
             for z in ("erect","stationary")}
    assert results["erect"]["root_nodding_posterior"]==0
    assert 0<=results["stationary"]["root_nodding_posterior"]<=1
    for obj in results.values():
        assert obj["fitted_symmetric_Mk_rate"]>0
        assert len(obj["top_edge_endpoint_U_to_D_probability_candidates"])==12
        for tip in obj["top_edge_endpoint_U_to_D_probability_candidates"]:
            assert 0<=tip["parent_erect_child_nodding_endpoint_probability"]<=1
            assert not tip["independent_U_to_D_gain_count_on_branch_confirmed"]
            assert not tip["historical_rain_at_node_observed"]


def test_unidentified_independent_global_rain_evolution_is_explicit_stop_rule():
    assert not m.FLAG_STOPS["individual_origin_events_identified"]
    assert not m.FLAG_STOPS["precipitation_dependent_U_to_D_rate_fitted"]
    assert not m.FLAG_STOPS["photos_or_exif_angles_used"]
    assert not m.FLAG_STOPS["full_251_source_taxa_in_one_validated_molecular_phylogeny"]
    assert m.NA63==63 and m.N_NODDING==5
