"""Model correctness for all-available S1 continuous phylogenetic history."""
from io import StringIO
from copy import deepcopy
import json
import numpy as np
import pytest
from Bio import Phylo
from analysis.v4 import s1_continuous_history as h

def make_tree():
    return Phylo.read(StringIO("((A:1,B:1):1,(C:0.5,D:0.5):1.5);"),"newick")

def make_tips():
    return {
      "A":np.array([42,0,1,2,3,4],float),
      "B":np.array([56,1,-1,3,5,2],float),
      "C":np.array([92,2,3,-1,5,1],float),
      "D":np.array([101,1,2,2,3,0],float),
    }

def test_pair_PIC_exact_and_root():
    tree=Phylo.read(StringIO("(A:1,B:1);"),"newick")
    raw={"A":np.array([10,3,2,1,4,2],float),
         "B":np.array([20,5,4,3,2,0],float)}
    info=h.ancestral_PIC(tree,raw)
    assert info["PIC"].shape==(1,6)
    assert abs(info["PIC"][0,0])==pytest.approx(10/np.sqrt(2))
    assert info["root"][0]==pytest.approx(15)
    assert info["root_variance_branch_units"]==pytest.approx(.5)
    assert info["original_zero_length_branches"]==0

def test_PIC_gls_equivalence_all_five_variables():
    """Independent contrasts must equal exact Brownian GLS with intercept."""
    rng=np.random.default_rng(144)
    nw="((A:1,B:1):1,((C:1,D:1):.8,(E:1,F:1):.8):.2);"
    tree=Phylo.read(StringIO(nw),"newick")
    tips=[t.name for t in tree.get_terminals()]
    arr=rng.normal(size=(len(tips),6))
    observed=dict(zip(tips,arr))
    info=h.ancestral_PIC(tree,observed)
    p=info["PIC"]
    beta=np.linalg.lstsq(p[:,1:],p[:,0],rcond=None)[0]
    # C_i,j = common shared Brownian root-to-MRCA length.
    nodes={t.name:t for t in tree.get_terminals()}
    dep=tree.depths()
    cov=np.array([[dep[tree.common_ancestor(nodes[a],nodes[b])]
                   for b in tips] for a in tips],float)
    X=np.column_stack((np.ones(len(tips)),arr[:,1:]))
    W=np.linalg.solve(cov,X)
    beta_gls=np.linalg.solve(X.T@W,X.T@np.linalg.solve(cov,arr[:,0]))[1:]
    assert beta==pytest.approx(beta_gls,abs=1e-8)
    assert info["n_independent_contrasts_within_scenario"]==5

def test_zero_branch_regularization_and_polytomy_are_explicit():
    tree=Phylo.read(StringIO("(A:0,B:0,C:1);"),"newick")
    raw={n:np.ones(6)*i for i,n in enumerate(("A","B","C"))}
    p=h.ancestral_PIC(tree,raw)
    assert p["PIC"].shape==(2,6)
    assert p["original_zero_length_branches"]==2
    assert np.isfinite(p["PIC"]).all()
    with pytest.raises(ValueError,match="tips and actual S1"):
        h.ancestral_PIC(tree,{"A":raw["A"]})

def test_identical_same_five_climates_and_mean_median_statistical_ends():
    rng=np.random.default_rng(2)
    t=rng.normal(size=(120,6))
    t[:,0]+=t[:,1]*3.0
    result=h.fit_PIC_multi(t)
    assert set(result["coefficients_per_1_cohort_species_SD_climate"])==set(h.EXPOSURES)
    assert len(result["within_scenario_five_variable_Holm_P"])==5
    assert result["df"]==115
    assert result["coefficients_per_1_cohort_species_SD_climate"]["chelsa_bio01"]>0
    assert result["not_an_evolutionary_transition_rate_or_origin_count"] is True
    with pytest.raises(ValueError,match="Five-exposure independent"):
        h.fit_PIC_multi(np.ones((50,6)))

def rows(n=60):
    rng=np.random.default_rng(22)
    out=[]
    for i in range(n):
        vals=rng.normal(size=5)
        vals*=np.array([4,60,200,10,2])
        vals+=np.array([15,1200,700,35,13])
        out.append({
          "accepted_key":str(i),"taxon_source_names":[f"Cirsium sp{i}"],
          "tree_tip":f"Cirsium_sp{i}","join_status":"exact_independent_source_name_tip_match",
          "S1_observations":50,
          "image_angle_mean_deg":float(60+4*vals[0]/15+0.01*vals[2]+rng.normal()),
          "image_angle_median_deg":float(61+4*vals[0]/15+0.01*vals[2]+rng.normal()),
          **dict(zip(h.EXPOSURES,vals)),
        })
    return out

def test_tree_scenario_all_taxa_primary_and_sampled_limit():
    records=rows(60)
    tree=Phylo.read(StringIO("("+",".join(f"{r['tree_tip']}:2" for r in records)+");"),"newick")
    result=h.analyze_scenario("synthetic_star",tree,records,min_obs=1)
    assert result["status"].startswith("conditional_Brownian")
    assert result["tree_tips_n"]==60
    assert result["PIC_contrasts_n"]==59
    assert "chelsa_bio15" in result["mean_image_angle"]["coefficients_per_1_cohort_species_SD_climate"]
    assert result["no_binary_downward_origin_count"] is True
    records[:30]=[{**x,"S1_observations":1} for x in records[:30]]
    status=h.analyze_scenario("synthetic",tree,records,min_obs=5)
    assert status["status"]=="not_enough_high_observation_taxa"

def test_main_roster_denominator_and_frozen_chain_are_never_synthetic():
    import pathlib
    receipt=json.loads((pathlib.Path("analysis/v4/results/20261021_unified_abcde")/
         "S1_T2_JOINT_244_194_VERIFIED_20261021.json").read_text())
    assert receipt["results"]["model_scope"]["orientation_eligible_accepted_keys"]==244
    assert receipt["results"]["model_scope"]["orientation_eligible_observations"]==170066
    assert receipt["results"]["model_scope"]["placement_tree_covered_keys"]==194
    assert h.TRAITS[0]=="image_angle_mean_deg"
    assert h.EXTRA_TRAIT=="image_angle_median_deg"
