"""Cross-tree, same-photo-angle continuous Pagel lambda tests (not binary U/D)."""
from io import StringIO

from Bio import Phylo
import numpy as np
import pytest

from analysis.v4 import s1_continuous_nuclear_vs_grafted_signal as mod
from analysis.v4.joint_s1_t2_phylogeny import EXPOSURES


def fake_tree():
    a=[f"Cirsium_t{i:02d}:0.4" for i in range(31)]
    b=[f"Cirsium_t{i:02d}:0.4" for i in range(31,63)]
    return Phylo.read(StringIO("(("+",".join(a)+"):1.0,("+
             ",".join(b)+"):1.0);"),"newick")


def fake_s1():
    rng=np.random.default_rng(91)
    rows=[]
    for i in range(42):
        climate=rng.normal(size=5)
        angle=68+2*climate[2]+rng.normal(scale=3)+int(i>=31)*12
        rows.append({
          "tip":f"Cirsium_t{i:02d}", "accepted_key":f"source_{i}",
          "source_concept_collision":False,
          "n_eligible_images":15 if i<30 else 7,
          "image_mean_deg":float(angle),"image_median_deg":float(angle+1),
          "S1_five_climate_means":dict(zip(EXPOSURES,climate.tolist())),
        })
    return rows


def test_genuinely_published_tree_and_S1_must_have_identical_tips():
    t=fake_tree()
    rows=fake_s1()
    matched=mod.match_original_S1_to_nuclear(rows,t)
    assert len(matched)==42
    rows[0]["source_concept_collision"]=True
    with pytest.raises(ValueError,match="collision"):
        mod.match_original_S1_to_nuclear(rows,t)


def test_source_matched_nuclear_lambda_and_52_placements_one_shared_phenotype():
    nuclear=fake_tree()
    rows=fake_s1()
    scenarios=[(f"placement_{i}",fake_tree()) for i in range(52)]
    result=mod.compare(rows,nuclear,scenarios,bootstrap_reps=20)
    status=result["summary"]
    assert status["same_NA_S1_image_measured_taxa"]==42
    assert status["global_S1_phylogeny_covered_continuous_taxa"]==194
    assert status["no_true_nodding_origin_count"] is True
    assert status["same_tip_comparison"] is True
    assert status["cohort_image_count_sizes"]=={"all42":42,"at_least5":42,"at_least10":30}
    assert status["n_historical_GBOTB_grafts"]==52
    for cohort,panel in status["nuclear_v_grafted_panels"].items():
        key="species_mean|all_five_original_S1_climates|nugget0.05"
        assert panel[key]["n_tip_cohort"] in (42,30)
        assert panel[key]["grafted_52"]["n_scenarios"]==52
        assert 0<=panel[key]["nuclear"]["lambda_ML"]<=1
        assert panel[key]["nuclear"]["n_taxa"]==panel[key]["n_tip_cohort"]
    assert len(result["full_53_tree_fit_details"])==len(status["cohort_image_count_sizes"])*53*8
    assert 0<status["nuclear_original_fixed_tree_lambda0_parametric_bootstrap"]["all42"][
       "all_five_original_S1_climates__species_mean__nugget0.05"][
          "boundary_calibrated_bootstrap_p"]<=1


def test_nulls_do_not_create_nonexistent_species_or_confident_origin_numbers():
    nuclear=fake_tree()
    rows=fake_s1()
    rows[0]["tip"]="Cirsium_unplaced"
    with pytest.raises(ValueError,match="matched-tip count"):
        mod.match_original_S1_to_nuclear(rows,nuclear)
    rows=fake_s1()
    tree=fake_tree()
    tree.get_terminals()[0].branch_length=-2
    with pytest.raises(ValueError,match="Missing/negative"):
        mod.fit_one_tree(tree,rows,root_label="bad",
                         source_type="NA",bootstrap_only=False)
