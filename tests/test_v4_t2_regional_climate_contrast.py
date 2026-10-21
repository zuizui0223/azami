"""T2 five-exposure regional source-stratified diagnostic safety tests."""
from io import StringIO
import numpy as np
import pytest
from Bio import Phylo

from analysis.v4 import t2_regional_climate_contrast as mod

def synthetic_rows():
    rows=[]
    for g,source in ((0,mod.NA_SOURCE),(1,mod.EA_SOURCE)):
        for i in range(8):
            values={
              "chelsa_bio01": 8 + 1.5*i + 2*g,
              "chelsa_bio04": 670 + 12*i*i + 15*g,
              "chelsa_bio12": 410+30*i+13*g+33*(i%3),
              "chelsa_bio15": 11 + 6*(i%4) + .5*i + 3*g,
              "chelsa_rsds_mean": 8 + 2*np.sin(i)+.2*g
            }
            rows.append({
                "source":source,
                "source_taxon":f"Cirsium t{g}_{i}",
                "state":"nodding" if i in (1,6) else "erect",
                "s1_five_climate_centroids":values,
                "s1_image_orientation_median":float(i*4+70),
                "s1_observations":100+i,
                "join_status":"exact_S1_and_tree",
                "s1_phylogeny_194_key_match":True,
                "accepted_key":f"{g}-{i}"
            })
    return rows

def test_source_blocking_has_one_observed_map_and_no_causal_p():
    c=mod.cohort(synthetic_rows())
    assert c["n"]==16 and c["n_nodding"]==4
    assert c["X"].shape==(16,7)
    labels,meta=mod.source_blocked_state_maps(c)
    assert meta["maps"]==28*28
    assert meta["observed_assignment_in_universe"]==1
    assert np.all(labels.sum(axis=1)==4)
    base=mod.exact_scores(c,labels)
    assert base["p_value_status"]=="NOT_CALIBRATED_FOR_PHYLOGENETIC_DEPENDENCE"
    assert len(base["coefficient_is_linear_probability_change_per_cohort_SD"])==2
    for p,q in zip(base["combinatorial_reference_rank"],
                   base["two_axis_maxT_reference_rank"]):
        assert 0 <= p <= q <= 1

def test_same_cohort_linear_probability_identity_and_brownian_nugget():
    c=mod.cohort(synthetic_rows())
    labels,_=mod.source_blocked_state_maps(c)
    b=mod.exact_scores(c,labels)
    expected=np.linalg.lstsq(c["X"],c["y"],rcond=None)[0]
    idx=[c["climate_start"]+list(mod.EXPOSURES).index(v) for v in mod.TEST_AXES]
    assert b["coefficient_is_linear_probability_change_per_cohort_SD"]==pytest.approx(
        expected[idx].tolist(),abs=1e-10)
    names=c["tips"]
    tr=Phylo.read(StringIO("("+",".join(f"{name}:0.5" for name in names)+");"),"newick")
    a=mod.pruned_covariance(tr,names)
    assert np.allclose(a,np.eye(len(names)))
    ctree=mod.exact_scores(c,labels,covariance=a)
    assert ctree["coefficient_is_linear_probability_change_per_cohort_SD"]==pytest.approx(
        b["coefficient_is_linear_probability_change_per_cohort_SD"],abs=1e-9)
    with pytest.raises(ValueError,match="Invalid frozen"):
        mod.weights(c,np.ones((3,3)))

def test_single_NA_source_drops_duplicate_EA_indicator():
    c=mod.cohort([r for r in synthetic_rows() if r["source"]==mod.NA_SOURCE])
    assert c["X"].shape==(8,6)
    assert c["climate_start"]==1
    labels,meta=mod.source_blocked_state_maps(c)
    assert meta["maps"]==28
    assert mod.exact_scores(c,labels)["n_nodding"]==2

def test_fail_closed_when_source_morphology_climate_or_tree_is_ambiguous():
    rows=synthetic_rows()
    rows[0]["s1_five_climate_centroids"].pop("chelsa_bio12")
    receipt={"workflow_run_id":39584744194,"workflow_artifact_digest":mod.RECEIPT_DIGEST,
             "results":{"botanical_state_inventory":{
                 "NA_exact_S1_194_tree_matches":42,
                 "EastAsia_exact_S1_194_tree_matches":4}},
             "source_provenance":{"native_sha256":"x","EAzami_fdt4_git_blob":"y",
                                  "NA_head_position_git_blob":"z"}}
    source={"status":"empirical_full_244_native_join_source_sha_verified",
            "fixed_covariates":list(mod.EXPOSURES),
            "model_scope":{"orientation_eligible_accepted_keys":244,
             "orientation_eligible_observations":170066,
             "placement_tree_covered_keys":194,
             "placement_tree_covered_observations":169912,
             "grafted_placement_scenarios":52},
            "provenance":{"S1_model_source_sha256":"x",
                          "EAzami_source_git_blob":"y",
                          "NA_source_git_blob":"z",
                          "native_decoded_source_authenticated":True},
            "T2_independent_state_inventory":{"taxon_source_audit":rows}}
    with pytest.raises(ValueError,match="coverage or taxonomy"):
        mod.validate_joint(source,receipt)
    with pytest.raises(KeyError):
        mod.cohort(rows)

def test_balanced_source_composition_not_mistaken_for_52_replicates():
    rows=synthetic_rows()
    c=mod.cohort(rows)
    source=c["source_counts"]
    assert source[mod.NA_SOURCE]["nodding"]==2
    assert source[mod.EA_SOURCE]["nodding"]==2
    assert set(mod.TEST_AXES)=={"chelsa_bio12","chelsa_bio15"}
