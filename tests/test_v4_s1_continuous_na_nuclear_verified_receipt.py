"""Trust gate: same photographed NA42 nuclear vs original grafted phylogeny."""
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/"analysis/v4/results/20261021_unified_abcde/S1_CONTINUOUS_NA42_NUCLEAR_VS_GRAFTED_SIGNAL_VERIFIED_20261021.json"
REG=ROOT/"analysis/v4/RAIN_ORIENTATION_RESEARCH_DECISIONS_20261021.json"


def test_real_NA42_independent_tree_and_52_placements_same_photo_angle():
    r=json.loads(RECEIPT.read_text(encoding="utf-8"))
    s=r["result"]
    assert r["source_run_id"]==39596484704
    assert r["output_artifact"]["id"]==12344505559
    assert r["output_artifact"]["sha256"].startswith("sha256:")
    assert s["status"]=="REAL_S1_SAME_NA42_CONTINUOUS_ANGLE_INDEPENDENT_NUCLEAR_VS_52_GRAFTED_PHYLOGENETIC_SIGNAL"
    assert s["same_original_photo_phenotype_nuclear_vs_placements_n_taxa"]==42
    assert s["S1_five_original_covar"]==[
       "chelsa_bio01","chelsa_bio04","chelsa_bio12","chelsa_bio15","chelsa_rsds_mean"]
    assert s["published_independent_NA_nuclear"]["five_climate_adjusted"]["lambda_ML"]==0
    assert s["published_independent_NA_nuclear"]["unadjusted"]["lambda_ML"]==0
    p=s["matched_52_placement_phenotype_fiveclimate_lambda"]
    assert p["median"]==0
    assert p["min"]==0
    assert p["max"]==0.44419565891063845
    assert p["n_deltaAIC_over_2"]==0
    assert s["NA42_is_not_global_Cirsium_independent_phylogeny"] is True
    assert s["nodding_gravity_independent_evolution_or_rain_selection_uninferred"] is True


def test_nuclear_bootstrap_not_mislabeled_as_GBOTB_or_true_biological_origin():
    r=json.loads(RECEIPT.read_text(encoding="utf-8"))
    mk=r["result"]["same_NA42_nuclear_bootstrap"]
    for label in ("nuclear_unadjusted","nuclear_S1_five_climate_adjusted"):
        row=mk[label]
        assert row["observed_lambda_ML"]==0
        assert row["observed_lrt"]==0
        assert row["boundary_calibrated_bootstrap_p"]==1
        assert row["requested_draws"]==399
        assert "independently published author nuclear" in row["scope"]
        assert "NOT 52 V.PhyloMaker2 placements" in row["scope"]
    assert r["limitations"]["same_42_tips_use_S1_image_median_not_gravity"]
    assert r["limitations"]["not_global_complete_194_nuclear_tree"]
    assert r["limitations"]["no_inferred_biological_nodding_independent_origins"]


def test_decision_registry_preserves_continuous_S1_primary_and_external_nuclear_scope():
    reg=json.loads(REG.read_text(encoding="utf-8"))
    x=reg["S1_T2_primary_continuous_head_history"]
    y=x["independently_published_NA42_nuclear_same_S1_angle_signal_20261021"]
    assert x["original_S1_image_angle_used_as_continuous_primary"] is True
    assert x["native_S1_taxon_keys_with_measured_image_orientation"]==244
    assert x["historical_placement_matched_taxon_keys"]==194
    assert y["NA_original_S1_image_median_taxa"]==42
    assert y["independent_author_nuclear_Pagel_lambda_original_angle_fiveclimate_adjusted"]==0
    assert y["causal_historical_rain_selection_not_tested"] is True
    assert reg["T3"]["cross_tribe_convergence_inferred"] is False
