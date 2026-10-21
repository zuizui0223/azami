"""Regression gate: the empirical PGLS52/source-nuclear result stays exploratory."""
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]
REL = "analysis/v4/results/20261021_unified_abcde/T2_PGLS_52_VS_NUCLEAR_VERIFIED_20261021.json"
REGISTRY = "analysis/v4/RAIN_ORIENTATION_RESEARCH_DECISIONS_20261021.json"


def test_source_verified_same_tip_tree_sensitivity_never_becomes_rain_adaptation():
    receipt = json.loads((ROOT / REL).read_text(encoding="utf-8"))
    reg = json.loads((ROOT / REGISTRY).read_text(encoding="utf-8"))
    result = receipt["empirical_aggregate_from_successful_actions_log"]
    nu = result["independent_nuclear_reference"]
    t2 = reg["T2"]["exploratory_real_tree_and_author_trait_evidence"]
    assert receipt["workflow_run_id"] == 39579878018
    assert receipt["status"] == "EXECUTED_EXPLORATORY_NOT_GLOBAL_CAUSAL_EVOLUTIONARY_INFERENCE"
    assert result["n_scenario_results"] == 52
    assert result["matched_taxa_n_range"] == [44, 44]
    assert result["matched_nodding_n_range"] == [3, 3]
    assert result["unknown_root_min_steps_scenario_range"] == [2, 3]
    assert result["erect_to_nodding_q025_median_q975_scenario_median_range"] == [3., 3.]
    assert result["author_evidence"]["n_taxa"] == 63
    assert result["author_evidence"]["n_nodding"] == 5
    assert len(result["author_evidence"]["source_sha256"]) == 64
    assert nu["common_tip_comparison"]["n_taxa"] == 44
    assert nu["common_tip_comparison"]["unknown_root_min_steps"] == 3
    assert nu["complete_published_63_tip_context"]["n_taxa"] == 63
    assert nu["complete_published_63_tip_context"]["unknown_root_min_steps"] == 4
    assert result["rain_temperature_evolutionary_test_executed"] is False
    assert result["empirical_transition_dates_estimated"] is False
    assert nu["not_a_confirmed_dated_origin_analysis"] is True
    assert result["not_independent_repeated_selection"] is True
    assert t2["pgls_trees_are_not_independent_molecular_species_trees"] is True
    assert t2["origin_dating_or_rain_temperature_model_completed"] is False
    assert reg["T2"]["existing_52_trees_sufficient_for_dated_origin_counts"] is False
    assert reg["T3"]["cross_tribe_convergence_inferred"] is False
