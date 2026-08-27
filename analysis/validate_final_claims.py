#!/usr/bin/env python3
"""Validate the reviewer-revised Chapter 1 manuscript claim registry."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--claims", default="manuscript/final_claims.json")
    return parser.parse_args()


def main() -> None:
    claims = json.loads(Path(parse_args().claims).read_text(encoding="utf-8"))
    assert claims["release_label"] == "Azami Chapter 1 v2"
    assert claims["freeze_tag"] == "azami-ch1-v2-2026-08-27"
    assert Path(claims["defense_document"]).is_file()
    datasets = claims["datasets"]
    nested = claims["nested_visible_variance"]
    legacy = claims["legacy_precision_audit"]
    revised = claims["revised_primary_results"]
    inference = claims["within_species_inference_levels"]
    bias_control = claims["bias_control_reanalysis"]
    involucre = claims["auxiliary_involucre"]
    history = claims["historical_sensitivity"]
    molecular = claims["molecular_database_coverage"]
    full27 = claims["canonical_full27_full_environment"]
    sampling = full27["sampling_composition_sensitivity"]

    atlas = datasets["balanced_image_comparison_atlas"]
    assert atlas["n_taxa"] == 216
    assert atlas["n_observations"] == 3725
    assert atlas["n_photos"] == 3725
    assert atlas["n_heads"] == 6626
    assert atlas["one_photo_per_observation"] is True
    assert datasets["exhaustive_detected_layer"]["n_observations"] == 406582
    assert datasets["exhaustive_spatially_thinned_primary"]["n_spatially_thinned_observations"] == 46276
    assert datasets["exhaustive_spatially_thinned_primary"]["n_input_taxa"] == 259
    assert datasets["revised_precision_aware_cohort"]["n_taxa"] == 101
    assert datasets["legacy_lability_cohort"]["n_taxa"] == 102

    assert full27["n_registered_endpoints"] == 27
    assert full27["n_measured_endpoints"] == 22
    assert full27["n_spatial_observations"] == 46276
    assert full27["n_source_assigned_taxa"] == 259
    assert full27["within_global_q_supported_rows"] == 26
    assert full27["among_min5_global_q_supported_rows"] == 10
    candidates = full27["adaptive_pattern_candidates_under_current_controls"]
    assert [row["endpoint"] for row in candidates] == [
        "corolla_lab_chroma",
        "orientation_image_vertical_angle",
    ]
    assert all(row["n_historical_placement_trees_passed"] == 52 for row in candidates)
    assert "not anthocyanin concentration" in candidates[0]["mechanism_boundary"]
    assert "not gravity" in candidates[1]["mechanism_boundary"]
    assert all(row["decisive_next_test"] for row in candidates)
    assert sampling["status"] == "retrospective_selected_row_sensitivity"
    assert sampling["n_entry_pairs"] == 36
    assert sampling["n_scenario_rows"] == 674
    assert sampling["n_within_pairs_stable_all_declared_scenarios"] == 16
    assert sampling["n_within_pairs_direction_unstable"] == 10
    assert sampling["n_among_pairs_stable_all_declared_scenarios"] == 10
    assert sampling["n_broad_space_passes_stable_all_declared_scenarios"] == 6
    assert sampling["n_broad_space_passes_direction_unstable"] == 1
    assert sampling["direction_unstable_broad_space_passes"] == [
        "bract_projection_roughness~chelsa_rsds_mean"
    ]
    assert sampling[
        "both_adaptive_pattern_candidates_direction_stable_all_declared_scenarios"
    ] is True
    assert full27["validation"]["sampling_composition_checks_passed"] == 7
    assert full27["validation"]["sampling_composition_checks_total"] == 7

    assert nested["n_endpoints"] == 9
    for key in (
        "within_assigned_species_fraction_range",
        "among_photographs_within_species_fraction_range",
        "among_heads_within_photo_fraction_range",
        "one_head_per_photo_within_fraction_range",
        "balanced_10_photo_median_within_fraction_range",
    ):
        low, high = nested[key]
        assert 0 <= low <= high <= 1
    assert nested["within_assigned_species_fraction_range"][0] > 0.5
    assert nested["one_head_per_photo_within_fraction_range"][0] > 0.5
    assert nested["balanced_10_photo_median_within_fraction_range"][0] > 0.5
    assert nested["species_cluster_bootstrap_repeats"] == 2000
    assert nested["balanced_subsample_repeats"] == 500
    assert "not genetic" in nested["allowed_claim"].lower()

    assert abs(legacy["legacy_score_vs_median_slope_n_rho"]) > 0.9
    assert abs(legacy["legacy_score_vs_median_slope_se_rho"]) > 0.9
    assert legacy["partial_p_value"] > 0.05
    assert "withdraw" in legacy["decision"].lower()

    axis = revised["noise_adjusted_axis"]
    assert -1 <= axis["spearman_rho"] <= 1
    assert axis["p_value"] > 0.05
    assert axis["species_bootstrap_ci95"][0] < 0 < axis["species_bootstrap_ci95"][1]
    assert abs(axis["association_energy_vs_median_slope_n_rho"]) < 0.1
    assert abs(axis["association_energy_vs_median_slope_se_rho"]) < 0.1

    meta = revised["hierarchical_variance_meta_regression"]
    assert meta["n_slope_estimates"] == 2828
    assert meta["n_trait_predictor_groups"] == 28
    assert meta["profile_ci95"][0] < 0 < meta["profile_ci95"][1]
    assert meta["likelihood_ratio_p_value"] > 0.05

    exhaustive = inference["exhaustive_spatially_thinned_primary_36_component_models"]
    balanced = inference["balanced_atlas_le_10km_sensitivity"]
    assert exhaustive["n_bh_fdr_component_rows"] == 8
    assert exhaustive["n_bh_fdr_non_circular_linear_rows"] == 4
    assert exhaustive["n_bh_fdr_hue_component_rows"] == 4
    assert balanced["n_bh_fdr_main_rows"] == 0
    assert "different datasets" in inference["interpretation_rule"].lower()

    primary_control = bias_control["primary_season_and_dominant_taxon"]
    assert primary_control["n_observations"] == 46276
    assert primary_control["n_non_circular_rows_strictly_retained"] == 4
    assert primary_control["n_non_circular_rows_tested"] == 4
    assert len(primary_control["rows"]) == 4
    assert all(
        row["omission_beta_range"][0] <= row["season_adjusted_beta"] <= row["omission_beta_range"][1]
        for row in primary_control["rows"]
    )
    hemisphere_control = bias_control["hemisphere_season_sensitivity"]
    assert hemisphere_control["n_southern_observations"] == 2356
    assert hemisphere_control["n_taxa_sampled_both_hemispheres"] == 3
    assert hemisphere_control["n_non_circular_rows_strictly_retained"] == 4
    assert len(hemisphere_control["rows"]) == 4
    assert all(
        row["phase_aligned_omission_beta_range"][0]
        <= row["phase_aligned_beta"]
        <= row["phase_aligned_omission_beta_range"][1]
        and row["hemisphere_curve_omission_beta_range"][0]
        <= row["hemisphere_curve_beta"]
        <= row["hemisphere_curve_omission_beta_range"][1]
        for row in hemisphere_control["rows"]
    )
    native_control = bias_control["native_range_sensitivity"]
    assert native_control["n_observations"] == 46276
    assert native_control["n_taxa"] == 259
    assert native_control["n_resolved_taxa"] == 245
    assert sum(native_control["native_status_counts"].values()) == 46276
    assert native_control["native_status_counts"]["native"] == 27066
    assert native_control["n_non_circular_rows_tested"] == 4
    assert native_control["n_non_circular_rows_strictly_retained"] == 2
    assert len(native_control["rows"]) == 4
    retained_native_pairs = {
        (row["endpoint"], row["predictor"])
        for row in native_control["rows"]
        if row["native_range_robust"]
    }
    assert retained_native_pairs == {
        ("orientation_angle", "BIO1"),
        ("shape_aspect_ratio", "BIO4"),
    }
    assert claims["current_native_range_withdrawals"] == [
        "corolla_chroma-BIO12",
        "shape_aspect_ratio-BIO12",
    ]
    assert bias_control["repeat_photo_cohort"]["n_repeat_observations"] == 20073
    assert bias_control["repeat_photo_cohort"]["n_additional_photos"] == 38675

    assert involucre["n_strictly_retained_rows"] == 0
    assert "withdrawn" in involucre["current_status"]
    assert len(involucre["withdrawn_original_unadjusted_rows"]) == 3

    assert history["n_input_taxa"] == 216
    assert history["n_direct_backbone_taxa"] == 54
    assert history["n_random_trees"] == 50
    assert history["n_signal_fits"] == 636
    assert history["n_failed_signal_fits"] == 0

    for count in molecular.values():
        assert isinstance(count, int) and 0 <= count <= 216
    assert molecular["complete_plastome"] < molecular["plastid"]

    print(
        json.dumps(
            {
                "status": "valid",
                "release": claims["analysis_release"],
                "nested_within_fraction_range": nested["within_assigned_species_fraction_range"],
                "revised_precision_aware_taxa": datasets["revised_precision_aware_cohort"]["n_taxa"],
                "corrected_axis_rho": axis["spearman_rho"],
                "meta_regression_p": meta["likelihood_ratio_p_value"],
                "exhaustive_linear_fdr_rows": exhaustive["n_bh_fdr_non_circular_linear_rows"],
                "bias_control_linear_rows_retained": primary_control["n_non_circular_rows_strictly_retained"],
                "hemisphere_control_linear_rows_retained": hemisphere_control["n_non_circular_rows_strictly_retained"],
                "native_range_linear_rows_retained": native_control["n_non_circular_rows_strictly_retained"],
                "involucre_rows_retained": involucre["n_strictly_retained_rows"],
                "phylogenetic_fits": history["n_signal_fits"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
