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
    datasets = claims["datasets"]
    legacy = claims["legacy_precision_audit"]
    revised = claims["revised_primary_results"]
    inference = claims["within_species_inference_levels"]
    history = claims["historical_sensitivity"]
    molecular = claims["molecular_database_coverage"]

    assert datasets["balanced_image_comparison_atlas"]["n_taxa"] == 216
    assert datasets["balanced_image_comparison_atlas"]["n_heads"] == 6626
    assert datasets["exhaustive_detected_layer"]["n_observations"] == 406582
    assert datasets["exhaustive_spatially_thinned_primary"]["n_spatially_thinned_observations"] == 46276
    assert datasets["exhaustive_spatially_thinned_primary"]["n_input_taxa"] == 259
    assert datasets["revised_precision_aware_cohort"]["n_taxa"] == 101
    assert datasets["legacy_lability_cohort"]["n_taxa"] == 102

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
                "revised_precision_aware_taxa": datasets["revised_precision_aware_cohort"]["n_taxa"],
                "corrected_axis_rho": axis["spearman_rho"],
                "meta_regression_p": meta["likelihood_ratio_p_value"],
                "exhaustive_linear_fdr_rows": exhaustive["n_bh_fdr_non_circular_linear_rows"],
                "phylogenetic_fits": history["n_signal_fits"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
