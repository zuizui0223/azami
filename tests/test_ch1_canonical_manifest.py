from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_prepublication_documents_are_not_tracked() -> None:
    assert not (ROOT / "manuscript").exists()
    assert not (ROOT / "submission").exists()


def test_frozen_reproducibility_outputs_exist() -> None:
    required = (
        "analysis_outputs/v2_full27_environment_atlas_2026-08-27/v2_full27_environment_within.csv",
        "analysis_outputs/v2_full27_environment_atlas_2026-08-27/v2_full27_environment_among.csv",
        "analysis_outputs/v2_full27_environment_atlas_2026-08-27/v2_full27_environment_cross_scale.csv",
        "analysis_outputs/v2_full27_environment_atlas_2026-08-27/v2_full27_variance_decomposition.csv",
        "analysis_outputs/v2_full27_environment_atlas_2026-08-27/spatial/v2_full27_spatial_within.csv",
        "analysis_outputs/v2_full27_environment_atlas_2026-08-27/spatial/v2_full27_spatial_among.csv",
        "analysis_outputs/v2_full27_environment_atlas_2026-08-27/historical/v2_full27_historical_placement_summary.csv",
        "analysis_outputs/v2_full27_environment_atlas_2026-08-27/sampling/v2_full27_sampling_composition_summary.csv",
        "analysis_outputs/literature_orientation_external_validation_2026-09-01/literature_orientation_external_validation_report.json",
    )
    for relative in required:
        assert (ROOT / relative).is_file(), relative


def test_analysis_contracts_and_registries_exist() -> None:
    required = (
        "analysis/ch1/v2_full27_environment_atlas_contract.json",
        "analysis/ch1/v2_full27_sampling_composition_sensitivity_contract.json",
        "analysis/ch1/image_to_trait_automated_technical_audit_summary.json",
        "analysis/ch1/literature_orientation_external_validation_contract.json",
        "ch1_global/v2/ontology/ch1_continuous_trait_contract.csv",
        "ch1_global/v2/ontology/ch1_environment_predictor_registry.json",
        "reproducibility/contracts/cirsium_ch1_image_trait_dictionary.csv",
    )
    for relative in required:
        assert (ROOT / relative).is_file(), relative


def test_reproducibility_figures_are_materialized() -> None:
    stems = (
        "Figure_1_v2_measurement_pipeline",
        "Figure_2_v2_geographic_sampling_domain",
        "Figure_3_v2_taxon_mean_information_loss",
        "Figure_4_v2_within_among_scale_atlas",
        "Figure_5_v2_candidate_robustness",
        "Figure_S1_v2_image_to_trait_technical_audit",
        "Figure_S2_v2_image_to_trait_perturbation_audit",
        "Figure_S3_v2_endpoint_measurement_support",
        "Figure_S4_v2_sampling_composition_audit",
        "Figure_S5_v2_spatial_diagnostic_surface",
        "Figure_S6_v2_historical_placement_stability",
        "Figure_S7_v2_whole_capitulum_secondary_synthesis",
    )
    for stem in stems:
        for suffix in (".png", ".pdf"):
            assert (ROOT / "reproducibility/figures" / f"{stem}{suffix}").is_file(), f"{stem}{suffix}"


def test_technical_audit_registry_points_to_reproducibility_tree() -> None:
    payload = json.loads(
        (ROOT / "analysis/ch1/image_to_trait_automated_technical_audit_summary.json").read_text(
            encoding="utf-8"
        )
    )
    for relative in payload["repository_files"]:
        assert not relative.startswith("manuscript/")
        assert (ROOT / relative).is_file(), relative
