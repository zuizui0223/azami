#!/usr/bin/env python3
"""Build submission-facing Main and Supporting tables from frozen v2 outputs.

This script reshapes and joins already-frozen result tables.  It does not fit a
model, calculate a new P value, or select an additional hypothesis.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "analysis_outputs" / "v2_full27_environment_atlas_2026-08-27"
DEFAULT_OUTPUT = ROOT / "manuscript" / "tables" / "v2_submission"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_table(frame: pd.DataFrame, output: Path, filename: str) -> Path:
    path = output / filename
    frame.to_csv(path, index=False, lineterminator="\n")
    return path


def build_s10_full_roadmap() -> pd.DataFrame:
    rows = [
        (
            "Photo localization",
            "Where is a visible capitulum in each public photograph?",
            "406,582 detector-positive observations in the exhaustive source stream",
            "Visible-capitulum bounding boxes",
            "YOLO localization only; deterministic functions score traits after cropping",
            "A box defines the image unit; independent detector validation remains an external gate",
            "Figure 1",
            "Tables S1-S2",
        ),
        (
            "Cohort and realized domain",
            "Which coordinate-bearing records enter the spatial atlas?",
            "Detector-positive stream filtered by coordinate, positional-accuracy and taxon-by-cell rules",
            "46,276 observations; 259 source-assigned taxa",
            "Maximum 10-km positional uncertainty; one observation per taxon by 0.25-degree cell",
            "Freeze the realized sampling domain before examining endpoint values",
            "Figure 2",
            "Table S2",
        ),
        (
            "Continuous measurement contract",
            "Which visible capitulum dimensions are numerically recoverable?",
            "Tight and context crops from the frozen detector-positive stream",
            "27 registered endpoints; 22 measured; five explicit unexecuted",
            "Deterministic orientation, visible-colour, outline and fine-geometry measurements",
            "Missing or unassessable measurement remains missing evidence",
            "Figure 1",
            "Table S1",
        ),
        (
            "Taxon-mean information loss",
            "How much repeated image-phenotype variation is compressed by one taxon mean?",
            "All finite measurements in the frozen 46,276-observation cohort",
            "22 measured endpoints",
            "Raw sums-of-squares partition plus 500 two-observations-per-taxon resamples",
            "Describe visible image-phenotype information loss, not genetic variance",
            "Figure 3",
            "Table S3",
        ),
        (
            "Within/among environmental atlas",
            "Does environmental alignment occur within taxa, among taxa, or both?",
            "Frozen cohort for within-taxon models; taxon medians with at least five observations for primary among-taxon models",
            "26 inferential units by nine predictors in six blocks; 234 canonical rows",
            "Taxon-demeaned clustered models and taxon-label permutations; separate global BH families",
            "Classify every row as both, within-only, among-only, neither or not-comparable",
            "Figure 4",
            "Tables S4-S6",
        ),
        (
            "Sampling-composition audit",
            "Are selected directions stable to declared taxon, region, native-range and weighting perturbations?",
            "26 within-taxon and ten among-taxon globally supported rows",
            "674 evaluable retrospective scenarios",
            "Direction and effect-magnitude retention only; no post-selection P values",
            "Annotate sensitivity without converting the sample into a probability sample",
            "Figure 5",
            "Table S7",
        ),
        (
            "Broad-space gate",
            "Does support persist after a second-order spherical-coordinate basis and residual-spatial check?",
            "The same 36 globally supported rows",
            "26 within-taxon and ten among-taxon rows",
            "Spatial-basis refit, permutation P and eight-nearest-neighbour residual Moran diagnostic",
            "Five within and two among rows pass the declared broad-space screen",
            "Figure 5",
            "Table S8",
        ),
        (
            "Historical-placement gate",
            "Are the two broad-space-retained among-taxon directions stable to audited placement uncertainty?",
            "Two among-taxon rows retained by the broad-space gate",
            "Two deterministic and 50 randomized within-genus placement trees",
            "Pagel-lambda PGLS sensitivity on 52 trees per row",
            "Two adaptive-pattern candidates under current controls; not adaptation, mechanism or convergence proof",
            "Figure 5",
            "Table S8",
        ),
    ]
    return pd.DataFrame(
        rows,
        columns=[
            "analytical_stage",
            "scientific_question",
            "input_cohort",
            "endpoint_predictor_scope",
            "statistical_operation",
            "decision_or_gate",
            "main_figure",
            "supplementary_table",
        ],
    )


def build_table_1() -> pd.DataFrame:
    """Compact portrait version of the full analytical roadmap."""
    full = build_s10_full_roadmap()
    return pd.DataFrame({
        "analytical_stage": full["analytical_stage"],
        "frozen_scope": full["endpoint_predictor_scope"],
        "operation_and_gate": full["statistical_operation"] + ". Gate: " + full["decision_or_gate"],
        "linked_display": full["main_figure"] + "; " + full["supplementary_table"],
    })


def build_table_2(
    spatial_among: pd.DataFrame,
    sampling: pd.DataFrame,
    historical: pd.DataFrame,
) -> pd.DataFrame:
    """Expose the exact frozen evidence chain for the two candidate rows."""
    ratio_columns = [
        "dominant_taxon_omission_minimum_effect_magnitude_ratio",
        "leave_one_broad_region_out_minimum_effect_magnitude_ratio",
        "native_only_minimum_effect_magnitude_ratio",
        "equal_taxon_weight_minimum_effect_magnitude_ratio",
    ]
    order = [
        ("corolla_lab_chroma", "chelsa_rsds_mean", "Lower camera-recorded corolla chroma with higher shortwave radiation"),
        ("orientation_image_vertical_angle", "chelsa_bio12", "Larger image-referenced orientation angle with higher annual precipitation"),
    ]
    rows: list[dict[str, object]] = []
    for unit_id, predictor, label in order:
        spatial_row = spatial_among.loc[
            (spatial_among["unit_id"] == unit_id)
            & (spatial_among["predictor"] == predictor)
            & spatial_among["broad_spatial_sensitivity_pass"].astype(str).str.lower().eq("true")
        ].iloc[0]
        sampling_row = sampling.loc[
            sampling["scale"].eq("among_taxon")
            & sampling["unit_id"].eq(unit_id)
            & sampling["predictor"].eq(predictor)
        ].iloc[0]
        historical_row = historical.loc[
            historical["unit_id"].eq(unit_id) & historical["predictor"].eq(predictor)
        ].iloc[0]
        ratios = pd.to_numeric(sampling_row[ratio_columns], errors="coerce").dropna()
        rows.append({
            "candidate_pattern": label,
            "global_among_taxon": (
                f"n = {int(spatial_row['n_taxa'])}; beta = {float(spatial_row['base_beta_std']):+.3f}; "
                f"global BH q = {float(spatial_row['base_q_fdr_bh_global_family']):.3f}"
            ),
            "sampling_composition": (
                f"all declared directions stable; minimum effect-magnitude ratio = {float(ratios.min()):.3f}"
            ),
            "broad_space_gate": (
                f"beta = {float(spatial_row['spatial_beta_std']):+.3f}; permutation P = {float(spatial_row['spatial_permutation_p_value']):.3f}; "
                f"residual Moran P = {float(spatial_row['residual_morans_p_value']):.3f}"
            ),
            "historical_placement": (
                f"P < 0.05 on {int(historical_row['n_placement_trees_p_lt_0_05'])}/{int(historical_row['n_successful_placement_trees'])} trees; "
                f"P range {float(historical_row['minimum_p_value']):.2e}-{float(historical_row['maximum_p_value']):.2e}; "
                f"lambda range {float(historical_row['lambda_min']):.3f}-{float(historical_row['lambda_max']):.3f}"
            ),
        })
    return pd.DataFrame(rows)


def build_s1(inventory: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "endpoint_id",
        "module",
        "biological_construct",
        "unit",
        "measurement_status",
        "validation_status",
        "n_observations_measured",
        "n_taxa_measured",
        "interpretation",
        "not_interpretable_as",
    ]
    return inventory.loc[:, columns].copy()


def build_s2() -> pd.DataFrame:
    return pd.DataFrame(
        [
            ("detector_positive", 406582, 286, "at least one YOLO-localized visible capitulum", "manuscript/02_methods.md; manuscript/final_claims.json"),
            ("coordinate_usable", 392989, 271, "coordinate usable for environmental attachment", "manuscript/02_methods.md; ch1_global/v2/CONTINUOUS_TRAIT_REDESIGN.md"),
            ("positional_accuracy_le_10km", 297293, 259, "reported positional uncertainty at most 10 km", "manuscript/02_methods.md; ch1_global/v2/CONTINUOUS_TRAIT_REDESIGN.md"),
            ("taxon_by_0.25_degree_cell_thinning", 46276, 259, "one observation per source-assigned taxon by 0.25-degree cell", "manuscript/final_claims.json; v2_full27_environment_report.json"),
            ("final_environment_cohort", 46276, 259, "nine frozen predictors attached; minimum coverage gate passed", "v2_full27_environment_report.json; manuscript/final_claims.json"),
        ],
        columns=["filtering_stage", "n_observations", "n_source_assigned_taxa", "rule", "current_v2_source"],
    )


def build_s3(variance: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "endpoint_id",
        "status",
        "n_observations",
        "n_taxa",
        "n_taxa_with_at_least_two",
        "fraction_below_taxon_means",
        "equal_replication_fraction_median",
        "equal_replication_fraction_low_95",
        "equal_replication_fraction_high_95",
        "variance_grain",
        "interpretation_boundary",
    ]
    return variance.loc[:, columns].copy()


def build_s7(summary: pd.DataFrame) -> pd.DataFrame:
    families = [
        "dominant_taxon_omission",
        "leave_one_broad_region_out",
        "native_only",
        "equal_taxon_weight",
    ]
    rows: list[dict[str, object]] = []
    for _, row in summary.iterrows():
        for family in families:
            raw_n_scenarios = row[f"{family}_n_scenarios"]
            n_scenarios = 0 if pd.isna(raw_n_scenarios) else int(raw_n_scenarios)
            if n_scenarios == 0:
                continue
            rows.append(
                {
                    "scale": row["scale"],
                    "unit_id": row["unit_id"],
                    "predictor": row["predictor"],
                    "inferential_unit": row["inferential_unit"],
                    "scenario_family": family,
                    "n_scenarios": n_scenarios,
                    "n_evaluable": int(row[f"{family}_n_evaluable"]),
                    "all_evaluable": row[f"{family}_all_evaluable"],
                    "all_directions_stable_where_evaluable": row[f"{family}_all_directions_stable_where_evaluable"],
                    "minimum_effect_magnitude_ratio": row[f"{family}_minimum_effect_magnitude_ratio"],
                    "overall_stability_class": row["sampling_composition_stability_class"],
                }
            )
    return pd.DataFrame(rows)


def build_s8(
    spatial_within: pd.DataFrame,
    spatial_among: pd.DataFrame,
    sampling: pd.DataFrame,
    historical: pd.DataFrame,
) -> pd.DataFrame:
    spatial = pd.concat([spatial_within, spatial_among], ignore_index=True, sort=False)
    ratio_columns = [
        "dominant_taxon_omission_minimum_effect_magnitude_ratio",
        "leave_one_broad_region_out_minimum_effect_magnitude_ratio",
        "native_only_minimum_effect_magnitude_ratio",
        "equal_taxon_weight_minimum_effect_magnitude_ratio",
    ]
    sampling = sampling.copy()
    sampling["sampling_minimum_effect_magnitude_ratio"] = sampling[ratio_columns].apply(
        pd.to_numeric, errors="coerce"
    ).min(axis=1, skipna=True)
    sampling_columns = [
        "scale",
        "unit_id",
        "predictor",
        "all_declared_scenarios_evaluable",
        "all_directions_stable_where_evaluable",
        "sampling_minimum_effect_magnitude_ratio",
        "sampling_composition_stability_class",
    ]
    merged = spatial.merge(
        sampling[sampling_columns],
        on=["scale", "unit_id", "predictor"],
        how="left",
        validate="one_to_one",
    )
    historical_columns = [
        "unit_id",
        "predictor",
        "n_successful_placement_trees",
        "n_placement_trees_p_lt_0_05",
        "minimum_p_value",
        "maximum_p_value",
        "lambda_min",
        "lambda_max",
        "linear_direction_stable_across_trees",
        "historical_placement_sensitivity_pass",
        "candidate_class",
    ]
    merged = merged.merge(
        historical[historical_columns],
        on=["unit_id", "predictor"],
        how="left",
        validate="many_to_one",
    )
    merged["candidate_class"] = merged["candidate_class"].fillna("not_final_candidate")
    output_columns = [
        "scale",
        "unit_id",
        "member_endpoint_ids",
        "inferential_unit",
        "predictor",
        "environment_block",
        "base_beta_std",
        "base_q_fdr_bh_global_family",
        "sampling_composition_stability_class",
        "all_declared_scenarios_evaluable",
        "all_directions_stable_where_evaluable",
        "sampling_minimum_effect_magnitude_ratio",
        "status",
        "spatial_beta_std",
        "spatial_permutation_p_value",
        "same_linear_direction_as_base",
        "residual_morans_i",
        "residual_morans_p_value",
        "broad_spatial_sensitivity_pass",
        "n_successful_placement_trees",
        "n_placement_trees_p_lt_0_05",
        "minimum_p_value",
        "maximum_p_value",
        "lambda_min",
        "lambda_max",
        "linear_direction_stable_across_trees",
        "historical_placement_sensitivity_pass",
        "candidate_class",
    ]
    return merged.loc[:, output_columns].copy()


def build_s9(final_claims: dict[str, object]) -> pd.DataFrame:
    """Expose the frozen secondary-synthesis headline without recomputation."""
    secondary = final_claims["secondary_whole_capitulum"]
    assert isinstance(secondary, dict)
    common = {
        "synthesis": "expanded_complete18_whole_capitulum",
        "n_taxa": int(secondary["n_taxa"]),
        "n_endpoints": int(secondary["n_endpoints"]),
        "role": str(secondary["role"]),
    }
    return pd.DataFrame([
        {**common, "metric": "PC1_to_PC3_cumulative_variance", "estimate": float(secondary["pc1_to_pc3_cumulative_variance"]), "interpretation_boundary": "multidimensional phenotypic summary; not endpoint selection"},
        {**common, "metric": "within_among_association_spearman", "estimate": float(secondary["within_among_association_spearman"]), "interpretation_boundary": "partial geometry alignment; not shared causal mechanism"},
    ])


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    inventory = pd.read_csv(args.input_dir / "v2_full27_endpoint_inventory.csv")
    variance = pd.read_csv(args.input_dir / "v2_full27_variance_decomposition.csv")
    within = pd.read_csv(args.input_dir / "v2_full27_environment_within.csv")
    among = pd.read_csv(args.input_dir / "v2_full27_environment_among.csv")
    cross = pd.read_csv(args.input_dir / "v2_full27_environment_cross_scale.csv")
    sampling = pd.read_csv(args.input_dir / "sampling" / "v2_full27_sampling_composition_summary.csv")
    scenarios = pd.read_csv(args.input_dir / "sampling" / "v2_full27_sampling_composition_scenarios.csv")
    spatial_within = pd.read_csv(args.input_dir / "spatial" / "v2_full27_spatial_within.csv")
    spatial_among = pd.read_csv(args.input_dir / "spatial" / "v2_full27_spatial_among.csv")
    historical = pd.read_csv(args.input_dir / "historical" / "v2_full27_historical_placement_summary.csv")
    final_claims_path = ROOT / "manuscript" / "final_claims.json"
    final_claims = json.loads(final_claims_path.read_text(encoding="utf-8"))

    tables = {
        "Table_1_v2_analytical_roadmap.csv": build_table_1(),
        "Table_2_v2_candidate_evidence_chain.csv": build_table_2(spatial_among, sampling, historical),
        "Table_S1_v2_complete_27_endpoint_contract.csv": build_s1(inventory),
        "Table_S2_v2_cohort_geographic_filtering.csv": build_s2(),
        "Table_S3_v2_complete_variance_decomposition.csv": build_s3(variance),
        "Table_S4_v2_complete_within_taxon_environment_atlas.csv": within,
        "Table_S5_v2_complete_among_taxon_environment_atlas.csv": among,
        "Table_S6_v2_complete_cross_scale_classification.csv": cross,
        "Table_S7_v2_sampling_composition_robustness.csv": build_s7(sampling),
        "Table_S8_v2_spatial_historical_sequential_gate.csv": build_s8(
            spatial_within, spatial_among, sampling, historical
        ),
        "Table_S9_v2_secondary_whole_capitulum_summary.csv": build_s9(final_claims),
        "Table_S10_v2_full_analytical_roadmap.csv": build_s10_full_roadmap(),
    }

    assert len(inventory) == 27
    assert inventory["measurement_status"].value_counts().to_dict() == {
        "measured": 22,
        "unexecuted_no_measurement": 5,
    }
    assert len(cross) == 234
    assert cross["cross_scale_class"].value_counts().to_dict() == {
        "neither": 152,
        "not_comparable": 54,
        "within_only": 18,
        "both_scales": 7,
        "among_only": 3,
    }
    assert len(scenarios) == 674
    assert int(tables["Table_S7_v2_sampling_composition_robustness.csv"]["n_scenarios"].sum()) == 674
    assert len(tables["Table_S8_v2_spatial_historical_sequential_gate.csv"]) == 36
    assert (historical["n_successful_placement_trees"] == 52).all()
    assert (historical["n_placement_trees_p_lt_0_05"] == 52).all()
    assert tables["Table_S9_v2_secondary_whole_capitulum_summary.csv"]["n_taxa"].eq(127).all()
    assert tables["Table_S9_v2_secondary_whole_capitulum_summary.csv"]["n_endpoints"].eq(18).all()
    assert len(tables["Table_1_v2_analytical_roadmap.csv"]) == 8
    assert len(tables["Table_2_v2_candidate_evidence_chain.csv"]) == 2
    assert len(tables["Table_S10_v2_full_analytical_roadmap.csv"]) == 8

    output_paths = [write_table(frame, args.output_dir, name) for name, frame in tables.items()]
    input_paths = [
        args.input_dir / "v2_full27_endpoint_inventory.csv",
        args.input_dir / "v2_full27_variance_decomposition.csv",
        args.input_dir / "v2_full27_environment_within.csv",
        args.input_dir / "v2_full27_environment_among.csv",
        args.input_dir / "v2_full27_environment_cross_scale.csv",
        args.input_dir / "sampling" / "v2_full27_sampling_composition_summary.csv",
        args.input_dir / "sampling" / "v2_full27_sampling_composition_scenarios.csv",
        args.input_dir / "spatial" / "v2_full27_spatial_within.csv",
        args.input_dir / "spatial" / "v2_full27_spatial_among.csv",
        args.input_dir / "historical" / "v2_full27_historical_placement_summary.csv",
        final_claims_path,
    ]
    report = {
        "status": "PASS",
        "source_lane": "full27_full_environment_only",
        "operation": "reshape_and_join_frozen_outputs_only",
        "new_statistics_fitted": False,
        "table_rows": {path.name: len(tables[path.name]) for path in output_paths},
        "input_sha256": {path.relative_to(ROOT).as_posix(): sha256(path) for path in input_paths},
        "output_sha256": {path.name: sha256(path) for path in output_paths},
        "claim_boundary": "submission tables for traceability; not a new analysis or hypothesis screen",
    }
    (args.output_dir / "table_build_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
