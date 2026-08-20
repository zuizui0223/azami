#!/usr/bin/env python3
"""Validate committed Supplement tables against the frozen Chapter 1 claim registry.

This is a consistency check only. It never regenerates coefficients or changes
multiplicity families; it prevents the Supplement from silently drifting away
from the accepted manuscript results.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--claims", type=Path, default=Path("manuscript/final_claims.json"))
    p.add_argument("--tables-dir", type=Path, default=Path("manuscript/supplement/tables"))
    p.add_argument("--supplement", type=Path, default=Path("manuscript/supplement/SUPPLEMENTARY_INFORMATION.md"))
    return p.parse_args()


def close(a: float, b: float, tol: float = 1e-6) -> bool:
    return math.isclose(float(a), float(b), rel_tol=0.0, abs_tol=tol)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(f"Supplement validation failed: {message}")


def main() -> None:
    a = parse_args()
    claims = json.loads(a.claims.read_text(encoding="utf-8"))
    tables = a.tables_dir
    checks: list[str] = []

    # S1 cohort accounting.
    s1 = pd.read_csv(tables / "S01_canonical_cohorts.csv")
    by = s1.set_index("cohort")
    datasets = claims["datasets"]
    expected = {
        "balanced_image_comparison_atlas": (datasets["balanced_image_comparison_atlas"]["n_observations"], datasets["balanced_image_comparison_atlas"]["n_taxa"]),
        "exhaustive_detector_positive": (datasets["exhaustive_detected_layer"]["n_observations"], datasets["exhaustive_detected_layer"]["n_input_taxa"]),
        "exhaustive_spatially_thinned_primary": (datasets["exhaustive_spatially_thinned_primary"]["n_spatially_thinned_observations"], datasets["exhaustive_spatially_thinned_primary"]["n_input_taxa"]),
        "high_resolution_involucre": (datasets["high_resolution_involucre"]["n_observations"], datasets["high_resolution_involucre"]["n_taxa"]),
        "independent_detector_audit": (datasets["independent_detector_audit"]["n_source_images"], datasets["independent_detector_audit"]["n_taxa"]),
    }
    for cohort, (n_obs, n_taxa) in expected.items():
        require(cohort in by.index, f"S1 missing cohort {cohort}")
        require(str(by.loc[cohort, "observations"]) == str(n_obs), f"S1 observation count mismatch for {cohort}")
        require(int(by.loc[cohort, "taxa"]) == int(n_taxa), f"S1 taxon count mismatch for {cohort}")
    require(int(by.loc["balanced_image_comparison_atlas", "heads"]) == datasets["balanced_image_comparison_atlas"]["n_heads"], "S1 atlas head count")
    checks.append("S1 cohort counts")

    # S3 nested variance.
    s3 = pd.read_csv(tables / "S03_nested_visible_variance.csv")
    require(len(s3) == 9, "S3 must contain nine primary endpoints")
    nv = claims["nested_visible_variance"]
    lo, hi = nv["within_assigned_species_fraction_range"]
    require(close(s3["below_taxon_fraction"].min(), lo), "S3 minimum below-taxon fraction")
    require(close(s3["below_taxon_fraction"].max(), hi), "S3 maximum below-taxon fraction")
    one_lo, one_hi = min(nv["one_head_per_photo_within_fraction_range"]), max(nv["one_head_per_photo_within_fraction_range"])
    require(close(s3["one_head_per_photo_below_taxon_fraction"].min(), one_lo), "S3 one-head minimum")
    require(close(s3["one_head_per_photo_below_taxon_fraction"].max(), one_hi), "S3 one-head maximum")
    checks.append("S3 nested variance")

    # S4 is deliberately a presentation table: 28 linear rows + four joint hue vectors.
    s4 = pd.read_csv(tables / "S04_full_primary_within_taxon_coefficients.csv")
    require(len(s4) == 32, "S4 must contain 32 presentation rows")
    linear = s4.loc[s4["endpoint_type"].eq("linear")]
    hue = s4.loc[s4["endpoint_type"].eq("circular_joint")]
    require(len(linear) == 28, "S4 must contain 28 non-circular linear rows")
    require(len(hue) == 4, "S4 must contain four joint hue-vector rows")
    linear_supported = linear["fdr_significant_0_05"].astype(str).str.strip().str.lower().isin({"true", "1", "yes"}).sum()
    require(int(linear_supported) == claims["within_species_inference_levels"]["exhaustive_spatially_thinned_primary_36_component_models"]["n_bh_fdr_non_circular_linear_rows"], "S4 non-circular BH-supported count")
    require(hue["q_fdr_bh"].isna().all(), "S4 joint hue rows must not invent scalar q values")
    require(hue[["hue_beta_sin", "hue_beta_cos"]].notna().all().all(), "S4 joint hue vectors incomplete")
    checks.append("S4 36-component-family presentation")

    # S5/S6 grouped SPDE.
    s5 = pd.read_csv(tables / "S05_spde_model_group_summary.csv")
    require(len(s5) == claims["grouped_spde"]["n_fits"], "S5 grouped SPDE fit count")
    require(int(s5["cpo_failures"].sum()) == claims["grouped_spde"]["cpo_failures"], "S5 CPO failures")
    require(s5["trait"].nunique() == 9 and s5["model_group"].nunique() == 4, "S5 endpoint/group dimensions")
    s6 = pd.read_csv(tables / "S06_spde_fixed_effects_full.csv")
    require(len(s6) > 300, "S6 full fixed-effect table unexpectedly small")
    s6a = pd.read_csv(tables / "S06a_spde_global_BH_supported_rows.csv")
    require(not s6a.empty, "S6a supported SPDE rows missing")
    checks.append("S5-S6 grouped SPDE")

    # S7 niche permutation and S8 historical sensitivity.
    s7 = pd.read_csv(tables / "S07_niche_permutation_summary.csv")
    normalized = set(s7["Endpoint"].str.lower().str.replace(" ", "_", regex=False))
    require({"orientation", "chroma", "hue_sine", "hue_cosine", "width-profile_cv"}.issubset(normalized), "S7 supported all-taxa niche traits")
    require(claims["niche_permutation"]["n_permutations"] == 10000, "claim registry permutation count")
    s8 = pd.read_csv(tables / "S08_randomized_PGLS_retained_rows.csv")
    require(len(s8) == len(claims["historical_sensitivity"]["randomized_tree_rows_supported_50_of_50"]), "S8 retained randomized-tree row count")
    require(s8["Trees with FDR support"].astype(str).eq("50 / 50").all(), "S8 must be 50/50 retained rows")
    checks.append("S7-S8 among-taxon/historical")

    # S9 spatial diagnostics.
    s9 = pd.read_csv(tables / "S09_spatial_robustness.csv")
    require(len(s9) == 9, "S9 must contain nine endpoints")
    moran_lo, moran_hi = claims["spatial_robustness"]["residual_moran_i_range"]
    require(s9["residual_morans_I"].min() >= moran_lo - 5e-5, "S9 Moran minimum outside frozen range")
    require(s9["residual_morans_I"].max() <= moran_hi + 5e-5, "S9 Moran maximum outside frozen range")
    require((s9["permutation_P"] >= 0.05).all(), "S9 unexpectedly contains P < 0.05")
    checks.append("S9 spatial robustness")

    # S10 taxonomy.
    s10 = pd.read_csv(tables / "S10_wcvp_synonym_candidates.csv")
    require(len(s10) == claims["taxonomic_robustness"]["wcvp_synonym_conflicts_collapsed"], "S10 WCVP synonym count")
    checks.append("S10 taxonomy")

    # S11 precision audit.
    s11 = pd.read_csv(tables / "S11_precision_confounding_audit_summary.csv")
    legacy = claims["legacy_precision_audit"]
    row = s11.loc[s11["comparison"].eq("legacy variation vs raw absolute-slope RMS")]
    require(len(row) == 1 and close(row.iloc[0]["estimate"], legacy["variation_vs_legacy_absolute_slope_rms_rho"]), "S11 legacy rho")
    row = s11.loc[s11["comparison"].eq("revised noise-adjusted axis")]
    require(len(row) == 1 and close(row.iloc[0]["estimate"], claims["revised_primary_results"]["noise_adjusted_axis"]["spearman_rho"]), "S11 revised rho")
    checks.append("S11 precision audit")

    # S12 gates and wording guard.
    s12 = pd.read_csv(tables / "S12_submission_gates.csv")
    pending = set(s12.loc[s12["status"].eq("pending_external"), "gate"])
    require(pending == {"Independent detector human boxes", "Independent orientation/colour/outline reference measurements"}, "S12 external blocker set")
    text = a.supplement.read_text(encoding="utf-8")
    require("Complete 36-row exhaustive within-taxon climate coefficient family" not in text, "Supplement must not describe 32-row joint-hue presentation as 36 rows")
    require("28 ordinary linear endpoint–predictor rows plus four joint hue-vector summaries" in text, "Supplement S4 joint-hue explanation missing")
    checks.append("S12 gates and wording")

    print(json.dumps({"status": "pass", "checks": checks}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
