#!/usr/bin/env python3
"""Build the GEB v2 multiscale environmental evidence atlas.

The atlas preserves the established manuscript flow: global continuous phenotype
pattern first, then within-taxon environment, spatial context, among-taxon
sorting and historical sensitivity. PR69 controls grade robustness instead of
acting as universal vetoes.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

PRIMARY_NAME_MAP = {
    "orientation_angle": "orientation_image_vertical_angle",
    "corolla_chroma": "corolla_lab_chroma",
    "hue_sin": "corolla_hue",
    "hue_cos": "corolla_hue",
    "shape_aspect_ratio": "capitulum_outline_aspect_ratio",
}
PREDICTOR_MAP = {
    "BIO1": "chelsa_bio01",
    "BIO4": "chelsa_bio04",
    "BIO12": "chelsa_bio12",
    "BIO15": "chelsa_bio15",
    "BIO1 annual mean temperature": "chelsa_bio01",
    "BIO4 temperature seasonality": "chelsa_bio04",
    "BIO12 annual precipitation": "chelsa_bio12",
    "BIO15 precipitation seasonality": "chelsa_bio15",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--contract", required=True)
    p.add_argument("--continuous-coefficients", required=True)
    p.add_argument("--candidate-quality", required=True)
    p.add_argument("--primary-supported", required=True)
    p.add_argument("--season", required=True)
    p.add_argument("--hemisphere", required=True)
    p.add_argument("--native", required=True)
    p.add_argument("--spde", required=True)
    p.add_argument("--niche", required=True)
    p.add_argument("--pgls", required=True)
    p.add_argument("--out-dir", required=True)
    return p.parse_args()


def direction(value: float) -> str:
    return "positive" if value > 0 else "negative" if value < 0 else "zero"


def module_for(endpoint_id: str, contract: pd.DataFrame) -> str:
    if endpoint_id == "corolla_hue":
        return "colour"
    hit = contract[contract["endpoint_id"].eq(endpoint_id)]
    return str(hit.iloc[0]["module"]) if not hit.empty else "unknown"


def main() -> None:
    args = parse_args()
    contract = pd.read_csv(args.contract, dtype=str, keep_default_na=False)
    continuous = pd.read_csv(args.continuous_coefficients, low_memory=False)
    quality = pd.read_csv(args.candidate_quality, low_memory=False)
    primary = pd.read_csv(args.primary_supported, low_memory=False)
    season = pd.read_csv(args.season, low_memory=False)
    hemisphere = pd.read_csv(args.hemisphere, low_memory=False)
    native = pd.read_csv(args.native, low_memory=False)
    spde = pd.read_csv(args.spde, low_memory=False)
    niche = pd.read_csv(args.niche, low_memory=False)
    pgls = pd.read_csv(args.pgls, low_memory=False)

    rows: list[dict[str, object]] = []

    # Frozen global primary results remain the primary environmental estimand.
    for record in primary.to_dict("records"):
        endpoint_raw = str(record["endpoint"])
        endpoint_id = PRIMARY_NAME_MAP.get(endpoint_raw, endpoint_raw)
        predictor = PREDICTOR_MAP.get(str(record["predictor"]), str(record["predictor"]))
        beta = float(record["standardized_beta"])
        q_text = str(record["bh_q"]).replace("<", "")
        q = float(q_text)
        grade = "B_supported_global"
        robustness = "frozen global primary BH support"
        if endpoint_raw in {"orientation_angle", "corolla_chroma", "shape_aspect_ratio"}:
            n = native[(native["endpoint"].eq(endpoint_raw)) & (native["predictor"].eq(str(record["predictor"])))]
            s = season[(season["endpoint"].eq(endpoint_raw)) & (season["predictor"].eq(str(record["predictor"])))]
            h = hemisphere[(hemisphere["endpoint"].eq(endpoint_raw)) & (hemisphere["predictor"].eq(str(record["predictor"])))]
            native_ok = bool(n.iloc[0]["native_range_robust"]) if not n.empty else False
            if native_ok:
                grade = "A_robust"
                robustness = "global + seasonal + dominant-taxon + hemisphere + native-range support"
            else:
                grade = "B_supported_range_sensitive"
                robustness = "global + seasonal/dominant-taxon/hemisphere support; native-only loses BH support with direction retained"
            if not s.empty and not bool(s.iloc[0]["strict_retention"]):
                robustness += "; seasonal audit not retained"
            if not h.empty and not bool(h.iloc[0]["strict_retention"]):
                robustness += "; hemisphere audit not retained"
        elif endpoint_raw.startswith("hue_"):
            grade = "B_supported_colour_uncalibrated"
            robustness = "global BH-supported circular component; interpret hue jointly; camera-colour calibration remains open"
        rows.append({
            "analysis_layer": "within_taxon_global_primary",
            "endpoint_id": endpoint_id,
            "module": module_for(endpoint_id, contract),
            "predictor": predictor,
            "effect_direction": direction(beta),
            "effect_size": beta,
            "q_value": q,
            "evidence_grade": grade,
            "robustness_context": robustness,
            "causal_scope": "spatial association; not plasticity/adaptation/selection",
        })

    # PR70 candidate continuous traits: report FDR-supported pattern, then grade image sensitivity.
    candidate_sig = continuous[
        continuous["analysis_tier"].eq("candidate")
        & continuous["status"].eq("ok")
        & pd.to_numeric(continuous["q_fdr_bh_within_tier"], errors="coerce").lt(0.05)
    ].copy()
    for record in candidate_sig.to_dict("records"):
        endpoint_id = str(record["endpoint_id"])
        predictor = str(record["predictor"])
        beta = float(record.get("beta_std_within", np.nan))
        q = float(record["q_fdr_bh_within_tier"])
        match = quality[(quality["endpoint_id"].eq(endpoint_id)) & (quality["predictor"].eq(predictor))]
        grade = "C_exploratory_continuous_signal"
        robustness = "candidate-tier FDR-supported continuous image phenotype; botanical calibration pending"
        if not match.empty:
            m = match.iloc[0]
            q2 = float(m["quality_adjusted_q"])
            same = bool(m["same_sign_all_successful_strata"])
            if q2 < 0.05 and same:
                grade = "C_exploratory_quality_robust"
                robustness = "candidate FDR support persists after resolution/sharpness adjustment with same-sign successful strata; botanical calibration pending"
            elif not same and int(m["successful_resolution_strata"]) > 0:
                grade = "C_exploratory_image_sensitive"
                robustness = "candidate FDR signal shows resolution-stratum sign instability; retain as image-sensitive pattern only"
            else:
                grade = "C_exploratory_quality_weakened"
                robustness = "candidate FDR signal weakens after resolution/sharpness adjustment; no automatic erasure of global image-phenotype pattern"
        rows.append({
            "analysis_layer": "within_taxon_pr70_candidate",
            "endpoint_id": endpoint_id,
            "module": module_for(endpoint_id, contract),
            "predictor": predictor,
            "effect_direction": direction(beta),
            "effect_size": beta,
            "q_value": q,
            "evidence_grade": grade,
            "robustness_context": robustness,
            "causal_scope": "exploratory image-derived association; botanical identity/function not established",
        })

    # Spatial-model context from the established main flow.
    for record in spde.to_dict("records"):
        if str(record.get("direction", "")) == "none":
            continue
        endpoint_raw = str(record["endpoint"])
        endpoint_id = PRIMARY_NAME_MAP.get(endpoint_raw, endpoint_raw)
        rows.append({
            "analysis_layer": "grouped_spde_context",
            "endpoint_id": endpoint_id,
            "module": module_for(endpoint_id, contract),
            "predictor": PREDICTOR_MAP.get(str(record["predictor"]), str(record["predictor"])),
            "effect_direction": str(record["direction"]),
            "effect_size": np.nan,
            "q_value": np.nan,
            "evidence_grade": "B_spatial_context",
            "robustness_context": str(record.get("notes", "stable grouped SPDE pattern")),
            "causal_scope": "spatially controlled association",
        })

    # Among-taxon niche sorting: one row per supported endpoint, predictor is multivariate environment.
    for record in niche.to_dict("records"):
        support = str(record.get("Supported pattern", ""))
        if not support:
            continue
        endpoint_name = str(record["Endpoint"])
        endpoint_lookup = {
            "Orientation": "orientation_image_vertical_angle",
            "Lightness": "corolla_lab_lightness",
            "Chroma": "corolla_lab_chroma",
            "Hue sine": "corolla_hue",
            "Hue cosine": "corolla_hue",
            "Width-profile CV": "capitulum_width_profile_cv",
        }
        endpoint_id = endpoint_lookup.get(endpoint_name, endpoint_name)
        rows.append({
            "analysis_layer": "among_taxon_environmental_sorting",
            "endpoint_id": endpoint_id,
            "module": module_for(endpoint_id, contract),
            "predictor": "multivariate_environmental_space",
            "effect_direction": "separation",
            "effect_size": np.nan,
            "q_value": np.nan,
            "evidence_grade": "B_among_taxon_pattern",
            "robustness_context": support,
            "causal_scope": "among-taxon environmental sorting; not within-taxon response",
        })

    # Historical layer stays sensitivity/context.
    for record in pgls.to_dict("records"):
        endpoint_lookup = {
            "Orientation angle": "orientation_image_vertical_angle",
            "Hue sine": "corolla_hue",
            "Solidity": "capitulum_outline_solidity",
            "Width-profile CV": "capitulum_width_profile_cv",
        }
        endpoint_id = endpoint_lookup.get(str(record["Endpoint"]), str(record["Endpoint"]))
        beta = float(str(record["Median β"]).replace("+", "").replace("−", "-"))
        rows.append({
            "analysis_layer": "historical_placement_sensitivity",
            "endpoint_id": endpoint_id,
            "module": module_for(endpoint_id, contract),
            "predictor": PREDICTOR_MAP.get(str(record["Predictor"]), str(record["Predictor"])),
            "effect_direction": direction(beta),
            "effect_size": beta,
            "q_value": np.nan,
            "evidence_grade": "D_historical_context_only",
            "robustness_context": f"{record['Trees with FDR support']} randomized trees; incomplete direct-backbone coverage",
            "causal_scope": "historical placement sensitivity, not resolved phylogenetic correction",
        })

    atlas = pd.DataFrame(rows)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    atlas.to_csv(out / "geb_v2_environmental_evidence_atlas.csv", index=False)

    contract_inferential = contract[contract["analysis_tier"].isin(["primary", "candidate"])]
    report = {
        "positioning": "global continuous within-taxon public-image phenomics",
        "main_flow": [
            "continuous trait measurement beyond species means/categories",
            "global phenotypic diversity and below-taxon variation",
            "within-taxon environmental structure",
            "spatial-model context",
            "among-taxon environmental sorting",
            "historical placement sensitivity",
            "robustness grading and mechanistic hypotheses",
        ],
        "n_contract_endpoints": int(len(contract)),
        "n_inferential_continuous_endpoints": int(len(contract_inferential)),
        "n_primary_endpoints": int((contract["analysis_tier"] == "primary").sum()),
        "n_candidate_endpoints": int((contract["analysis_tier"] == "candidate").sum()),
        "n_environmental_evidence_rows": int(len(atlas)),
        "n_pr70_candidate_fdr_signals": int(len(candidate_sig)),
        "evidence_grade_counts": atlas["evidence_grade"].value_counts().to_dict(),
        "claim_rule": "PR69 sensitivities grade evidence strength; they are not universal vetoes on the frozen global estimand.",
        "causal_boundary": "Pattern and environment-trait structure are the target; plasticity, adaptation and causal mechanisms are not claimed.",
    }
    (out / "geb_v2_evidence_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
