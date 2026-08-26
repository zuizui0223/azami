#!/usr/bin/env python3
"""Compare registry-wide candidate signals with locked PR #69 bias controls.

The registry-wide climate pass is deliberately a screen.  This module makes
the handoff to the frozen reviewer controls explicit so a small generic-screen
q value cannot be promoted directly into a manuscript claim.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--screen-coefficients", required=True)
    parser.add_argument("--adjusted-models", required=True)
    parser.add_argument("--headline-audit", required=True)
    parser.add_argument("--out-dir", required=True)
    return parser.parse_args()


def _as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def compare(
    screen: pd.DataFrame,
    adjusted: pd.DataFrame,
    headline: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    required_screen = {
        "endpoint_id",
        "analysis_tier",
        "predictor",
        "beta_std_within",
        "p_value",
        "q_fdr_bh_within_tier",
        "fdr_significant_0_05",
        "n_observations",
        "n_taxa",
    }
    required_adjusted = {
        "endpoint",
        "cohort",
        "predictor",
        "beta_standardized",
        "p_value",
        "q_fdr_bh_climate_12",
        "n_observations",
        "n_taxa",
    }
    required_headline = {
        "endpoint",
        "adjusted_bio4_q",
        "positive_in_all_successful_strata",
        "strict_resolution_control_retained",
    }
    for label, frame, required in (
        ("screen", screen, required_screen),
        ("adjusted", adjusted, required_adjusted),
        ("headline", headline, required_headline),
    ):
        missing = sorted(required.difference(frame.columns))
        if missing:
            raise ValueError(f"{label} table is missing columns: {missing}")

    signals = screen.loc[
        screen["analysis_tier"].eq("candidate")
        & _as_bool(screen["fdr_significant_0_05"])
    ].copy()
    if signals.empty:
        columns = [
            "endpoint_id", "predictor", "screen_beta", "screen_q",
            "quality_adjusted_beta", "quality_adjusted_q", "retention_status",
        ]
        output = pd.DataFrame(columns=columns)
        report = {
            "status": "complete_no_candidate_screen_signals",
            "n_candidate_screen_signals": 0,
            "n_with_matching_pr69_quality_model": 0,
            "n_submission_eligible": 0,
        }
        return output, report

    quality = adjusted.loc[adjusted["cohort"].eq("all_resolution_adjusted")].copy()
    quality = quality.rename(columns={
        "endpoint": "endpoint_id",
        "beta_standardized": "quality_adjusted_beta",
        "p_value": "quality_adjusted_p",
        "q_fdr_bh_climate_12": "quality_adjusted_q",
        "n_observations": "quality_adjusted_n_observations",
        "n_taxa": "quality_adjusted_n_taxa",
    })
    quality = quality[[
        "endpoint_id", "predictor", "quality_adjusted_beta", "quality_adjusted_p",
        "quality_adjusted_q", "quality_adjusted_n_observations", "quality_adjusted_n_taxa",
    ]]
    if quality.duplicated(["endpoint_id", "predictor"]).any():
        raise ValueError("Adjusted model table is not unique by endpoint/predictor")

    output = signals.rename(columns={
        "beta_std_within": "screen_beta",
        "p_value": "screen_p",
        "q_fdr_bh_within_tier": "screen_q",
        "n_observations": "screen_n_observations",
        "n_taxa": "screen_n_taxa",
    }).merge(
        quality,
        on=["endpoint_id", "predictor"],
        how="left",
        validate="one_to_one",
        indicator="quality_gate_join",
    )

    headline_payload = headline[[
        "endpoint", "positive_in_all_successful_strata",
        "strict_resolution_control_retained",
    ]].rename(columns={"endpoint": "endpoint_id"})
    if headline_payload["endpoint_id"].duplicated().any():
        raise ValueError("Headline audit is not unique by endpoint")
    output = output.merge(headline_payload, on="endpoint_id", how="left", validate="many_to_one")

    reasons: list[str] = []
    for row in output.to_dict("records"):
        row_reasons: list[str] = []
        adjusted_q = pd.to_numeric(pd.Series([row.get("quality_adjusted_q")]), errors="coerce").iloc[0]
        if row.get("quality_gate_join") != "both":
            row_reasons.append("missing_endpoint_specific_pr69_quality_model")
        elif not np.isfinite(adjusted_q) or float(adjusted_q) >= 0.05:
            row_reasons.append("fails_locked_quality_adjusted_bh")
        if row.get("predictor") == "env_chelsa_bio04_native":
            all_positive = str(row.get("positive_in_all_successful_strata")).strip().lower() in {
                "true", "1", "yes"
            }
            if not all_positive:
                row_reasons.append("fails_resolution_stratum_sign_consistency")
        row_reasons.append("independent_botanical_reference_validation_pending")
        reasons.append(";".join(row_reasons))
    output["failure_reasons"] = reasons
    output["retention_status"] = "not_retained_after_pr69_gates"
    output["submission_claim_eligible"] = False

    ordered = [
        "endpoint_id", "predictor", "screen_beta", "screen_p", "screen_q",
        "screen_n_observations", "screen_n_taxa", "quality_adjusted_beta",
        "quality_adjusted_p", "quality_adjusted_q", "quality_adjusted_n_observations",
        "quality_adjusted_n_taxa", "positive_in_all_successful_strata",
        "strict_resolution_control_retained", "failure_reasons", "retention_status",
        "submission_claim_eligible",
    ]
    output = output[ordered].sort_values(["screen_q", "endpoint_id", "predictor"]).reset_index(drop=True)
    report = {
        "status": "complete",
        "n_candidate_screen_signals": int(len(output)),
        "n_with_matching_pr69_quality_model": int(
            output["quality_adjusted_q"].notna().sum()
        ),
        "n_passing_quality_adjusted_bh": int(
            pd.to_numeric(output["quality_adjusted_q"], errors="coerce").lt(0.05).sum()
        ),
        "n_passing_full_pr69_retention_rule": int(
            _as_bool(output["strict_resolution_control_retained"]).sum()
        ),
        "n_submission_eligible": int(_as_bool(output["submission_claim_eligible"]).sum()),
        "interpretation": (
            "Generic candidate-screen q values are discovery diagnostics. PR69 endpoint-specific "
            "quality, resolution-stratum and independent botanical-validation gates control promotion."
        ),
    }
    return output, report


def main() -> None:
    args = parse_args()
    screen = pd.read_csv(args.screen_coefficients, low_memory=False)
    adjusted = pd.read_csv(args.adjusted_models, low_memory=False)
    headline = pd.read_csv(args.headline_audit, low_memory=False)
    output, report = compare(screen, adjusted, headline)
    destination = Path(args.out_dir)
    destination.mkdir(parents=True, exist_ok=True)
    output.to_csv(destination / "candidate_screen_vs_pr69_gates.csv", index=False, encoding="utf-8-sig")
    (destination / "candidate_screen_vs_pr69_gates_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
