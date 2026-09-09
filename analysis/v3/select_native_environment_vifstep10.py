#!/usr/bin/env python3
"""Deterministic VIF-step selection for the frozen nine GEB-v2 predictors.

This is a user-requested comparative reanalysis.  GEB v2 originally treated VIF
as a retrospective diagnostic because its primary atlas fitted one predictor at
a time.  Here VIFstep is used only to decide which of the same nine marginal
predictor slopes enter the repeated atlas.  It does not turn the atlas into a
joint multivariable model.

Selection is scale-specific because GEB v2 explicitly distinguishes within-
taxon and among-taxon environmental structure:
  * within: observation-level predictors demeaned within taxon;
  * among: taxon medians.
At each step the predictor with the largest VIF is removed until max(VIF) <=
the requested threshold.  The complete-case matrix is frozen before stepping.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

PREDICTORS = [
    "chelsa_bio01",
    "chelsa_bio04",
    "chelsa_bio12",
    "chelsa_bio15",
    "chelsa_rsds_mean",
    "chelsa_vpd_mean",
    "chelsa_sfcwind_mean",
    "chelsa_gsp",
    "chelsa_npp",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--environment", required=True, type=Path)
    parser.add_argument("--out-json", required=True, type=Path)
    parser.add_argument("--out-csv", required=True, type=Path)
    parser.add_argument("--threshold", type=float, default=10.0)
    return parser.parse_args()


def vif_table(values: pd.DataFrame, predictors: list[str]) -> list[dict]:
    if not predictors:
        return []
    z = values[predictors].astype(float)
    z = (z - z.mean()) / z.std(ddof=0)
    rows: list[dict] = []
    for predictor in predictors:
        response = z[predictor].to_numpy(float)
        others = [name for name in predictors if name != predictor]
        if not others:
            r_squared = 0.0
            vif = 1.0
        else:
            design = np.column_stack([np.ones(len(z)), z[others].to_numpy(float)])
            coefficients = np.linalg.lstsq(design, response, rcond=None)[0]
            residual = response - design @ coefficients
            total_ss = float(response @ response)
            residual_ss = float(residual @ residual)
            r_squared = 1.0 - residual_ss / total_ss
            vif = float("inf") if r_squared >= 1.0 else 1.0 / (1.0 - r_squared)
        rows.append(
            {
                "predictor": predictor,
                "vif": float(vif),
                "r_squared_against_remaining": float(r_squared),
            }
        )
    return rows


def step_scope(scope: str, values: pd.DataFrame, threshold: float) -> tuple[dict, list[dict]]:
    complete = values[PREDICTORS].dropna().copy()
    if len(complete) < 20:
        raise ValueError(f"{scope}: too few complete units for VIFstep: {len(complete)}")
    zero = [name for name in PREDICTORS if not np.isfinite(complete[name].std(ddof=0)) or complete[name].std(ddof=0) <= 0]
    if zero:
        raise ValueError(f"{scope}: zero-variance predictors: {zero}")

    remaining = list(PREDICTORS)
    history: list[dict] = []
    step = 0
    while True:
        current = vif_table(complete, remaining)
        for row in current:
            history.append(
                {
                    "scope": scope,
                    "step": step,
                    "n_units": int(len(complete)),
                    "predictor": row["predictor"],
                    "vif": row["vif"],
                    "r_squared_against_remaining": row["r_squared_against_remaining"],
                    "status_at_step": "candidate",
                }
            )
        worst = max(current, key=lambda row: (row["vif"], -PREDICTORS.index(row["predictor"])))
        if worst["vif"] <= threshold:
            break
        dropped = worst["predictor"]
        for row in reversed(history):
            if row["scope"] == scope and row["step"] == step and row["predictor"] == dropped:
                row["status_at_step"] = "dropped_max_vif"
                break
        remaining.remove(dropped)
        step += 1
        if len(remaining) <= 1:
            break

    final = vif_table(complete, remaining)
    final_map = {row["predictor"]: row["vif"] for row in final}
    result = {
        "scope": scope,
        "threshold": threshold,
        "n_complete_units": int(len(complete)),
        "selected_predictors": remaining,
        "dropped_predictors": [name for name in PREDICTORS if name not in remaining],
        "final_vif": final_map,
        "final_max_vif": float(max(final_map.values())) if final_map else None,
        "n_steps": step,
    }
    return result, history


def main() -> int:
    args = parse_args()
    frame = pd.read_csv(args.environment, usecols=["taxon_name", *PREDICTORS], low_memory=False)
    for predictor in PREDICTORS:
        frame[predictor] = pd.to_numeric(frame[predictor], errors="coerce")

    raw = frame[PREDICTORS]
    within = frame[PREDICTORS] - frame.groupby("taxon_name")[PREDICTORS].transform("mean")
    among = frame.groupby("taxon_name")[PREDICTORS].median(numeric_only=True)

    results: dict[str, dict] = {}
    history: list[dict] = []
    for scope, values in [
        ("observation_raw", raw),
        ("within_taxon_demeaned", within),
        ("among_taxon_median", among),
    ]:
        result, rows = step_scope(scope, values, args.threshold)
        results[scope] = result
        history.extend(rows)

    report = {
        "analysis_id": "native_vifstep10_scale_specific_v1",
        "status": "user_requested_comparative_predictor_filter",
        "candidate_predictors": PREDICTORS,
        "threshold_rule": "iteratively remove the current maximum VIF while max(VIF) > 10; retain when max(VIF) <= 10",
        "selection_role": "filter marginal predictor tests only; downstream v2 models remain separate univariate standardized slopes",
        "scales": results,
        "downstream_predictor_sets": {
            "within_taxon": results["within_taxon_demeaned"]["selected_predictors"],
            "among_taxon": results["among_taxon_median"]["selected_predictors"],
        },
    }
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    pd.DataFrame(history).to_csv(args.out_csv, index=False, float_format="%.10g")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
