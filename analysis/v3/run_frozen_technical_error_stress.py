#!/usr/bin/env python3
"""Stress-test the frozen orientation headline using only frozen technical-audit summaries.

Only orientation has a frozen value-scale perturbation summary suitable for a
quantitative stress envelope. The frozen audit does not preserve per-image error
records, and visible colour has no frozen numeric error distribution. Therefore
this script does NOT call the simulations an empirical measurement-error model.
It uses two explicit zero-mean independent Gaussian stress conventions whose
absolute-error p95 values equal the frozen orientation summaries:

- 4.67 degrees: post-QC horizontal-mirror discrepancy p95;
- 54.1 degrees: discrepancy p95 under a deliberate 5% bounding-box shift.

The second is intentionally severe. Neither convention addresses camera roll,
systematic/environment-dependent error, detector selection, or gravity accuracy.
The frozen v2 inference is never overwritten.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

ORIENTATION = "orientation_image_vertical_angle"
PREDICTOR = "chelsa_bio12"
FROZEN_BETA = 0.30435928589775146
NORMAL_ABS_P95 = 1.959963984540054


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--traits", type=Path, required=True)
    p.add_argument("--environment", type=Path, required=True)
    p.add_argument("--technical-audit-summary", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--replicates", type=int, default=2000)
    p.add_argument("--seed", type=int, default=20260910)
    return p.parse_args()


def as_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def z(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    sd = float(np.std(x, ddof=0))
    if not np.isfinite(x).all() or sd <= 0:
        raise ValueError("nonfinite or constant vector")
    return (x - float(np.mean(x))) / sd


def baseline_data(traits_path: Path, env_path: Path) -> tuple[list[str], list[np.ndarray], np.ndarray, float, int]:
    t = pd.read_csv(traits_path, low_memory=False)
    e = pd.read_csv(env_path, low_memory=False)
    t["obs_id"] = t.obs_id.astype(str)
    t["taxon_name"] = t.taxon_name.astype(str)
    e["taxon_name"] = e.taxon_name.astype(str)
    t = t[t.endpoint_id.eq(ORIENTATION)].copy()
    t["measurement_available"] = as_bool(t.measurement_available)
    t["value"] = pd.to_numeric(t.value, errors="coerce")
    t = t[t.measurement_available & t.value.notna()].copy()
    counts = t.groupby("taxon_name").size()
    eligible = sorted(counts[counts >= 5].index)
    t = t[t.taxon_name.isin(eligible)].copy()
    env_med = e.groupby("taxon_name")[PREDICTOR].median().reindex(eligible)
    if env_med.isna().any():
        raise ValueError("missing environment median in eligible orientation taxa")
    groups = [t.loc[t.taxon_name.eq(taxon), "value"].to_numpy(float) for taxon in eligible]
    y = np.array([np.median(g) for g in groups], dtype=float)
    xz = z(env_med.to_numpy(float))
    yz = z(y)
    beta = float(np.dot(xz, yz) / np.dot(xz, xz))
    return eligible, groups, xz, beta, int(len(t))


def simulate(groups: list[np.ndarray], xz: np.ndarray, sigma: float, reps: int, seed: int, batch: int = 100) -> np.ndarray:
    rng = np.random.default_rng(seed)
    out = np.empty(reps, dtype=float)
    den = float(np.dot(xz, xz))
    start = 0
    while start < reps:
        size = min(batch, reps - start)
        med = np.empty((size, len(groups)), dtype=float)
        for j, values in enumerate(groups):
            noise = rng.normal(0.0, sigma, size=(size, len(values)))
            perturbed = np.clip(values[None, :] + noise, 0.0, 180.0)
            med[:, j] = np.median(perturbed, axis=1)
        means = med.mean(axis=1, keepdims=True)
        sds = med.std(axis=1, ddof=0, keepdims=True)
        yz = (med - means) / sds
        out[start : start + size] = (yz @ xz) / den
        start += size
    return out


def summarize(name: str, p95: float, beta: np.ndarray, baseline: float) -> dict[str, object]:
    return {
        "scenario": name,
        "calibration_absolute_error_p95_degrees": p95,
        "simulation_sigma_degrees": p95 / NORMAL_ABS_P95,
        "replicates": int(len(beta)),
        "beta_median": float(np.median(beta)),
        "beta_low95": float(np.quantile(beta, 0.025)),
        "beta_high95": float(np.quantile(beta, 0.975)),
        "sign_retention_fraction": float(np.mean(np.sign(beta) == np.sign(baseline))),
        "positive_fraction": float(np.mean(beta > 0)),
        "median_effect_ratio_to_unperturbed": float(np.median(beta / baseline)),
        "interpretation": "synthetic independent symmetric stress envelope, not an empirical error distribution",
    }


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    audit = json.loads(args.technical_audit_summary.read_text(encoding="utf-8"))
    mirror_p95 = float(audit["mirror_audit"]["orientation_discrepancy_p95_degrees_after_qc"])
    bbox_p95 = float(audit["perturbation_audit"]["orientation_bbox_shift_5pct_error_p95_degrees"])
    colour_summary = str(audit["perturbation_audit"].get("colour_summary", ""))

    taxa, groups, xz, baseline, n_available = baseline_data(args.traits, args.environment)
    if abs(baseline - FROZEN_BETA) > 1e-10:
        raise ValueError(f"orientation-BIO12 baseline mismatch: {baseline} vs {FROZEN_BETA}")

    routine = simulate(groups, xz, mirror_p95 / NORMAL_ABS_P95, args.replicates, args.seed + 1)
    severe = simulate(groups, xz, bbox_p95 / NORMAL_ABS_P95, args.replicates, args.seed + 2)
    results = pd.DataFrame([
        summarize("mirror_p95_calibrated_random_stress", mirror_p95, routine, baseline),
        summarize("five_percent_bbox_shift_p95_calibrated_severe_stress", bbox_p95, severe, baseline),
    ])
    results.to_csv(args.out_dir / "orientation_bio12_technical_stress.csv", index=False)

    coverage = pd.DataFrame([
        {
            "trait_or_module": "orientation",
            "quantitative_propagation_status": "stress_test_available",
            "frozen_numeric_information": f"mirror discrepancy p95={mirror_p95} deg; 5% bbox-shift discrepancy p95={bbox_p95} deg",
            "limitation": "no per-image frozen error distribution; synthetic independent symmetric stress only",
        },
        {
            "trait_or_module": "visible_colour",
            "quantitative_propagation_status": "unsupported_by_frozen_summary",
            "frozen_numeric_information": colour_summary or "generally stable; no single numeric error summary frozen",
            "limitation": "cannot invent a chroma error distribution; independent calibrated colour validation remains required",
        },
        {
            "trait_or_module": "gross_outline",
            "quantitative_propagation_status": "unsupported_on_value_scale",
            "frozen_numeric_information": f"half-resolution rank agreement approximately {audit['perturbation_audit']['outline_half_resolution_spearman_rho_approx']}",
            "limitation": "rank correlation does not identify endpoint-specific value-scale error distributions",
        },
    ])
    coverage.to_csv(args.out_dir / "technical_error_propagation_coverage.csv", index=False)

    report = {
        "analysis_id": "ch1_v3_frozen_technical_error_stress_20260910",
        "headline": "presentation_angle_x_annual_precipitation",
        "endpoint": ORIENTATION,
        "predictor": PREDICTOR,
        "eligible_taxa": len(taxa),
        "available_orientation_observations_in_eligible_taxa": n_available,
        "unperturbed_beta": baseline,
        "frozen_beta_reference": FROZEN_BETA,
        "stress_scenarios": results.to_dict("records"),
        "colour_chroma_status": "quantitative error propagation not executed because the frozen audit preserves no numeric colour error distribution",
        "claim_boundary": "post-hoc synthetic technical stress; does not address systematic camera roll, environment-dependent measurement error, detector-bbox selection, gravity accuracy, or biological validation; frozen v2 inference unchanged",
    }
    (args.out_dir / "frozen_technical_error_stress_report.json").write_text(json.dumps(report, indent=2, allow_nan=False), encoding="utf-8")
    print(json.dumps(report, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
