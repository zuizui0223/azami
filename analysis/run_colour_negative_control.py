#!/usr/bin/env python3
"""Measure non-flower Lab controls and test paired region-by-BIO12 slopes."""
from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import cv2
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests


CONTROL_KINDS = ("border", "green")
METRICS = ("lightness", "chroma")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--packet-root", default=".")
    parser.add_argument("--environment", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--min-control-pixels", type=int, default=200)
    parser.add_argument("--min-taxon-n", type=int, default=5)
    return parser.parse_args()


def resolve_path(value: object, root: Path) -> Path:
    path = Path(str(value))
    return path if path.is_absolute() else root / path


def lab_summary(lab: np.ndarray, mask: np.ndarray) -> tuple[float, float, int]:
    n_pixels = int(mask.sum())
    if n_pixels == 0:
        return float("nan"), float("nan"), 0
    lightness = lab[..., 0].astype(np.float32) * 100.0 / 255.0
    a_axis = lab[..., 1].astype(np.float32) - 128.0
    b_axis = lab[..., 2].astype(np.float32) - 128.0
    chroma = np.sqrt(a_axis ** 2 + b_axis ** 2)
    return (
        float(np.median(lightness[mask])),
        float(np.median(chroma[mask])),
        n_pixels,
    )


def measure_controls(image: np.ndarray, minimum: int = 200) -> dict[str, object]:
    if image is None or min(image.shape[:2]) < 40:
        return {"status": "image_too_small"}
    height, width = image.shape[:2]
    yy, xx = np.mgrid[0:height, 0:width]
    radial = np.sqrt(
        ((xx - (width - 1) / 2) / max(width / 2, 1)) ** 2
        + ((yy - (height - 1) / 2) / max(height / 2, 1)) ** 2
    )
    # The central detector region is excluded; the outer context annulus is the
    # paired imaging control. Thresholds are fixed before outcome execution.
    border = (radial >= 0.58) & (radial <= 1.0)
    hsv = cv2.cvtColor(image, cv2.COLOR_BGR2HSV)
    hue, saturation, value = [hsv[..., index] for index in range(3)]
    green = (
        border
        & (hue >= 25) & (hue <= 95)
        & (saturation >= 35)
        & (value >= 30) & (value <= 245)
    )
    lab = cv2.cvtColor(image, cv2.COLOR_BGR2LAB)
    border_l, border_c, border_n = lab_summary(lab, border)
    green_l, green_c, green_n = lab_summary(lab, green)
    return {
        "status": "usable" if border_n >= minimum else "insufficient_border_pixels",
        "border_lab_lightness": border_l if border_n >= minimum else float("nan"),
        "border_lab_chroma": border_c if border_n >= minimum else float("nan"),
        "border_n_pixels": border_n,
        "green_lab_lightness": green_l if green_n >= minimum else float("nan"),
        "green_lab_chroma": green_c if green_n >= minimum else float("nan"),
        "green_n_pixels": green_n,
        "green_status": "usable" if green_n >= minimum else "insufficient_green_pixels",
    }


def measure_manifest(manifest: pd.DataFrame, root: Path, minimum: int) -> pd.DataFrame:
    required = {"annotation_unit_id", "obs_id", "taxon_name", "context_crop_path"}
    missing = required.difference(manifest.columns)
    if missing:
        raise ValueError(f"Manifest is missing: {sorted(missing)}")
    rows: list[dict[str, object]] = []
    for record in manifest.to_dict("records"):
        path = resolve_path(record["context_crop_path"], root)
        image = cv2.imread(str(path), cv2.IMREAD_COLOR)
        result = measure_controls(image, minimum)
        rows.append({
            "annotation_unit_id": record["annotation_unit_id"],
            "obs_id": str(record["obs_id"]),
            "taxon_name": record["taxon_name"],
            "context_crop_path": str(record["context_crop_path"]),
            **result,
        })
    return pd.DataFrame(rows)


def aggregate_observations(measurements: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "border_lab_lightness", "border_lab_chroma",
        "green_lab_lightness", "green_lab_chroma",
    ]
    usable = measurements.copy()
    for column in columns:
        usable[column] = pd.to_numeric(usable[column], errors="coerce")
    return (
        usable.groupby(["obs_id", "taxon_name"], as_index=False)
        .agg(
            n_context_crops=("annotation_unit_id", "nunique"),
            **{f"{column}_median": (column, "median") for column in columns},
        )
    )


def paired_region_fit(
    frame: pd.DataFrame,
    metric: str,
    control_kind: str,
    minimum_taxon_n: int,
) -> dict[str, object]:
    flower = f"corolla_lab_{metric}_median"
    background = f"{control_kind}_lab_{metric}_median"
    columns = ["obs_id", "taxon_name", "chelsa_bio12", flower, background]
    work = frame[columns].copy()
    for column in ["chelsa_bio12", flower, background]:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    work = work.dropna()
    counts = work.groupby("taxon_name").size()
    work = work.loc[work["taxon_name"].isin(counts.loc[counts.ge(minimum_taxon_n)].index)].copy()
    if len(work) < 50 or work["taxon_name"].nunique() < 10:
        return {"status": "insufficient"}

    climate = work["chelsa_bio12"] - work.groupby("taxon_name")["chelsa_bio12"].transform("mean")
    climate_sd = float(climate.std(ddof=1))
    if not math.isfinite(climate_sd) or climate_sd <= 0:
        return {"status": "no_climate_variance"}
    climate = climate / climate_sd
    long = pd.concat([
        pd.DataFrame({
            "obs_id": work["obs_id"], "taxon_name": work["taxon_name"],
            "region": "flower", "response": work[flower], "climate": climate,
        }),
        pd.DataFrame({
            "obs_id": work["obs_id"], "taxon_name": work["taxon_name"],
            "region": "background", "response": work[background], "climate": climate,
        }),
    ], ignore_index=True)
    long["response"] = long["response"] - long.groupby(
        ["taxon_name", "region"]
    )["response"].transform("mean")
    long["background"] = long["region"].eq("background").astype(float)
    design = pd.DataFrame({
        "constant": 1.0,
        "climate": long["climate"],
        "background": long["background"],
        "climate_x_background": long["climate"] * long["background"],
    })
    fit = sm.OLS(long["response"], design).fit(
        cov_type="cluster", cov_kwds={"groups": long["taxon_name"]}
    )
    flower_slope = float(fit.params["climate"])
    interaction = float(fit.params["climate_x_background"])
    background_slope = flower_slope + interaction
    covariance = fit.cov_params()
    background_variance = float(
        covariance.loc["climate", "climate"]
        + covariance.loc["climate_x_background", "climate_x_background"]
        + 2 * covariance.loc["climate", "climate_x_background"]
    )
    return {
        "status": "ok",
        "metric": metric,
        "control_kind": control_kind,
        "n_paired_observations": int(len(work)),
        "n_taxa": int(work["taxon_name"].nunique()),
        "flower_slope_lab_units": flower_slope,
        "background_slope_lab_units": background_slope,
        "background_slope_se": math.sqrt(max(background_variance, 0.0)),
        "background_minus_flower": interaction,
        "interaction_se": float(fit.bse["climate_x_background"]),
        "interaction_p_value": float(fit.pvalues["climate_x_background"]),
    }


def main() -> None:
    args = parse_args()
    output = Path(args.out_dir)
    output.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(args.manifest, dtype=str, keep_default_na=False)
    measurements = measure_manifest(manifest, Path(args.packet_root), args.min_control_pixels)
    observations = aggregate_observations(measurements)
    environment = pd.read_csv(args.environment, low_memory=False)
    environment["obs_id"] = environment["obs_id"].astype(str)
    required = {
        "obs_id", "taxon_name", "chelsa_bio12",
        "corolla_lab_lightness_median", "corolla_lab_chroma_median",
    }
    missing = required.difference(environment.columns)
    if missing:
        raise ValueError(f"Environment table is missing: {sorted(missing)}")
    joined = observations.merge(
        environment[[
            "obs_id", "chelsa_bio12", "corolla_lab_lightness_median",
            "corolla_lab_chroma_median",
        ]].drop_duplicates("obs_id"),
        on="obs_id", how="inner", validate="one_to_one",
    )
    tests = pd.DataFrame([
        paired_region_fit(joined, metric, control, args.min_taxon_n)
        for metric in METRICS for control in CONTROL_KINDS
    ])
    tests["interaction_q_fdr_bh_4"] = np.nan
    successful = tests["status"].eq("ok")
    if int(successful.sum()) == 4:
        _, q_values, _, _ = multipletests(
            tests.loc[successful, "interaction_p_value"].to_numpy(float),
            alpha=0.05, method="fdr_bh",
        )
        tests.loc[successful, "interaction_q_fdr_bh_4"] = q_values

    measurements.to_csv(output / "colour_negative_control_head_measurements.csv", index=False)
    observations.to_csv(output / "colour_negative_control_observations.csv", index=False)
    tests.to_csv(output / "colour_negative_control_paired_tests.csv", index=False)
    report = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "complete" if int(successful.sum()) == 4 else "incomplete",
        "n_manifest_heads": int(len(manifest)),
        "n_measured_heads": int(measurements["status"].eq("usable").sum()),
        "n_joined_observations": int(len(joined)),
        "test_family": "four paired region-by-BIO12 interactions: two Lab metrics by two background definitions",
        "interpretation_rule": (
            "Specificity requires a flower-versus-background interaction; a nonsignificant background slope alone is insufficient."
        ),
    }
    (output / "colour_negative_control_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
