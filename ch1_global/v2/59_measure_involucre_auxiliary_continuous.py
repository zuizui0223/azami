#!/usr/bin/env python3
"""Measure exploratory continuous involucre/bract proxies on close head crops.

These variables are image-geometry proxies, not categorical botanical truth.
They are retained only when resolution, sharpness, segmentation and horizontal-
flip repeatability pass explicit QC. Failure never means biological absence.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import cv2
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", required=True)
    p.add_argument("--packet-root", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--min-dimension", type=int, default=150)
    p.add_argument("--min-sharpness", type=float, default=55.0)
    p.add_argument("--min-mask-quality", type=float, default=0.55)
    return p.parse_args()


def fast_segment(image: np.ndarray) -> tuple[np.ndarray, np.ndarray | None, dict[str, Any]]:
    original_height, original_width = image.shape[:2]
    scale = min(1.0, 320.0 / max(original_height, original_width))
    resized = cv2.resize(image, None, fx=scale, fy=scale, interpolation=cv2.INTER_AREA) if scale < 1 else image.copy()
    height, width = resized.shape[:2]
    lab = cv2.cvtColor(resized, cv2.COLOR_BGR2LAB).astype(np.float32)
    strip_h, strip_w = max(2, height // 20), max(2, width // 20)
    border = np.concatenate([
        lab[:strip_h].reshape(-1, 3), lab[-strip_h:].reshape(-1, 3),
        lab[:, :strip_w].reshape(-1, 3), lab[:, -strip_w:].reshape(-1, 3),
    ])
    background = np.median(border, axis=0)
    distance = np.linalg.norm(lab - background, axis=2)
    distance8 = np.uint8(np.clip(distance / max(np.percentile(distance, 99), 1) * 255, 0, 255))
    _, mask = cv2.threshold(distance8, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    yy, xx = np.mgrid[0:height, 0:width]
    central = (((xx - (width - 1) / 2) / (0.54 * width)) ** 2 + ((yy - (height - 1) / 2) / (0.54 * height)) ** 2) <= 1
    mask = ((mask > 0) & central).astype(np.uint8)
    kernel_size = max(3, (int(round(min(height, width) * 0.025)) // 2) * 2 + 1)
    mask = cv2.morphologyEx(mask, cv2.MORPH_CLOSE, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (kernel_size, kernel_size)))
    mask = cv2.morphologyEx(mask, cv2.MORPH_OPEN, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3)))
    contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    centre = np.array([width / 2, height / 2])
    candidates: list[tuple[float, np.ndarray]] = []
    for contour in contours:
        area = cv2.contourArea(contour)
        if area < 0.02 * height * width:
            continue
        moments = cv2.moments(contour)
        if moments["m00"] <= 0:
            continue
        centroid = np.array([moments["m10"] / moments["m00"], moments["m01"] / moments["m00"]])
        offset = np.linalg.norm((centroid - centre) / np.array([width, height]))
        candidates.append((area / (1 + 6 * offset), contour))
    if not candidates:
        return resized, None, {}
    contour = max(candidates, key=lambda item: item[0])[1]
    foreground = np.zeros((height, width), np.uint8)
    cv2.drawContours(foreground, [contour], -1, 1, -1)
    hull = cv2.convexHull(contour)
    solidity = cv2.contourArea(contour) / max(cv2.contourArea(hull), 1)
    moments = cv2.moments(contour)
    centroid = np.array([moments["m10"] / moments["m00"], moments["m01"] / moments["m00"]])
    centre_offset = float(np.linalg.norm((centroid - centre) / np.array([width / 2, height / 2])))
    points = contour[:, 0, :]
    border_touch = float(np.mean((points[:, 0] <= 1) | (points[:, 0] >= width - 2) | (points[:, 1] <= 1) | (points[:, 1] >= height - 2)))
    area_fraction = float(foreground.mean())
    quality = float(np.clip(0.45 * solidity + 0.25 * (1 - centre_offset) + 0.20 * (1 - border_touch) + 0.10 * np.clip(area_fraction / 0.3, 0, 1), 0, 1))
    return resized, foreground, {
        "contour": contour,
        "mask_quality": quality,
        "mask_solidity": solidity,
        "mask_centre_offset": centre_offset,
        "mask_border_touch": border_touch,
        "mask_area_fraction": area_fraction,
    }


def floral_mask(image: np.ndarray, foreground: np.ndarray) -> np.ndarray:
    hsv = cv2.cvtColor(image, cv2.COLOR_BGR2HSV)
    gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    hue, saturation, value = [hsv[..., index] for index in range(3)]
    base = foreground.astype(bool)
    purple = base & (value >= 80) & (saturation >= 35) & (hue >= 105) & (hue <= 164)
    redmagenta = base & (value >= 80) & (saturation >= 25) & ((hue >= 165) | (hue <= 10))
    yellow = base & (value >= 115) & (saturation >= 45) & (hue >= 11) & (hue <= 38)
    occupied = purple | redmagenta | yellow
    laplacian = np.abs(cv2.Laplacian(gray, cv2.CV_32F))
    edge_threshold = np.percentile(laplacian[base], 55) if base.any() else 0
    white = base & ~occupied & (value >= 170) & (saturation <= 58) & (laplacian >= edge_threshold)
    return (purple | redmagenta | yellow | white).astype(np.uint8)


def circular_smooth(values: np.ndarray, width: int = 11) -> np.ndarray:
    if width % 2 == 0:
        width += 1
    padding = width // 2
    extended = np.r_[values[-padding:], values, values[:padding]]
    return np.convolve(extended, np.ones(width) / width, mode="valid")


def measure_once(image: np.ndarray) -> dict[str, Any]:
    original_height, original_width = image.shape[:2]
    sharpness = float(cv2.Laplacian(cv2.cvtColor(image, cv2.COLOR_BGR2GRAY), cv2.CV_64F).var())
    resized, foreground, quality = fast_segment(image)
    if foreground is None:
        return {"measurement_status": "mask_failed", "sharpness": sharpness, "min_dimension": min(original_height, original_width)}
    contour = quality.pop("contour")
    floral = floral_mask(resized, foreground)
    rows, columns = np.nonzero(foreground)
    floral_rows, floral_columns = np.nonzero(floral)
    floral_fraction = float(len(floral_columns) / max(len(columns), 1))
    base = {
        "sharpness": sharpness,
        "min_dimension": min(original_height, original_width),
        "floral_fraction": floral_fraction,
        **quality,
    }
    if len(columns) < 100:
        return {"measurement_status": "mask_failed", **base}
    if len(floral_columns) < max(12, int(0.02 * len(columns))):
        return {"measurement_status": "floral_end_not_recovered", **base}

    coordinates = np.column_stack([columns, rows]).astype(np.float32)
    mean, eigenvectors, _ = cv2.PCACompute2(coordinates, mean=None, maxComponents=2)
    centre = mean[0].astype(float)
    axis = eigenvectors[0].astype(float)
    axis /= max(np.linalg.norm(axis), 1e-9)
    floral_centroid = np.array([np.median(floral_columns), np.median(floral_rows)])
    if np.dot(floral_centroid - centre, axis) < 0:
        axis = -axis

    contour_points = contour[:, 0, :].astype(float)
    projection = (contour_points - centre) @ axis
    involucre = contour_points[projection <= np.quantile(projection, 0.72)]
    if len(involucre) < 30:
        return {"measurement_status": "involucre_contour_not_recovered", **base}
    vectors = involucre - centre
    angles = (np.arctan2(vectors[:, 1], vectors[:, 0]) + 2 * np.pi) % (2 * np.pi)
    radii = np.linalg.norm(vectors, axis=1)
    bins = 120
    identifiers = np.floor(angles / (2 * np.pi) * bins).astype(int) % bins
    profile = np.full(bins, np.nan)
    for identifier in np.unique(identifiers):
        profile[identifier] = radii[identifiers == identifier].max()
    valid = np.isfinite(profile)
    if valid.sum() < 25:
        return {"measurement_status": "involucre_contour_not_recovered", **base}
    positions = np.arange(bins)
    valid_positions = positions[valid]
    valid_values = profile[valid]
    interpolated = np.interp(
        positions,
        np.r_[valid_positions - bins, valid_positions, valid_positions + bins],
        np.r_[valid_values, valid_values, valid_values],
    )
    residual = interpolated - circular_smooth(interpolated)
    equivalent_radius = math.sqrt(float(foreground.sum()) / math.pi)
    positive = np.maximum(residual[valid], 0)
    threshold = max(0.045 * equivalent_radius, float(np.percentile(positive, 85)))
    peak_count = sum(
        1 for index in range(bins)
        if valid[index]
        and residual[index] > threshold
        and residual[index] >= residual[(index - 1) % bins]
        and residual[index] >= residual[(index + 1) % bins]
    )
    return {
        "measurement_status": "measured",
        **base,
        "involucre_projection_roughness": float(np.std(residual[valid]) / max(equivalent_radius, 1)),
        "involucre_projection_p95": float(np.percentile(positive, 95) / max(equivalent_radius, 1)),
        "involucre_projection_max": float(np.max(positive) / max(equivalent_radius, 1)),
        "involucre_spread_fraction": float(np.mean(positive > 0.04 * equivalent_radius)),
        "spine_peak_count_proxy": int(peak_count),
        "spine_relative_length_max_proxy": float(np.max(positive) / max(equivalent_radius, 1)),
        "involucre_contour_bins": int(valid.sum()),
    }


def combine(original: dict[str, Any], mirrored: dict[str, Any], args: argparse.Namespace) -> dict[str, Any]:
    output = {f"original_{key}": value for key, value in original.items() if key not in {"measurement_status"}}
    output.update({f"mirror_{key}": value for key, value in mirrored.items() if key not in {"measurement_status"}})
    if original.get("measurement_status") != "measured" or mirrored.get("measurement_status") != "measured":
        output["involucre_status"] = original.get("measurement_status", "failed")
        output["mirror_measurement_status"] = mirrored.get("measurement_status", "failed")
        return output
    reasons: list[str] = []
    if original["min_dimension"] < args.min_dimension:
        reasons.append("low_resolution")
    if original["sharpness"] < args.min_sharpness:
        reasons.append("low_sharpness")
    if original["mask_quality"] < args.min_mask_quality or original["mask_border_touch"] > 0.18:
        reasons.append("low_mask_quality")
    if original["floral_fraction"] < 0.025:
        reasons.append("low_floral_evidence")

    metrics = (
        "involucre_projection_roughness",
        "involucre_projection_p95",
        "involucre_projection_max",
        "involucre_spread_fraction",
        "spine_relative_length_max_proxy",
    )
    thresholds = {
        "involucre_projection_roughness": 0.04,
        "involucre_projection_p95": 0.05,
        "involucre_projection_max": 0.07,
        "involucre_spread_fraction": 0.12,
        "spine_relative_length_max_proxy": 0.07,
    }
    for metric in metrics:
        delta = abs(float(original[metric]) - float(mirrored[metric]))
        output[f"{metric}_flip_delta"] = delta
        output[metric] = (float(original[metric]) + float(mirrored[metric])) / 2
        if delta > thresholds[metric]:
            reasons.append("unstable_under_horizontal_flip")
    output["spine_peak_count_proxy"] = int(round((original["spine_peak_count_proxy"] + mirrored["spine_peak_count_proxy"]) / 2))
    output["involucre_contour_bins"] = int(round((original["involucre_contour_bins"] + mirrored["involucre_contour_bins"]) / 2))
    output["involucre_status"] = "usable" if not reasons else ";".join(sorted(set(reasons)))
    output["involucre_confidence"] = float(np.clip(
        0.30 * min(1, original["min_dimension"] / 250)
        + 0.25 * np.clip(original["sharpness"] / 180, 0, 1)
        + 0.35 * original["mask_quality"]
        + 0.10 * np.clip(original["floral_fraction"] / 0.12, 0, 1),
        0, 1,
    ))
    return output


def aggregate(table: pd.DataFrame, group: str, prefix: str) -> pd.DataFrame:
    usable = table.loc[table["involucre_status"].eq("usable")].copy()
    metrics = [
        "involucre_projection_roughness", "involucre_projection_p95",
        "involucre_projection_max", "involucre_spread_fraction",
        "spine_peak_count_proxy", "spine_relative_length_max_proxy",
    ]
    rows: list[dict[str, Any]] = []
    for key, part in table.groupby(group, sort=True):
        accepted = usable.loc[usable[group].eq(key)]
        row: dict[str, Any] = {group: key, f"n_heads_{prefix}": int(len(part)), f"n_usable_heads_{prefix}": int(len(accepted))}
        for metric in metrics:
            values = pd.to_numeric(accepted[metric], errors="coerce").dropna()
            row[f"{metric}_{prefix}_median"] = values.median() if len(values) else np.nan
            row[f"{metric}_{prefix}_mad"] = (values - values.median()).abs().median() if len(values) else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    manifest = pd.read_csv(args.manifest, dtype=str, keep_default_na=False)
    metadata = pd.read_csv(args.metadata, dtype=str, keep_default_na=False)
    required_manifest = {"annotation_unit_id", "crop_path"}
    if required_manifest.difference(manifest.columns):
        raise ValueError("Manifest lacks annotation_unit_id/crop_path")
    metadata["obs_id"] = metadata["obs_id"].astype(str)
    joined = manifest.merge(
        metadata[["annotation_unit_id", "obs_id", "photo_id", "taxon_name", "latitude", "longitude"]],
        on="annotation_unit_id", how="left", validate="one_to_one",
    )
    if joined["taxon_name"].eq("").any() or joined["taxon_name"].isna().any():
        raise ValueError("Manifest rows lack metadata")

    root = Path(args.packet_root)
    rows: list[dict[str, Any]] = []
    for record in joined.to_dict("records"):
        path = root / record["crop_path"]
        image = cv2.imread(str(path))
        if image is None:
            result = {"involucre_status": "image_unreadable"}
        else:
            result = combine(measure_once(image), measure_once(cv2.flip(image, 1)), args)
        rows.append({**record, **result})
    table = pd.DataFrame(rows)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    table.to_csv(out / "involucre_auxiliary_head_level.csv", index=False, encoding="utf-8-sig")
    observation = aggregate(table, "obs_id", "observation")
    observation = observation.merge(
        joined[["obs_id", "taxon_name", "latitude", "longitude"]].drop_duplicates("obs_id"),
        on="obs_id", how="left", validate="one_to_one",
    )
    observation.to_csv(out / "involucre_auxiliary_observation_level.csv", index=False, encoding="utf-8-sig")
    species = aggregate(table, "taxon_name", "species")
    species.to_csv(out / "involucre_auxiliary_species_level.csv", index=False, encoding="utf-8-sig")
    summary = table["involucre_status"].value_counts(dropna=False).rename_axis("status").reset_index(name="n_heads")
    summary["fraction"] = summary["n_heads"] / len(table)
    summary.to_csv(out / "involucre_auxiliary_status_summary.csv", index=False)
    report = {
        "n_heads": int(len(table)),
        "n_usable_heads": int(table["involucre_status"].eq("usable").sum()),
        "n_observations_with_usable": int(table.loc[table["involucre_status"].eq("usable"), "obs_id"].nunique()),
        "n_species_with_usable": int(table.loc[table["involucre_status"].eq("usable"), "taxon_name"].nunique()),
        "status_counts": summary.to_dict("records"),
        "semantic_status": (
            "Exploratory high-resolution continuous contour proxies for involucre projection and "
            "spine-like protrusion. They do not identify botanical bract states, hair or mucilage."
        ),
    }
    (out / "involucre_auxiliary_provenance.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
