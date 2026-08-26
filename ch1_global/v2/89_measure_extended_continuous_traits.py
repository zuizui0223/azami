#!/usr/bin/env python3
"""Measure category-free extended capitulum traits from high-resolution crops.

This module adds continuous display, involucral architecture, armature and
surface-texture measurements without pretending that generic photographs reveal
botanical states.  Architecture and surface traits have separate assessability
gates.  Every measurement is repeated after horizontal mirroring; unstable
measurements remain unavailable rather than becoming a biological zero.
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


ARCHITECTURE_METRICS = (
    "involucre_length_width_ratio",
    "involucre_apical_taper_ratio",
    "involucre_basal_taper_ratio",
    "involucre_projection_roughness",
    "involucre_projection_p95",
    "involucre_projection_max",
    "involucre_spread_fraction",
    "bract_projection_peak_density",
    "bract_projection_asymmetry",
)
SURFACE_METRICS = (
    "involucre_surface_edge_density",
    "involucre_surface_lbp_entropy",
    "involucre_surface_high_frequency_energy",
    "involucre_surface_specular_fraction_proxy",
)
FLIP_TOLERANCE = {
    "involucre_length_width_ratio": 0.15,
    "involucre_apical_taper_ratio": 0.16,
    "involucre_basal_taper_ratio": 0.16,
    "involucre_projection_roughness": 0.04,
    "involucre_projection_p95": 0.05,
    "involucre_projection_max": 0.07,
    "involucre_spread_fraction": 0.12,
    "bract_projection_peak_density": 0.08,
    "bract_projection_asymmetry": 0.10,
    "involucre_surface_edge_density": 0.10,
    "involucre_surface_lbp_entropy": 0.12,
    "involucre_surface_high_frequency_energy": 0.60,
    "involucre_surface_specular_fraction_proxy": 0.08,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--metadata")
    parser.add_argument("--image-root", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--min-architecture-dimension", type=int, default=150)
    parser.add_argument("--min-surface-dimension", type=int, default=300)
    parser.add_argument("--min-architecture-sharpness", type=float, default=45.0)
    parser.add_argument("--min-surface-sharpness", type=float, default=80.0)
    return parser.parse_args()


def _largest_central_component(binary: np.ndarray) -> np.ndarray:
    count, labels, stats, centroids = cv2.connectedComponentsWithStats(binary.astype(np.uint8), 8)
    if count <= 1:
        return np.zeros_like(binary, dtype=bool)
    height, width = binary.shape
    centre = np.array([(width - 1) / 2, (height - 1) / 2], dtype=float)
    central = labels[int(round(centre[1])), int(round(centre[0]))]
    if central > 0:
        return labels == central
    best_label = 0
    best_score = -np.inf
    for label in range(1, count):
        area = float(stats[label, cv2.CC_STAT_AREA])
        distance = np.linalg.norm((centroids[label] - centre) / np.array([max(width, 1), max(height, 1)]))
        score = math.log1p(area) - 5.0 * distance
        if score > best_score:
            best_score = score
            best_label = label
    return labels == best_label


def segment_head(image: np.ndarray) -> tuple[np.ndarray | None, dict[str, float]]:
    height, width = image.shape[:2]
    scale = min(1.0, 640.0 / max(height, width))
    resized = cv2.resize(image, None, fx=scale, fy=scale, interpolation=cv2.INTER_AREA) if scale < 1 else image.copy()
    lab = cv2.cvtColor(resized, cv2.COLOR_BGR2LAB).astype(np.float32)
    border = np.concatenate((lab[0], lab[-1], lab[:, 0], lab[:, -1]), axis=0)
    reference = np.median(border, axis=0)
    distance = np.linalg.norm(lab - reference, axis=2)
    normalized = np.clip(distance / max(np.percentile(distance, 99), 1.0) * 255, 0, 255).astype(np.uint8)
    threshold, binary = cv2.threshold(normalized, 0, 1, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    kernel_size = max(3, int(round(min(resized.shape[:2]) * 0.018)) | 1)
    kernel = np.ones((kernel_size, kernel_size), np.uint8)
    binary = cv2.morphologyEx(binary, cv2.MORPH_CLOSE, kernel)
    binary = cv2.morphologyEx(binary, cv2.MORPH_OPEN, np.ones((3, 3), np.uint8))
    foreground = _largest_central_component(binary.astype(bool))
    area_fraction = float(foreground.mean())
    if foreground.sum() < 100 or not 0.06 <= area_fraction <= 0.92:
        return None, {"mask_quality": 0.0, "mask_area_fraction": area_fraction, "segmentation_threshold": float(threshold)}
    border_touch = float(np.concatenate((foreground[0], foreground[-1], foreground[:, 0], foreground[:, -1])).mean())
    centre_window = foreground[
        int(foreground.shape[0] * 0.25):int(foreground.shape[0] * 0.75),
        int(foreground.shape[1] * 0.25):int(foreground.shape[1] * 0.75),
    ]
    central_coverage = float(centre_window.mean()) if centre_window.size else 0.0
    quality = float(np.clip(0.55 * central_coverage + 0.30 * (1 - border_touch) + 0.15 * min(area_fraction / 0.25, 1), 0, 1))
    return foreground, {
        "mask_quality": quality,
        "mask_area_fraction": area_fraction,
        "mask_border_touch": border_touch,
        "segmentation_threshold": float(threshold),
        "resize_scale": float(scale),
    }


def floral_mask(image: np.ndarray, foreground: np.ndarray) -> np.ndarray:
    if image.shape[:2] != foreground.shape:
        image = cv2.resize(image, (foreground.shape[1], foreground.shape[0]), interpolation=cv2.INTER_AREA)
    hsv = cv2.cvtColor(image, cv2.COLOR_BGR2HSV)
    hue, saturation, value = [hsv[..., index] for index in range(3)]
    purple = foreground & (value >= 80) & (saturation >= 35) & (hue >= 105) & (hue <= 164)
    redmagenta = foreground & (value >= 80) & (saturation >= 25) & ((hue >= 165) | (hue <= 10))
    yellow = foreground & (value >= 115) & (saturation >= 45) & (hue >= 11) & (hue <= 38)
    occupied = purple | redmagenta | yellow
    gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    laplacian = np.abs(cv2.Laplacian(gray, cv2.CV_32F))
    edge_threshold = float(np.percentile(laplacian[foreground], 45)) if foreground.any() else 0.0
    white = foreground & ~occupied & (value >= 170) & (saturation <= 58) & (laplacian >= edge_threshold)
    return purple | redmagenta | yellow | white


def _profile_metrics(
    contour: np.ndarray,
    centre: np.ndarray,
    axis: np.ndarray,
    equivalent_radius: float,
    bins: int = 120,
) -> dict[str, float] | None:
    points = contour[:, 0, :].astype(float)
    projection = (points - centre) @ axis
    involucre = points[projection <= np.quantile(projection, 0.72)]
    if len(involucre) < 30:
        return None
    vectors = involucre - centre
    angles = (np.arctan2(vectors[:, 1], vectors[:, 0]) + 2 * np.pi) % (2 * np.pi)
    radii = np.linalg.norm(vectors, axis=1)
    identifiers = np.floor(angles / (2 * np.pi) * bins).astype(int) % bins
    profile = np.full(bins, np.nan)
    for identifier in np.unique(identifiers):
        profile[identifier] = radii[identifiers == identifier].max()
    valid = np.isfinite(profile)
    if valid.sum() < 25:
        return None
    positions = np.arange(bins)
    vp = positions[valid]
    vv = profile[valid]
    interpolated = np.interp(positions, np.r_[vp - bins, vp, vp + bins], np.r_[vv, vv, vv])
    kernel = np.ones(11, dtype=float) / 11
    smooth = np.convolve(np.r_[interpolated[-5:], interpolated, interpolated[:5]], kernel, mode="valid")
    residual = interpolated - smooth
    positive = np.maximum(residual, 0)
    threshold = max(0.045 * equivalent_radius, float(np.percentile(positive[valid], 85)))
    peaks = [
        index
        for index in range(bins)
        if valid[index]
        and residual[index] > threshold
        and residual[index] >= residual[(index - 1) % bins]
        and residual[index] >= residual[(index + 1) % bins]
    ]
    bin_angles = (positions + 0.5) / bins * 2 * np.pi
    bin_vectors = np.column_stack((np.cos(bin_angles), np.sin(bin_angles)))
    perpendicular = np.array([-axis[1], axis[0]])
    left = valid & ((bin_vectors @ perpendicular) >= 0)
    right = valid & ~left
    left_mean = float(positive[left].mean()) if left.any() else 0.0
    right_mean = float(positive[right].mean()) if right.any() else 0.0
    return {
        "involucre_projection_roughness": float(np.std(residual[valid]) / max(equivalent_radius, 1.0)),
        "involucre_projection_p95": float(np.percentile(positive[valid], 95) / max(equivalent_radius, 1.0)),
        "involucre_projection_max": float(np.max(positive[valid]) / max(equivalent_radius, 1.0)),
        "involucre_spread_fraction": float(np.mean(positive[valid] > 0.04 * equivalent_radius)),
        "bract_projection_peak_density": float(len(peaks) / valid.sum()),
        "bract_projection_asymmetry": float(abs(left_mean - right_mean) / max(equivalent_radius, 1.0)),
        "involucre_contour_valid_bins": int(valid.sum()),
    }


def _band_width(longitudinal: np.ndarray, transverse: np.ndarray, low: float, high: float) -> float:
    lo, hi = np.quantile(longitudinal, [low, high])
    selected = transverse[(longitudinal >= lo) & (longitudinal <= hi)]
    if len(selected) < 10:
        return float("nan")
    return float(np.percentile(selected, 95) - np.percentile(selected, 5))


def _lbp_entropy(gray: np.ndarray, mask: np.ndarray) -> float:
    if gray.shape[0] < 3 or gray.shape[1] < 3:
        return float("nan")
    centre = gray[1:-1, 1:-1]
    neighbours = (
        gray[:-2, :-2], gray[:-2, 1:-1], gray[:-2, 2:], gray[1:-1, 2:],
        gray[2:, 2:], gray[2:, 1:-1], gray[2:, :-2], gray[1:-1, :-2],
    )
    code = np.zeros_like(centre, dtype=np.uint8)
    for bit, neighbour in enumerate(neighbours):
        code |= ((neighbour >= centre).astype(np.uint8) << bit)
    selected = code[mask[1:-1, 1:-1]]
    if len(selected) < 50:
        return float("nan")
    counts = np.bincount(selected, minlength=256).astype(float)
    probability = counts[counts > 0] / counts.sum()
    return float(-(probability * np.log(probability)).sum() / np.log(256))


def measure_once(image: np.ndarray, args: argparse.Namespace) -> dict[str, Any]:
    original_height, original_width = image.shape[:2]
    min_dimension = int(min(original_height, original_width))
    gray_original = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    sharpness = float(cv2.Laplacian(gray_original, cv2.CV_64F).var())
    foreground, quality = segment_head(image)
    base: dict[str, Any] = {"min_dimension": min_dimension, "sharpness": sharpness, **quality}
    if foreground is None:
        return {**base, "architecture_status": "mask_failed", "surface_status": "mask_failed"}

    resized = cv2.resize(image, (foreground.shape[1], foreground.shape[0]), interpolation=cv2.INTER_AREA)
    floral = floral_mask(resized, foreground)
    rows, columns = np.nonzero(foreground)
    floral_rows, floral_columns = np.nonzero(floral)
    floral_fraction = float(floral.sum() / max(foreground.sum(), 1))
    base["visible_floret_fraction_extended"] = floral_fraction
    if len(columns) < 100 or len(floral_columns) < max(12, int(0.02 * len(columns))):
        return {**base, "architecture_status": "floral_end_not_recovered", "surface_status": "floral_end_not_recovered"}

    coordinates = np.column_stack((columns, rows)).astype(np.float32)
    mean, eigenvectors, _ = cv2.PCACompute2(coordinates, mean=None, maxComponents=2)
    centre = mean[0].astype(float)
    axis = eigenvectors[0].astype(float)
    axis /= max(np.linalg.norm(axis), 1e-9)
    floral_centroid = np.array([np.median(floral_columns), np.median(floral_rows)])
    if np.dot(floral_centroid - centre, axis) < 0:
        axis = -axis
    perpendicular = np.array([-axis[1], axis[0]])
    longitudinal_all = (coordinates - centre) @ axis
    cutoff = np.quantile(longitudinal_all, 0.72)
    involucre_mask = np.zeros_like(foreground, dtype=bool)
    keep = longitudinal_all <= cutoff
    involucre_mask[rows[keep], columns[keep]] = True
    involucre_mask &= ~cv2.dilate(floral.astype(np.uint8), np.ones((3, 3), np.uint8), iterations=1).astype(bool)
    if involucre_mask.sum() < 100:
        return {**base, "architecture_status": "involucre_region_not_recovered", "surface_status": "involucre_region_not_recovered"}

    irows, icols = np.nonzero(involucre_mask)
    icoordinates = np.column_stack((icols, irows)).astype(float)
    longitudinal = (icoordinates - centre) @ axis
    transverse = (icoordinates - centre) @ perpendicular
    length = float(np.percentile(longitudinal, 97.5) - np.percentile(longitudinal, 2.5))
    width = float(np.percentile(transverse, 97.5) - np.percentile(transverse, 2.5))
    basal_width = _band_width(longitudinal, transverse, 0.05, 0.25)
    middle_width = _band_width(longitudinal, transverse, 0.35, 0.65)
    apical_width = _band_width(longitudinal, transverse, 0.75, 0.95)

    contours, _ = cv2.findContours(involucre_mask.astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    if not contours:
        return {**base, "architecture_status": "involucre_contour_not_recovered", "surface_status": "involucre_contour_not_recovered"}
    contour = max(contours, key=cv2.contourArea)
    equivalent_radius = math.sqrt(float(foreground.sum()) / math.pi)
    profile = _profile_metrics(contour, centre, axis, equivalent_radius)
    if profile is None:
        return {**base, "architecture_status": "involucre_contour_not_recovered", "surface_status": "involucre_contour_not_recovered"}

    architecture_reasons: list[str] = []
    if min_dimension < args.min_architecture_dimension:
        architecture_reasons.append("low_resolution")
    if sharpness < args.min_architecture_sharpness:
        architecture_reasons.append("low_sharpness")
    if quality["mask_quality"] < 0.30 or quality["mask_border_touch"] > 0.20:
        architecture_reasons.append("low_mask_quality")

    result: dict[str, Any] = {
        **base,
        **profile,
        "involucre_length_width_ratio": length / max(width, 1.0),
        "involucre_apical_taper_ratio": apical_width / max(middle_width, 1.0),
        "involucre_basal_taper_ratio": basal_width / max(middle_width, 1.0),
        "architecture_status": "usable" if not architecture_reasons else ";".join(sorted(set(architecture_reasons))),
    }

    surface_reasons = list(architecture_reasons)
    if min_dimension < args.min_surface_dimension:
        surface_reasons.append("low_surface_resolution")
    if sharpness < args.min_surface_sharpness:
        surface_reasons.append("low_surface_sharpness")
    erosion_size = max(3, int(round(min(foreground.shape) * 0.025)) | 1)
    surface = cv2.erode(involucre_mask.astype(np.uint8), np.ones((erosion_size, erosion_size), np.uint8)).astype(bool)
    if surface.sum() < 100:
        surface_reasons.append("surface_region_too_small")
        result["surface_status"] = ";".join(sorted(set(surface_reasons)))
        return result

    gray = cv2.cvtColor(resized, cv2.COLOR_BGR2GRAY)
    median = float(np.median(gray[surface]))
    low = int(max(0, 0.66 * median))
    high = int(min(255, max(low + 1, 1.33 * median)))
    edges = cv2.Canny(gray, low, high) > 0
    laplacian = np.abs(cv2.Laplacian(gray, cv2.CV_32F))
    q25, q75 = np.percentile(gray[surface], [25, 75])
    hsv = cv2.cvtColor(resized, cv2.COLOR_BGR2HSV)
    saturation, value = hsv[..., 1], hsv[..., 2]
    value_threshold = max(210.0, float(np.percentile(value[surface], 95)))
    result.update({
        "involucre_surface_edge_density": float(edges[surface].mean()),
        "involucre_surface_lbp_entropy": _lbp_entropy(gray, surface),
        "involucre_surface_high_frequency_energy": float(np.median(laplacian[surface]) / max(q75 - q25, 5.0)),
        "involucre_surface_specular_fraction_proxy": float(((saturation <= 45) & (value >= value_threshold) & surface).sum() / surface.sum()),
        "surface_status": "usable" if not surface_reasons else ";".join(sorted(set(surface_reasons))),
    })
    return result


def combine(original: dict[str, Any], mirrored: dict[str, Any]) -> dict[str, Any]:
    output: dict[str, Any] = {
        "min_dimension": original.get("min_dimension"),
        "sharpness": original.get("sharpness"),
        "mask_quality": original.get("mask_quality"),
        "visible_floret_fraction_extended": np.nanmean([
            original.get("visible_floret_fraction_extended", np.nan),
            mirrored.get("visible_floret_fraction_extended", np.nan),
        ]),
    }
    for group, metrics in (("architecture", ARCHITECTURE_METRICS), ("surface", SURFACE_METRICS)):
        group_reasons: list[str] = []
        for side in (original, mirrored):
            if side.get(f"{group}_status") != "usable":
                group_reasons.extend(str(side.get(f"{group}_status", "failed")).split(";"))
        metric_statuses: list[str] = []
        for metric in metrics:
            reasons = list(group_reasons)
            left = float(original.get(metric, np.nan))
            right = float(mirrored.get(metric, np.nan))
            if np.isfinite(left) and np.isfinite(right):
                output[metric] = (left + right) / 2
                output[f"{metric}_flip_delta"] = abs(left - right)
                if abs(left - right) > FLIP_TOLERANCE[metric]:
                    reasons.append("unstable_under_horizontal_flip")
            else:
                output[metric] = np.nan
                output[f"{metric}_flip_delta"] = np.nan
                reasons.append("measurement_missing")
            metric_status = "usable" if not reasons else ";".join(sorted(set(filter(None, reasons))))
            output[f"{metric}_status"] = metric_status
            metric_statuses.append(metric_status)
        output[f"{group}_status"] = (
            "usable" if all(status == "usable" for status in metric_statuses)
            else "one_or_more_endpoint_failures"
        )
    return output


def _mad(values: pd.Series) -> float:
    numeric = pd.to_numeric(values, errors="coerce").dropna().to_numpy(float)
    if not len(numeric):
        return float("nan")
    median = np.median(numeric)
    return float(np.median(np.abs(numeric - median)))


def aggregate(table: pd.DataFrame, group_column: str, suffix: str) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for key, group in table.groupby(group_column, dropna=False, sort=True):
        row: dict[str, Any] = {group_column: key}
        if group_column == "obs_id" and "taxon_name" in group:
            taxa = group["taxon_name"].dropna().astype(str).unique()
            row["taxon_name"] = taxa[0] if len(taxa) == 1 else ""
        for _, metrics in (("architecture", ARCHITECTURE_METRICS), ("surface", SURFACE_METRICS)):
            for metric in metrics:
                usable = group[group[f"{metric}_status"].eq("usable")]
                row[f"n_usable_heads_{metric}"] = int(usable["annotation_unit_id"].nunique())
                values = pd.to_numeric(usable[metric], errors="coerce")
                row[f"{metric}_{suffix}_median"] = float(values.median()) if values.notna().any() else np.nan
                row[f"{metric}_{suffix}_mad"] = _mad(values)
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    manifest = pd.read_csv(args.manifest, dtype=str, keep_default_na=False)
    if args.metadata:
        metadata = pd.read_csv(args.metadata, dtype=str, keep_default_na=False)
        if "annotation_unit_id" not in metadata:
            required_metadata = {"queue_id", "det_index"}
            missing_metadata = sorted(required_metadata.difference(metadata.columns))
            if missing_metadata:
                raise ValueError(f"Metadata cannot derive annotation units; missing={missing_metadata}")
            metadata["annotation_unit_id"] = [
                f"{queue}_head_{int(float(index)) + 1:02d}"
                for queue, index in zip(metadata["queue_id"], metadata["det_index"])
            ]
        payload = metadata[[
            column for column in ("annotation_unit_id", "obs_id", "taxon_name", "photo_id", "latitude", "longitude", "positional_accuracy")
            if column in metadata
        ]].drop_duplicates("annotation_unit_id")
        manifest = manifest.merge(payload, on="annotation_unit_id", how="left", validate="one_to_one")
    required = {"annotation_unit_id", "crop_path", "obs_id", "taxon_name"}
    missing = sorted(required.difference(manifest.columns))
    if missing:
        raise ValueError(f"Manifest is missing required columns: {missing}")
    if manifest["annotation_unit_id"].duplicated().any():
        raise ValueError("Manifest annotation_unit_id must be unique")
    root = Path(args.image_root)
    rows: list[dict[str, Any]] = []
    for record in manifest.to_dict("records"):
        image = cv2.imread(str(root / record["crop_path"]))
        if image is None:
            measured = {
                "architecture_status": "image_unreadable",
                "surface_status": "image_unreadable",
                **{
                    metric: np.nan
                    for metric in (*ARCHITECTURE_METRICS, *SURFACE_METRICS)
                },
                **{
                    f"{metric}_status": "image_unreadable"
                    for metric in (*ARCHITECTURE_METRICS, *SURFACE_METRICS)
                },
            }
        else:
            measured = combine(measure_once(image, args), measure_once(cv2.flip(image, 1), args))
        rows.append({**record, **measured})
    head = pd.DataFrame(rows).sort_values("annotation_unit_id").reset_index(drop=True)
    observation = aggregate(head, "obs_id", "observation")
    species = aggregate(head, "taxon_name", "species")
    output = Path(args.out_dir)
    output.mkdir(parents=True, exist_ok=True)
    head.to_csv(output / "extended_continuous_head_level.csv", index=False, encoding="utf-8-sig")
    observation.to_csv(output / "extended_continuous_observation_level.csv", index=False, encoding="utf-8-sig")
    species.to_csv(output / "extended_continuous_species_level.csv", index=False, encoding="utf-8-sig")
    report = {
        "n_heads": int(len(head)),
        "n_observations": int(head["obs_id"].nunique()),
        "n_taxa": int(head["taxon_name"].nunique()),
        "n_architecture_usable": int(head["architecture_status"].eq("usable").sum()),
        "n_surface_usable": int(head["surface_status"].eq("usable").sum()),
        "endpoint_usable_heads": {
            metric: int(head[f"{metric}_status"].eq("usable").sum())
            for metric in (*ARCHITECTURE_METRICS, *SURFACE_METRICS)
        },
        "semantic_status": (
            "Category-free continuous image measurements. Architecture metrics are contour-derived; "
            "surface metrics are validation-only proxies until botanical reference annotations exist."
        ),
        "forbidden_claims": [
            "projection maximum is botanical spine length",
            "edge density is hair density",
            "specular fraction is mucilage presence or secretion amount",
        ],
    }
    (output / "extended_continuous_provenance.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
