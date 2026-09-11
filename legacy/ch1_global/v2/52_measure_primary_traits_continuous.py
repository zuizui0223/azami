#!/usr/bin/env python3
"""Continuous automated measurements for three whole-capitulum traits.

The script measures colour, outline geometry and head-peduncle angle directly.
It does not train on human category labels and does not force each image into a
morphological class. Horizontal-flip repeatability and visual-quality fields are
exported so downstream analyses can retain only stable measurements.
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

NON_SCOREABLE = "unassessable"
PRIMARY_TRAITS = ("capitulum_orientation", "corolla_colour_class", "capitulum_shape")


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Measure continuous whole-capitulum traits from detected-head images.")
    p.add_argument("--packet-manifest", required=True)
    p.add_argument("--packet-root", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--colour-confidence-floor", type=float, default=0.55)
    p.add_argument("--shape-confidence-floor", type=float, default=0.50)
    p.add_argument("--orientation-confidence-floor", type=float, default=0.65)
    p.add_argument("--flip-angle-tolerance", type=float, default=20.0)
    return p.parse_args()


def packet_path(value: Any, root: Path) -> Path:
    relative = Path(text(value))
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError(f"Unsafe packet path: {relative}")
    return root / relative


def read_bgr(path: Path) -> np.ndarray:
    image = cv2.imread(str(path), cv2.IMREAD_COLOR)
    if image is None:
        raise ValueError(f"Unreadable image: {path}")
    return image


def largest_central_component(mask: np.ndarray) -> np.ndarray:
    binary = (mask > 0).astype(np.uint8)
    n_labels, labels, stats, centroids = cv2.connectedComponentsWithStats(binary, 8)
    if n_labels <= 1:
        return np.zeros_like(binary)
    height, width = binary.shape
    centre = np.array([width / 2, height / 2], dtype=float)
    best_label = None
    best_score = -1e18
    central_window = labels[
        int(height * 0.25):int(height * 0.75),
        int(width * 0.25):int(width * 0.75),
    ]
    for label in range(1, n_labels):
        area = float(stats[label, cv2.CC_STAT_AREA])
        centroid = np.array(centroids[label], dtype=float)
        distance = np.linalg.norm(
            (centroid - centre) / np.array([max(width, 1), max(height, 1)])
        )
        central_overlap = float(np.count_nonzero(central_window == label))
        score = math.log1p(area) + 1.5 * math.log1p(central_overlap) - 4.0 * distance
        if score > best_score:
            best_score = score
            best_label = label
    return (labels == best_label).astype(np.uint8)


def central_foreground(image: np.ndarray) -> tuple[np.ndarray, dict[str, float]]:
    """Fast centre-biased foreground estimate from border-colour contrast."""
    height, width = image.shape[:2]
    if min(height, width) < 20:
        return np.zeros((height, width), np.uint8), {"mask_quality": 0.0}
    lab = cv2.cvtColor(image, cv2.COLOR_BGR2LAB).astype(np.float32)
    border_width = max(2, int(round(min(height, width) * 0.08)))
    border_mask = np.zeros((height, width), np.uint8)
    border_mask[:border_width, :] = 1
    border_mask[-border_width:, :] = 1
    border_mask[:, :border_width] = 1
    border_mask[:, -border_width:] = 1
    border_pixels = lab[border_mask.astype(bool)]
    background = np.median(border_pixels, axis=0)
    distance = np.linalg.norm(lab - background, axis=2)
    scaled = cv2.normalize(distance, None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
    otsu_threshold, binary = cv2.threshold(
        scaled, 0, 1, cv2.THRESH_BINARY + cv2.THRESH_OTSU
    )

    yy, xx = np.mgrid[0:height, 0:width]
    centre_x, centre_y = (width - 1) / 2, (height - 1) / 2
    radius_x, radius_y = max(width * 0.46, 1), max(height * 0.46, 1)
    centre = (
        ((xx - centre_x) / radius_x) ** 2
        + ((yy - centre_y) / radius_y) ** 2
        <= 1
    )
    lower = max(6.0, float(np.percentile(distance[centre], 42)))
    binary = ((binary > 0) | (centre & (distance >= lower))).astype(np.uint8)

    kernel_size = max(3, int(round(min(height, width) * 0.028)) | 1)
    kernel = cv2.getStructuringElement(
        cv2.MORPH_ELLIPSE, (kernel_size, kernel_size)
    )
    binary = cv2.morphologyEx(binary, cv2.MORPH_CLOSE, kernel)
    binary = cv2.morphologyEx(binary, cv2.MORPH_OPEN, kernel)
    foreground = largest_central_component(binary)

    area_fraction = float(foreground.mean())
    border = np.concatenate([
        foreground[0], foreground[-1], foreground[:, 0], foreground[:, -1]
    ])
    border_fraction = float(border.mean())
    area_score = max(0.0, 1.0 - abs(area_fraction - 0.42) / 0.42)
    border_score = max(0.0, 1.0 - border_fraction * 2.5)
    if foreground.any():
        contrast_score = float(np.clip(
            (
                np.median(distance[foreground.astype(bool)])
                - np.median(distance[border_mask.astype(bool)])
            ) / 35,
            0,
            1,
        ))
    else:
        contrast_score = 0.0
    mask_quality = float(np.clip(
        0.45 * area_score + 0.35 * border_score + 0.20 * contrast_score,
        0,
        1,
    ))
    return foreground, {
        "mask_quality": mask_quality,
        "foreground_area_fraction": area_fraction,
        "foreground_border_fraction": border_fraction,
        "foreground_otsu_threshold": float(otsu_threshold),
    }


def colour_measurement(image: np.ndarray) -> dict[str, Any]:
    foreground, quality = central_foreground(image)
    hsv = cv2.cvtColor(image, cv2.COLOR_BGR2HSV)
    lab = cv2.cvtColor(image, cv2.COLOR_BGR2LAB)
    hue, saturation, value = [hsv[..., index] for index in range(3)]
    lightness, a_axis, b_axis = [
        lab[..., index].astype(np.float32) for index in range(3)
    ]
    yy, xx = np.mgrid[0:image.shape[0], 0:image.shape[1]]
    centre_x, centre_y = (image.shape[1] - 1) / 2, (image.shape[0] - 1) / 2
    radius_x = max(image.shape[1] * 0.52, 1)
    radius_y = max(image.shape[0] * 0.52, 1)
    centre = (
        ((xx - centre_x) / radius_x) ** 2
        + ((yy - centre_y) / radius_y) ** 2
        <= 1
    ).astype(np.uint8)
    base = (foreground & centre).astype(bool)
    if base.sum() < 20:
        return {"state": NON_SCOREABLE, "confidence": 0.0, **quality}

    bright = value >= 80
    purple = base & bright & (saturation >= 35) & (hue >= 105) & (hue <= 169)
    redmagenta = base & bright & (saturation >= 25) & ((hue >= 165) | (hue <= 10))
    yellow = base & (value >= 115) & (saturation >= 45) & (hue >= 11) & (hue <= 38)
    gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    laplacian = np.abs(cv2.Laplacian(gray, cv2.CV_32F))
    edge_threshold = float(np.percentile(laplacian[base], 45)) if base.any() else 0.0
    white = base & (value >= 170) & (saturation <= 58) & (laplacian >= edge_threshold)

    raw_counts = {
        "purple_blue": int(purple.sum()),
        "redmagenta": int(redmagenta.sum()),
        "yellow_or_other": int(yellow.sum()),
        "white_or_cream": int(white.sum()),
    }
    total = sum(raw_counts.values())
    coverage = total / max(int(base.sum()), 1)
    if total < max(12, int(base.sum() * 0.025)) or coverage < 0.035:
        return {
            "state": NON_SCOREABLE,
            "confidence": float(np.clip(
                quality["mask_quality"] * coverage / 0.06, 0, 0.45
            )),
            "floral_pixel_fraction": coverage,
            "colour_counts_json": json.dumps(raw_counts, sort_keys=True),
            **quality,
        }

    proportions = {name: count / total for name, count in raw_counts.items()}
    floral = purple | redmagenta | yellow | white
    chromatic = floral & (saturation >= 25)
    median_lightness = (
        float(np.median(lightness[floral]) * 100.0 / 255.0)
        if floral.any() else float("nan")
    )
    if chromatic.any():
        hue_radians = np.deg2rad(
            hue[chromatic].astype(np.float32) * 2.0
        )
        mean_hue = float((
            np.rad2deg(np.arctan2(
                np.sin(hue_radians).mean(), np.cos(hue_radians).mean()
            )) + 360.0
        ) % 360.0)
    else:
        mean_hue = float("nan")

    ordered = sorted(proportions.items(), key=lambda item: (-item[1], item[0]))
    top_name, top_fraction = ordered[0]
    second_fraction = ordered[1][1]
    selected = redmagenta if top_name == "redmagenta" else base
    median_saturation = float(np.median(saturation[selected]))
    median_chroma = float(np.median(np.sqrt(
        (a_axis[selected] - 128) ** 2 + (b_axis[selected] - 128) ** 2
    )))
    dominance = max(0.0, top_fraction - second_fraction)
    exposure = float(np.clip((np.median(value[base]) - 55) / 140, 0, 1))
    coverage_score = float(np.clip((coverage - 0.025) / 0.18, 0, 1))
    confidence = float(np.clip(
        0.32 * quality["mask_quality"]
        + 0.30 * dominance
        + 0.23 * coverage_score
        + 0.15 * exposure,
        0,
        1,
    ))
    return {
        "state": top_name,
        "confidence": confidence,
        "floral_pixel_fraction": coverage,
        "dominant_colour_fraction": top_fraction,
        "second_colour_fraction": second_fraction,
        "median_saturation": median_saturation,
        "median_lab_chroma": median_chroma,
        "median_lab_lightness": median_lightness,
        "mean_hue_degrees": mean_hue,
        "white_fraction": proportions["white_or_cream"],
        "redmagenta_fraction": proportions["redmagenta"],
        "purple_fraction": proportions["purple_blue"],
        "yellow_fraction": proportions["yellow_or_other"],
        "colour_counts_json": json.dumps(raw_counts, sort_keys=True),
        **quality,
    }


def rotate_mask_to_major_axis(mask: np.ndarray) -> tuple[np.ndarray, float]:
    rows, columns = np.nonzero(mask)
    if len(columns) < 20:
        return mask, 0.0
    coordinates = np.column_stack([columns, rows]).astype(np.float32)
    _, eigenvectors = cv2.PCACompute(coordinates, mean=None, maxComponents=2)
    vector = eigenvectors[0]
    angle = math.degrees(math.atan2(float(vector[1]), float(vector[0])))
    rotation = 90.0 - angle
    height, width = mask.shape
    matrix = cv2.getRotationMatrix2D((width / 2, height / 2), rotation, 1.0)
    rotated = cv2.warpAffine(
        mask.astype(np.uint8), matrix, (width, height), flags=cv2.INTER_NEAREST
    )
    return rotated, rotation


def shape_measurement(image: np.ndarray) -> dict[str, Any]:
    foreground, quality = central_foreground(image)
    height, width = foreground.shape
    kernel_size = max(3, int(round(min(height, width) * 0.055)) | 1)
    kernel = cv2.getStructuringElement(
        cv2.MORPH_ELLIPSE, (kernel_size, kernel_size)
    )
    body = cv2.morphologyEx(foreground, cv2.MORPH_OPEN, kernel)
    body = cv2.morphologyEx(body, cv2.MORPH_CLOSE, kernel)
    body = largest_central_component(body)
    if body.sum() < 40:
        return {"state": NON_SCOREABLE, "confidence": 0.0, **quality}

    rotated, rotation = rotate_mask_to_major_axis(body)
    rows, columns = np.nonzero(rotated)
    x0, x1, y0, y1 = columns.min(), columns.max(), rows.min(), rows.max()
    crop = rotated[y0:y1 + 1, x0:x1 + 1]
    crop_height, crop_width = crop.shape
    if crop_height < 8 or crop_width < 8:
        return {"state": NON_SCOREABLE, "confidence": 0.0, **quality}
    aspect_ratio = crop_height / crop_width
    widths = crop.sum(axis=1).astype(float)
    widths /= max(widths.max(), 1.0)
    n_rows = len(widths)

    def band(low: float, high: float) -> float:
        start = max(0, int(round(low * n_rows)))
        end = min(n_rows, max(start + 1, int(round(high * n_rows))))
        return float(np.median(widths[start:end]))

    end_a = band(0.05, 0.20)
    middle = band(0.38, 0.62)
    end_b = band(0.80, 0.95)
    narrow_end = min(end_a, end_b)
    width_values = widths[widths > 0.05]
    width_cv = float(
        np.std(width_values) / max(np.mean(width_values), 1e-6)
    )
    contours, _ = cv2.findContours(
        (crop * 255).astype(np.uint8),
        cv2.RETR_EXTERNAL,
        cv2.CHAIN_APPROX_SIMPLE,
    )
    contour = max(contours, key=cv2.contourArea)
    area = float(cv2.contourArea(contour))
    perimeter = float(cv2.arcLength(contour, True))
    circularity = float(4 * math.pi * area / max(perimeter * perimeter, 1e-6))
    hull = cv2.convexHull(contour)
    solidity = float(area / max(cv2.contourArea(hull), 1e-6))
    confidence = float(np.clip(
        0.45 * quality["mask_quality"]
        + 0.30 * np.clip(solidity, 0, 1)
        + 0.25 * np.clip(circularity, 0, 1),
        0,
        1,
    ))
    return {
        "state": "measured",
        "confidence": confidence,
        "aspect_ratio": float(aspect_ratio),
        "circularity": circularity,
        "solidity": solidity,
        "width_cv": width_cv,
        "end_width_a": end_a,
        "middle_width": middle,
        "end_width_b": end_b,
        "narrow_end_width": narrow_end,
        "rotation_degrees": float(rotation),
        **quality,
    }


def template_head_box(
    head: np.ndarray, context: np.ndarray
) -> tuple[tuple[int, int, int, int] | None, float]:
    if head.shape[0] > context.shape[0] or head.shape[1] > context.shape[1]:
        return None, 0.0
    result = cv2.matchTemplate(context, head, cv2.TM_CCOEFF_NORMED)
    _, score, _, location = cv2.minMaxLoc(result)
    x_coordinate, y_coordinate = location
    return (
        x_coordinate,
        y_coordinate,
        x_coordinate + head.shape[1],
        y_coordinate + head.shape[0],
    ), float(score)


def point_distance_to_box(
    point: tuple[float, float], box: tuple[int, int, int, int]
) -> float:
    x_coordinate, y_coordinate = point
    x0, y0, x1, y1 = box
    delta_x = max(x0 - x_coordinate, 0, x_coordinate - x1)
    delta_y = max(y0 - y_coordinate, 0, y_coordinate - y1)
    return math.hypot(delta_x, delta_y)


def orientation_measurement(
    head: np.ndarray, context: np.ndarray
) -> dict[str, Any]:
    max_side = max(context.shape[:2])
    if max_side > 320:
        scale = 320.0 / max_side
        context = cv2.resize(
            context, None, fx=scale, fy=scale, interpolation=cv2.INTER_AREA
        )
        head = cv2.resize(
            head, None, fx=scale, fy=scale, interpolation=cv2.INTER_AREA
        )
    box, template_score = template_head_box(head, context)
    if box is None or template_score < 0.35:
        return {
            "state": NON_SCOREABLE,
            "confidence": 0.0,
            "template_score": template_score,
        }
    x0, y0, x1, y1 = box
    centre = np.array([(x0 + x1) / 2, (y0 + y1) / 2], dtype=float)
    gray = cv2.cvtColor(context, cv2.COLOR_BGR2GRAY)
    blurred = cv2.GaussianBlur(gray, (5, 5), 0)
    edges = cv2.Canny(blurred, 40, 120)
    edges[max(0, y0):min(edges.shape[0], y1), max(0, x0):min(edges.shape[1], x1)] = 0
    minimum_length = max(12, int(round(min(head.shape[:2]) * 0.35)))
    lines = cv2.HoughLinesP(
        edges,
        1,
        np.pi / 180,
        threshold=max(12, minimum_length // 2),
        minLineLength=minimum_length,
        maxLineGap=max(5, minimum_length // 3),
    )
    if lines is None:
        return {
            "state": NON_SCOREABLE,
            "confidence": 0.0,
            "template_score": template_score,
        }

    best = None
    best_score = -1e9
    head_diagonal = math.hypot(x1 - x0, y1 - y0)
    for line in lines[:, 0]:
        point_a = np.array([line[0], line[1]], dtype=float)
        point_b = np.array([line[2], line[3]], dtype=float)
        length = float(np.linalg.norm(point_b - point_a))
        distance_a = point_distance_to_box(tuple(point_a), box)
        distance_b = point_distance_to_box(tuple(point_b), box)
        near, far = (
            (point_a, point_b) if distance_a <= distance_b else (point_b, point_a)
        )
        near_distance = min(distance_a, distance_b)
        far_distance = max(distance_a, distance_b)
        score = (
            1.2 * length / head_diagonal
            - 2.0 * near_distance / head_diagonal
            + 0.35 * far_distance / head_diagonal
        )
        direction = near - far
        to_centre = centre - near
        if np.linalg.norm(direction) and np.linalg.norm(to_centre):
            alignment = float(np.dot(
                direction / np.linalg.norm(direction),
                to_centre / np.linalg.norm(to_centre),
            ))
            score += 0.8 * alignment
        if score > best_score:
            best_score = score
            best = (near, far, length)
    if best is None:
        return {
            "state": NON_SCOREABLE,
            "confidence": 0.0,
            "template_score": template_score,
        }

    near, far, line_length = best
    stem_up = near - far
    if np.linalg.norm(stem_up) < 1:
        return {
            "state": NON_SCOREABLE,
            "confidence": 0.0,
            "template_score": template_score,
        }
    foreground, head_quality = central_foreground(head)
    rows, columns = np.nonzero(foreground)
    if len(columns) < 30:
        return {
            "state": NON_SCOREABLE,
            "confidence": 0.0,
            "template_score": template_score,
        }
    coordinates = np.column_stack([columns + x0, rows + y0]).astype(np.float32)
    mean, eigenvectors = cv2.PCACompute(coordinates, mean=None, maxComponents=2)
    axis = eigenvectors[0].astype(float)
    axis /= max(np.linalg.norm(axis), 1e-9)
    head_centre = mean[0].astype(float)
    endpoint_a = head_centre + axis * head_diagonal * 0.5
    endpoint_b = head_centre - axis * head_diagonal * 0.5
    floral_end = (
        endpoint_a
        if np.linalg.norm(endpoint_a - near) >= np.linalg.norm(endpoint_b - near)
        else endpoint_b
    )
    head_axis = floral_end - near
    if np.linalg.norm(head_axis) < 1:
        return {
            "state": NON_SCOREABLE,
            "confidence": 0.0,
            "template_score": template_score,
        }
    stem_up /= np.linalg.norm(stem_up)
    head_axis /= np.linalg.norm(head_axis)
    cosine = float(np.clip(np.dot(stem_up, head_axis), -1, 1))
    angle = math.degrees(math.acos(cosine))
    boundaries = [22, 55, 105, 150]
    nearest_boundary = min(abs(angle - boundary) for boundary in boundaries)
    boundary_score = float(np.clip(nearest_boundary / 28, 0, 1))
    line_score = float(np.clip((best_score + 0.5) / 2.5, 0, 1))
    confidence = float(np.clip(
        0.25 * template_score
        + 0.30 * line_score
        + 0.25 * head_quality["mask_quality"]
        + 0.20 * boundary_score,
        0,
        1,
    ))
    return {
        "state": "measured",
        "confidence": confidence,
        "angle_degrees": float(angle),
        "template_score": template_score,
        "stem_line_score": float(best_score),
        "stem_line_length": float(line_length),
        **{f"head_{key}": value for key, value in head_quality.items()},
    }


def finite_mean(*values: Any) -> float:
    numbers = [
        float(value)
        for value in values
        if value is not None and np.isfinite(float(value))
    ]
    return float(np.mean(numbers)) if numbers else float("nan")


def circular_difference_degrees(a_value: Any, b_value: Any) -> float:
    try:
        a_value = float(a_value)
        b_value = float(b_value)
    except (TypeError, ValueError):
        return float("nan")
    if not np.isfinite(a_value) or not np.isfinite(b_value):
        return float("nan")
    return float(abs((a_value - b_value + 180.0) % 360.0 - 180.0))


def numeric(result: dict[str, Any], key: str) -> float:
    value = result.get(key, float("nan"))
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def status_colour(
    original: dict[str, Any], flipped: dict[str, Any], confidence_floor: float
) -> tuple[str, float]:
    confidence = finite_mean(original.get("confidence"), flipped.get("confidence"))
    delta_chroma = abs(
        numeric(original, "median_lab_chroma")
        - numeric(flipped, "median_lab_chroma")
    )
    delta_lightness = abs(
        numeric(original, "median_lab_lightness")
        - numeric(flipped, "median_lab_lightness")
    )
    delta_hue = circular_difference_degrees(
        original.get("mean_hue_degrees"), flipped.get("mean_hue_degrees")
    )
    coverage = finite_mean(
        original.get("floral_pixel_fraction"),
        flipped.get("floral_pixel_fraction"),
    )
    if not np.isfinite(coverage) or coverage < 0.035:
        return "low_visible_corolla_fraction", confidence
    if confidence < confidence_floor:
        return "low_colour_quality", confidence
    if (
        delta_chroma > 12
        or delta_lightness > 15
        or (np.isfinite(delta_hue) and delta_hue > 30)
    ):
        return "unstable_under_horizontal_flip", confidence
    return "usable", confidence


def status_shape(
    original: dict[str, Any], flipped: dict[str, Any], confidence_floor: float
) -> tuple[str, float]:
    confidence = finite_mean(original.get("confidence"), flipped.get("confidence"))
    delta_aspect = abs(
        numeric(original, "aspect_ratio") - numeric(flipped, "aspect_ratio")
    )
    delta_circularity = abs(
        numeric(original, "circularity") - numeric(flipped, "circularity")
    )
    delta_solidity = abs(
        numeric(original, "solidity") - numeric(flipped, "solidity")
    )
    if confidence < confidence_floor:
        return "low_outline_quality", confidence
    if delta_aspect > 0.14 or delta_circularity > 0.14 or delta_solidity > 0.12:
        return "unstable_under_horizontal_flip", confidence
    return "usable", confidence


def status_orientation(
    original: dict[str, Any],
    flipped: dict[str, Any],
    confidence_floor: float,
    tolerance: float,
) -> tuple[str, float]:
    confidence = finite_mean(original.get("confidence"), flipped.get("confidence"))
    angle_a = numeric(original, "angle_degrees")
    angle_b = numeric(flipped, "angle_degrees")
    if not np.isfinite(angle_a) or not np.isfinite(angle_b):
        return "head_peduncle_axis_not_recovered", confidence
    if confidence < confidence_floor:
        return "low_axis_quality", confidence
    if abs(angle_a - angle_b) > tolerance:
        return "unstable_under_horizontal_flip", confidence
    return "usable", confidence


def main() -> None:
    args = parse_args()
    root = Path(args.packet_root)
    manifest = pd.read_csv(args.packet_manifest, dtype=str, keep_default_na=False)
    required = {"annotation_unit_id", "crop_path", "context_crop_path"}
    missing = required.difference(manifest.columns)
    if missing or manifest.empty or manifest["annotation_unit_id"].duplicated().any():
        raise ValueError(
            f"Packet manifest must contain unique heads; missing={sorted(missing)}"
        )

    records: list[dict[str, Any]] = []
    for index, row in enumerate(manifest.to_dict("records"), start=1):
        unit = text(row["annotation_unit_id"])
        head = read_bgr(packet_path(row["crop_path"], root))
        context = read_bgr(packet_path(row["context_crop_path"], root))
        head_flip = cv2.flip(head, 1)
        context_flip = cv2.flip(context, 1)

        colour_a = colour_measurement(head)
        colour_b = colour_measurement(head_flip)
        shape_a = shape_measurement(head)
        shape_b = shape_measurement(head_flip)
        orientation_a = orientation_measurement(head, context)
        orientation_b = orientation_measurement(head_flip, context_flip)

        colour_status, colour_confidence = status_colour(
            colour_a, colour_b, args.colour_confidence_floor
        )
        shape_status, shape_confidence = status_shape(
            shape_a, shape_b, args.shape_confidence_floor
        )
        orientation_status, orientation_confidence = status_orientation(
            orientation_a,
            orientation_b,
            args.orientation_confidence_floor,
            args.flip_angle_tolerance,
        )

        hue_a = numeric(colour_a, "mean_hue_degrees")
        hue_b = numeric(colour_b, "mean_hue_degrees")
        if np.isfinite(hue_a) and np.isfinite(hue_b):
            hue_radians = np.deg2rad([hue_a, hue_b])
            hue_mean = float((
                np.rad2deg(np.arctan2(
                    np.sin(hue_radians).mean(), np.cos(hue_radians).mean()
                )) + 360.0
            ) % 360.0)
        else:
            hue_mean = finite_mean(hue_a, hue_b)

        records.append({
            "annotation_unit_id": unit,
            "colour_status": colour_status,
            "colour_confidence": colour_confidence,
            "corolla_visible_fraction": finite_mean(
                colour_a.get("floral_pixel_fraction"),
                colour_b.get("floral_pixel_fraction"),
            ),
            "corolla_lab_lightness": finite_mean(
                colour_a.get("median_lab_lightness"),
                colour_b.get("median_lab_lightness"),
            ),
            "corolla_lab_chroma": finite_mean(
                colour_a.get("median_lab_chroma"),
                colour_b.get("median_lab_chroma"),
            ),
            "corolla_hue_degrees": hue_mean,
            "corolla_hue_sin": (
                float(np.sin(np.deg2rad(hue_mean)))
                if np.isfinite(hue_mean) else float("nan")
            ),
            "corolla_hue_cos": (
                float(np.cos(np.deg2rad(hue_mean)))
                if np.isfinite(hue_mean) else float("nan")
            ),
            "corolla_white_fraction": finite_mean(
                colour_a.get("white_fraction"), colour_b.get("white_fraction")
            ),
            "corolla_redmagenta_fraction": finite_mean(
                colour_a.get("redmagenta_fraction"),
                colour_b.get("redmagenta_fraction"),
            ),
            "corolla_purple_fraction": finite_mean(
                colour_a.get("purple_fraction"), colour_b.get("purple_fraction")
            ),
            "corolla_yellow_fraction": finite_mean(
                colour_a.get("yellow_fraction"), colour_b.get("yellow_fraction")
            ),
            "colour_flip_delta_chroma": abs(
                numeric(colour_a, "median_lab_chroma")
                - numeric(colour_b, "median_lab_chroma")
            ),
            "colour_flip_delta_lightness": abs(
                numeric(colour_a, "median_lab_lightness")
                - numeric(colour_b, "median_lab_lightness")
            ),
            "colour_flip_delta_hue": circular_difference_degrees(hue_a, hue_b),
            "shape_status": shape_status,
            "shape_confidence": shape_confidence,
            "shape_aspect_ratio": finite_mean(
                shape_a.get("aspect_ratio"), shape_b.get("aspect_ratio")
            ),
            "shape_circularity": finite_mean(
                shape_a.get("circularity"), shape_b.get("circularity")
            ),
            "shape_solidity": finite_mean(
                shape_a.get("solidity"), shape_b.get("solidity")
            ),
            "shape_width_cv": finite_mean(
                shape_a.get("width_cv"), shape_b.get("width_cv")
            ),
            "shape_narrow_end_width": finite_mean(
                shape_a.get("narrow_end_width"),
                shape_b.get("narrow_end_width"),
            ),
            "shape_middle_width": finite_mean(
                shape_a.get("middle_width"), shape_b.get("middle_width")
            ),
            "shape_end_width_a": finite_mean(
                shape_a.get("end_width_a"), shape_b.get("end_width_a")
            ),
            "shape_end_width_b": finite_mean(
                shape_a.get("end_width_b"), shape_b.get("end_width_b")
            ),
            "shape_flip_delta_aspect": abs(
                numeric(shape_a, "aspect_ratio")
                - numeric(shape_b, "aspect_ratio")
            ),
            "shape_flip_delta_circularity": abs(
                numeric(shape_a, "circularity")
                - numeric(shape_b, "circularity")
            ),
            "shape_flip_delta_solidity": abs(
                numeric(shape_a, "solidity")
                - numeric(shape_b, "solidity")
            ),
            "orientation_status": orientation_status,
            "orientation_confidence": orientation_confidence,
            "orientation_angle_degrees": finite_mean(
                orientation_a.get("angle_degrees"),
                orientation_b.get("angle_degrees"),
            ),
            "orientation_flip_delta_degrees": abs(
                numeric(orientation_a, "angle_degrees")
                - numeric(orientation_b, "angle_degrees")
            ),
            "orientation_template_score": finite_mean(
                orientation_a.get("template_score"),
                orientation_b.get("template_score"),
            ),
            "orientation_stem_line_score": finite_mean(
                orientation_a.get("stem_line_score"),
                orientation_b.get("stem_line_score"),
            ),
            "orientation_head_mask_quality": finite_mean(
                orientation_a.get("head_mask_quality"),
                orientation_b.get("head_mask_quality"),
            ),
            "measurement_basis": (
                "detector_crop_plus_border_contrast_mask_plus_outline_geometry_"
                "plus_head_peduncle_axis_with_horizontal_flip_repeatability"
            ),
        })
        if index % 50 == 0 or index == len(manifest):
            print(f"[INFO] measured {index}/{len(manifest)} heads")

    output = Path(args.out_dir)
    output.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(records).sort_values("annotation_unit_id").reset_index(drop=True)
    frame.to_csv(
        output / "primary_trait_continuous_head_measurements.csv",
        index=False,
        encoding="utf-8-sig",
    )
    quality_rows = []
    for trait in ("colour", "shape", "orientation"):
        counts = frame[f"{trait}_status"].value_counts(dropna=False)
        for status, count in counts.items():
            quality_rows.append({
                "trait": trait,
                "status": status,
                "n_heads": int(count),
                "fraction": float(count / len(frame)),
            })
    quality = pd.DataFrame(quality_rows)
    quality.to_csv(
        output / "primary_trait_continuous_quality_summary.csv",
        index=False,
        encoding="utf-8-sig",
    )
    report = {
        "created_at_utc": pd.Timestamp.utcnow().isoformat(),
        "semantic_status": (
            "Fully automated continuous image measurements. Human category "
            "labels are not used for fitting or threshold selection."
        ),
        "n_heads": int(len(frame)),
        "quality_thresholds": {
            "colour_confidence_floor": args.colour_confidence_floor,
            "shape_confidence_floor": args.shape_confidence_floor,
            "orientation_confidence_floor": args.orientation_confidence_floor,
            "flip_angle_tolerance": args.flip_angle_tolerance,
        },
        "usable_fraction": {
            trait: float(frame[f"{trait}_status"].eq("usable").mean())
            for trait in ("colour", "shape", "orientation")
        },
        "analysis_guidance": {
            "orientation": (
                "Use orientation_angle_degrees as a continuous 0-180 degree measurement."
            ),
            "colour": (
                "Use Lab lightness/chroma and hue sine/cosine; retain colour_status=usable."
            ),
            "shape": (
                "Use aspect ratio, circularity, solidity and width-profile variables; retain shape_status=usable."
            ),
            "categories": (
                "Any categorical bins are downstream display conventions, not learned ground truth."
            ),
        },
    }
    (output / "primary_trait_continuous_provenance.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
