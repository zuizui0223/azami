#!/usr/bin/env python3
"""Run the corrected v2 continuous primary-trait measurements.

V2 fixes two issues detected by the first complete run:

1. Colour-family masks are mutually exclusive, so visible-corolla coverage
   cannot exceed one and the four colour fractions sum to one.
2. Orientation is the signed longitudinal head-axis angle relative to the
   EXIF-oriented image vertical (0 degrees up, 90 horizontal, 180 down).
   The sign of the PCA axis is chosen from the visible floral-pixel centroid.

Human category labels are not used by either correction.
"""
from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path
from typing import Any

import cv2
import numpy as np
import pandas as pd


_COMPAT_PATH = Path(__file__).with_name("55_run_primary_traits_continuous.py")
_SPEC = importlib.util.spec_from_file_location("ch1_continuous_v1_compat", _COMPAT_PATH)
if _SPEC is None or _SPEC.loader is None:
    raise RuntimeError(f"Unable to load compatibility entrypoint: {_COMPAT_PATH}")
COMPAT = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(COMPAT)
MEASURE = COMPAT.MEASURE


def _exclusive_colour_masks(
    image: np.ndarray,
    foreground: np.ndarray,
) -> tuple[np.ndarray, dict[str, np.ndarray], dict[str, np.ndarray]]:
    """Return mutually exclusive floral colour masks within the head foreground."""
    hsv = cv2.cvtColor(image, cv2.COLOR_BGR2HSV)
    lab = cv2.cvtColor(image, cv2.COLOR_BGR2LAB)
    hue, saturation, value = [hsv[..., index] for index in range(3)]
    lightness, a_axis, b_axis = [
        lab[..., index].astype(np.float32) for index in range(3)
    ]

    height, width = foreground.shape
    yy, xx = np.mgrid[0:height, 0:width]
    centre_x, centre_y = (width - 1) / 2, (height - 1) / 2
    radius_x, radius_y = max(width * 0.52, 1), max(height * 0.52, 1)
    central = (
        ((xx - centre_x) / radius_x) ** 2
        + ((yy - centre_y) / radius_y) ** 2
        <= 1
    )
    base = foreground.astype(bool) & central

    bright = value >= 80
    # OpenCV hue is [0, 179]. Boundaries are intentionally non-overlapping.
    purple = base & bright & (saturation >= 35) & (hue >= 105) & (hue <= 164)
    redmagenta = base & bright & (saturation >= 25) & (
        (hue >= 165) | (hue <= 10)
    )
    yellow = base & (value >= 115) & (saturation >= 45) & (
        (hue >= 11) & (hue <= 38)
    )
    occupied = purple | redmagenta | yellow

    gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    laplacian = np.abs(cv2.Laplacian(gray, cv2.CV_32F))
    edge_threshold = (
        float(np.percentile(laplacian[base], 45)) if base.any() else 0.0
    )
    white = (
        base
        & ~occupied
        & (value >= 170)
        & (saturation <= 58)
        & (laplacian >= edge_threshold)
    )
    floral = purple | redmagenta | yellow | white
    masks = {
        "purple_blue": purple,
        "redmagenta": redmagenta,
        "yellow_or_other": yellow,
        "white_or_cream": white,
    }
    channels = {
        "hue": hue,
        "saturation": saturation,
        "value": value,
        "lightness": lightness,
        "a_axis": a_axis,
        "b_axis": b_axis,
        "base": base,
    }
    return floral, masks, channels


def colour_measurement_v2(image: np.ndarray) -> dict[str, Any]:
    foreground, quality = MEASURE.central_foreground(image)
    floral, masks, channels = _exclusive_colour_masks(image, foreground)
    base = channels["base"]
    if int(base.sum()) < 20:
        return {"state": MEASURE.NON_SCOREABLE, "confidence": 0.0, **quality}

    counts = {name: int(mask.sum()) for name, mask in masks.items()}
    total = int(floral.sum())
    coverage = total / max(int(base.sum()), 1)
    if total < max(12, int(base.sum() * 0.025)) or coverage < 0.035:
        return {
            "state": MEASURE.NON_SCOREABLE,
            "confidence": float(np.clip(
                quality["mask_quality"] * coverage / 0.06, 0, 0.45
            )),
            "floral_pixel_fraction": coverage,
            "colour_counts_json": json.dumps(counts, sort_keys=True),
            **quality,
        }

    proportions = {name: count / total for name, count in counts.items()}
    if not math.isclose(sum(proportions.values()), 1.0, abs_tol=1e-9):
        raise RuntimeError("Exclusive colour fractions do not sum to one")
    if coverage > 1.0 + 1e-9:
        raise RuntimeError("Visible-corolla coverage exceeds one")

    hue = channels["hue"]
    saturation = channels["saturation"]
    value = channels["value"]
    lightness = channels["lightness"]
    a_axis = channels["a_axis"]
    b_axis = channels["b_axis"]
    chromatic = floral & (saturation >= 25)
    median_lightness = float(np.median(lightness[floral]) * 100.0 / 255.0)
    if chromatic.any():
        hue_radians = np.deg2rad(hue[chromatic].astype(np.float32) * 2.0)
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
    selected = masks["redmagenta"] if top_name == "redmagenta" else floral
    median_saturation = float(np.median(saturation[selected]))
    median_chroma = float(np.median(np.sqrt(
        (a_axis[selected] - 128) ** 2 + (b_axis[selected] - 128) ** 2
    )))

    state = top_name
    if top_name == "redmagenta":
        state = (
            "pale_pink"
            if median_saturation < 72 or median_chroma < 24
            else "pink"
        )
    if top_fraction < 0.68 and second_fraction >= 0.24 and coverage >= 0.08:
        state = "mixed"

    dominance = max(0.0, top_fraction - second_fraction)
    exposure = float(np.clip((np.median(value[floral]) - 55) / 140, 0, 1))
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
        "state": state,
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
        "colour_counts_json": json.dumps(counts, sort_keys=True),
        **quality,
    }


def orientation_measurement_v2(
    head: np.ndarray,
    context: np.ndarray,
) -> dict[str, Any]:
    """Measure signed head-axis angle relative to the oriented image vertical."""
    max_side = max(context.shape[:2])
    if max_side > 320:
        scale = 320.0 / max_side
        context = cv2.resize(
            context, None, fx=scale, fy=scale, interpolation=cv2.INTER_AREA
        )
        head = cv2.resize(
            head, None, fx=scale, fy=scale, interpolation=cv2.INTER_AREA
        )

    box, template_score = MEASURE.template_head_box(head, context)
    if box is None or template_score < 0.35:
        return {
            "state": MEASURE.NON_SCOREABLE,
            "confidence": 0.0,
            "template_score": template_score,
        }

    foreground, head_quality = MEASURE.central_foreground(head)
    rows, columns = np.nonzero(foreground)
    if len(columns) < 30:
        return {
            "state": MEASURE.NON_SCOREABLE,
            "confidence": 0.0,
            "template_score": template_score,
            **{f"head_{key}": value for key, value in head_quality.items()},
        }

    coordinates = np.column_stack([columns, rows]).astype(np.float32)
    mean, eigenvectors, eigenvalues = cv2.PCACompute2(
        coordinates, mean=None, maxComponents=2
    )
    axis = eigenvectors[0].astype(float)
    axis /= max(np.linalg.norm(axis), 1e-9)

    floral, _, _ = _exclusive_colour_masks(head, foreground)
    floral_rows, floral_columns = np.nonzero(floral)
    minimum_floral_pixels = max(12, int(len(columns) * 0.025))
    if len(floral_columns) < minimum_floral_pixels:
        return {
            "state": MEASURE.NON_SCOREABLE,
            "confidence": 0.0,
            "template_score": template_score,
            **{f"head_{key}": value for key, value in head_quality.items()},
        }

    centre = mean[0].astype(float)
    floral_coordinates = np.column_stack([
        floral_columns, floral_rows
    ]).astype(float)
    projections = (floral_coordinates - centre) @ axis
    signed_projection = float(np.median(projections))
    diagonal = math.hypot(head.shape[0], head.shape[1])
    if abs(signed_projection) < 0.004 * diagonal:
        return {
            "state": MEASURE.NON_SCOREABLE,
            "confidence": 0.0,
            "template_score": template_score,
            **{f"head_{key}": value for key, value in head_quality.items()},
        }

    signed_axis = axis if signed_projection >= 0 else -axis
    image_up = np.array([0.0, -1.0])
    angle = math.degrees(math.acos(float(np.clip(
        np.dot(signed_axis, image_up), -1, 1
    ))))

    anisotropy = float(
        eigenvalues[0, 0] / max(eigenvalues[1, 0], 1e-9)
    )
    anisotropy_score = float(np.clip((anisotropy - 1) / 1.5, 0, 1))
    floral_coverage = float(floral.sum() / max(foreground.sum(), 1))
    coverage_score = float(np.clip((floral_coverage - 0.025) / 0.25, 0, 1))
    sign_score = float(np.clip(
        abs(signed_projection) / max(0.035 * diagonal, 1e-9), 0, 1
    ))
    confidence = float(np.clip(
        0.20 * template_score
        + 0.30 * head_quality["mask_quality"]
        + 0.20 * anisotropy_score
        + 0.20 * sign_score
        + 0.10 * coverage_score,
        0,
        1,
    ))
    return {
        "state": "measured",
        "confidence": confidence,
        "angle_degrees": float(angle),
        "template_score": template_score,
        # Retain the legacy output slot as a generic signed-axis evidence score.
        "stem_line_score": sign_score,
        "stem_line_length": float("nan"),
        **{f"head_{key}": value for key, value in head_quality.items()},
    }


def _postprocess_output(output_directory: Path) -> None:
    measurement_path = (
        output_directory / "primary_trait_continuous_head_measurements.csv"
    )
    table = pd.read_csv(measurement_path)
    table["continuous_trait_version"] = "v2"
    table["orientation_reference"] = "exif_oriented_image_vertical"
    table["measurement_basis"] = (
        "v2_exclusive_colour_masks_plus_signed_pca_floral_axis_relative_to_"
        "image_vertical_plus_continuous_outline_geometry"
    )

    visible = pd.to_numeric(table["corolla_visible_fraction"], errors="coerce")
    if visible.dropna().gt(1.0 + 1e-9).any():
        raise RuntimeError("V2 output contains corolla-visible fractions above one")
    fractions = table[[
        "corolla_white_fraction",
        "corolla_redmagenta_fraction",
        "corolla_purple_fraction",
        "corolla_yellow_fraction",
    ]].apply(pd.to_numeric, errors="coerce")
    scoreable = fractions.notna().all(axis=1)
    if scoreable.any() and not np.allclose(
        fractions.loc[scoreable].sum(axis=1), 1.0, atol=1e-6
    ):
        raise RuntimeError("V2 output colour-family fractions do not sum to one")
    table.to_csv(measurement_path, index=False, encoding="utf-8-sig")

    provenance_path = output_directory / "primary_trait_continuous_provenance.json"
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    provenance["continuous_trait_version"] = "v2"
    provenance["corrections"] = {
        "colour": (
            "Colour-family masks are mutually exclusive; visible-corolla "
            "coverage is the union of the four masks and is bounded by one."
        ),
        "orientation": (
            "Angle is signed by the visible floral-pixel centroid and measured "
            "relative to EXIF-oriented image vertical: 0 up, 90 horizontal, "
            "180 down. It is not the v1 head-peduncle bend angle."
        ),
    }
    provenance["analysis_guidance"]["orientation"] = (
        "Use orientation_angle_degrees as an image-referenced signed continuous "
        "head-axis angle (0 up, 90 horizontal, 180 down), retaining only "
        "orientation_status=usable."
    )
    provenance_path.write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2), encoding="utf-8"
    )


def main() -> None:
    args = MEASURE.parse_args()
    MEASURE.colour_measurement = colour_measurement_v2
    MEASURE.orientation_measurement = orientation_measurement_v2
    MEASURE.main()
    _postprocess_output(Path(args.out_dir))


if __name__ == "__main__":
    main()
