"""Versioned, location-blind adapters for all 27 registered image endpoints.

Reuse the corrected primary engine (56, not the superseded orientation in 52)
and extended engine 89. Raw original/mirror outputs are retained even when a
particular endpoint is ineligible. These are automated image features, not
botanical labels or independently calibrated measurements.
"""
from __future__ import annotations

import argparse
import csv
from functools import lru_cache
import importlib.util
import math
from pathlib import Path

import cv2
import numpy as np

from .workflow import ROOT, text_digest

LEGACY = ROOT / "ch1_global/v2"
REGISTRY = LEGACY / "ontology/ch1_continuous_trait_contract.csv"
SOURCES = [LEGACY / name for name in ("52_measure_primary_traits_continuous.py",
           "55_run_primary_traits_continuous.py", "56_run_primary_traits_continuous_v2.py",
           "89_measure_extended_continuous_traits.py")]
PARAMETERS = {"version": "v3_cached_27_image_features_v1", "primary_min_dimension": 96,
              "colour_confidence_floor": .55, "shape_confidence_floor": .50,
              "orientation_confidence_floor": .65, "flip_angle_tolerance": 20.,
              "min_architecture_dimension": 150, "min_surface_dimension": 300,
              "min_architecture_sharpness": 45., "min_surface_sharpness": 80.,
              "context_green_hue_opencv": [39, 104], "context_green_min_saturation": 35,
              "context_green_min_value": 40, "context_min_pixels": 100,
              "context_min_total_area_fraction": .01,
              "composition_tolerance": 1e-6,
              "context_exclusion": "union of every detector head-padded box, intersected with this context",
              "paired_colour": "same OpenCV uint8 Lab/hue statistics on union floral pixels, all non-head context and green non-head context; no dominant-colour switching",
              "legacy_chroma": "unchanged engine 56: redmagenta-only when dominant redmagenta, otherwise union floral pixels",
              "primary_qc": "legacy mirror quality plus both values finite, registry bounds and minimum pixel dimension; retain all raw values",
              "extended_qc": "engine 89 endpoint-specific mirror, geometry, resolution and sharpness flags plus registry bounds",
              "new_diagnostics_are_not_endpoints": True}


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@lru_cache(maxsize=1)
def engines():
    # Load the compatibility wrapper only once; it patches OpenCV HoughLinesP.
    return load(SOURCES[2], "v3_primary56"), load(SOURCES[3], "v3_extended89")


def registry():
    with REGISTRY.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 27 or len({r["endpoint_id"] for r in rows}) != 27:
        raise ValueError("The original 27-endpoint registry changed")
    return rows


def specification():
    return {"parameters": PARAMETERS, "source_sha256_text_lf": {
        path.relative_to(ROOT).as_posix(): text_digest(path)
        for path in [*SOURCES, REGISTRY, Path(__file__)]}}


def clean(value):
    """Strict JSON: numeric missingness is null, never an invented zero."""
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [clean(v) for v in value]
    if isinstance(value, np.generic):
        return clean(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def number(value):
    try:
        return float(value) if math.isfinite(float(value)) else None
    except (ValueError, TypeError):
        return None


def mean_pair(a, b):
    finite = [v for v in (number(a), number(b)) if v is not None]
    return sum(finite) / len(finite) if finite else None


def bound_reasons(values, row):
    reasons = []
    if any(number(v) is None for v in values):
        reasons.append("measurement_missing")
    for v in values:
        v = number(v)
        if v is not None and ((row["lower_bound"] and v < float(row["lower_bound"])) or
                              (row["upper_bound"] and v > float(row["upper_bound"]))):
            reasons.append("outside_registry_bounds")
    return reasons


def colour_summary(image, mask):
    mask = mask.astype(bool)
    count = int(mask.sum())
    result = {"n_pixels": count, "lab_lightness": None, "lab_chroma": None,
              "hue_sin": None, "hue_cos": None, "hue_resultant": None}
    if not count:
        return result
    lab = cv2.cvtColor(image, cv2.COLOR_BGR2LAB).astype(float)
    hsv = cv2.cvtColor(image, cv2.COLOR_BGR2HSV)
    result.update(lab_lightness=float(np.median(lab[..., 0][mask])*100/255),
                  lab_chroma=float(np.median(np.hypot(lab[..., 1][mask]-128, lab[..., 2][mask]-128))))
    chromatic = mask & (hsv[..., 1] >= 25)
    if chromatic.any():
        radians = np.deg2rad(hsv[..., 0][chromatic].astype(float)*2)
        sine, cosine = float(np.sin(radians).mean()), float(np.cos(radians).mean())
        result.update(hue_sin=sine, hue_cos=cosine, hue_resultant=math.hypot(sine, cosine))
    return result


def context_masks(context, global_context_box, all_head_boxes):
    x1, y1, x2, y2 = global_context_box
    if context.shape[:2] != (y2-y1, x2-x1):
        raise ValueError("Context image and global crop coordinates differ")
    available = np.ones(context.shape[:2], dtype=bool)
    for bx1, by1, bx2, by2 in all_head_boxes:
        left, top, right, bottom = max(x1, bx1), max(y1, by1), min(x2, bx2), min(y2, by2)
        if right > left and bottom > top:
            available[top-y1:bottom-y1, left-x1:right-x1] = False
    hsv = cv2.cvtColor(context, cv2.COLOR_BGR2HSV)
    lo, hi = PARAMETERS["context_green_hue_opencv"]
    green = available & (hsv[..., 0] >= lo) & (hsv[..., 0] <= hi) & (hsv[..., 1] >= PARAMETERS["context_green_min_saturation"]) & (hsv[..., 2] >= PARAMETERS["context_green_min_value"])
    return available, green


def measure(head, context, global_context_box, all_head_boxes):
    primary, extended = engines()
    base = primary.MEASURE
    raw = {}
    errors = {}
    for side, h, c in (("original", head, context), ("mirror", cv2.flip(head, 1), cv2.flip(context, 1))):
        raw[side] = {}
        for group, function in (("colour", lambda: primary.colour_measurement_v2(h)),
                                ("shape", lambda: base.shape_measurement(h)),
                                ("orientation", lambda: primary.orientation_measurement_v2(h, c)),
                                ("extended", lambda: extended.measure_once(h, argparse.Namespace(**PARAMETERS)))):
            try:
                raw[side][group] = function()
            except (ValueError, RuntimeError, cv2.error) as error:
                raw[side][group] = {}
                errors[f"{side}:{group}"] = {"type": type(error).__name__, "message": str(error)}
    a, b = raw["original"], raw["mirror"]
    statuses = {}
    for group, fn, args in (("colour", base.status_colour, (PARAMETERS["colour_confidence_floor"],)),
                            ("shape", base.status_shape, (PARAMETERS["shape_confidence_floor"],)),
                            ("orientation", base.status_orientation, (PARAMETERS["orientation_confidence_floor"], PARAMETERS["flip_angle_tolerance"]))):
        status, confidence = fn(a[group], b[group], *args)
        if status == "head_peduncle_axis_not_recovered":
            status = "signed_head_axis_not_recovered"
        statuses[group] = {"status": status, "confidence": confidence}
    combined_extended = extended.combine(a["extended"], b["extended"])
    mapping = {"orientation_angle_degrees": ("orientation", "angle_degrees"),
               "corolla_visible_fraction": ("colour", "floral_pixel_fraction"),
               "corolla_lab_lightness": ("colour", "median_lab_lightness"),
               "corolla_lab_chroma": ("colour", "median_lab_chroma"),
               **{f"corolla_{k}_fraction": ("colour", f"{k}_fraction") for k in ("white", "redmagenta", "purple", "yellow")},
               **{f"shape_{k}": ("shape", k) for k in ("aspect_ratio", "circularity", "solidity", "width_cv")}}
    endpoints = []
    for row in registry():
        variable = row["measurement_variable"]
        if variable.startswith("corolla_hue_"):
            group = "colour"
            angles = [number(v[group].get("mean_hue_degrees")) for v in (a, b)]
            trig = math.sin if variable.endswith("sin") else math.cos
            values = [trig(math.radians(v)) if v is not None else None for v in angles]
            finite_angles = [math.radians(v) for v in angles if v is not None]
            s = sum(math.sin(v) for v in finite_angles)
            c = sum(math.cos(v) for v in finite_angles)
            value = trig(math.atan2(s, c)) if finite_angles and math.hypot(s, c) > 1e-12 else None
        elif variable in mapping:
            group, key = mapping[variable]
            values = [v[group].get(key) for v in (a, b)]
            value = mean_pair(*values)
        else:
            group = "extended"
            values = [v[group].get(variable) for v in (a, b)]
            value = combined_extended[variable]
        reasons = bound_reasons([*values, value], row)
        if group == "extended":
            status = combined_extended[f"{variable}_status"]
        else:
            status = statuses[group]["status"]
            if number(statuses[group]["confidence"]) is None:
                reasons.append("quality_missing")
            if min(head.shape[:2]) < PARAMETERS["primary_min_dimension"]:
                reasons.append("low_resolution")
        if status != "usable":
            reasons.extend(status.split(";"))
        if any(f"{side}:{group}" in errors for side in ("original", "mirror")):
            reasons.append("engine_error")
        if row["compositional_group"]:
            for side in (a, b):
                fractions = [number(side["colour"].get(f"{k}_fraction")) for k in ("white", "redmagenta", "purple", "yellow")]
                if any(v is None for v in fractions) or abs(sum(fractions)-1) > PARAMETERS["composition_tolerance"]:
                    reasons.append("incomplete_or_unclosed_composition")
        endpoints.append({"endpoint_id": row["endpoint_id"], "unit": row["unit"], "value": value,
                          "original": values[0], "mirror": values[1],
                          "mirror_abs_difference": abs(values[0]-values[1]) if all(number(v) is not None for v in values) else None,
                          "status": ";".join(sorted(set(reasons))) if reasons else "usable"})
    # These masks/statistics are new diagnostic outputs, separate from the 27.
    foreground, quality = base.central_foreground(head)
    floral, _, _ = primary._exclusive_colour_masks(head, foreground)
    available, green = context_masks(context, global_context_box, all_head_boxes)
    paired = {"floral_union": colour_summary(head, floral),
              "non_head_context": colour_summary(context, available), "green_non_head_context": colour_summary(context, green)}
    for kind in ("non_head_context", "green_non_head_context"):
        n = paired[kind]["n_pixels"]
        paired[kind]["support_status"] = "available" if n >= max(PARAMETERS["context_min_pixels"], PARAMETERS["context_min_total_area_fraction"]*context.shape[0]*context.shape[1]) else "insufficient_pixels"
    diagnostics = {"head_width_px": head.shape[1], "head_height_px": head.shape[0],
                   "head_min_dimension_px": min(head.shape[:2]),
                   "head_laplacian_variance": float(cv2.Laplacian(cv2.cvtColor(head, cv2.COLOR_BGR2GRAY), cv2.CV_64F).var()),
                   "primary_foreground_quality": quality, "paired_colour": paired,
                   "context_is_verified_leaf_tissue": False}
    return clean({"endpoints": endpoints, "raw": raw, "legacy_qc": statuses,
                  "extended_combined": combined_extended, "engine_errors": errors, "diagnostics": diagnostics}), {
                      "head_foreground": foreground, "head_floral_union": floral,
                      "context_non_head": available, "context_green_non_head": green}
