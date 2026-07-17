#!/usr/bin/env python3
"""Build inspectable Figure 1 measurement panels from frozen v2 inputs.

The input selection CSV is deliberately manual and licence-gated. Each row names one
head and its source, tight and context images. The script reuses the production v2
measurement functions to draw overlays, preventing figure-only reimplementation.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
from pathlib import Path
from typing import Any

import cv2
import numpy as np
import pandas as pd

REQUIRED = {
    "panel_id", "annotation_unit_id", "panel_role", "source_image_path",
    "tight_crop_path", "context_crop_path", "bbox_x1", "bbox_y1", "bbox_x2", "bbox_y2",
    "photo_license_code", "photo_attribution", "licence_verified",
}
OPEN_FIGURE_LICENSES = {"cc0", "cc-by", "cc-by-4.0", "cc-by-3.0", "cc-by-sa", "cc-by-sa-4.0", "cc-by-sa-3.0"}


def load_v2_module() -> Any:
    path = Path(__file__).with_name("56_run_primary_traits_continuous_v2.py")
    spec = importlib.util.spec_from_file_location("ch1_primary_traits_v2_for_overlay", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import production measurement module: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create Figure 1 measurement audit overlays.")
    parser.add_argument("--selection", required=True, help="Licence-audited Figure 1 selection CSV.")
    parser.add_argument("--measurements", required=True, help="primary_trait_continuous_head_measurements.csv")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--allow-nonopen-licence", action="store_true", help="Bypass open-licence gate only after publisher/legal review.")
    return parser.parse_args()


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def as_bool(value: Any) -> bool:
    return text(value).lower() in {"1", "true", "yes", "y"}


def read_image(path_value: Any) -> np.ndarray:
    path = Path(text(path_value))
    image = cv2.imread(str(path), cv2.IMREAD_COLOR)
    if image is None:
        raise ValueError(f"Unreadable image: {path}")
    return image


def finite(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def source_bbox_overlay(image: np.ndarray, row: pd.Series) -> np.ndarray:
    out = image.copy()
    coordinates = [int(round(float(row[name]))) for name in ("bbox_x1", "bbox_y1", "bbox_x2", "bbox_y2")]
    thickness = max(2, round(min(out.shape[:2]) / 250))
    cv2.rectangle(out, (coordinates[0], coordinates[1]), (coordinates[2], coordinates[3]), (0, 0, 255), thickness)
    cv2.putText(out, text(row["panel_id"]), (coordinates[0], max(24, coordinates[1] - 8)), cv2.FONT_HERSHEY_SIMPLEX, 0.7, (0, 0, 255), 2, cv2.LINE_AA)
    return out


def colour_overlay(v2: Any, tight: np.ndarray) -> tuple[np.ndarray, dict[str, Any]]:
    foreground, _ = v2.MEASURE.central_foreground(tight)
    floral, masks, _ = v2._exclusive_colour_masks(tight, foreground)
    out = tight.copy()
    colours = {"redmagenta": (0, 0, 255), "purple_blue": (255, 0, 120), "white_or_cream": (0, 255, 0), "yellow_or_other": (0, 255, 255)}
    for name, mask in masks.items():
        selected = mask.astype(bool)
        tint = np.zeros_like(out)
        tint[:] = colours[name]
        out[selected] = cv2.addWeighted(out[selected], 0.55, tint[selected], 0.45, 0)
    contours, _ = cv2.findContours((floral.astype(np.uint8) * 255), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    cv2.drawContours(out, contours, -1, (255, 255, 255), 1)
    return out, v2.colour_measurement_v2(tight)


def shape_overlay(v2: Any, tight: np.ndarray) -> tuple[np.ndarray, dict[str, Any]]:
    foreground, _ = v2.MEASURE.central_foreground(tight)
    height, width = foreground.shape
    kernel_size = max(3, int(round(min(height, width) * 0.055)) | 1)
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (kernel_size, kernel_size))
    body = cv2.morphologyEx(foreground, cv2.MORPH_OPEN, kernel)
    body = cv2.morphologyEx(body, cv2.MORPH_CLOSE, kernel)
    body = v2.MEASURE.largest_central_component(body)
    out = tight.copy()
    contours, _ = cv2.findContours((body * 255).astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    if contours:
        contour = max(contours, key=cv2.contourArea)
        hull = cv2.convexHull(contour)
        cv2.drawContours(out, [hull], -1, (0, 255, 255), 2)
        cv2.drawContours(out, [contour], -1, (0, 0, 255), 2)
        moments = cv2.moments(contour)
        if moments["m00"]:
            centre = np.array([moments["m10"] / moments["m00"], moments["m01"] / moments["m00"]])
            rows, columns = np.nonzero(body)
            if len(columns) >= 20:
                coords = np.column_stack([columns, rows]).astype(np.float32)
                _, eigenvectors = cv2.PCACompute(coords, mean=None, maxComponents=2)
                axis = eigenvectors[0].astype(float)
                length = 0.45 * math.hypot(height, width)
                a = tuple(np.round(centre - axis * length).astype(int))
                b = tuple(np.round(centre + axis * length).astype(int))
                cv2.line(out, a, b, (255, 0, 0), 2, cv2.LINE_AA)
    return out, v2.MEASURE.shape_measurement(tight)


def orientation_overlay(v2: Any, tight: np.ndarray, context: np.ndarray) -> tuple[np.ndarray, dict[str, Any]]:
    result = v2.orientation_measurement_v2(tight, context)
    out = context.copy()
    box, _ = v2.MEASURE.template_head_box(tight, context)
    if box is not None:
        x0, y0, x1, y1 = box
        cv2.rectangle(out, (x0, y0), (x1, y1), (0, 0, 255), 2)
        centre = np.array([(x0 + x1) / 2, (y0 + y1) / 2], dtype=float)
        angle = finite(result.get("angle_degrees"))
        scale = 0.45 * max(x1 - x0, y1 - y0)
        cv2.line(out, tuple(np.round(centre).astype(int)), tuple(np.round(centre + np.array([0, -scale])).astype(int)), (0, 255, 0), 2, cv2.LINE_AA)
        if angle is not None:
            radians = math.radians(angle)
            endpoint = centre + np.array([math.sin(radians), -math.cos(radians)]) * scale
            cv2.arrowedLine(out, tuple(np.round(centre).astype(int)), tuple(np.round(endpoint).astype(int)), (255, 0, 0), 2, cv2.LINE_AA, tipLength=0.15)
            cv2.putText(out, f"{angle:.1f} deg", (max(5, x0), max(24, y0 - 8)), cv2.FONT_HERSHEY_SIMPLEX, 0.65, (255, 0, 0), 2, cv2.LINE_AA)
    return out, result


def save(path: Path, image: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not cv2.imwrite(str(path), image):
        raise OSError(f"Failed to write {path}")


def main() -> None:
    args = parse_args()
    selection = pd.read_csv(args.selection, dtype=str, keep_default_na=False)
    missing = REQUIRED.difference(selection.columns)
    if missing:
        raise ValueError(f"Selection missing columns: {sorted(missing)}")
    if selection["panel_id"].duplicated().any() or selection["annotation_unit_id"].duplicated().any():
        raise ValueError("panel_id and annotation_unit_id must be unique")
    selection["licence_verified_bool"] = selection["licence_verified"].map(as_bool)
    if not selection["licence_verified_bool"].all():
        raise ValueError("Every selected image must have licence_verified=true")
    normalized_licenses = selection["photo_license_code"].map(text).str.lower()
    unsupported = sorted(set(normalized_licenses).difference(OPEN_FIGURE_LICENSES))
    if unsupported and not args.allow_nonopen_licence:
        raise ValueError(f"Non-open or unrecognized Figure licences: {unsupported}")
    measurements = pd.read_csv(args.measurements, dtype=str, keep_default_na=False, low_memory=False)
    if "annotation_unit_id" not in measurements.columns:
        raise ValueError("Measurements missing annotation_unit_id")
    merged = selection.merge(measurements, on="annotation_unit_id", how="left", validate="one_to_one", suffixes=("", "_measurement"))
    out_dir = Path(args.out_dir)
    v2 = load_v2_module()
    audit_rows: list[dict[str, Any]] = []
    for _, row in merged.iterrows():
        panel_id = text(row["panel_id"])
        panel_dir = out_dir / "panels" / panel_id
        try:
            source = read_image(row["source_image_path"])
            tight = read_image(row["tight_crop_path"])
            context = read_image(row["context_crop_path"])
            colour_image, colour_result = colour_overlay(v2, tight)
            shape_image, shape_result = shape_overlay(v2, tight)
            orientation_image, orientation_result = orientation_overlay(v2, tight, context)
            save(panel_dir / "01_source_bbox.png", source_bbox_overlay(source, row))
            save(panel_dir / "02_tight_crop.png", tight)
            save(panel_dir / "03_context_crop.png", context)
            save(panel_dir / "04_orientation_overlay.png", orientation_image)
            save(panel_dir / "05_colour_overlay.png", colour_image)
            save(panel_dir / "06_shape_overlay.png", shape_image)
            status, message = "success", ""
        except Exception as error:
            colour_result = shape_result = orientation_result = {}
            status, message = "failed", str(error)
        audit_rows.append({"panel_id": panel_id, "annotation_unit_id": text(row["annotation_unit_id"]), "panel_role": text(row["panel_role"]), "status": status, "message": message, "photo_license_code": text(row["photo_license_code"]), "photo_attribution": text(row["photo_attribution"]), "orientation_state_recomputed": text(orientation_result.get("state", "")), "orientation_angle_recomputed": orientation_result.get("angle_degrees", ""), "colour_state_recomputed": text(colour_result.get("state", "")), "shape_state_recomputed": text(shape_result.get("state", ""))})
    audit = pd.DataFrame(audit_rows)
    out_dir.mkdir(parents=True, exist_ok=True)
    audit.to_csv(out_dir / "figure1_panel_audit.csv", index=False, encoding="utf-8-sig")
    provenance = {"n_selected": int(len(selection)), "n_success": int(audit["status"].eq("success").sum()), "n_failed": int(audit["status"].eq("failed").sum()), "open_figure_licences": sorted(OPEN_FIGURE_LICENSES), "overlay_basis": "production v2 measurement functions imported from 56_run_primary_traits_continuous_v2.py", "interpretation": "Overlays are face-validity displays; frozen numerical outputs remain the source of manuscript values."}
    (out_dir / "figure1_panel_provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    print(json.dumps(provenance, indent=2))


if __name__ == "__main__":
    main()
