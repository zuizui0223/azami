#!/usr/bin/env python3
"""Build an actual-photo Figure 1 from executed Chapter 1 artifacts.

The detector boxes and continuous trait measurements are read from frozen workflow
artifacts. Source photographs are re-downloaded from their recorded iNaturalist URLs.
The production v2 measurement functions are imported to draw colour, shape, and
orientation overlays. No detector or ecological model is refit.
"""
from __future__ import annotations

import argparse
import importlib.util
import io
import json
import math
from pathlib import Path
from typing import Any

import cv2
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from PIL import Image

SELECTED_UNITS = [
    "ch1atlas_0002870_head_01",  # purple, upright, CC BY
    "ch1atlas_0002331_head_01",  # pale/white, inclined, CC0
    "ch1atlas_0003663_head_02",  # pink, inclined, CC0
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--global-artifact-dir", required=True)
    p.add_argument("--continuous-artifact-dir", required=True)
    p.add_argument("--metadata-artifact-dir", required=True)
    p.add_argument("--repo-root", default=".")
    p.add_argument("--out-dir", required=True)
    return p.parse_args()


def load_module(path: Path, name: str) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def download_image(url: str) -> Image.Image:
    response = requests.get(url, timeout=90, headers={"User-Agent": "azami-figure1/1.0"})
    response.raise_for_status()
    with Image.open(io.BytesIO(response.content)) as img:
        return img.convert("RGB")


def clamp_box(width: int, height: int, row: pd.Series) -> tuple[int, int, int, int]:
    left = max(0, min(width - 1, math.floor(float(row.bbox_x1))))
    top = max(0, min(height - 1, math.floor(float(row.bbox_y1))))
    right = max(left + 1, min(width, math.ceil(float(row.bbox_x2))))
    bottom = max(top + 1, min(height, math.ceil(float(row.bbox_y2))))
    return left, top, right, bottom


def crop_pair(image: Image.Image, row: pd.Series, pad_ratio: float = 1.5) -> tuple[Image.Image, Image.Image]:
    width, height = image.size
    left, top, right, bottom = clamp_box(width, height, row)
    head = image.crop((left, top, right, bottom)).convert("RGB")
    bw, bh = right - left, bottom - top
    context = image.crop((
        max(0, math.floor(left - pad_ratio * bw)),
        max(0, math.floor(top - pad_ratio * bh)),
        min(width, math.ceil(right + pad_ratio * bw)),
        min(height, math.ceil(bottom + pad_ratio * bh)),
    )).convert("RGB")
    return head, context


def bgr(path: Path) -> np.ndarray:
    image = cv2.imread(str(path), cv2.IMREAD_COLOR)
    if image is None:
        raise RuntimeError(f"Cannot read {path}")
    return image


def rgb(path: Path) -> np.ndarray:
    return cv2.cvtColor(bgr(path), cv2.COLOR_BGR2RGB)


def panel_label(ax: plt.Axes, label: str, title: str) -> None:
    ax.text(0.01, 0.99, label, transform=ax.transAxes, ha="left", va="top",
            fontsize=13, fontweight="bold", bbox=dict(facecolor="white", alpha=0.78, edgecolor="none", pad=2))
    ax.set_title(title, fontsize=10, pad=4)
    ax.axis("off")


def main() -> None:
    a = parse_args()
    repo = Path(a.repo_root)
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    work = out / "working"
    work.mkdir(exist_ok=True)

    global_dir = Path(a.global_artifact_dir)
    cont_dir = Path(a.continuous_artifact_dir)
    meta_dir = Path(a.metadata_artifact_dir)

    yolo = pd.read_csv(global_dir / "merged_shards/global_ai_merged_yolo_crop_metadata.csv", low_memory=False)
    yolo["annotation_unit_id"] = yolo["queue_id"].astype(str) + "_head_" + (pd.to_numeric(yolo["det_index"]).astype(int) + 1).astype(str).str.zfill(2)
    yolo = yolo[yolo.annotation_unit_id.isin(SELECTED_UNITS)].copy()
    if set(yolo.annotation_unit_id) != set(SELECTED_UNITS):
        raise RuntimeError("Selected units not all present in frozen YOLO metadata")

    measurements = pd.read_csv(cont_dir / "continuous/primary_trait_continuous_head_measurements.csv", low_memory=False)
    measurements = measurements[measurements.annotation_unit_id.isin(SELECTED_UNITS)].copy()
    metadata = pd.read_csv(
        meta_dir / "photo_metadata_merged.csv",
        usecols=["photo_id", "photo_license_code", "photo_attribution", "user_login"],
        low_memory=False,
    ).drop_duplicates("photo_id")
    selected = yolo.merge(measurements, on="annotation_unit_id", how="left", validate="one_to_one", suffixes=("", "_measurement"))
    selected = selected.merge(metadata, on="photo_id", how="left", validate="many_to_one")
    selected["display_order"] = selected.annotation_unit_id.map({v: i for i, v in enumerate(SELECTED_UNITS)})
    selected = selected.sort_values("display_order")

    selection_rows = []
    for _, row in selected.iterrows():
        unit = row.annotation_unit_id
        source = download_image(str(row.medium_image_url))
        head, context = crop_pair(source, row)
        source_path = work / f"{unit}_source.jpg"
        head_path = work / f"{unit}_head.jpg"
        context_path = work / f"{unit}_context.jpg"
        source.save(source_path, quality=95)
        head.save(head_path, quality=95)
        context.save(context_path, quality=95)
        selection_rows.append({
            "panel_id": f"YOLO {float(row.yolo_conf):.2f}",
            "annotation_unit_id": unit,
            "panel_role": ["purple_upright", "pale_inclined", "pink_inclined"][int(row.display_order)],
            "source_image_path": str(source_path),
            "tight_crop_path": str(head_path),
            "context_crop_path": str(context_path),
            "bbox_x1": row.bbox_x1,
            "bbox_y1": row.bbox_y1,
            "bbox_x2": row.bbox_x2,
            "bbox_y2": row.bbox_y2,
            "photo_license_code": row.photo_license_code,
            "photo_attribution": row.photo_attribution,
            "licence_verified": "true",
        })
    selection = pd.DataFrame(selection_rows)
    selection_path = work / "figure1_selection.csv"
    selection.to_csv(selection_path, index=False)

    audit_module = load_module(repo / "ch1_global/v2/78_build_figure1_measurement_audit.py", "figure1_audit")
    v2 = audit_module.load_v2_module()
    panel_dirs: dict[str, Path] = {}
    audit_rows = []
    for _, sel in selection.iterrows():
        unit = sel.annotation_unit_id
        panel_dir = work / "panels" / unit
        panel_dir.mkdir(parents=True, exist_ok=True)
        source_bgr = bgr(Path(sel.source_image_path))
        tight_bgr = bgr(Path(sel.tight_crop_path))
        context_bgr = bgr(Path(sel.context_crop_path))
        colour_image, colour_result = audit_module.colour_overlay(v2, tight_bgr)
        shape_image, shape_result = audit_module.shape_overlay(v2, tight_bgr)
        orientation_image, orientation_result = audit_module.orientation_overlay(v2, tight_bgr, context_bgr)
        audit_module.save(panel_dir / "01_source_bbox.png", audit_module.source_bbox_overlay(source_bgr, sel))
        audit_module.save(panel_dir / "02_tight_crop.png", tight_bgr)
        audit_module.save(panel_dir / "03_context_crop.png", context_bgr)
        audit_module.save(panel_dir / "04_orientation_overlay.png", orientation_image)
        audit_module.save(panel_dir / "05_colour_overlay.png", colour_image)
        audit_module.save(panel_dir / "06_shape_overlay.png", shape_image)
        panel_dirs[unit] = panel_dir
        audit_rows.append({
            "annotation_unit_id": unit,
            "orientation_recomputed": orientation_result.get("angle_degrees"),
            "colour_state_recomputed": colour_result.get("state"),
            "shape_state_recomputed": shape_result.get("state"),
        })

    main_unit = SELECTED_UNITS[0]
    main_row = selected[selected.annotation_unit_id.eq(main_unit)].iloc[0]
    fig = plt.figure(figsize=(13.2, 10.6))
    gs = fig.add_gridspec(3, 4, height_ratios=[1.05, 1.0, 1.0], hspace=0.25, wspace=0.08)

    for col, unit in enumerate(SELECTED_UNITS):
        row = selected[selected.annotation_unit_id.eq(unit)].iloc[0]
        ax = fig.add_subplot(gs[0, col])
        ax.imshow(rgb(panel_dirs[unit] / "01_source_bbox.png"))
        descriptor = ["purple / upright", "pale / inclined", "pink / inclined"][col]
        panel_label(ax, chr(ord('A') + col), f"{row.taxon_name}\n{descriptor}; YOLO p={row.yolo_conf:.2f}")
    ax = fig.add_subplot(gs[0, 3])
    ax.axis("off")
    ax.text(0.04, 0.92, "Actual production records", fontsize=14, fontweight="bold", va="top")
    ax.text(0.04, 0.78,
            "• source photos used in the executed atlas\n"
            "• stored YOLO visible-capitulum boxes\n"
            "• open photo licences verified\n"
            "• final v2 continuous measurements",
            fontsize=11, va="top", linespacing=1.5)
    ax.text(0.04, 0.28,
            "YOLO defines the observation unit.\nColour, outline and orientation are then\nmeasured by deterministic image processing.",
            fontsize=11, va="top", fontweight="bold")

    bottom_specs = [
        (gs[1, 0], "D", "Tight head crop", "02_tight_crop.png"),
        (gs[1, 1], "E", "Context crop", "03_context_crop.png"),
        (gs[1, 2], "F", "Corolla colour mask", "05_colour_overlay.png"),
        (gs[1, 3], "G", "Outline + convex hull", "06_shape_overlay.png"),
        (gs[2, 0:2], "H", "Orientation relative to image vertical", "04_orientation_overlay.png"),
    ]
    for spec, label, title, filename in bottom_specs:
        ax = fig.add_subplot(spec)
        ax.imshow(rgb(panel_dirs[main_unit] / filename))
        panel_label(ax, label, title)

    ax = fig.add_subplot(gs[2, 2:4])
    ax.axis("off")
    values = [
        ("Orientation", f"{main_row.orientation_angle_degrees:.1f}°", main_row.orientation_status),
        ("Lightness", f"{main_row.corolla_lab_lightness:.1f}", main_row.colour_status),
        ("Chroma", f"{main_row.corolla_lab_chroma:.1f}", main_row.colour_status),
        ("Hue", f"{main_row.corolla_hue_degrees:.1f}°", main_row.colour_status),
        ("Aspect ratio", f"{main_row.shape_aspect_ratio:.2f}", main_row.shape_status),
        ("Circularity", f"{main_row.shape_circularity:.2f}", main_row.shape_status),
        ("Solidity", f"{main_row.shape_solidity:.2f}", main_row.shape_status),
    ]
    ax.text(0.02, 0.96, "I   Continuous outputs for the detailed example", fontsize=13, fontweight="bold", va="top")
    y = 0.82
    for name, value, status in values:
        ax.text(0.05, y, name, fontsize=11, fontweight="bold")
        ax.text(0.48, y, value, fontsize=11)
        ax.text(0.70, y, f"status: {status}", fontsize=10)
        y -= 0.105
    source_note = (
        f"Detailed example: {main_row.taxon_name}; observation {int(main_row.obs_id)}; "
        f"photo {int(main_row.photo_id)}; {main_row.photo_attribution}; "
        f"licence {str(main_row.photo_license_code).upper()}"
    )
    fig.text(0.5, 0.012, source_note, ha="center", va="bottom", fontsize=8.5)

    fig.suptitle("Figure 1. Actual photographs, YOLO capitulum detection and continuous trait extraction", fontsize=17, fontweight="bold", y=0.985)
    fig.savefig(out / "figure1_real_photo_yolo_pipeline.png", dpi=300, bbox_inches="tight")
    fig.savefig(out / "figure1_real_photo_yolo_pipeline.svg", bbox_inches="tight")
    plt.close(fig)

    provenance_cols = [
        "annotation_unit_id", "taxon_name", "obs_id", "photo_id", "medium_image_url",
        "photo_license_code", "photo_attribution", "yolo_conf", "bbox_x1", "bbox_y1", "bbox_x2", "bbox_y2",
        "colour_status", "shape_status", "orientation_status", "orientation_angle_degrees",
        "corolla_lab_lightness", "corolla_lab_chroma", "corolla_hue_degrees",
        "shape_aspect_ratio", "shape_circularity", "shape_solidity",
    ]
    selected[provenance_cols].to_csv(out / "figure1_real_photo_provenance.csv", index=False)
    (out / "figure1_real_photo_build.json").write_text(json.dumps({
        "global_artifact_id": 8099953404,
        "continuous_artifact_id": 8225059018,
        "metadata_artifact_id": 8066010557,
        "selected_units": SELECTED_UNITS,
        "detector_role": "YOLO visible-capitulum detection only",
        "measurement_role": "production v2 deterministic continuous colour, outline, and orientation measurements",
        "detector_accuracy_note": "This figure demonstrates production outputs; independent precision/recall remains a separate manual audit gate.",
    }, indent=2), encoding="utf-8")
    pd.DataFrame(audit_rows).to_csv(out / "figure1_overlay_recomputation_audit.csv", index=False)
    print(json.dumps({"status": "success", "out_dir": str(out), "selected": SELECTED_UNITS}, indent=2))


if __name__ == "__main__":
    main()
