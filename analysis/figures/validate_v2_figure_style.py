#!/usr/bin/env python3
"""Validate canonical physical properties of externally rendered figures."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from PIL import Image

import azami_figstyle as fs

MAIN_STEMS = [
    "Figure_1_v2_measurement_pipeline",
    "Figure_2_v2_geographic_sampling_domain",
    "Figure_3_v2_taxon_mean_information_loss",
    "Figure_4_v2_within_among_scale_atlas",
    "Figure_5_v2_candidate_robustness",
]
SUPPORTING_STEMS = [
    "Figure_S1_v2_image_to_trait_technical_audit",
    "Figure_S2_v2_image_to_trait_perturbation_audit",
    "Figure_S3_v2_endpoint_measurement_support",
    "Figure_S4_v2_sampling_composition_audit",
    "Figure_S5_v2_spatial_diagnostic_surface",
    "Figure_S6_v2_historical_placement_stability",
    "Figure_S7_v2_whole_capitulum_secondary_synthesis",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--figure-dir", type=Path, required=True)
    return parser.parse_args()


def validate_source_contract() -> None:
    here = Path(__file__).resolve().parent
    builders = [here / "build_v2_figures.py", here / "run_image_to_trait_perturbation_audit.py"]
    for path in builders:
        text = path.read_text(encoding="utf-8")
        if "bbox_inches=\"tight\"" in text or "bbox_inches='tight'" in text:
            raise AssertionError(f"tight bounding-box export reintroduced in {path.name}")
        if "dpi=300" in text:
            raise AssertionError(f"independent 300-dpi export reintroduced in {path.name}")
        if "azami_figstyle" not in text:
            raise AssertionError(f"shared figure style is not imported by {path.name}")

    if len(fs.order(split_hue=False)) != 26:
        raise AssertionError("analysis-unit registry is not 26")
    if len(fs.order(split_hue=True, measured_only=True)) != 22:
        raise AssertionError("measured-endpoint registry is not 22")
    if len(fs.PCA_ENDPOINT_KEYS) != 18:
        raise AssertionError("PCA endpoint registry is not 18")
    if len(fs.S7_MATRIX_KEYS) != 17:
        raise AssertionError("S7 matrix registry is not 17")
    if not fs.SCALE_CLASS["not_comparable"]["hatch"]:
        raise AssertionError("not-comparable atlas cells lost their hatch")
    if not fs.STABILITY["unstable"]["hatch"]:
        raise AssertionError("direction-unstable bars lost their hatch")


def validate_rendered_figures(figure_dir: Path) -> list[dict]:
    expected_width_px = round(fs.WIDTH["double"] * fs.DPI_RASTER)
    maximum_height_px = round(fs.MAX_HEIGHT_IN * fs.DPI_RASTER)
    rows: list[dict] = []

    for stem in MAIN_STEMS + SUPPORTING_STEMS:
        png = figure_dir / f"{stem}.png"
        pdf = figure_dir / f"{stem}.pdf"
        if not png.is_file() or not pdf.is_file():
            raise FileNotFoundError(f"missing canonical PNG/PDF pair for {stem}")
        with Image.open(png) as image:
            width_px, height_px = image.size
            dpi = image.info.get("dpi")

        if abs(width_px - expected_width_px) > 2:
            raise AssertionError(
                f"{stem}: width {width_px}px differs from canonical {expected_width_px}px"
            )
        if height_px > maximum_height_px:
            raise AssertionError(
                f"{stem}: height {height_px}px exceeds canonical maximum {maximum_height_px}px"
            )
        if dpi is not None:
            xdpi, ydpi = dpi
            if not (590 <= xdpi <= 610 and 590 <= ydpi <= 610):
                raise AssertionError(f"{stem}: PNG dpi metadata {dpi} is not approximately 600")

        rows.append({
            "stem": stem,
            "width_px": width_px,
            "height_px": height_px,
            "width_in": width_px / fs.DPI_RASTER,
            "height_in": height_px / fs.DPI_RASTER,
            "dpi_metadata": dpi,
        })
    return rows


def main() -> None:
    args = parse_args()
    validate_source_contract()
    rows = validate_rendered_figures(args.figure_dir)
    report = {
        "status": "PASS",
        "target_width_inches": fs.WIDTH["double"],
        "raster_dpi": fs.DPI_RASTER,
        "max_height_inches": fs.MAX_HEIGHT_IN,
        "analysis_units": len(fs.order(split_hue=False)),
        "measured_endpoints": len(fs.order(split_hue=True, measured_only=True)),
        "pca_endpoints": len(fs.PCA_ENDPOINT_KEYS),
        "s7_units": len(fs.S7_MATRIX_KEYS),
        "figures": rows,
    }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
