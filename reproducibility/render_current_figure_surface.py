#!/usr/bin/env python3
"""Build the current GEB Chapter 1 figure surface and checksum manifest.

This is a presentation/release assembler. It does not refit ecological models.
It combines frozen v2 display pieces with current v3 renderers and applies the
frozen GEB Main/SI role map.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import tempfile

from reproducibility.render_construct_environment import render as render_construct_environment
from reproducibility.render_layout_revisions import render as render_layout_revisions
from reproducibility.render_scale_integration import render as render_scale_integration

ROOT = Path(__file__).resolve().parents[1]
FROZEN = ROOT / "reproducibility" / "figures"
AXES = ROOT / "reproducibility" / "current_reference" / "axes"

# manuscript_label, canonical output stem, source kind, source stem
FIGURES = (
    ("Figure 1", "Figure_1_measurement_construct_workflow", "layout", "Figure_1_v2_measurement_pipeline"),
    ("Figure 2", "Figure_2_global_sampling_domain", "frozen", "Figure_2_v2_geographic_sampling_domain"),
    ("Figure 3", "Figure_3_scale_integration", "scale", "Figure_v3_scale_integration"),
    ("Figure 4", "Figure_4_construct_environment", "construct", "Figure_4_construct_environment"),
    ("Figure 5", "Figure_5_ecological_anchor_robustness", "frozen", "Figure_5_v2_candidate_robustness"),
    ("Figure S1.1", "Figure_S1_1_image_to_trait_technical_audit", "frozen", "Figure_S1_v2_image_to_trait_technical_audit"),
    ("Figure S1.2", "Figure_S1_2_image_to_trait_perturbation_audit", "frozen", "Figure_S2_v2_image_to_trait_perturbation_audit"),
    ("Figure S1.3", "Figure_S1_3_endpoint_measurement_support", "frozen", "Figure_S3_v2_endpoint_measurement_support"),
    ("Figure S1.4", "Figure_S1_4_sampling_composition_audit", "frozen", "Figure_S4_v2_sampling_composition_audit"),
    ("Figure S1.5", "Figure_S1_5_spatial_diagnostic_surface", "layout", "Figure_S5_v2_spatial_diagnostic_surface"),
    ("Figure S1.6", "Figure_S1_6_historical_placement_stability", "frozen", "Figure_S6_v2_historical_placement_stability"),
    ("Figure S1.7", "Figure_S1_7_whole_capitulum_secondary_synthesis", "frozen", "Figure_S7_v2_whole_capitulum_secondary_synthesis"),
    ("Figure S2.1", "Figure_S2_1_taxon_mean_information_loss", "frozen", "Figure_3_v2_taxon_mean_information_loss"),
)

PROVENANCE_FILES = (
    ("layout/layout_receipt.json", "layout_receipt.json"),
    ("scale/Figure_v3_scale_integration_provenance.json", "Figure_3_scale_integration_provenance.json"),
    ("construct/construct_environment_receipt.json", "Figure_4_construct_environment_receipt.json"),
    ("frozen/figure_build_report.json", "frozen_figure_build_report.json"),
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def copy(src: Path, dst: Path) -> None:
    if not src.is_file():
        raise FileNotFoundError(src)
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def _source(kind: str, stem: str, ext: str, roots: dict[str, Path]) -> Path:
    if kind == "frozen":
        return FROZEN / f"{stem}.{ext}"
    if kind == "layout":
        return roots["layout"] / f"{stem}.{ext}"
    if kind == "scale":
        return roots["scale"] / f"{stem}.{ext}"
    if kind == "construct":
        return roots["construct"] / f"{stem}.{ext}"
    raise ValueError(kind)


def build(output: Path) -> dict:
    output = output.resolve()
    if output == FROZEN or FROZEN in output.parents:
        raise ValueError("Current figure surface must not overwrite the frozen v2 archive")
    output.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="azami-current-figures-") as td:
        td = Path(td)
        roots = {
            "layout": td / "layout",
            "scale": td / "scale",
            "construct": td / "construct",
            "frozen": FROZEN,
        }
        render_layout_revisions(roots["layout"])
        render_scale_integration(roots["scale"])
        render_construct_environment(AXES, roots["construct"])

        figure_rows = []
        for label, out_stem, kind, source_stem in FIGURES:
            extensions = ("png", "pdf")
            for ext in extensions:
                src = _source(kind, source_stem, ext, roots)
                dst = output / "figures" / f"{out_stem}.{ext}"
                copy(src, dst)
                figure_rows.append({
                    "path": dst.relative_to(output).as_posix(),
                    "sha256": sha256(dst),
                    "manuscript_label": label,
                    "role": "main" if label.startswith("Figure ") and "S" not in label.split()[1] else "supporting",
                    "source_kind": kind,
                    "source_stem": source_stem,
                    "format": ext,
                })

        provenance_map = {
            "layout/layout_receipt.json": roots["layout"] / "layout_receipt.json",
            "scale/Figure_v3_scale_integration_provenance.json":
                roots["scale"] / "Figure_v3_scale_integration_provenance.json",
            "construct/construct_environment_receipt.json":
                roots["construct"] / "construct_environment_receipt.json",
            "frozen/figure_build_report.json": FROZEN / "figure_build_report.json",
        }
        provenance_rows = []
        for key, out_name in PROVENANCE_FILES:
            src = provenance_map[key]
            dst = output / "provenance" / out_name
            copy(src, dst)
            provenance_rows.append({
                "path": dst.relative_to(output).as_posix(),
                "sha256": sha256(dst),
                "role": "figure_provenance",
            })

    all_rows = figure_rows + provenance_rows
    manifest = {
        "schema_version": 1,
        "surface_id": "azami_ch1_geb_current_figure_surface_20260926",
        "scientific_outputs_changed": False,
        "main_figure_count": 5,
        "supporting_figure_count": 8,
        "main_labels": [f"Figure {i}" for i in range(1, 6)],
        "supporting_labels": [f"Figure S1.{i}" for i in range(1, 8)] + ["Figure S2.1"],
        "claim_boundary": (
            "Current GEB display surface assembled from frozen/reproducible numerical outputs; "
            "does not refit models or establish physical-trait accuracy."
        ),
        "document_pagination_validated": False,
        "files": all_rows,
    }
    manifest_path = output / "final_figure_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    receipt = {
        "status": "PASS",
        "manifest": manifest_path.name,
        "manifest_sha256": sha256(manifest_path),
        "figure_files": len(figure_rows),
        "provenance_files": len(provenance_rows),
        "main_figures": 5,
        "supporting_figures": 8,
    }
    (output / "figure_surface_receipt.json").write_text(
        json.dumps(receipt, indent=2) + "\n", encoding="utf-8"
    )
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "work/current-figure-surface")
    args = parser.parse_args()
    print(json.dumps(build(args.output), indent=2))


if __name__ == "__main__":
    main()
