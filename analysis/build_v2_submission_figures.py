#!/usr/bin/env python3
"""Build submission figures only from the frozen full-27 v2 result tables."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch, FancyBboxPatch
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "analysis_outputs" / "v2_full27_environment_atlas_2026-08-27"
DEFAULT_OUTPUT = ROOT / "manuscript" / "figures" / "v2_submission"

ENDPOINT_LABELS = {
    "orientation_image_vertical_angle": "orientation angle",
    "corolla_lab_lightness": "corolla lightness",
    "corolla_lab_chroma": "corolla chroma",
    "corolla_hue_sin": "hue sine",
    "corolla_hue_cos": "hue cosine",
    "capitulum_outline_aspect_ratio": "outline aspect ratio",
    "capitulum_outline_circularity": "outline circularity",
    "capitulum_outline_solidity": "outline solidity",
    "capitulum_width_profile_cv": "width-profile CV",
    "involucre_length_width_ratio": "involucre L/W",
    "involucre_apical_taper_ratio": "apical taper",
    "involucre_basal_taper_ratio": "basal taper",
    "bract_projection_roughness": "projection roughness",
    "bract_projection_p95": "projection p95",
    "bract_projection_maximum": "projection maximum",
    "bract_spread_fraction": "spread fraction",
    "bract_projection_peak_density": "projection-peak density",
    "bract_projection_asymmetry": "projection asymmetry",
    "involucre_surface_edge_density": "surface edge density",
    "involucre_surface_lbp_entropy": "surface LBP entropy",
    "involucre_surface_high_frequency_energy": "surface high-frequency energy",
    "involucre_surface_specular_fraction": "surface specular fraction",
    "corolla_hue": "joint corolla hue",
    "corolla_purple_pixel_fraction": "purple pixel fraction",
    "corolla_redmagenta_pixel_fraction": "red-magenta pixel fraction",
    "corolla_white_pixel_fraction": "white pixel fraction",
    "corolla_yellow_pixel_fraction": "yellow pixel fraction",
    "visible_floret_fraction": "visible floret fraction",
}

PREDICTOR_LABELS = {
    "chelsa_bio01": "BIO1",
    "chelsa_bio04": "BIO4",
    "chelsa_bio12": "BIO12",
    "chelsa_bio15": "BIO15",
    "chelsa_rsds_mean": "radiation",
    "chelsa_vpd_mean": "VPD",
    "chelsa_sfcwind_mean": "wind",
    "chelsa_gsp": "GSP",
    "chelsa_npp": "NPP",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def save(fig: plt.Figure, output: Path, stem: str) -> None:
    fig.savefig(output / f"{stem}.png", dpi=300, bbox_inches="tight", facecolor="white")
    fixed_time = datetime(2026, 8, 27, tzinfo=timezone.utc)
    fig.savefig(
        output / f"{stem}.pdf",
        bbox_inches="tight",
        facecolor="white",
        metadata={
            "Title": stem,
            "Author": "",
            "Creator": "Azami Chapter 1 v2 frozen figure builder",
            "CreationDate": fixed_time,
            "ModDate": fixed_time,
        },
    )
    plt.close(fig)


def figure_method(output: Path) -> None:
    stages = [
        ("Public photographs", "repeated + coordinates"),
        ("27 endpoints", "continuous contract"),
        ("22 measured", "5 explicit unexecuted"),
        ("Present geography", "within / among taxa"),
        ("Robustness gates", "sampling + space + history"),
        ("Two candidates", "mechanisms still open"),
    ]
    fig, ax = plt.subplots(figsize=(13.5, 3.4))
    ax.set_xlim(0, 13.5)
    ax.set_ylim(0, 3)
    ax.axis("off")
    colors = ["#e8f1f8", "#dcefe6", "#f2ecd9", "#e9e4f3", "#f5e2dc", "#dbe8d4"]
    for index, ((title, subtitle), color) in enumerate(zip(stages, colors)):
        x = 0.2 + index * 2.2
        box = FancyBboxPatch(
            (x, 0.85), 1.92, 1.25,
            boxstyle="round,pad=0.08,rounding_size=0.08",
            linewidth=1.2, edgecolor="#314b5f", facecolor=color,
        )
        ax.add_patch(box)
        ax.text(x + 0.96, 1.62, title, ha="center", va="center", fontsize=9.2, weight="bold")
        ax.text(x + 0.96, 1.22, subtitle, ha="center", va="center", fontsize=7.6, color="#374151")
        if index < len(stages) - 1:
            ax.annotate("", xy=(x + 2.17, 1.48), xytext=(x + 1.94, 1.48),
                        arrowprops={"arrowstyle": "->", "lw": 1.4, "color": "#314b5f"})
    ax.text(0.2, 2.55, "From citizen photographs to a gated spatial phenotype atlas", fontsize=14, weight="bold")
    ax.text(0.2, 0.28, "Missing measurement remains missing evidence; candidate patterns are not causal mechanisms.", fontsize=9.5, color="#4b5563")
    save(fig, output, "Figure_1_v2_method_flow")


def figure_information_loss(input_dir: Path, output: Path) -> None:
    frame = pd.read_csv(input_dir / "v2_full27_variance_decomposition.csv")
    frame = frame.loc[frame["status"].eq("ok")].copy()
    frame = frame.sort_values("equal_replication_fraction_median")
    labels = [ENDPOINT_LABELS.get(value, value) for value in frame["endpoint_id"]]
    y = np.arange(len(frame))
    fig, ax = plt.subplots(figsize=(9.2, 8.5))
    ax.barh(y + 0.18, frame["fraction_below_taxon_means"] * 100, height=0.34,
            color="#8fb9d7", label="raw observed sample")
    ax.barh(y - 0.18, frame["equal_replication_fraction_median"] * 100, height=0.34,
            color="#d9a66f", label="two observations per taxon")
    ax.set_yticks(y, labels, fontsize=8.4)
    ax.set_xlim(0, 100)
    ax.set_xlabel("Visible variation below source-assigned taxon means (%)")
    ax.set_title("Taxon means compress repeated image-phenotype variation", loc="left", weight="bold")
    ax.grid(axis="x", color="#d1d5db", linewidth=0.6)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, loc="lower right")
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    save(fig, output, "Figure_2_v2_taxon_mean_information_loss")


def figure_scale_atlas(input_dir: Path, output: Path) -> None:
    frame = pd.read_csv(input_dir / "v2_full27_environment_cross_scale.csv")
    units = list(dict.fromkeys(frame["unit_id"]))
    predictors = list(PREDICTOR_LABELS)
    code = {"among_only": 0, "neither": 1, "within_only": 2, "both_scales": 3, "not_comparable": 4}
    matrix = np.full((len(units), len(predictors)), 4, dtype=int)
    for _, row in frame.iterrows():
        if row["unit_id"] in units and row["predictor"] in predictors:
            matrix[units.index(row["unit_id"]), predictors.index(row["predictor"])] = code[row["cross_scale_class"]]
    cmap = ListedColormap(["#b66a75", "#ececec", "#5f95bd", "#59488f", "#ffffff"])
    fig, ax = plt.subplots(figsize=(10.5, 9.5))
    ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=-0.5, vmax=4.5)
    ax.set_xticks(np.arange(len(predictors)), [PREDICTOR_LABELS[p] for p in predictors], rotation=40, ha="right")
    ax.set_yticks(np.arange(len(units)), [ENDPOINT_LABELS.get(u, u) for u in units], fontsize=8.2)
    ax.set_title("Environmental alignment differs within and among taxa", loc="left", weight="bold")
    ax.set_xlabel("Environmental gradient")
    ax.set_ylabel("Continuous endpoint or joint unit")
    ax.set_xticks(np.arange(-0.5, len(predictors), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(units), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.8)
    ax.tick_params(which="minor", bottom=False, left=False)
    legend = [
        Patch(facecolor="#59488f", label="both scales"),
        Patch(facecolor="#5f95bd", label="within only"),
        Patch(facecolor="#b66a75", label="among only"),
        Patch(facecolor="#ececec", label="neither"),
        Patch(facecolor="white", edgecolor="#bbbbbb", label="not comparable"),
    ]
    ax.legend(handles=legend, frameon=False, ncol=3, bbox_to_anchor=(0, -0.13), loc="upper left")
    save(fig, output, "Figure_3_v2_within_among_scale_atlas")


def figure_candidates(input_dir: Path, output: Path) -> None:
    among = pd.read_csv(input_dir / "spatial" / "v2_full27_spatial_among.csv")
    selected = among.loc[among["broad_spatial_sensitivity_pass"].astype(str).str.lower().eq("true")].copy()
    order = ["corolla_lab_chroma", "orientation_image_vertical_angle"]
    selected["order"] = selected["unit_id"].map({value: index for index, value in enumerate(order)})
    selected = selected.sort_values("order")
    labels = ["chroma ~\nradiation", "orientation ~\nprecipitation"]
    y = np.arange(2)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.2), gridspec_kw={"width_ratios": [1.35, 1]})
    ax1.axvline(0, color="#6b7280", linewidth=0.8)
    ax1.barh(y + 0.16, selected["base_beta_std"], height=0.28, color="#8fb9d7", label="global among-taxon")
    ax1.barh(y - 0.16, selected["spatial_beta_std"], height=0.28, color="#d9a66f", label="broad-space")
    ax1.set_yticks(y, labels)
    ax1.invert_yaxis()
    ax1.set_xlabel("Standardized coefficient")
    ax1.set_title("A  Effect direction survives broad-space control", loc="left", weight="bold", fontsize=12)
    ax1.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.54, -0.16), ncol=2)
    ax1.grid(axis="x", color="#d1d5db", linewidth=0.6)
    ax1.set_axisbelow(True)

    ratios = [0.614876, 0.780783]
    ax2.barh(y, ratios, color=["#7e6aa2", "#5d9471"], height=0.48)
    ax2.set_yticks(y, labels)
    ax2.invert_yaxis()
    ax2.set_xlim(0, 1.05)
    ax2.set_xlabel("Minimum effect-magnitude ratio")
    ax2.set_title("B  Sampling direction stable", loc="left", weight="bold", fontsize=12)
    for index, value in enumerate(ratios):
        ax2.text(value + 0.025, index, f"{value:.3f}\n52/52 trees", va="center", fontsize=9)
    ax2.grid(axis="x", color="#d1d5db", linewidth=0.6)
    ax2.set_axisbelow(True)
    fig.suptitle("Two adaptive-pattern candidates under current sequential controls", x=0.06, y=0.97, ha="left", fontsize=14, weight="bold")
    fig.subplots_adjust(left=0.22, right=0.98, bottom=0.24, top=0.80, wspace=0.48)
    fig.text(0.06, 0.055, "Stability across tested sampling, broad-space and placement uncertainty does not demonstrate adaptation or mechanism.", fontsize=9, color="#4b5563")
    save(fig, output, "Figure_4_v2_candidate_robustness")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    figure_method(args.output_dir)
    figure_information_loss(args.input_dir, args.output_dir)
    figure_scale_atlas(args.input_dir, args.output_dir)
    figure_candidates(args.input_dir, args.output_dir)
    inputs = [
        args.input_dir / "v2_full27_variance_decomposition.csv",
        args.input_dir / "v2_full27_environment_cross_scale.csv",
        args.input_dir / "spatial" / "v2_full27_spatial_among.csv",
        args.input_dir / "historical" / "v2_full27_historical_placement_summary.csv",
        args.input_dir / "sampling" / "v2_full27_sampling_composition_summary.csv",
    ]
    report = {
        "status": "ok",
        "source_lane": "full27_full_environment_only",
        "input_sha256": {path.relative_to(ROOT).as_posix(): sha256(path) for path in inputs},
        "figures": sorted(path.name for path in args.output_dir.glob("Figure_*.png")),
        "claim_boundary": "submission visualization of frozen exploratory v2 results; not a new analysis",
    }
    (args.output_dir / "figure_build_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
