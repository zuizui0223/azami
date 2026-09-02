#!/usr/bin/env python3
"""Build submission figures from frozen Chapter 1 v2 result tables.

Only frozen outputs and presentation/source-provenance assets are reshaped.
No model, P value, or additional hypothesis is fitted or selected here.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LogNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, FancyBboxPatch
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "analysis_outputs" / "v2_full27_environment_atlas_2026-08-27"
DEFAULT_OUTPUT = ROOT / "manuscript" / "figures" / "v2_submission"
SOURCE = DEFAULT_OUTPUT / "source"

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
    "chelsa_bio01": "BIO1", "chelsa_bio04": "BIO4", "chelsa_bio12": "BIO12",
    "chelsa_bio15": "BIO15", "chelsa_rsds_mean": "radiation", "chelsa_vpd_mean": "VPD",
    "chelsa_sfcwind_mean": "wind", "chelsa_gsp": "GSP", "chelsa_npp": "NPP",
}
MAIN_STEMS = [
    "Figure_1_v2_measurement_pipeline",
    "Figure_2_v2_geographic_sampling_domain",
    "Figure_3_v2_taxon_mean_information_loss",
    "Figure_4_v2_within_among_scale_atlas",
    "Figure_5_v2_candidate_robustness",
]
SUPP_STEMS = [
    "Figure_S3_v2_endpoint_measurement_support",
    "Figure_S4_v2_sampling_composition_audit",
    "Figure_S5_v2_spatial_diagnostic_surface",
    "Figure_S6_v2_historical_placement_stability",
    "Figure_S7_v2_whole_capitulum_secondary_synthesis",
]
OBSOLETE_STEMS = [
    "Figure_1_v2_method_flow", "Figure_2_v2_taxon_mean_information_loss",
    "Figure_3_v2_within_among_scale_atlas", "Figure_4_v2_candidate_robustness",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def save(fig: plt.Figure, output: Path, stem: str) -> None:
    fig.savefig(output / f"{stem}.png", dpi=300, bbox_inches="tight", facecolor="white")
    fixed_time = datetime(2026, 8, 28, tzinfo=timezone.utc)
    fig.savefig(output / f"{stem}.pdf", bbox_inches="tight", facecolor="white", metadata={
        "Title": stem, "Author": "", "Creator": "Azami Chapter 1 v2 frozen figure builder",
        "CreationDate": fixed_time, "ModDate": fixed_time,
    })
    plt.close(fig)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def box(ax: plt.Axes, xy: tuple[float, float], width: float, height: float,
        title: str, subtitle: str, color: str, fontsize: float = 9.0,
        text_color: str = "#111111", pad: float = 0.025,
        rounding_size: float = 0.025) -> None:
    x, y = xy
    ax.add_patch(FancyBboxPatch(
        (x, y), width, height, boxstyle=f"round,pad={pad},rounding_size={rounding_size}",
        linewidth=1.0, edgecolor="#314b5f", facecolor=color,
    ))
    ax.text(x + width / 2, y + height * 0.63, title, ha="center", va="center", weight="bold", fontsize=fontsize, color=text_color)
    ax.text(x + width / 2, y + height * 0.27, subtitle, ha="center", va="center", color=text_color, fontsize=fontsize - 1.3)


def figure_measurement(output: Path) -> None:
    image_path = SOURCE / "Figure_1_real_photo_measurement_provenance.png"
    provenance_path = SOURCE / "Figure_1_real_photo_measurement_provenance.json"
    assert json.loads(provenance_path.read_text(encoding="utf-8"))
    source_image = plt.imread(image_path)
    fig = plt.figure(figsize=(6.6, 8.6))
    grid = fig.add_gridspec(4, 1, height_ratios=[2.25, 2.15, 2.2, 0.78], hspace=0.34)

    photo_grid = grid[0].subgridspec(1, 4, wspace=0.09)
    photos = [
        (source_image[390:1100, 90:650], "(a) $\\it{Cirsium\\ pannonicum}$\npurple; upright"),
        (source_image[390:1100, 760:1450], "(b) $\\it{Cirsium\\ mohavense}$\npale; inclined"),
        (source_image[505:1015, 1510:2200], "(c) $\\it{Cirsium\\ subcoriaceum}$\npink; inclined"),
    ]
    for index, (photo, title) in enumerate(photos):
        ax = fig.add_subplot(photo_grid[index])
        ax.imshow(photo)
        ax.set_title(title, fontsize=7.4, pad=3)
        ax.axis("off")
    records = fig.add_subplot(photo_grid[3])
    records.axis("off")
    records.text(0, 0.96, "(d) Why localize?", fontsize=8.2, weight="bold", va="top")
    records.text(0, 0.79, "A full scene may contain\nleaves, stems, background\nand several heads.", fontsize=7.4, va="top", linespacing=1.35)
    records.text(0, 0.30, "YOLO fixes one capitulum\nas the measurement unit.\nIt does not score traits.", fontsize=7.2, weight="bold", va="top", linespacing=1.28)

    crop_grid = grid[1].subgridspec(1, 4, wspace=0.09)
    crops = [
        (source_image[1300:2010, 100:610], "(e) Tight crop"),
        (source_image[1300:2010, 860:1370], "(f) Context crop"),
        (source_image[1300:2010, 1580:2100], "(g) Corolla mask"),
        (source_image[1300:2010, 2340:2860], "(h) Outline + hull"),
    ]
    for index, (photo, title) in enumerate(crops):
        ax = fig.add_subplot(crop_grid[index])
        ax.imshow(photo)
        ax.set_title(title, fontsize=7.6, pad=3)
        ax.axis("off")

    output_grid = grid[2].subgridspec(1, 2, width_ratios=[0.82, 1.18], wspace=0.10)
    orientation = fig.add_subplot(output_grid[0])
    orientation.imshow(source_image[2190:2860, 480:1000])
    orientation.set_title("(i) Orientation relative to image vertical", fontsize=7.6, pad=3)
    orientation.axis("off")
    values = fig.add_subplot(output_grid[1])
    values.axis("off")
    values.text(0, 0.98, "(j) Continuous outputs", fontsize=8.3, weight="bold", va="top")
    rows = [
        ("Orientation", "0.7°"), ("Lightness", "56.5"), ("Chroma", "33.3"),
        ("Hue", "300.2°"), ("Aspect ratio", "1.59"), ("Circularity", "0.77"), ("Solidity", "0.95"),
    ]
    for index, (label, value) in enumerate(rows):
        y = 0.83 - index * 0.105
        values.text(0.02, y, label, fontsize=7.8, weight="bold")
        values.text(0.56, y, value, fontsize=7.8)
        values.text(0.78, y, "usable", fontsize=7.4, color="#374151")
    values.text(0.02, 0.045, "Detailed example: $\\it{Cirsium\\ pannonicum}$; observation 295734197; photo 532420148; CC BY.", fontsize=6.7, color="#4b5563", wrap=True)

    footer = fig.add_subplot(grid[3])
    footer.set(xlim=(0, 1), ylim=(0, 1))
    footer.axis("off")
    cards = [
        ("27 registered", "continuous endpoints", "#e8f1f8"),
        ("22 measured", "finite v2 outputs", "#dcefe6"),
        ("5 unexecuted", "remain missing evidence", "#f2ecd9"),
        ("Detector = crop", "does not score traits", "#e9e4f3"),
    ]
    for i, (title, subtitle, color) in enumerate(cards):
        box(footer, (0.015 + i * 0.247, 0.31), 0.225, 0.54, title, subtitle, color, 7.3,
            pad=0.006, rounding_size=0.014)
    footer.text(0.015, 0.08,
        "Three open-licensed photographs are presentation/source provenance only; the v2 cohort and results come exclusively from the prespecified full-27 lane.",
        fontsize=6.5, color="#4b5563")
    save(fig, output, MAIN_STEMS[0])


def draw_world(ax: plt.Axes) -> None:
    payload = json.loads((SOURCE / "natural_earth_110m_land_outline.json").read_text(encoding="utf-8"))
    for coords in payload["parts"]:
        values = np.asarray(coords, dtype=float)
        longitude = np.deg2rad(values[:, 0])
        latitude = np.deg2rad(values[:, 1])
        ax.fill(longitude, latitude, facecolor="#f1f1ed", edgecolor="#b6b8b5", linewidth=0.35, zorder=0)
    ax.set_longitude_grid(60)
    ax.set_latitude_grid(30)
    ax.grid(color="#e0e3e5", linewidth=0.4, zorder=-1)
    scale_x0, scale_x1, scale_y = 0.08, 0.13, 0.12
    ax.plot([scale_x0, scale_x1], [scale_y, scale_y], transform=ax.transAxes, color="#30343b", linewidth=2.2, zorder=4)
    ax.plot([scale_x0, scale_x0], [scale_y - 0.012, scale_y + 0.012], transform=ax.transAxes, color="#30343b", linewidth=1.1, zorder=4)
    ax.plot([scale_x1, scale_x1], [scale_y - 0.012, scale_y + 0.012], transform=ax.transAxes, color="#30343b", linewidth=1.1, zorder=4)
    ax.text((scale_x0 + scale_x1) / 2, scale_y + 0.027, "2,000 km at Equator", transform=ax.transAxes,
        ha="center", va="bottom", fontsize=7.5, color="#30343b")


def figure_geographic_domain(output: Path) -> None:
    cells = pd.read_csv(SOURCE / "Figure_2_v2_realized_sampling_cells.csv")
    assert int(cells["n_observations"].sum()) == 46_276
    fig = plt.figure(figsize=(6.6, 8.1))
    # A dedicated colour-bar column keeps all three content panels on the
    # same left and right boundaries. Automatic colour bars previously
    # narrowed panels (a) and (c) by different amounts.
    grid = fig.add_gridspec(
        3, 2, height_ratios=[3.1, 1.28, 3.1], width_ratios=[24, 0.72],
        hspace=0.43, wspace=0.08, left=0.105, right=0.93, top=0.985, bottom=0.085,
    )
    ax_map = fig.add_subplot(grid[0, 0], projection="mollweide")
    draw_world(ax_map)
    points = ax_map.scatter(np.deg2rad(cells["cell_lon_center"]), np.deg2rad(cells["cell_lat_center"]), c=cells["n_observations"],
        cmap="viridis", norm=LogNorm(vmin=1, vmax=max(10, cells["n_observations"].max())),
        s=13, marker="s", linewidth=0, alpha=0.92, zorder=2)
    cb = fig.colorbar(points, cax=fig.add_subplot(grid[0, 1]))
    cb.set_label("Observations per 2° cell")
    ax_map.set_title("(a) Realized geographic sampling domain", loc="left", weight="bold", fontsize=9.5)

    ax_flow = fig.add_subplot(grid[1, 0])
    ax_flow.set(xlim=(0, 1), ylim=(0, 1))
    ax_flow.axis("off")
    steps = [
        ("Detector-positive", "406,582 obs / 286 taxa"), ("Coordinates", "392,989 / 271"),
        ("Accuracy ≤10 km", "297,293 / 259"), ("Taxon × 0.25°", "46,276 / 259"),
    ]
    box_width = 0.19
    box_gap = (0.98 - len(steps) * box_width) / (len(steps) - 1)
    box_x = [0.01 + i * (box_width + box_gap) for i in range(len(steps))]
    for i, (x, (title, subtitle)) in enumerate(zip(box_x, steps)):
        box(ax_flow, (x, 0.32), box_width, 0.42, title, subtitle,
            ["#e8f1f8", "#dcefe6", "#f2ecd9", "#e9e4f3"][i], 7.2,
            pad=0.006, rounding_size=0.014)
    for x, next_x in zip(box_x[:-1], box_x[1:]):
        ax_flow.annotate("", xy=(next_x - 0.012, 0.53), xytext=(x + box_width + 0.012, 0.53),
            arrowprops={"arrowstyle": "-|>", "lw": 1.25, "color": "#314b5f", "mutation_scale": 10}, zorder=6)
    ax_flow.text(0.01, 0.92, "(b) Prespecified cohort flow", fontsize=9.5, weight="bold")
    ax_flow.text(0.01, 0.08, "Europe + North America: 92.26% of retained observations", fontsize=8, color="#4b5563")

    ax_env = fig.add_subplot(grid[2, 0])
    env = ax_env.scatter(cells["bio1_degrees_c_mean"], cells["bio12_mm_mean"],
        s=np.clip(np.sqrt(cells["n_observations"].clip(lower=1)) * 2.2, 5, 70),
        c=cells["cell_lat_center"], cmap="coolwarm", alpha=0.55, linewidth=0)
    cb2 = fig.colorbar(env, cax=fig.add_subplot(grid[2, 1]))
    cb2.set_label("Cell latitude")
    ax_env.set(xlabel="Annual mean temperature, BIO1 (°C)", ylabel="Annual precipitation, BIO12 (mm)")
    ax_env.set_title("(c) Realized climatic domain", loc="left", weight="bold", fontsize=9.5)
    ax_env.grid(color="#e5e7eb", linewidth=0.5)
    ax_env.set_axisbelow(True)
    # The Mollweide footprint is slightly narrower than its allocated grid
    # cell. Inset the rectangular panels by the same amount so (a)-(c) have
    # equal visible content widths without distorting the map projection.
    for aligned_ax in (ax_flow, ax_env):
        position = aligned_ax.get_position()
        inset = position.width * 0.038
        aligned_ax.set_position([position.x0 + inset, position.y0, position.width - 2 * inset, position.height])
    fig.text(0.105, 0.016, "The panels show where retained records occur; they do not establish representativeness.", fontsize=7.5, color="#4b5563")
    save(fig, output, MAIN_STEMS[1])


def figure_information_loss(input_dir: Path, output: Path) -> None:
    frame = pd.read_csv(input_dir / "v2_full27_variance_decomposition.csv")
    frame = frame.loc[frame["status"].eq("ok")].copy().sort_values("equal_replication_fraction_median")
    assert len(frame) == 22
    labels = [ENDPOINT_LABELS.get(value, value) for value in frame["endpoint_id"]]
    y = np.arange(len(frame))
    raw_within = frame["fraction_below_taxon_means"].to_numpy() * 100
    balanced_within = frame["equal_replication_fraction_median"].to_numpy() * 100
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 6.7), sharey=True)
    panel_data = [
        (axes[0], raw_within, "(a) Raw observed sample"),
        (axes[1], balanced_within, "(b) Two observations per taxon"),
    ]
    for ax, within, title in panel_data:
        among = 100 - within
        ax.barh(y, within, height=0.68, color="#6f9fc1", edgecolor="white", linewidth=0.25)
        ax.barh(y, among, left=within, height=0.68, color="#dedede", edgecolor="white", linewidth=0.25)
        ax.set(xlim=(0, 100), xticks=[0, 25, 50, 75, 100])
        # Reserve a blank row below the lowest bar for the range annotation.
        # Keeping the box in data space prevents it from covering the hue-sine
        # bar when figure dimensions or tight bounding boxes change.
        ax.set_ylim(-2.35, len(frame) - 0.3)
        ax.set_title(title, loc="left", weight="bold", fontsize=8.7, pad=8)
        ax.grid(axis="x", color="#d1d5db", linewidth=0.55)
        ax.set_axisbelow(True)
        ax.text(2.0, -1.35, f"Within-taxon range\n{within.min():.1f}–{within.max():.1f}%",
            ha="left", va="center", fontsize=7.0, color="#374151",
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "#c7cbd1", "alpha": 0.96})
        for spine in ("top", "right", "left"):
            ax.spines[spine].set_visible(False)
    axes[0].set_yticks(y, labels, fontsize=7.4)
    axes[1].tick_params(axis="y", left=False, labelleft=False)
    fig.supxlabel("Share of total observed image-phenotype variation (%)", y=0.085, fontsize=9)
    fig.legend(
        handles=[Patch(facecolor="#6f9fc1", label="within source-assigned taxa"), Patch(facecolor="#dedede", label="among source-assigned taxa")],
        loc="lower center", bbox_to_anchor=(0.61, 0.035), frameon=False, ncol=2, fontsize=8.0,
    )
    fig.subplots_adjust(left=0.33, right=0.985, top=0.95, bottom=0.16, wspace=0.16)
    fig.text(0.33, 0.012, "Each bar sums to 100%. Measurement, photography and biological variation are not separated.", fontsize=7.0, color="#4b5563")
    save(fig, output, MAIN_STEMS[2])


def figure_scale_atlas(input_dir: Path, output: Path) -> None:
    frame = pd.read_csv(input_dir / "v2_full27_environment_cross_scale.csv")
    counts = frame["cross_scale_class"].value_counts().to_dict()
    expected = {"both_scales": 7, "within_only": 18, "among_only": 3, "neither": 152, "not_comparable": 54}
    assert len(frame) == 234 and counts == expected
    units = list(dict.fromkeys(frame["unit_id"]))
    predictors = list(PREDICTOR_LABELS)
    code = {"among_only": 0, "neither": 1, "within_only": 2, "both_scales": 3, "not_comparable": 4}
    matrix = np.full((len(units), len(predictors)), 4, dtype=int)
    for _, row in frame.iterrows():
        matrix[units.index(row["unit_id"]), predictors.index(row["predictor"])] = code[row["cross_scale_class"]]
    cmap = ListedColormap(["#b66a75", "#ececec", "#5f95bd", "#59488f", "#ffffff"])
    fig, ax = plt.subplots(figsize=(6.6, 6.9))
    ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=-0.5, vmax=4.5)
    ax.set_xticks(np.arange(len(predictors)), [PREDICTOR_LABELS[p] for p in predictors], rotation=40, ha="right", fontsize=7.5)
    ax.set_yticks(np.arange(len(units)), [ENDPOINT_LABELS.get(u, u) for u in units], fontsize=7.2)
    ax.set(xlabel="Environmental gradient", ylabel="Continuous endpoint or joint circular unit")
    ax.set_xticks(np.arange(-0.5, len(predictors), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(units), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.8)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.legend(handles=[
        Patch(facecolor="#59488f", label=f"both scales ({counts['both_scales']})"),
        Patch(facecolor="#5f95bd", label=f"within only ({counts['within_only']})"),
        Patch(facecolor="#b66a75", label=f"among only ({counts['among_only']})"),
        Patch(facecolor="#ececec", label=f"neither ({counts['neither']})"),
        Patch(facecolor="white", edgecolor="#bbbbbb", label=f"not comparable ({counts['not_comparable']})"),
    ], frameon=False, ncol=2, bbox_to_anchor=(0, -0.14), loc="upper left", fontsize=8.0,
       handlelength=1.7, columnspacing=1.1, labelspacing=0.5)
    fig.subplots_adjust(left=0.29, right=0.985, top=0.985, bottom=0.18)
    save(fig, output, MAIN_STEMS[3])


def figure_candidates(input_dir: Path, output: Path) -> None:
    spatial = pd.read_csv(input_dir / "spatial" / "v2_full27_spatial_among.csv")
    selected = spatial.loc[spatial["broad_spatial_sensitivity_pass"].astype(str).str.lower().eq("true")].copy()
    sampling = pd.read_csv(input_dir / "sampling" / "v2_full27_sampling_composition_summary.csv")
    historical = pd.read_csv(input_dir / "historical" / "v2_full27_historical_placement_summary.csv")
    order = [("corolla_lab_chroma", "chelsa_rsds_mean"), ("orientation_image_vertical_angle", "chelsa_bio12")]
    assert len(selected) == 2 and len(historical) == 2
    rows = []
    for unit, predictor in order:
        s = selected.loc[(selected["unit_id"] == unit) & (selected["predictor"] == predictor)].iloc[0]
        r = sampling.loc[(sampling["scale"] == "among_taxon") & (sampling["unit_id"] == unit) & (sampling["predictor"] == predictor)].iloc[0]
        h = historical.loc[(historical["unit_id"] == unit) & (historical["predictor"] == predictor)].iloc[0]
        ratios = [float(r[c]) for c in r.index if c.endswith("minimum_effect_magnitude_ratio") and pd.notna(r[c]) and str(r[c]) != ""]
        rows.append((s, min(ratios), h))

    fig = plt.figure(figsize=(6.6, 8.1))
    grid = fig.add_gridspec(3, 1, height_ratios=[2.0, 2.35, 2.6], hspace=0.38)
    flow = fig.add_subplot(grid[0])
    flow.set(xlim=(0, 1), ylim=(0, 1))
    flow.axis("off")
    steps = [("Canonical grid", "234 rows"), ("Among-taxon BH", "10 supported"),
        ("Sampling audit", "10/10 direction-stable"), ("Broad-space gate", "2 pass"),
        ("Historical gate", "2 × 52/52 trees"), ("Candidate ceiling", "2 patterns")]
    colors = ["#e8f1f8", "#e8f1f8", "#dcefe6", "#f2ecd9", "#e9e4f3", "#dbe8d4"]
    positions = [(0.02, 0.57), (0.36, 0.57), (0.70, 0.57), (0.70, 0.16), (0.36, 0.16), (0.02, 0.16)]
    box_width = 0.28
    box_height = 0.23
    for (x, y0), ((title, subtitle), color) in zip(positions, zip(steps, colors)):
        box(flow, (x, y0), box_width, box_height, title, subtitle, color, 7.4,
            pad=0.006, rounding_size=0.014)
    for start, end in [(0, 1), (1, 2), (3, 4), (4, 5)]:
        x0, y0 = positions[start]
        x1, y1 = positions[end]
        if y0 == y1 and x1 > x0:
            xytext, xy = (x0 + box_width + 0.012, y0 + box_height / 2), (x1 - 0.012, y1 + box_height / 2)
        else:
            xytext, xy = (x0 - 0.012, y0 + box_height / 2), (x1 + box_width + 0.012, y1 + box_height / 2)
        flow.annotate("", xy=xy, xytext=xytext,
            arrowprops={"arrowstyle": "-|>", "lw": 1.2, "color": "#314b5f", "mutation_scale": 10}, zorder=6)
    x2, y2 = positions[2]
    x3, y3 = positions[3]
    flow.annotate("", xy=(x3 + box_width / 2, y3 + box_height + 0.012), xytext=(x2 + box_width / 2, y2 - 0.012),
        arrowprops={"arrowstyle": "-|>", "lw": 1.2, "color": "#314b5f", "mutation_scale": 10}, zorder=6)
    flow.text(0.005, 0.93, "(a) Sequential attrition: why only two among-taxon rows remain", fontsize=9.5, weight="bold")
    flow.text(0.005, 0.02, "Sampling stability annotates rows; promotion requires broad-space and historical-placement gates.", fontsize=7.2, color="#4b5563")

    bars = fig.add_subplot(grid[1])
    y = np.arange(2)
    labels = ["lower chroma ~\nhigher radiation", "larger orientation angle ~\nhigher precipitation"]
    bars.axvline(0, color="#6b7280", linewidth=0.8)
    bars.barh(y + 0.16, [float(r[0]["base_beta_std"]) for r in rows], height=0.28, color="#8fb9d7", label="global among-taxon")
    bars.barh(y - 0.16, [float(r[0]["spatial_beta_std"]) for r in rows], height=0.28, color="#d9a66f", label="broad-space")
    bars.set_yticks(y, labels)
    bars.invert_yaxis()
    bars.set_xlabel("Standardized coefficient")
    bars.set_title("(b) Direction after broad-space control", loc="left", weight="bold", fontsize=9.5, pad=25)
    bars.grid(axis="x", color="#d1d5db", linewidth=0.6)
    bars.set_axisbelow(True)
    bars.legend(frameon=False, loc="lower center", bbox_to_anchor=(0.5, 1.005), ncol=2, fontsize=8)

    gates = fig.add_subplot(grid[2])
    gates.set(xlim=(0, 1), ylim=(-0.6, 1.85))
    gates.axis("off")
    gates.text(0, 1.72, "(c) Gate evidence retained in the fixed ledger", fontsize=9.5, weight="bold")
    for i, ((s, ratio, h), label) in enumerate(zip(rows, labels)):
        y0 = 1.05 - i * 1.05
        gates.add_patch(FancyBboxPatch((0.0, y0 - 0.30), 0.98, 0.76, boxstyle="round,pad=0.02", facecolor="#f7f8f8", edgecolor="#c8cdd1"))
        gates.text(0.03, y0 + 0.31, label.replace("\n", " "), fontsize=8.2, weight="bold")
        gates.text(0.03, y0 + 0.08, f"global β = {float(s['base_beta_std']):+.3f}; q = {float(s['base_q_fdr_bh_global_family']):.3f}", fontsize=7.8)
        gates.text(0.03, y0 - 0.10, f"sampling minimum ratio = {ratio:.3f}; broad-space β = {float(s['spatial_beta_std']):+.3f}", fontsize=7.8)
        gates.text(0.03, y0 - 0.27, f"spatial P = {float(s['spatial_permutation_p_value']):.3f}; residual Moran P = {float(s['residual_morans_p_value']):.3f}; trees = {int(h['n_placement_trees_p_lt_0_05'])}/{int(h['n_successful_placement_trees'])}", fontsize=7.6)
    fig.text(0.07, 0.012, "Stability under these gates does not demonstrate adaptation, convergence, historical cause or mechanism.", fontsize=7.2, color="#4b5563")
    save(fig, output, MAIN_STEMS[4])


def figure_s1_endpoint_support(input_dir: Path, output: Path) -> None:
    frame = pd.read_csv(input_dir / "v2_full27_endpoint_inventory.csv")
    frame = frame.loc[frame["measurement_status"].eq("measured")].copy().sort_values("n_observations_measured")
    labels = [ENDPOINT_LABELS.get(value, value) for value in frame["endpoint_id"]]
    y = np.arange(len(frame))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.6, 7.0), gridspec_kw={"width_ratios": [1.35, 1]})
    ax1.barh(y, frame["n_observations_measured"], height=0.62, color="#7ea6c4")
    ax1.set_yticks(y, labels, fontsize=7.2)
    ax1.set_xscale("log")
    ax1.set(xlabel="Observations (log scale)")
    ax1.set_title("(a) Observation support", loc="left", weight="bold")
    ax1.grid(axis="x", color="#d1d5db", linewidth=0.6)
    ax1.set_axisbelow(True)
    ax2.barh(y, frame["n_taxa_measured"], height=0.62, color="#b18a6a")
    ax2.set_yticks(y, [""] * len(y))
    ax2.set(xlabel="Taxa measured")
    ax2.set_title("(b) Taxonomic support", loc="left", weight="bold")
    ax2.grid(axis="x", color="#d1d5db", linewidth=0.6)
    ax2.set_axisbelow(True)
    fig.subplots_adjust(left=0.37, right=0.98, bottom=0.09, top=0.96, wspace=0.08)
    save(fig, output, SUPP_STEMS[0])


def figure_s2_sampling_audit(input_dir: Path, output: Path) -> None:
    frame = pd.read_csv(input_dir / "sampling" / "v2_full27_sampling_composition_summary.csv")
    ratio_cols = ["dominant_taxon_omission_minimum_effect_magnitude_ratio",
        "leave_one_broad_region_out_minimum_effect_magnitude_ratio", "native_only_minimum_effect_magnitude_ratio",
        "equal_taxon_weight_minimum_effect_magnitude_ratio"]
    for col in ratio_cols:
        frame[col] = pd.to_numeric(frame[col], errors="coerce")
    frame["minimum_ratio"] = frame[ratio_cols].min(axis=1, skipna=True)
    frame["stable"] = frame["all_directions_stable_where_evaluable"].astype(str).str.lower().eq("true")
    frame["label"] = frame.apply(lambda row: f"{'A' if row['scale']=='among_taxon' else 'W'}: {ENDPOINT_LABELS.get(row['unit_id'], row['unit_id'])} ~ {PREDICTOR_LABELS.get(row['predictor'], row['predictor'])}", axis=1)
    fig, axes = plt.subplots(
        2, 1, figsize=(6.6, 7.6), sharex=True,
        gridspec_kw={"height_ratios": [26, 10], "hspace": 0.12},
    )
    panels = [("within_taxon", "(a) Within taxa"), ("among_taxon", "(b) Among taxa")]
    for ax, (scale, title) in zip(axes, panels):
        subset = frame.loc[frame["scale"].eq(scale)].copy()
        subset = subset.sort_values(["stable", "minimum_ratio"], ascending=[True, True]).reset_index(drop=True)
        y = np.arange(len(subset))
        colors = np.where(subset["stable"], "#668f78", "#b86d6d")
        ax.barh(y, subset["minimum_ratio"], color=colors, height=0.64)
        ax.axvline(1.0, color="#6b7280", linewidth=0.8, linestyle="--")
        ax.set_yticks(y, subset["label"].str.replace(r"^[AW]: ", "", regex=True), fontsize=7.2)
        ax.set_title(title, loc="left", weight="bold", fontsize=9.5)
        ax.set_xlim(0, 1.055)
        ax.grid(axis="x", color="#d1d5db", linewidth=0.6)
        ax.set_axisbelow(True)
        for spine in ("top", "right", "left"):
            ax.spines[spine].set_visible(False)
    axes[-1].set_xlabel("Minimum effect ratio across perturbations")
    axes[-1].legend(
        handles=[Patch(facecolor="#668f78", label="direction stable"), Patch(facecolor="#b86d6d", label="direction unstable")],
        frameon=False, loc="lower right", fontsize=7.5,
    )
    fig.subplots_adjust(left=0.43, right=0.98, bottom=0.08, top=0.98)
    save(fig, output, SUPP_STEMS[1])


def figure_s3_spatial_diagnostics(input_dir: Path, output: Path) -> None:
    frame = pd.concat([pd.read_csv(input_dir / "spatial" / "v2_full27_spatial_within.csv"),
        pd.read_csv(input_dir / "spatial" / "v2_full27_spatial_among.csv")], ignore_index=True)
    frame = frame.loc[frame["status"].eq("ok")].copy()
    frame["spatial_permutation_p_value"] = pd.to_numeric(frame["spatial_permutation_p_value"], errors="coerce")
    frame["residual_morans_p_value"] = pd.to_numeric(frame["residual_morans_p_value"], errors="coerce")
    frame = frame.dropna(subset=["spatial_permutation_p_value", "residual_morans_p_value"])
    frame["neglog10_spatial_p"] = -np.log10(frame["spatial_permutation_p_value"].clip(lower=1e-6))
    frame["pass"] = frame["broad_spatial_sensitivity_pass"].astype(str).str.lower().eq("true")
    frame["is_among"] = frame["scale"].eq("among_taxon")
    fig, ax = plt.subplots(figsize=(6.6, 4.9))
    groups = [(False, False, "within: fail", "o", "#9aa7b2"), (False, True, "within: pass", "o", "#507fa3"),
        (True, False, "among: fail", "s", "#c59b7a"), (True, True, "among: pass", "s", "#925c6a")]
    for is_among, passed, label, marker, color in groups:
        subset = frame.loc[(frame["is_among"] == is_among) & (frame["pass"] == passed)]
        ax.scatter(subset["neglog10_spatial_p"], subset["residual_morans_p_value"],
            s=48 if passed else 30, marker=marker, color=color, alpha=0.9, label=label)
    ax.axvline(-np.log10(0.05), color="#6b7280", linestyle="--", linewidth=0.9)
    ax.axhline(0.05, color="#6b7280", linestyle="--", linewidth=0.9)
    ax.set(xlabel="Spatial association support, -log10(permutation P)", ylabel="Residual Moran P")
    ax.set_xlim(-0.09, 3.16)
    ax.grid(color="#e5e7eb", linewidth=0.5)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, loc="upper left", fontsize=7.5)
    for _, row in frame.loc[frame["pass"]].iterrows():
        label = f"{ENDPOINT_LABELS.get(row['unit_id'], row['unit_id'])} – {PREDICTOR_LABELS.get(row['predictor'], row['predictor'])}"
        align_right = row["neglog10_spatial_p"] >= 2.5
        ax.annotate(
            label,
            (row["neglog10_spatial_p"], row["residual_morans_p_value"]),
            xytext=(-5 if align_right else 5, 5),
            textcoords="offset points",
            ha="right" if align_right else "left",
            fontsize=7.0,
        )
    fig.subplots_adjust(left=0.10, right=0.985, bottom=0.12, top=0.98)
    save(fig, output, SUPP_STEMS[2])


def figure_s4_historical_placement(input_dir: Path, output: Path) -> None:
    """Display the frozen 52-tree placement sensitivity without refitting models."""
    frame = pd.read_csv(input_dir / "historical" / "v2_full27_historical_placement_models.csv")
    assert len(frame) == 104
    assert frame.groupby(["unit_id", "predictor"]).size().eq(52).all()
    order = [
        ("corolla_lab_chroma", "chelsa_rsds_mean", "chroma ~ radiation", "#4f7fa4"),
        ("orientation_image_vertical_angle", "chelsa_bio12", "orientation ~ precipitation", "#a96f45"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(6.6, 3.6), sharey=True)
    panels = [
        ("pgls_beta_std", "PGLS coefficient", lambda values: values),
        ("p_value", "Placement P (-log10)", lambda values: -np.log10(values)),
        ("lambda", "Pagel lambda", lambda values: values),
    ]
    for panel_index, (column, xlabel, transform) in enumerate(panels):
        ax = axes[panel_index]
        for row_index, (unit, predictor, label, color) in enumerate(order):
            subset = frame.loc[(frame["unit_id"] == unit) & (frame["predictor"] == predictor)].copy()
            random_rows = subset.loc[subset["scenario"].eq("S2_random")].sort_values("replicate")
            deterministic = subset.loc[~subset["scenario"].eq("S2_random")].sort_values("scenario")
            random_y = row_index + np.linspace(-0.16, 0.16, len(random_rows))
            ax.scatter(
                transform(pd.to_numeric(random_rows[column], errors="raise")), random_y,
                s=21, facecolors="white", edgecolors=color, linewidths=0.8, alpha=0.78,
            )
            deterministic_y = row_index + np.array([-0.055, 0.055])
            ax.scatter(
                transform(pd.to_numeric(deterministic[column], errors="raise")), deterministic_y,
                s=52, marker="D", color=color, edgecolors="#30343b", linewidths=0.5, zorder=4,
            )
        ax.grid(axis="x", color="#d8dde1", linewidth=0.55)
        ax.set_axisbelow(True)
        ax.set_xlabel(xlabel)
        ax.set_title(f"({chr(97 + panel_index)})", loc="left", weight="bold", fontsize=9.5)
        for spine in ("top", "right", "left"):
            ax.spines[spine].set_visible(False)
    axes[0].axvline(0, color="#5f6368", linestyle="--", linewidth=0.8)
    axes[1].axvline(-np.log10(0.05), color="#5f6368", linestyle="--", linewidth=0.8)
    axes[0].set_yticks([0, 1], [order[0][2], order[1][2]], fontsize=7.2)
    axes[0].invert_yaxis()
    axes[1].tick_params(left=False, labelleft=False)
    axes[2].tick_params(left=False, labelleft=False)
    axes[2].set_xlim(left=-0.002)
    fig.legend(handles=[
        Line2D([0], [0], marker="D", linestyle="none", markerfacecolor="#6b7280", markeredgecolor="#30343b", label="two deterministic placements"),
        Line2D([0], [0], marker="o", linestyle="none", markerfacecolor="white", markeredgecolor="#6b7280", label="50 randomized placements"),
    ], frameon=False, loc="lower center", bbox_to_anchor=(0.5, 0.01), ncol=2, fontsize=7)
    fig.subplots_adjust(left=0.16, right=0.985, bottom=0.27, top=0.94, wspace=0.16)
    save(fig, output, SUPP_STEMS[3])


def figure_s5_whole_capitulum(input_dir: Path, output: Path) -> None:
    """Show the two frozen secondary syntheses without conflating their cohorts."""
    secondary = input_dir / "secondary"
    scores = pd.read_csv(secondary / "pca" / "geb_v2_taxon_morphospace_scores.csv")
    scores = scores.loc[scores["scope"].eq("all_inferential_endpoints")].copy()
    assert len(scores) == 127
    matrix_rows = pd.read_csv(
        secondary / "multilevel" / "capitulum_space_unit_strength_matrices.csv"
    )
    matrix_rows = matrix_rows.loc[matrix_rows["scope"].eq("complete18_min5")].copy()
    order = [
        "orientation_image_vertical_angle",
        "corolla_lab_lightness", "corolla_lab_chroma", "corolla_hue",
        "capitulum_outline_aspect_ratio", "capitulum_outline_circularity",
        "capitulum_outline_solidity", "capitulum_width_profile_cv",
        "involucre_length_width_ratio", "involucre_apical_taper_ratio",
        "involucre_basal_taper_ratio", "bract_projection_asymmetry",
        "bract_projection_p95", "bract_projection_roughness", "bract_spread_fraction",
        "bract_projection_maximum", "bract_projection_peak_density",
    ]
    short = [
        "orient.", "L*", "C*", "hue", "aspect", "circular.", "solidity", "width CV",
        "inv. L/W", "apical", "basal", "asym.", "proj. p95", "rough.", "spread", "proj. max", "peaks",
    ]

    def matrix(scale: str) -> np.ndarray:
        values = np.eye(len(order), dtype=float)
        index = {name: position for position, name in enumerate(order)}
        subset = matrix_rows.loc[matrix_rows["scale"].eq(scale)]
        assert len(subset) == len(order) * (len(order) - 1) // 2
        for row in subset.itertuples(index=False):
            left, right = index[row.left], index[row.right]
            values[left, right] = values[right, left] = float(row.value)
        return values

    fig = plt.figure(figsize=(6.7, 7.4))
    grid = fig.add_gridspec(2, 2, height_ratios=(0.78, 1.35), hspace=0.54, wspace=0.14)
    ax_pca = fig.add_subplot(grid[0, :])
    ax_pca.scatter(scores["PC1"], scores["PC2"], s=20, color="#3f7f93", alpha=0.70, linewidths=0)
    ax_pca.axhline(0, color="#c7cdd1", linewidth=0.6)
    ax_pca.axvline(0, color="#c7cdd1", linewidth=0.6)
    ax_pca.set_xlabel("PC1 (18.49%)")
    ax_pca.set_ylabel("PC2 (12.01%)")
    ax_pca.set_title("(a) Taxon-median morphospace: 127 taxa, 18 endpoints", loc="left", weight="bold", fontsize=9.5)
    ax_pca.text(
        0.99, 0.03, "PC1–PC3 cumulative variance: 42.28%", transform=ax_pca.transAxes,
        ha="right", va="bottom", fontsize=7.5, color="#30343b",
    )
    for spine in ("top", "right"):
        ax_pca.spines[spine].set_visible(False)

    heat_axes = [fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1])]
    images = []
    for panel_index, (ax, scale, title) in enumerate(zip(
        heat_axes,
        ("within_taxon", "among_taxon"),
        ("(b) Within taxa: 1,734 observations, 42 taxa", "(c) Among taxon medians: 42 taxa"),
    )):
        image = ax.imshow(matrix(scale), vmin=0, vmax=1, cmap="viridis", interpolation="nearest")
        images.append(image)
        ax.set_xticks(range(len(short)), short, rotation=68, ha="right", fontsize=5.4)
        ax.set_yticks(range(len(short)), short if panel_index == 0 else [], fontsize=5.4)
        ax.set_title(title, loc="left", weight="bold", fontsize=7.8)
        for boundary in (0.5, 3.5, 7.5, 14.5):
            ax.axhline(boundary, color="white", linewidth=0.8)
            ax.axvline(boundary, color="white", linewidth=0.8)
        ax.tick_params(length=0)
    colourbar_axis = fig.add_axes([0.34, 0.035, 0.32, 0.018])
    colourbar = fig.colorbar(images[-1], cax=colourbar_axis, orientation="horizontal")
    colourbar.set_label("association strength", fontsize=7.5)
    colourbar.ax.tick_params(labelsize=6.5)
    fig.subplots_adjust(left=0.12, right=0.98, bottom=0.16, top=0.98)
    save(fig, output, SUPP_STEMS[4])


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    figure_measurement(args.output_dir)
    figure_geographic_domain(args.output_dir)
    figure_information_loss(args.input_dir, args.output_dir)
    figure_scale_atlas(args.input_dir, args.output_dir)
    figure_candidates(args.input_dir, args.output_dir)
    figure_s1_endpoint_support(args.input_dir, args.output_dir)
    figure_s2_sampling_audit(args.input_dir, args.output_dir)
    figure_s3_spatial_diagnostics(args.input_dir, args.output_dir)
    figure_s4_historical_placement(args.input_dir, args.output_dir)
    figure_s5_whole_capitulum(args.input_dir, args.output_dir)
    for stem in OBSOLETE_STEMS:
        if stem not in MAIN_STEMS:
            for suffix in (".png", ".pdf"):
                (args.output_dir / f"{stem}{suffix}").unlink(missing_ok=True)
    inputs = [args.input_dir / "v2_full27_endpoint_inventory.csv",
        args.input_dir / "v2_full27_variance_decomposition.csv", args.input_dir / "v2_full27_environment_cross_scale.csv",
        args.input_dir / "spatial" / "v2_full27_spatial_within.csv", args.input_dir / "spatial" / "v2_full27_spatial_among.csv",
        args.input_dir / "historical" / "v2_full27_historical_placement_summary.csv",
        args.input_dir / "historical" / "v2_full27_historical_placement_models.csv",
        args.input_dir / "sampling" / "v2_full27_sampling_composition_summary.csv",
        args.input_dir / "secondary" / "source_provenance.json",
        args.input_dir / "secondary" / "pca" / "geb_v2_taxon_morphospace_scores.csv",
        args.input_dir / "secondary" / "pca" / "geb_v2_taxon_morphospace_loadings.csv",
        args.input_dir / "secondary" / "pca" / "geb_v2_trait_geography_report.json",
        args.input_dir / "secondary" / "multilevel" / "capitulum_space_scope_summary.csv",
        args.input_dir / "secondary" / "multilevel" / "capitulum_space_module_integration.csv",
        args.input_dir / "secondary" / "multilevel" / "capitulum_space_unit_strength_matrices.csv",
        args.input_dir / "secondary" / "multilevel" / "capitulum_space_eigenspectra.csv",
        args.input_dir / "secondary" / "multilevel" / "capitulum_space_report.json",
        SOURCE / "Figure_1_real_photo_measurement_provenance.csv", SOURCE / "Figure_1_real_photo_measurement_provenance.json",
        SOURCE / "Figure_1_real_photo_measurement_provenance.png", SOURCE / "Figure_2_v2_realized_sampling_cells.csv",
        SOURCE / "Figure_2_v2_realized_sampling_cells_provenance.json", SOURCE / "natural_earth_110m_land_outline.json"]
    produced = sorted([p.name for stem in MAIN_STEMS + SUPP_STEMS for p in [args.output_dir / f"{stem}.png"] if p.exists()])
    assert len(produced) == 10
    report = {"status": "ok", "source_lane": "full27_full_environment_only",
        "main_figure_stems": MAIN_STEMS, "supporting_figure_stems": SUPP_STEMS,
        "input_sha256": {path.relative_to(ROOT).as_posix(): sha256(path) for path in inputs}, "figures": produced,
        "presentation_source_provenance_firewall": "Figure 1 photographs and detector/crop examples document the source measurement pipeline only; they do not contribute v1 result counts or effect estimates.",
        "claim_boundary": "submission visualization of frozen exploratory v2 results; not a new analysis"}
    (args.output_dir / "figure_build_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
