#!/usr/bin/env python3
"""Build canonical Chapter 1 figures from external frozen v2 artifacts.

Prepublication inputs, generated figures and result tables are intentionally not
tracked in this repository.  The caller must supply explicit input, figure-source
and output directories.

This script is presentation-only: it reshapes frozen outputs but does not fit a
model, calculate a P value, alter multiplicity, or select a hypothesis.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Patch, Rectangle
import numpy as np
import pandas as pd

import azami_figstyle as fs

PREDICTOR_ORDER = [
    "chelsa_bio01", "chelsa_bio04", "chelsa_bio12", "chelsa_bio15",
    "chelsa_gsp", "chelsa_rsds_mean", "chelsa_vpd_mean",
    "chelsa_sfcwind_mean", "chelsa_npp",
]
PREDICTOR_LABELS = {
    "chelsa_bio01": "BIO1",
    "chelsa_bio04": "BIO4",
    "chelsa_bio12": "BIO12",
    "chelsa_bio15": "BIO15",
    "chelsa_gsp": "GSP",
    "chelsa_rsds_mean": "radiation",
    "chelsa_vpd_mean": "VPD",
    "chelsa_sfcwind_mean": "wind",
    "chelsa_npp": "NPP",
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input-dir", type=Path, required=True,
        help="External frozen v2 analysis-output directory.",
    )
    parser.add_argument(
        "--source-dir", type=Path, required=True,
        help="External presentation/source-provenance assets used by Figures 1-2.",
    )
    parser.add_argument(
        "--output-dir", type=Path, required=True,
        help="External destination for generated PNG/PDF figures.",
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(fig: plt.Figure, output: Path, stem: str) -> None:
    fs.savefig(fig, stem, width="double", outdir=output)
    plt.close(fig)


def endpoint_label(key: str, *, short: bool = False) -> str:
    return fs.labels([key], short=short)[0]


def clean_axis(ax: plt.Axes, *, grid_axis: str = "x") -> None:
    ax.grid(axis=grid_axis, linewidth=0.4, color="#D9D9D9")
    ax.set_axisbelow(True)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)


def panel_box(
    ax: plt.Axes,
    xy: tuple[float, float],
    width: float,
    height: float,
    title: str,
    subtitle: str,
    color: str,
) -> None:
    x, y = xy
    ax.add_patch(FancyBboxPatch(
        (x, y), width, height,
        boxstyle="round,pad=0.006,rounding_size=0.014",
        linewidth=0.6, edgecolor="#4b5563", facecolor=color,
    ))
    ax.text(
        x + width / 2, y + height * 0.62, title,
        ha="center", va="center", weight="bold", fontsize=fs.FONT["annot"],
    )
    ax.text(
        x + width / 2, y + height * 0.25, subtitle,
        ha="center", va="center", fontsize=fs.FONT["footnote"], color="#374151",
    )


def _trim_legacy_corner(image: np.ndarray, fraction: float = 0.045) -> np.ndarray:
    """Trim old baked panel-letter margins from a legacy provenance composite."""
    height, width = image.shape[:2]
    return image[int(height * fraction):, int(width * fraction):]


def _pad_square(image: np.ndarray) -> np.ndarray:
    """Pad an image to a square without distorting its aspect ratio."""
    height, width = image.shape[:2]
    side = max(height, width)
    if image.ndim == 3:
        channels = image.shape[2]
        fill = 1.0 if np.issubdtype(image.dtype, np.floating) else 255
        canvas = np.full((side, side, channels), fill, dtype=image.dtype)
    else:
        fill = 1.0 if np.issubdtype(image.dtype, np.floating) else 255
        canvas = np.full((side, side), fill, dtype=image.dtype)
    y0 = (side - height) // 2
    x0 = (side - width) // 2
    canvas[y0:y0 + height, x0:x0 + width] = image
    return canvas


def _species_title(name: str) -> str:
    parts = str(name).split()
    if len(parts) >= 2:
        return rf"$\it{{{parts[0]}\ {parts[1]}}}$"
    return str(name)


# ---------------------------------------------------------------------------
# Main figures
# ---------------------------------------------------------------------------


def figure_measurement(source_dir: Path, output: Path) -> None:
    """Figure 1: photographs and measurement operations, without prose panels."""
    image_path = source_dir / "Figure_1_real_photo_measurement_provenance.png"
    provenance_path = source_dir / "Figure_1_real_photo_measurement_provenance.csv"
    assert image_path.is_file() and provenance_path.is_file()
    source_image = plt.imread(image_path)
    provenance = pd.read_csv(provenance_path)
    if len(provenance) < 3:
        raise ValueError("Figure 1 provenance table must contain at least three examples")

    # Coordinates refer to the preserved external provenance composite.  The
    # old in-bitmap panel-letter margins are trimmed, then examples are padded
    # rather than stretched so all source photographs share one display size.
    source_regions = [
        source_image[390:1100, 90:650],
        source_image[390:1100, 760:1450],
        source_image[505:1015, 1510:2200],
    ]
    operation_regions = [
        source_image[1300:2010, 100:610],
        source_image[1300:2010, 860:1370],
        source_image[1300:2010, 1580:2100],
        source_image[1300:2010, 2340:2860],
    ]
    orientation_region = source_image[2190:2860, 480:1000]

    fig = fs.figure(width="double", height=5.80)
    grid = fig.add_gridspec(3, 1, height_ratios=[1.55, 1.48, 1.25], hspace=0.34)

    top = grid[0].subgridspec(1, 3, wspace=0.08)
    for index, (region, row) in enumerate(zip(source_regions, provenance.iloc[:3].itertuples(index=False))):
        ax = fig.add_subplot(top[index])
        ax.imshow(_pad_square(_trim_legacy_corner(region)))
        angle = float(getattr(row, "orientation_angle_degrees"))
        chroma = float(getattr(row, "corolla_lab_chroma"))
        ax.set_title(
            f"({chr(97 + index)}) {_species_title(getattr(row, 'taxon_name'))}\n"
            f"angle {angle:.1f}°; chroma {chroma:.1f}",
            loc="left", fontsize=fs.FONT["annot"], fontweight="bold", pad=3,
        )
        ax.axis("off")

    middle = grid[1].subgridspec(1, 4, wspace=0.08)
    operation_titles = [
        ("d", "Tight crop"),
        ("e", "Context crop"),
        ("f", "Corolla mask"),
        ("g", "Outline + hull"),
    ]
    for index, (region, (letter, title)) in enumerate(zip(operation_regions, operation_titles)):
        ax = fig.add_subplot(middle[index])
        ax.imshow(_pad_square(_trim_legacy_corner(region)))
        fs.panel(ax, letter, title)
        ax.axis("off")

    bottom = grid[2].subgridspec(1, 3, width_ratios=[1.0, 1.0, 1.0], wspace=0.16)
    orientation = fig.add_subplot(bottom[0:2])
    orientation.imshow(_trim_legacy_corner(orientation_region))
    fs.panel(orientation, "h", "Orientation relative to image vertical")
    orientation.axis("off")

    # A compact non-tabular glyph summarizes what the image operations return.
    # It deliberately avoids recreating the former prose/table panels.
    outputs = fig.add_subplot(bottom[2])
    outputs.set_xlim(0, 1)
    outputs.set_ylim(0, 1)
    outputs.axis("off")
    fs.panel(outputs, "i", "Continuous outputs")
    example = provenance.iloc[0]
    glyphs = [
        (0.80, fs.MODULE_COLOUR["orientation"], "angle", f"{float(example['orientation_angle_degrees']):.1f}°"),
        (0.53, fs.MODULE_COLOUR["colour"], "chroma", f"{float(example['corolla_lab_chroma']):.1f}"),
        (0.26, fs.MODULE_COLOUR["outline"], "aspect", f"{float(example['shape_aspect_ratio']):.2f}"),
    ]
    for y, colour, label, value in glyphs:
        outputs.scatter([0.12], [y], s=72, color=colour, edgecolor="#333333", linewidth=0.35)
        outputs.text(0.24, y + 0.055, label, va="center", fontsize=fs.FONT["footnote"], color="#4b5563")
        outputs.text(0.24, y - 0.045, value, va="center", fontsize=fs.FONT["axis"], fontweight="bold")

    fig.subplots_adjust(left=0.035, right=0.985, top=0.98, bottom=0.03)
    save(fig, output, MAIN_STEMS[0])


def draw_world(ax: plt.Axes, source_dir: Path) -> None:
    payload = json.loads((source_dir / "natural_earth_110m_land_outline.json").read_text(encoding="utf-8"))
    for coords in payload["parts"]:
        values = np.asarray(coords, dtype=float)
        ax.fill(
            np.deg2rad(values[:, 0]), np.deg2rad(values[:, 1]),
            facecolor="#F3F3EF", edgecolor="#B8BBB8", linewidth=0.3, zorder=0,
        )
    ax.set_longitude_grid(60)
    ax.set_latitude_grid(30)
    ax.grid(color="#D9D9D9", linewidth=0.35, zorder=-1)
    x0, x1, y = 0.08, 0.13, 0.11
    ax.plot([x0, x1], [y, y], transform=ax.transAxes, color="#30343b", linewidth=1.8, zorder=4)
    ax.plot([x0, x0], [y - 0.010, y + 0.010], transform=ax.transAxes, color="#30343b", linewidth=0.8, zorder=4)
    ax.plot([x1, x1], [y - 0.010, y + 0.010], transform=ax.transAxes, color="#30343b", linewidth=0.8, zorder=4)
    ax.text(
        (x0 + x1) / 2, y + 0.024, "2,000 km at Equator",
        transform=ax.transAxes, ha="center", va="bottom", fontsize=fs.FONT["footnote"],
    )


def figure_geographic_domain(source_dir: Path, output: Path) -> None:
    cells = pd.read_csv(source_dir / "Figure_2_v2_realized_sampling_cells.csv")
    if int(cells["n_observations"].sum()) != 46_276:
        raise ValueError("Figure 2 source cells do not match the frozen cohort count")

    fig = fs.figure(width="double", height=7.20)
    grid = fig.add_gridspec(
        3, 2, height_ratios=[3.05, 1.05, 3.0], width_ratios=[24, 0.72],
        hspace=0.40, wspace=0.08, left=0.10, right=0.93, top=0.98, bottom=0.08,
    )

    ax_map = fig.add_subplot(grid[0, 0], projection="mollweide")
    draw_world(ax_map, source_dir)
    points = ax_map.scatter(
        np.deg2rad(cells["cell_lon_center"]), np.deg2rad(cells["cell_lat_center"]),
        c=cells["n_observations"], cmap=fs.CMAP_SEQUENTIAL,
        norm=LogNorm(vmin=1, vmax=max(10, cells["n_observations"].max())),
        s=11, marker="s", linewidth=0, alpha=0.90, zorder=2,
    )
    colourbar = fig.colorbar(points, cax=fig.add_subplot(grid[0, 1]))
    colourbar.set_label("Observations per 2° cell (display only)")
    fs.panel(ax_map, "a", "Realized geographic sampling domain")

    ax_flow = fig.add_subplot(grid[1, 0])
    ax_flow.set(xlim=(0, 1), ylim=(0, 1))
    ax_flow.axis("off")
    steps = [
        ("Detector-positive", "406,582 / 286 taxa"),
        ("Coordinates", "392,989 / 271"),
        ("Accuracy ≤10 km", "297,293 / 259"),
        ("Taxon × 0.25°", "46,276 / 259"),
    ]
    colours = [fs.C["skyblue"], "#DCEFE6", "#F2ECD9", "#E9E4F3"]
    width = 0.19
    gap = (0.98 - len(steps) * width) / (len(steps) - 1)
    xs = [0.01 + i * (width + gap) for i in range(len(steps))]
    for x, (title, subtitle), colour in zip(xs, steps, colours):
        panel_box(ax_flow, (x, 0.30), width, 0.42, title, subtitle, colour)
    for x, next_x in zip(xs[:-1], xs[1:]):
        ax_flow.annotate(
            "", xy=(next_x - 0.010, 0.51), xytext=(x + width + 0.010, 0.51),
            arrowprops={"arrowstyle": "-|>", "lw": 0.9, "color": "#4b5563", "mutation_scale": 8},
        )
    ax_flow.text(0.01, 0.92, "(b) Prespecified cohort flow", fontsize=fs.FONT["panel"], weight="bold")

    ax_env = fig.add_subplot(grid[2, 0])
    env = ax_env.scatter(
        cells["bio1_degrees_c_mean"], cells["bio12_mm_mean"],
        s=np.clip(np.sqrt(cells["n_observations"].clip(lower=1)) * 1.8, 4, 55),
        c=cells["cell_lat_center"], cmap=fs.CMAP_DIVERGING, alpha=0.50, linewidth=0,
    )
    cb2 = fig.colorbar(env, cax=fig.add_subplot(grid[2, 1]))
    cb2.set_label("Cell latitude")
    ax_env.set_xlabel("Annual mean temperature, BIO1 (°C)")
    ax_env.set_ylabel("Annual precipitation, BIO12 (mm)")
    ax_env.grid(color="#D9D9D9", linewidth=0.4)
    fs.panel(ax_env, "c", "Realized climatic domain")

    for aligned_ax in (ax_flow, ax_env):
        position = aligned_ax.get_position()
        inset = position.width * 0.038
        aligned_ax.set_position([position.x0 + inset, position.y0, position.width - 2 * inset, position.height])

    save(fig, output, MAIN_STEMS[1])


def figure_information_loss(input_dir: Path, output: Path) -> None:
    frame = pd.read_csv(input_dir / "v2_full27_variance_decomposition.csv")
    frame = frame.loc[frame["status"].eq("ok")].copy().set_index("endpoint_id")
    keys = fs.order(split_hue=True, measured_only=True)
    missing = [key for key in keys if key not in frame.index]
    if missing:
        raise KeyError(f"variance-decomposition endpoints missing from canonical registry: {missing}")
    frame = frame.loc[keys].reset_index()

    y = np.arange(len(frame))
    raw = frame["fraction_below_taxon_means"].to_numpy() * 100
    balanced = frame["equal_replication_fraction_median"].to_numpy() * 100

    fig, axes = plt.subplots(1, 2, figsize=(fs.WIDTH["double"], 6.55), sharey=True)
    for ax, within, letter, title in (
        (axes[0], raw, "a", "Raw observed sample"),
        (axes[1], balanced, "b", "Two observations per taxon"),
    ):
        among = 100 - within
        ax.barh(y, within, height=0.66, color=fs.C["blue"], edgecolor="white", linewidth=0.25)
        ax.barh(y, among, left=within, height=0.66, color=fs.C["lightgrey"], edgecolor="white", linewidth=0.25)
        ax.set_xlim(0, 100)
        ax.set_xticks([0, 25, 50, 75, 100])
        ax.grid(axis="x", color="#D9D9D9", linewidth=0.4)
        ax.invert_yaxis()
        fs.panel(ax, letter, title)
    axes[0].set_yticks(y, fs.labels(keys), fontsize=fs.FONT["tick"])
    axes[1].tick_params(axis="y", left=False, labelleft=False)
    fig.supxlabel("Share of total observed image-phenotype variation (%)", y=0.07, fontsize=fs.FONT["axis"])
    fig.legend(
        handles=[
            Patch(facecolor=fs.C["blue"], label="within source-assigned taxa"),
            Patch(facecolor=fs.C["lightgrey"], label="among source-assigned taxa"),
        ],
        loc="lower center", bbox_to_anchor=(0.61, 0.015), frameon=False, ncol=2,
    )
    fig.subplots_adjust(left=0.33, right=0.985, top=0.97, bottom=0.13, wspace=0.16)
    save(fig, output, MAIN_STEMS[2])


def figure_scale_atlas(input_dir: Path, output: Path) -> None:
    frame = pd.read_csv(input_dir / "v2_full27_environment_cross_scale.csv")
    counts = frame["cross_scale_class"].value_counts().to_dict()
    expected = {"both_scales": 7, "within_only": 18, "among_only": 3, "neither": 152, "not_comparable": 54}
    if len(frame) != 234 or counts != expected:
        raise ValueError(f"cross-scale frozen ledger drifted: rows={len(frame)}, counts={counts}")

    units = fs.order(split_hue=False, measured_only=False)
    lookup = frame.set_index(["unit_id", "predictor"])["cross_scale_class"]
    missing = [(unit, predictor) for unit in units for predictor in PREDICTOR_ORDER if (unit, predictor) not in lookup.index]
    if missing:
        raise KeyError(f"cross-scale rows missing: {missing[:10]}")

    fig = fs.figure(width="double", height=6.65)
    ax = fig.add_axes([0.30, 0.16, 0.68, 0.81])
    for row_index, unit in enumerate(units):
        for column_index, predictor in enumerate(PREDICTOR_ORDER):
            classification = lookup.loc[(unit, predictor)]
            role = fs.SCALE_CLASS[classification]
            ax.add_patch(Rectangle(
                (column_index - 0.5, row_index - 0.5), 1, 1,
                facecolor=role["color"], hatch=role["hatch"],
                edgecolor="#8A8A8A" if classification == "not_comparable" else "white",
                linewidth=0.35,
            ))
    ax.set_xlim(-0.5, len(PREDICTOR_ORDER) - 0.5)
    ax.set_ylim(len(units) - 0.5, -0.5)
    ax.set_xticks(
        range(len(PREDICTOR_ORDER)),
        [PREDICTOR_LABELS[predictor] for predictor in PREDICTOR_ORDER],
        rotation=38, ha="right",
    )
    ax.set_yticks(range(len(units)), fs.labels(units))
    ax.set_xlabel("Environmental gradient")
    ax.set_ylabel("Continuous endpoint or joint circular unit")
    ax.grid(False)
    ax.legend(
        handles=[
            Patch(
                facecolor=fs.SCALE_CLASS[key]["color"],
                hatch=fs.SCALE_CLASS[key]["hatch"],
                edgecolor="#808080" if key == "not_comparable" else "none",
                label=f"{fs.SCALE_CLASS[key]['label']} ({counts[key]})",
            )
            for key in ("both_scales", "within_only", "among_only", "neither", "not_comparable")
        ],
        frameon=False, ncol=2, bbox_to_anchor=(0, -0.13), loc="upper left",
    )
    save(fig, output, MAIN_STEMS[3])


def _candidate_rows(input_dir: Path):
    spatial = pd.read_csv(input_dir / "spatial" / "v2_full27_spatial_among.csv")
    sampling = pd.read_csv(input_dir / "sampling" / "v2_full27_sampling_composition_summary.csv")
    historical = pd.read_csv(input_dir / "historical" / "v2_full27_historical_placement_summary.csv")
    candidates = [
        ("corolla_lab_chroma", "chelsa_rsds_mean", "lower chroma ~ higher radiation"),
        ("orientation_image_vertical_angle", "chelsa_bio12", "larger orientation angle ~ higher precipitation"),
    ]
    rows = []
    for unit, predictor, label in candidates:
        spatial_row = spatial.loc[
            spatial["unit_id"].eq(unit)
            & spatial["predictor"].eq(predictor)
            & spatial["broad_spatial_sensitivity_pass"].astype(str).str.lower().eq("true")
        ].iloc[0]
        sampling_row = sampling.loc[
            sampling["scale"].eq("among_taxon")
            & sampling["unit_id"].eq(unit)
            & sampling["predictor"].eq(predictor)
        ].iloc[0]
        historical_row = historical.loc[
            historical["unit_id"].eq(unit) & historical["predictor"].eq(predictor)
        ].iloc[0]
        ratio_columns = [column for column in sampling_row.index if column.endswith("minimum_effect_magnitude_ratio")]
        ratios = pd.to_numeric(sampling_row[ratio_columns], errors="coerce").dropna()
        rows.append((label, spatial_row, float(ratios.min()), historical_row))
    if len(rows) != 2:
        raise ValueError("candidate ledger no longer contains exactly two retained among-taxon rows")
    return rows


def figure_candidates(input_dir: Path, output: Path) -> None:
    rows = _candidate_rows(input_dir)
    fig = fs.figure(width="double", height=5.25)
    grid = fig.add_gridspec(2, 1, height_ratios=[1.15, 2.6], hspace=0.36)

    flow = fig.add_subplot(grid[0])
    flow.set(xlim=(0, 1), ylim=(0, 1))
    flow.axis("off")
    steps = [
        ("Canonical grid", "234 rows"),
        ("Among BH", "10 supported"),
        ("Sampling audit", "directions annotated"),
        ("Broad-space", "2 pass"),
        ("52 trees", "2/2 retained"),
    ]
    width = 0.17
    gap = (0.98 - len(steps) * width) / (len(steps) - 1)
    xs = [0.01 + index * (width + gap) for index in range(len(steps))]
    colours = ["#E8F1F8", "#E8F1F8", "#DCEFE6", "#F2ECD9", "#E9E4F3"]
    for x, (title, subtitle), colour in zip(xs, steps, colours):
        panel_box(flow, (x, 0.29), width, 0.43, title, subtitle, colour)
    for x, next_x in zip(xs[:-1], xs[1:]):
        flow.annotate(
            "", xy=(next_x - 0.008, 0.505), xytext=(x + width + 0.008, 0.505),
            arrowprops={"arrowstyle": "-|>", "lw": 0.9, "color": "#4b5563", "mutation_scale": 8},
        )
    flow.text(0.0, 0.94, "(a) Sequential attrition", fontsize=fs.FONT["panel"], weight="bold")

    lower = grid[1].subgridspec(1, 2, width_ratios=[1.45, 1.0], wspace=0.20)
    coefficient = fig.add_subplot(lower[0])
    y = np.arange(2)
    global_beta = [float(row[1]["base_beta_std"]) for row in rows]
    spatial_beta = [float(row[1]["spatial_beta_std"]) for row in rows]
    coefficient.axvline(0, color="#808080", linewidth=0.7)
    coefficient.scatter(global_beta, y - 0.11, marker="o", s=36, color=fs.C["blue"], label="global among-taxon")
    coefficient.scatter(spatial_beta, y + 0.11, marker="s", s=36, color=fs.C["orange"], label="broad-space")
    coefficient.set_yticks(y, [row[0] for row in rows])
    coefficient.invert_yaxis()
    coefficient.set_xlabel("Standardized coefficient")
    coefficient.grid(axis="x", color="#D9D9D9", linewidth=0.4)
    fs.panel(coefficient, "b", "Retained coefficient direction")
    coefficient.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, frameon=False)

    gates = fig.add_subplot(lower[1])
    gates.set_xlim(-0.5, 2.5)
    gates.set_ylim(1.55, -0.55)
    gates.set_xticks([0, 1, 2], ["sampling", "spatial", "placements"])
    gates.set_yticks(y, ["chroma", "orientation"])
    for index, (_, spatial_row, ratio, historical_row) in enumerate(rows):
        gates.scatter(0, index, s=45, marker="o", color=fs.C["blue"])
        gates.text(0, index + 0.23, f"ratio {ratio:.2f}", ha="center", va="top", fontsize=fs.FONT["footnote"])
        gates.scatter(1, index, s=45, marker="s", color=fs.C["orange"])
        gates.text(
            1, index + 0.23,
            f"P={float(spatial_row['spatial_permutation_p_value']):.3f}\n"
            f"Moran P={float(spatial_row['residual_morans_p_value']):.3f}",
            ha="center", va="top", fontsize=fs.FONT["footnote"],
        )
        gates.scatter(2, index, s=55, marker="D", color=fs.C["green"])
        gates.text(
            2, index + 0.23,
            f"{int(historical_row['n_placement_trees_p_lt_0_05'])}/"
            f"{int(historical_row['n_successful_placement_trees'])}",
            ha="center", va="top", fontsize=fs.FONT["footnote"],
        )
    gates.grid(axis="x", color="#EEEEEE", linewidth=0.4)
    gates.tick_params(axis="x", rotation=25)
    fs.panel(gates, "c", "Gate evidence")
    save(fig, output, MAIN_STEMS[4])


# ---------------------------------------------------------------------------
# Supporting figures S3-S7
# ---------------------------------------------------------------------------


def figure_s3_endpoint_support(input_dir: Path, output: Path) -> None:
    frame = pd.read_csv(input_dir / "v2_full27_endpoint_inventory.csv")
    frame = frame.loc[frame["measurement_status"].eq("measured")].copy().set_index("endpoint_id")
    keys = fs.order(split_hue=True, measured_only=True)
    missing = [key for key in keys if key not in frame.index]
    if missing:
        raise KeyError(f"endpoint inventory is missing measured endpoints: {missing}")
    frame = frame.loc[keys]
    colours = [fs.MODULE_COLOUR[fs.module(key)] for key in keys]
    y = np.arange(len(keys))

    fig, axes = plt.subplots(
        1, 2, figsize=(fs.WIDTH["double"], 6.55),
        gridspec_kw={"width_ratios": [1.3, 1.0]}, sharey=True,
    )
    axes[0].barh(y, frame["n_observations_measured"], height=0.62, color=colours)
    axes[0].set_xscale("log")
    axes[0].set_xlabel("Observations (log scale)")
    axes[0].set_yticks(y, fs.labels(keys))
    axes[0].invert_yaxis()
    axes[0].grid(axis="x", color="#D9D9D9", linewidth=0.4)
    fs.panel(axes[0], "a", "Observation support")

    axes[1].barh(y, frame["n_taxa_measured"], height=0.62, color=colours)
    axes[1].set_xlabel("Taxa measured")
    axes[1].tick_params(axis="y", left=False, labelleft=False)
    axes[1].grid(axis="x", color="#D9D9D9", linewidth=0.4)
    fs.panel(axes[1], "b", "Taxonomic support")
    fig.subplots_adjust(left=0.36, right=0.98, bottom=0.08, top=0.97, wspace=0.08)
    save(fig, output, SUPP_STEMS[0])


def _sampling_panel_order(subset: pd.DataFrame) -> pd.DataFrame:
    unit_rank = {unit: index for index, unit in enumerate(fs.order(split_hue=False, measured_only=False))}
    predictor_rank = {predictor: index for index, predictor in enumerate(PREDICTOR_ORDER)}
    subset = subset.copy()
    subset["unit_rank"] = subset["unit_id"].map(unit_rank)
    subset["predictor_rank"] = subset["predictor"].map(predictor_rank)
    if subset[["unit_rank", "predictor_rank"]].isna().any().any():
        missing_units = sorted(subset.loc[subset["unit_rank"].isna(), "unit_id"].unique())
        missing_predictors = sorted(subset.loc[subset["predictor_rank"].isna(), "predictor"].unique())
        raise KeyError(f"sampling audit contains unknown units={missing_units}, predictors={missing_predictors}")
    return subset.sort_values(["unit_rank", "predictor_rank"]).reset_index(drop=True)


def figure_s4_sampling_audit(input_dir: Path, output: Path) -> None:
    frame = pd.read_csv(input_dir / "sampling" / "v2_full27_sampling_composition_summary.csv")
    ratio_columns = [
        "dominant_taxon_omission_minimum_effect_magnitude_ratio",
        "leave_one_broad_region_out_minimum_effect_magnitude_ratio",
        "native_only_minimum_effect_magnitude_ratio",
        "equal_taxon_weight_minimum_effect_magnitude_ratio",
    ]
    for column in ratio_columns:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame["minimum_ratio"] = frame[ratio_columns].min(axis=1, skipna=True)
    frame["stable"] = frame["all_directions_stable_where_evaluable"].astype(str).str.lower().eq("true")

    fig, axes = plt.subplots(
        1, 2, figsize=(fs.WIDTH["double"], 7.35),
        gridspec_kw={"width_ratios": [1.18, 1.0]},
    )
    for ax, scale, letter, title in (
        (axes[0], "within_taxon", "a", "Within taxa"),
        (axes[1], "among_taxon", "b", "Among taxa"),
    ):
        subset = _sampling_panel_order(frame.loc[frame["scale"].eq(scale)])
        y = np.arange(len(subset))
        for row_index, row in subset.iterrows():
            role = fs.STABILITY["stable" if row["stable"] else "unstable"]
            ax.barh(
                row_index, row["minimum_ratio"], height=0.62,
                color=role["color"], hatch=role["hatch"],
                edgecolor="#666666" if role["hatch"] else "none", linewidth=0.35,
            )
        labels = [
            f"{endpoint_label(row.unit_id, short=True)} ~ {PREDICTOR_LABELS.get(row.predictor, row.predictor)}"
            for row in subset.itertuples(index=False)
        ]
        ax.set_yticks(y, labels, fontsize=6.2)
        ax.invert_yaxis()
        ax.axvline(1.0, color="#808080", linewidth=0.7, linestyle="--")
        ax.set_xlim(0, 1.055)
        ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.set_xlabel("Minimum effect ratio")
        ax.grid(axis="x", color="#D9D9D9", linewidth=0.4)
        fs.panel(ax, letter, title)
    axes[1].legend(
        handles=[
            Patch(facecolor=fs.STABILITY["stable"]["color"], label="direction stable"),
            Patch(
                facecolor=fs.STABILITY["unstable"]["color"],
                hatch=fs.STABILITY["unstable"]["hatch"], edgecolor="#666666",
                label="direction unstable",
            ),
        ],
        loc="lower right", frameon=False,
    )
    fig.subplots_adjust(left=0.30, right=0.985, bottom=0.07, top=0.97, wspace=0.58)
    save(fig, output, SUPP_STEMS[1])


def figure_s5_spatial_diagnostics(input_dir: Path, output: Path) -> None:
    frame = pd.concat([
        pd.read_csv(input_dir / "spatial" / "v2_full27_spatial_within.csv"),
        pd.read_csv(input_dir / "spatial" / "v2_full27_spatial_among.csv"),
    ], ignore_index=True)
    frame = frame.loc[frame["status"].eq("ok")].copy()
    frame["spatial_permutation_p_value"] = pd.to_numeric(frame["spatial_permutation_p_value"], errors="coerce")
    frame["residual_morans_p_value"] = pd.to_numeric(frame["residual_morans_p_value"], errors="coerce")
    frame = frame.dropna(subset=["spatial_permutation_p_value", "residual_morans_p_value"])
    frame["neglog10_spatial_p"] = -np.log10(frame["spatial_permutation_p_value"].clip(lower=1e-6))
    frame["pass"] = frame["broad_spatial_sensitivity_pass"].astype(str).str.lower().eq("true")
    frame["is_among"] = frame["scale"].eq("among_taxon")

    fig = fs.figure(width="double", height=4.60)
    ax = fig.add_axes([0.11, 0.14, 0.87, 0.82])
    groups = [
        (False, False, "within: fail", "o", fs.C["grey"]),
        (False, True, "within: pass", "o", fs.C["blue"]),
        (True, False, "among: fail", "s", fs.C["orange"]),
        (True, True, "among: pass", "s", fs.C["purple"]),
    ]
    for is_among, passed, label, marker, colour in groups:
        subset = frame.loc[(frame["is_among"] == is_among) & (frame["pass"] == passed)]
        ax.scatter(
            subset["neglog10_spatial_p"], subset["residual_morans_p_value"],
            s=42 if passed else 26, marker=marker, color=colour, alpha=0.9, label=label,
        )

    x_threshold = -np.log10(0.05)
    ax.axvline(x_threshold, color="#808080", linestyle="--", linewidth=0.8)
    ax.axhline(0.05, color="#808080", linestyle="--", linewidth=0.8)
    ax.text(
        x_threshold + 0.03, 0.97, "permutation P = 0.05",
        transform=ax.get_xaxis_transform(), rotation=90, va="top", ha="left",
        fontsize=fs.FONT["footnote"],
    )
    ax.text(
        0.99, 0.052, "residual Moran P = 0.05",
        transform=ax.get_yaxis_transform(), va="bottom", ha="right",
        fontsize=fs.FONT["footnote"],
    )
    ax.set_xlabel("Spatial association support, -log10(permutation P)")
    ax.set_ylabel("Residual Moran P")
    ax.grid(color="#E5E7EB", linewidth=0.4)
    ax.legend(frameon=False, loc="upper left")

    for row in frame.loc[frame["pass"]].itertuples(index=False):
        label = f"{endpoint_label(row.unit_id)} – {PREDICTOR_LABELS.get(row.predictor, row.predictor)}"
        align_right = row.neglog10_spatial_p >= 2.5
        ax.annotate(
            label, (row.neglog10_spatial_p, row.residual_morans_p_value),
            xytext=(-6 if align_right else 6, 6), textcoords="offset points",
            ha="right" if align_right else "left", fontsize=fs.FONT["footnote"],
            arrowprops={"arrowstyle": "-", "lw": 0.45, "color": "#808080"},
        )
    save(fig, output, SUPP_STEMS[2])


def figure_s6_historical_placement(input_dir: Path, output: Path) -> None:
    frame = pd.read_csv(input_dir / "historical" / "v2_full27_historical_placement_models.csv")
    if len(frame) != 104 or not frame.groupby(["unit_id", "predictor"]).size().eq(52).all():
        raise ValueError("historical placement model table is not the frozen 2 × 52 tree ledger")

    candidates = [
        ("corolla_lab_chroma", "chelsa_rsds_mean", "chroma ~ radiation", fs.C["blue"]),
        ("orientation_image_vertical_angle", "chelsa_bio12", "orientation ~ precipitation", fs.C["orange"]),
    ]
    metrics = [
        ("pgls_beta_std", "PGLS coefficient", lambda values: values),
        ("p_value", "Placement P (-log10)", lambda values: -np.log10(np.clip(values, 1e-12, None))),
        ("lambda", "Pagel λ", lambda values: values),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(fs.WIDTH["double"], 3.35), sharey=True)
    for panel_index, (column, xlabel, transform) in enumerate(metrics):
        ax = axes[panel_index]
        for row_index, (unit, predictor, _label, colour) in enumerate(candidates):
            subset = frame.loc[frame["unit_id"].eq(unit) & frame["predictor"].eq(predictor)]
            random_values = transform(
                pd.to_numeric(subset.loc[subset["scenario"].eq("S2_random"), column], errors="raise").to_numpy()
            )
            deterministic = transform(
                pd.to_numeric(subset.loc[~subset["scenario"].eq("S2_random"), column], errors="raise").to_numpy()
            )
            parts = ax.violinplot(
                [random_values], positions=[row_index], vert=False, widths=0.55,
                showmeans=False, showextrema=False, showmedians=True,
            )
            for body in parts["bodies"]:
                body.set_facecolor(colour)
                body.set_edgecolor(colour)
                body.set_alpha(0.28)
            parts["cmedians"].set_color(colour)
            parts["cmedians"].set_linewidth(1.0)
            ax.scatter(
                deterministic, [row_index - 0.10, row_index + 0.10],
                s=30, marker="D", color=colour, edgecolor="#30343b", linewidth=0.35, zorder=4,
            )
        ax.grid(axis="x", color="#D9D9D9", linewidth=0.4)
        ax.set_xlabel(xlabel)
        fs.panel(ax, chr(97 + panel_index))
    axes[0].axvline(0, color="#808080", linestyle="--", linewidth=0.7)
    axes[1].axvline(-np.log10(0.05), color="#808080", linestyle="--", linewidth=0.7)
    axes[0].set_yticks([0, 1], [candidates[0][2], candidates[1][2]])
    axes[0].invert_yaxis()
    axes[1].tick_params(left=False, labelleft=False)
    axes[2].tick_params(left=False, labelleft=False)
    fig.legend(
        handles=[
            Patch(facecolor=fs.C["grey"], alpha=0.35, label="50 randomized placements"),
            Line2D(
                [0], [0], marker="D", linestyle="none", markerfacecolor="#808080",
                markeredgecolor="#30343b", label="two deterministic placements",
            ),
        ],
        frameon=False, loc="lower center", bbox_to_anchor=(0.5, 0.01), ncol=2,
    )
    fig.subplots_adjust(left=0.18, right=0.985, bottom=0.25, top=0.94, wspace=0.17)
    save(fig, output, SUPP_STEMS[3])


def _strength_matrix(matrix_rows: pd.DataFrame, scale: str) -> np.ndarray:
    keys = list(fs.S7_MATRIX_KEYS)
    values = np.eye(len(keys), dtype=float)
    index = {name: position for position, name in enumerate(keys)}
    subset = matrix_rows.loc[matrix_rows["scale"].eq(scale)]
    expected = len(keys) * (len(keys) - 1) // 2
    if len(subset) != expected:
        raise ValueError(f"{scale} strength matrix has {len(subset)} rows; expected {expected}")
    for row in subset.itertuples(index=False):
        if row.left not in index or row.right not in index:
            raise KeyError(f"unknown S7 inferential unit pair: {row.left}, {row.right}")
        left, right = index[row.left], index[row.right]
        values[left, right] = values[right, left] = float(row.value)
    return values


def figure_s7_whole_capitulum(input_dir: Path, output: Path) -> None:
    secondary = input_dir / "secondary"
    scores = pd.read_csv(secondary / "pca" / "geb_v2_taxon_morphospace_scores.csv")
    scores = scores.loc[scores["scope"].eq("all_inferential_endpoints")].copy()
    loadings = pd.read_csv(secondary / "pca" / "geb_v2_taxon_morphospace_loadings.csv")
    loadings = loadings.loc[loadings["scope"].eq("all_inferential_endpoints")].copy().set_index("endpoint_id")
    matrix_rows = pd.read_csv(secondary / "multilevel" / "capitulum_space_unit_strength_matrices.csv")
    matrix_rows = matrix_rows.loc[matrix_rows["scope"].eq("complete18_min5")].copy()

    if len(scores) != 127:
        raise ValueError(f"PCA score table has {len(scores)} rows; expected the frozen 127-taxon scope")
    missing_loadings = [key for key in fs.PCA_ENDPOINT_KEYS if key not in loadings.index]
    if missing_loadings:
        raise KeyError(f"PCA loadings missing canonical endpoints: {missing_loadings}")

    fig = fs.figure(width="double", height=7.10)
    grid = fig.add_gridspec(2, 2, height_ratios=[0.82, 1.28], hspace=0.46, wspace=0.14)
    ax_pca = fig.add_subplot(grid[0, :])
    ax_pca.scatter(scores["PC1"], scores["PC2"], s=15, color=fs.C["blue"], alpha=0.55, linewidths=0)
    ax_pca.axhline(0, color="#D0D0D0", linewidth=0.5)
    ax_pca.axvline(0, color="#D0D0D0", linewidth=0.5)
    ax_pca.set_xlabel("PC1 (18.49%)")
    ax_pca.set_ylabel("PC2 (12.01%)")
    fs.panel(ax_pca, "a", "Taxon-median morphospace: 127 taxa, 18 endpoints")

    loadings = loadings.loc[list(fs.PCA_ENDPOINT_KEYS)].copy()
    magnitudes = np.hypot(loadings["PC1_loading"], loadings["PC2_loading"])
    label_keys = set(magnitudes.nlargest(8).index)
    score_radius = 0.28 * min(float(np.ptp(scores["PC1"])), float(np.ptp(scores["PC2"])))
    loading_radius = float(magnitudes.max())
    scale = score_radius / loading_radius if loading_radius > 0 else 1.0
    for key, row in loadings.iterrows():
        dx = float(row["PC1_loading"]) * scale
        dy = float(row["PC2_loading"]) * scale
        colour = fs.MODULE_COLOUR[fs.module(key)]
        ax_pca.arrow(
            0, 0, dx, dy, width=0.005, head_width=0.06, head_length=0.08,
            length_includes_head=True, color=colour, alpha=0.60,
        )
        if key in label_keys:
            ax_pca.text(
                dx * 1.06, dy * 1.06, endpoint_label(key, short=True),
                fontsize=fs.FONT["footnote"], color=colour, ha="center", va="center",
            )

    heat_axes = [fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1])]
    images = []
    short_labels = fs.labels(fs.S7_MATRIX_KEYS, short=True)
    for panel_index, (ax, scale_name, title) in enumerate(zip(
        heat_axes,
        ("within_taxon", "among_taxon"),
        ("Within taxa: 1,734 observations, 42 taxa", "Among taxon medians: 42 taxa"),
    )):
        image = ax.imshow(
            _strength_matrix(matrix_rows, scale_name), vmin=0, vmax=1,
            cmap=fs.CMAP_SEQUENTIAL, interpolation="nearest",
        )
        images.append(image)
        ax.set_xticks(range(len(short_labels)), short_labels, rotation=68, ha="right", fontsize=5.5)
        ax.set_yticks(range(len(short_labels)), short_labels if panel_index == 0 else [], fontsize=5.5)
        fs.panel(ax, "b" if panel_index == 0 else "c", title)
        for boundary in (0.5, 3.5, 7.5, 10.5):
            ax.axhline(boundary, color="white", linewidth=0.7)
            ax.axvline(boundary, color="white", linewidth=0.7)
        ax.tick_params(length=0)
    colourbar_axis = fig.add_axes([0.34, 0.035, 0.32, 0.016])
    colourbar = fig.colorbar(images[-1], cax=colourbar_axis, orientation="horizontal")
    colourbar.set_label("Association strength")
    fig.subplots_adjust(left=0.12, right=0.98, bottom=0.15, top=0.98)
    save(fig, output, SUPP_STEMS[4])


# ---------------------------------------------------------------------------
# Build report
# ---------------------------------------------------------------------------


def main() -> None:
    args = parse_args()
    fs.use(grid=False)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    figure_measurement(args.source_dir, args.output_dir)
    figure_geographic_domain(args.source_dir, args.output_dir)
    figure_information_loss(args.input_dir, args.output_dir)
    figure_scale_atlas(args.input_dir, args.output_dir)
    figure_candidates(args.input_dir, args.output_dir)
    figure_s3_endpoint_support(args.input_dir, args.output_dir)
    figure_s4_sampling_audit(args.input_dir, args.output_dir)
    figure_s5_spatial_diagnostics(args.input_dir, args.output_dir)
    figure_s6_historical_placement(args.input_dir, args.output_dir)
    figure_s7_whole_capitulum(args.input_dir, args.output_dir)

    produced = [
        args.output_dir / f"{stem}{suffix}"
        for stem in MAIN_STEMS + SUPP_STEMS
        for suffix in (".png", ".pdf")
    ]
    missing = [str(path) for path in produced if not path.exists()]
    if missing:
        raise FileNotFoundError(f"figure builder did not produce expected outputs: {missing}")

    inputs = [
        args.input_dir / "v2_full27_endpoint_inventory.csv",
        args.input_dir / "v2_full27_variance_decomposition.csv",
        args.input_dir / "v2_full27_environment_cross_scale.csv",
        args.input_dir / "spatial" / "v2_full27_spatial_within.csv",
        args.input_dir / "spatial" / "v2_full27_spatial_among.csv",
        args.input_dir / "historical" / "v2_full27_historical_placement_summary.csv",
        args.input_dir / "historical" / "v2_full27_historical_placement_models.csv",
        args.input_dir / "sampling" / "v2_full27_sampling_composition_summary.csv",
        args.input_dir / "secondary" / "pca" / "geb_v2_taxon_morphospace_scores.csv",
        args.input_dir / "secondary" / "pca" / "geb_v2_taxon_morphospace_loadings.csv",
        args.input_dir / "secondary" / "multilevel" / "capitulum_space_unit_strength_matrices.csv",
        args.source_dir / "Figure_1_real_photo_measurement_provenance.png",
        args.source_dir / "Figure_1_real_photo_measurement_provenance.csv",
        args.source_dir / "Figure_2_v2_realized_sampling_cells.csv",
        args.source_dir / "natural_earth_110m_land_outline.json",
    ]
    unavailable = [str(path) for path in inputs if not path.is_file()]
    if unavailable:
        raise FileNotFoundError(f"required external figure inputs are missing: {unavailable}")

    report = {
        "status": "ok",
        "style_module": "analysis/figures/azami_figstyle.py",
        "physical_width_inches": fs.WIDTH["double"],
        "raster_dpi": fs.DPI_RASTER,
        "analysis_units": len(fs.order(split_hue=False)),
        "measured_endpoints": len(fs.order(split_hue=True, measured_only=True)),
        "pca_endpoint_count": len(fs.PCA_ENDPOINT_KEYS),
        "s7_matrix_unit_count": len(fs.S7_MATRIX_KEYS),
        "main_figure_stems": MAIN_STEMS,
        "supporting_figure_stems": SUPP_STEMS,
        "input_sha256": {str(path): sha256(path) for path in inputs},
        "claim_boundary": "presentation-only rebuild of external frozen v2 outputs; no new inference or selection",
    }
    (args.output_dir / "figure_build_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
