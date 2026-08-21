#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize

plt.rcParams["svg.hashsalt"] = "azami-main-figure2"

EXPECTED = {
    "all_detected": (406582, 286),
    "coordinate_usable": (392989, 271),
    "strict_10km": (297293, 259),
    "primary": (46276, 259),
}

PRIMARY_TRAITS = [
    "orientation_angle_degrees_median",
    "corolla_lab_lightness_median",
    "corolla_lab_chroma_median",
    "corolla_hue_sin_median",
    "corolla_hue_cos_median",
    "shape_aspect_ratio_median",
    "shape_circularity_median",
    "shape_solidity_median",
    "shape_width_cv_median",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--all-detected", required=True)
    p.add_argument("--coordinate-usable", required=True)
    p.add_argument("--strict-10km", required=True)
    p.add_argument("--primary", required=True)
    p.add_argument("--environment", required=True)
    p.add_argument("--countries", required=True)
    p.add_argument("--out-dir", required=True)
    return p.parse_args()


def load(path: str) -> pd.DataFrame:
    return pd.read_csv(path, low_memory=False)


def validate(name: str, df: pd.DataFrame) -> None:
    n, taxa = EXPECTED[name]
    if len(df) != n:
        raise ValueError(f"{name}: expected {n} rows, got {len(df)}")
    if df["taxon_name"].nunique() != taxa:
        raise ValueError(f"{name}: expected {taxa} taxa, got {df['taxon_name'].nunique()}")
    for c in ("latitude", "longitude"):
        if c not in df:
            raise ValueError(f"{name}: missing {c}")


def valid_coords(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["latitude"] = pd.to_numeric(out["latitude"], errors="coerce")
    out["longitude"] = pd.to_numeric(out["longitude"], errors="coerce")
    out = out.loc[
        out["latitude"].between(-90, 90)
        & out["longitude"].between(-180, 180)
    ].copy()
    return out


def density_grid(df: pd.DataFrame, deg: float, stage: str) -> pd.DataFrame:
    x = valid_coords(df)
    x["lon_cell"] = np.floor((x["longitude"] + 180.0) / deg) * deg - 180.0 + deg / 2
    x["lat_cell"] = np.floor((x["latitude"] + 90.0) / deg) * deg - 90.0 + deg / 2
    g = x.groupby(["lon_cell", "lat_cell"], as_index=False).size().rename(columns={"size": "n"})
    g.insert(0, "stage", stage)
    return g


def assessability_grid(df: pd.DataFrame, deg: float = 2.0, min_n: int = 20) -> pd.DataFrame:
    x = valid_coords(df)
    missing = [c for c in PRIMARY_TRAITS if c not in x]
    if missing:
        raise ValueError(f"Assessability source missing trait columns: {missing}")
    x["lon_cell"] = np.floor((x["longitude"] + 180.0) / deg) * deg - 180.0 + deg / 2
    x["lat_cell"] = np.floor((x["latitude"] + 90.0) / deg) * deg - 90.0 + deg / 2
    x["orientation_usable"] = pd.to_numeric(x[PRIMARY_TRAITS[0]], errors="coerce").notna()
    x["colour_usable"] = (
        pd.to_numeric(x[PRIMARY_TRAITS[1]], errors="coerce").notna()
        & pd.to_numeric(x[PRIMARY_TRAITS[2]], errors="coerce").notna()
        & pd.to_numeric(x[PRIMARY_TRAITS[3]], errors="coerce").notna()
        & pd.to_numeric(x[PRIMARY_TRAITS[4]], errors="coerce").notna()
    )
    x["outline_usable"] = (
        pd.to_numeric(x[PRIMARY_TRAITS[5]], errors="coerce").notna()
        & pd.to_numeric(x[PRIMARY_TRAITS[6]], errors="coerce").notna()
        & pd.to_numeric(x[PRIMARY_TRAITS[7]], errors="coerce").notna()
        & pd.to_numeric(x[PRIMARY_TRAITS[8]], errors="coerce").notna()
    )
    trait_ok = np.column_stack([
        pd.to_numeric(x[c], errors="coerce").notna().to_numpy(float)
        for c in PRIMARY_TRAITS
    ])
    x["mean_endpoint_usable"] = trait_ok.mean(axis=1)
    g = x.groupby(["lon_cell", "lat_cell"], as_index=False).agg(
        n=("obs_id", "size"),
        orientation_usable=("orientation_usable", "mean"),
        colour_usable=("colour_usable", "mean"),
        outline_usable=("outline_usable", "mean"),
        mean_endpoint_usable=("mean_endpoint_usable", "mean"),
    )
    return g.loc[g["n"] >= min_n].copy()


def world_layer(path: str) -> gpd.GeoDataFrame:
    world = gpd.read_file(path)
    if "CONTINENT" in world.columns:
        world = world.loc[world["CONTINENT"] != "Antarctica"].copy()
    return world


def base_map(ax, world: gpd.GeoDataFrame) -> None:
    world.plot(ax=ax, facecolor="#f2f2f2", edgecolor="#777777", linewidth=0.35, zorder=0)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-60, 85)
    ax.set_xticks([-120, -60, 0, 60, 120])
    ax.set_yticks([-30, 0, 30, 60])
    ax.tick_params(labelsize=7, length=2)
    ax.set_xlabel("Longitude", fontsize=8)
    ax.set_ylabel("Latitude", fontsize=8)


def panel_label(ax, label: str) -> None:
    ax.text(0.01, 0.98, label, transform=ax.transAxes, va="top", ha="left", fontsize=12, fontweight="bold")


def save(fig, out: Path, stem: str) -> None:
    fig.savefig(out / f"{stem}.svg", bbox_inches="tight", metadata={"Date": None})
    fig.savefig(out / f"{stem}.png", dpi=300, bbox_inches="tight", metadata={"Software": "azami"})
    fig.savefig(out / f"{stem}.pdf", bbox_inches="tight", metadata={"Creator": "azami", "CreationDate": None, "ModDate": None})
    plt.close(fig)


def plot_main(world, densities, primary, env, out: Path) -> None:
    fig = plt.figure(figsize=(12.5, 8.4))
    gs = fig.add_gridspec(2, 2, hspace=0.30, wspace=0.20)

    ax = fig.add_subplot(gs[0, 0])
    base_map(ax, world)
    a = densities.loc[densities["stage"] == "all_detected"]
    sc = ax.scatter(a["lon_cell"], a["lat_cell"], c=np.log10(a["n"] + 1), s=8, marker="s", linewidths=0, cmap="viridis", zorder=2)
    cb = fig.colorbar(sc, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label("log10(observations per 1° cell + 1)", fontsize=7)
    cb.ax.tick_params(labelsize=6)
    ax.set_title("Detector-positive observations before spatial filtering\n406,582 observations · 286 taxa", fontsize=9)
    panel_label(ax, "A")

    ax = fig.add_subplot(gs[0, 1])
    base_map(ax, world)
    p = valid_coords(primary)
    ax.scatter(p["longitude"], p["latitude"], s=1.7, alpha=0.20, linewidths=0, rasterized=True, zorder=2)
    ax.set_title("Primary spatially thinned analytical cohort\n46,276 observations · 259 taxa · one taxon × 0.25° cell", fontsize=9)
    panel_label(ax, "B")

    ax = fig.add_subplot(gs[1, 0])
    ax.axis("off")
    panel_label(ax, "C")
    ax.set_title("Cohort flow and analysis roles", fontsize=9, pad=4)
    boxes = [
        (0.06, 0.83, "Detector-positive\n406,582 obs · 286 taxa"),
        (0.06, 0.61, "Coordinate usable\n392,989 obs · 271 taxa"),
        (0.06, 0.39, "≤10 km positional accuracy\n297,293 obs · 259 taxa"),
        (0.06, 0.17, "Primary thinned\n46,276 obs · 259 taxa"),
    ]
    for x, y, txt in boxes:
        ax.text(x, y, txt, transform=ax.transAxes, ha="left", va="center", fontsize=8,
                bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor="#555555", linewidth=0.8))
    for y1, y2 in [(0.76, 0.68), (0.54, 0.46), (0.32, 0.24)]:
        ax.annotate("", xy=(0.17, y2), xytext=(0.17, y1), xycoords=ax.transAxes,
                    arrowprops=dict(arrowstyle="->", linewidth=0.9, color="#555555"))
    ax.text(0.55, 0.77, "Balanced image atlas\n3,725 observations\n6,626 heads · 216 taxa\nvariance + PCA + historical", transform=ax.transAxes,
            ha="center", va="center", fontsize=8,
            bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor="#555555", linewidth=0.8))
    ax.text(0.55, 0.45, "Grouped SPDE-INLA\n31,666–34,472 obs/endpoint\n139–141 taxa\nspatial environmental models", transform=ax.transAxes,
            ha="center", va="center", fontsize=8,
            bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor="#555555", linewidth=0.8))
    ax.text(0.55, 0.18, "High-resolution involucre\n1,443 heads · 1,292 obs\n210 taxa\n904 obs · 165 taxa at ≤10 km", transform=ax.transAxes,
            ha="center", va="center", fontsize=8,
            bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor="#555555", linewidth=0.8))
    ax.annotate("", xy=(0.43, 0.77), xytext=(0.30, 0.83), xycoords=ax.transAxes,
                arrowprops=dict(arrowstyle="->", linewidth=0.8, color="#777777"))
    ax.annotate("", xy=(0.43, 0.45), xytext=(0.30, 0.17), xycoords=ax.transAxes,
                arrowprops=dict(arrowstyle="->", linewidth=0.8, color="#777777"))
    ax.annotate("", xy=(0.43, 0.18), xytext=(0.30, 0.39), xycoords=ax.transAxes,
                arrowprops=dict(arrowstyle="->", linewidth=0.8, color="#777777"))

    ax = fig.add_subplot(gs[1, 1])
    bio1_c = pd.to_numeric(env["chelsa_bio01"], errors="coerce") * 0.1 - 273.15
    bio12 = pd.to_numeric(env["chelsa_bio12"], errors="coerce")
    ok = np.isfinite(bio1_c) & np.isfinite(bio12)
    hb = ax.hexbin(bio1_c[ok], bio12[ok], gridsize=48, bins="log", mincnt=1, cmap="viridis")
    cb = fig.colorbar(hb, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label("log10(observations per hexbin)", fontsize=7)
    cb.ax.tick_params(labelsize=6)
    ax.set_xlabel("CHELSA BIO1 annual mean temperature (°C)", fontsize=8)
    ax.set_ylabel("CHELSA BIO12 annual precipitation (mm)", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.set_title("Environmental domain of the primary cohort\n46,276 spatially thinned observations", fontsize=9)
    panel_label(ax, "D")

    fig.suptitle("Geographic sampling and analytical domain", fontsize=13, y=0.995)
    save(fig, out, "Figure_2_geographic_sampling_and_analysis_domain")


def plot_s6(world, densities, out: Path) -> None:
    stages = [
        ("all_detected", "A  Detector-positive", 406582),
        ("coordinate_usable", "B  Coordinate usable", 392989),
        ("strict_10km", "C  ≤10 km accuracy", 297293),
        ("primary", "D  Primary thinned", 46276),
    ]
    vmax = float(np.log10(densities["n"].max() + 1))
    norm = Normalize(0, vmax)
    fig, axes = plt.subplots(2, 2, figsize=(12, 7.5))
    for ax, (stage, title, n) in zip(axes.flat, stages):
        base_map(ax, world)
        x = densities.loc[densities["stage"] == stage]
        sc = ax.scatter(x["lon_cell"], x["lat_cell"], c=np.log10(x["n"] + 1), s=7, marker="s",
                        linewidths=0, cmap="viridis", norm=norm, zorder=2)
        ax.set_title(f"{title}\n{n:,} observations", fontsize=9)
    cb = fig.colorbar(sc, ax=axes.ravel().tolist(), fraction=0.018, pad=0.015)
    cb.set_label("log10(observations per 1° cell + 1)", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    fig.suptitle("Sampling geography across filtering stages", fontsize=13)
    save(fig, out, "Figure_S6_sampling_geography_across_filters")


def plot_s7(world, assess, out: Path) -> None:
    cols = [
        ("orientation_usable", "A  Orientation usable"),
        ("colour_usable", "B  Colour usable"),
        ("outline_usable", "C  Outline usable"),
        ("mean_endpoint_usable", "D  Mean usable fraction across 9 endpoints"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(12, 7.5))
    norm = Normalize(0, 1)
    for ax, (col, title) in zip(axes.flat, cols):
        base_map(ax, world)
        sc = ax.scatter(assess["lon_cell"], assess["lat_cell"], c=assess[col], s=16, marker="s",
                        linewidths=0, cmap="viridis", norm=norm, zorder=2)
        ax.set_title(title + "\n2° cells with ≥20 coordinate-usable observations", fontsize=9)
    cb = fig.colorbar(sc, ax=axes.ravel().tolist(), fraction=0.018, pad=0.015)
    cb.set_label("Usable fraction", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    fig.suptitle("Geographic variation in image-trait assessability", fontsize=13)
    save(fig, out, "Figure_S7_geographic_trait_assessability")


def main() -> None:
    args = parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    all_detected = load(args.all_detected)
    coordinate = load(args.coordinate_usable)
    strict = load(args.strict_10km)
    primary = load(args.primary)
    env = load(args.environment)
    for name, frame in [
        ("all_detected", all_detected),
        ("coordinate_usable", coordinate),
        ("strict_10km", strict),
        ("primary", primary),
    ]:
        validate(name, frame)
    if len(env) != EXPECTED["primary"][0] or env["taxon_name"].nunique() != EXPECTED["primary"][1]:
        raise ValueError("Enriched environment does not match frozen primary cohort")
    if set(env["obs_id"].astype(str)) != set(primary["obs_id"].astype(str)):
        raise ValueError("Environment and primary cohort obs_id sets differ")

    densities = pd.concat([
        density_grid(all_detected, 1.0, "all_detected"),
        density_grid(coordinate, 1.0, "coordinate_usable"),
        density_grid(strict, 1.0, "strict_10km"),
        density_grid(primary, 1.0, "primary"),
    ], ignore_index=True)
    assess = assessability_grid(coordinate)

    densities.to_csv(out / "Figure2_S6_sampling_density_1deg.csv", index=False)
    assess.to_csv(out / "FigureS7_trait_assessability_2deg.csv", index=False)
    pd.DataFrame({
        "obs_id": env["obs_id"].astype(str),
        "bio1_c": pd.to_numeric(env["chelsa_bio01"], errors="coerce") * 0.1 - 273.15,
        "bio12_mm": pd.to_numeric(env["chelsa_bio12"], errors="coerce"),
    }).to_csv(out / "Figure2D_primary_environment_points.csv", index=False)

    world = world_layer(args.countries)
    plot_main(world, densities, primary, env, out)
    plot_s6(world, densities, out)
    plot_s7(world, assess, out)

    summary = {
        "frozen_counts": {k: {"observations": v[0], "taxa": v[1]} for k, v in EXPECTED.items()},
        "coordinate_valid_all_detected": int(len(valid_coords(all_detected))),
        "sampling_density_cells_1deg": {s: int((densities["stage"] == s).sum()) for s in densities["stage"].unique()},
        "assessability_cells_2deg_n_ge_20": int(len(assess)),
        "bio1_c_range": [float((env["chelsa_bio01"] * 0.1 - 273.15).min()), float((env["chelsa_bio01"] * 0.1 - 273.15).max())],
        "bio12_mm_range": [float(env["chelsa_bio12"].min()), float(env["chelsa_bio12"].max())],
        "interpretation": "Visualization only; no new inferential family or biological claim.",
    }
    (out / "figure2_geography_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    expected_files = 3 * 3 + 4
    if len(list(out.iterdir())) != expected_files:
        raise RuntimeError(f"Unexpected release file count: {len(list(out.iterdir()))} != {expected_files}")


if __name__ == "__main__":
    main()
