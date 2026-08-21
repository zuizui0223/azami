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
    for col in ("obs_id", "taxon_name", "latitude", "longitude"):
        if col not in df:
            raise ValueError(f"{name}: missing {col}")


def valid_coords(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["latitude"] = pd.to_numeric(out["latitude"], errors="coerce")
    out["longitude"] = pd.to_numeric(out["longitude"], errors="coerce")
    return out.loc[
        out["latitude"].between(-90, 90)
        & out["longitude"].between(-180, 180)
    ].copy()


def density_grid(df: pd.DataFrame, stage: str, deg: float = 1.0) -> pd.DataFrame:
    x = valid_coords(df)
    x["lon_cell"] = np.floor((x["longitude"] + 180.0) / deg) * deg - 180.0 + deg / 2
    x["lat_cell"] = np.floor((x["latitude"] + 90.0) / deg) * deg - 90.0 + deg / 2
    out = x.groupby(["lon_cell", "lat_cell"], as_index=False).size().rename(columns={"size": "n"})
    out.insert(0, "stage", stage)
    return out


def assessability_grid(df: pd.DataFrame, deg: float = 2.0, min_n: int = 20) -> pd.DataFrame:
    x = valid_coords(df)
    missing = [c for c in PRIMARY_TRAITS if c not in x]
    if missing:
        raise ValueError(f"Assessability source missing trait columns: {missing}")
    x["lon_cell"] = np.floor((x["longitude"] + 180.0) / deg) * deg - 180.0 + deg / 2
    x["lat_cell"] = np.floor((x["latitude"] + 90.0) / deg) * deg - 90.0 + deg / 2
    numeric = {c: pd.to_numeric(x[c], errors="coerce") for c in PRIMARY_TRAITS}
    x["orientation_usable"] = numeric[PRIMARY_TRAITS[0]].notna()
    x["colour_usable"] = np.logical_and.reduce([numeric[c].notna() for c in PRIMARY_TRAITS[1:5]])
    x["outline_usable"] = np.logical_and.reduce([numeric[c].notna() for c in PRIMARY_TRAITS[5:9]])
    x["mean_endpoint_usable"] = np.column_stack([numeric[c].notna().to_numpy(float) for c in PRIMARY_TRAITS]).mean(axis=1)
    out = x.groupby(["lon_cell", "lat_cell"], as_index=False).agg(
        n=("obs_id", "size"),
        orientation_usable=("orientation_usable", "mean"),
        colour_usable=("colour_usable", "mean"),
        outline_usable=("outline_usable", "mean"),
        mean_endpoint_usable=("mean_endpoint_usable", "mean"),
    )
    return out.loc[out["n"] >= min_n].copy()


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


def flow_box(ax, x: float, y: float, text: str, fontsize: float = 7.25) -> None:
    ax.text(
        x, y, text, transform=ax.transAxes, ha="center", va="center", fontsize=fontsize,
        bbox=dict(boxstyle="round,pad=0.28", facecolor="white", edgecolor="#555555", linewidth=0.75),
    )


def arrow(ax, start: tuple[float, float], end: tuple[float, float]) -> None:
    ax.annotate("", xy=end, xytext=start, xycoords=ax.transAxes,
                arrowprops=dict(arrowstyle="->", linewidth=0.8, color="#666666"))


def plot_main(world, densities, primary, env, out: Path, mappable_counts: dict[str, int]) -> None:
    fig = plt.figure(figsize=(12.5, 8.4))
    gs = fig.add_gridspec(2, 2, hspace=0.30, wspace=0.20)

    ax = fig.add_subplot(gs[0, 0])
    base_map(ax, world)
    a = densities.loc[densities["stage"] == "all_detected"]
    sc = ax.scatter(a["lon_cell"], a["lat_cell"], c=np.log10(a["n"] + 1), s=8, marker="s", linewidths=0, cmap="viridis", zorder=2)
    cb = fig.colorbar(sc, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label("log10(observations per 1° cell + 1)", fontsize=7)
    cb.ax.tick_params(labelsize=6)
    ax.set_title(
        "Detector-positive observation density before spatial filtering\n"
        f"{mappable_counts['all_detected']:,} mappable of 406,582 total · 286 taxa",
        fontsize=9,
    )
    panel_label(ax, "A")

    ax = fig.add_subplot(gs[0, 1])
    base_map(ax, world)
    p = valid_coords(primary)
    ax.scatter(p["longitude"], p["latitude"], s=1.7, alpha=0.22, linewidths=0, rasterized=True, zorder=2)
    ax.set_title("Primary spatially thinned analytical cohort\n46,276 observations · 259 taxa · one taxon × 0.25° cell", fontsize=9)
    panel_label(ax, "B")

    ax = fig.add_subplot(gs[1, 0])
    ax.axis("off")
    panel_label(ax, "C")
    ax.set_title("Two analysis streams and derived cohorts", fontsize=9, pad=4)
    flow_box(ax, 0.50, 0.93, "Public biodiversity photographs", fontsize=7.4)

    flow_box(ax, 0.23, 0.75, "Balanced image atlas\n3,725 obs · 6,626 heads · 216 taxa")
    flow_box(ax, 0.14, 0.49, "Nested variance · PCA · PGLS")
    flow_box(ax, 0.32, 0.24, "High-resolution involucre\n1,443 heads · 1,292 obs · 210 taxa\n≤10 km: 904 obs · 165 taxa", fontsize=6.9)
    arrow(ax, (0.46, 0.89), (0.27, 0.80))
    arrow(ax, (0.20, 0.68), (0.15, 0.56))
    arrow(ax, (0.28, 0.68), (0.31, 0.33))

    flow_box(ax, 0.76, 0.77, "Detector-positive\n406,582 obs · 286 taxa")
    flow_box(ax, 0.76, 0.62, "Coordinate usable\n392,989 obs · 271 taxa")
    flow_box(ax, 0.76, 0.47, "≤10 km accuracy\n297,293 obs · 259 taxa")
    flow_box(ax, 0.76, 0.32, "Primary thinned\n46,276 obs · 259 taxa")
    flow_box(ax, 0.76, 0.15, "Grouped SPDE-INLA\n31,666–34,472 obs/endpoint · 139–141 taxa", fontsize=6.9)
    arrow(ax, (0.54, 0.89), (0.72, 0.82))
    arrow(ax, (0.76, 0.71), (0.76, 0.67))
    arrow(ax, (0.76, 0.56), (0.76, 0.52))
    arrow(ax, (0.76, 0.41), (0.76, 0.37))
    arrow(ax, (0.76, 0.26), (0.76, 0.21))

    ax.text(0.02, 0.03, "Separate streams answer different scales; cohorts are not pooled under one FDR family.",
            transform=ax.transAxes, fontsize=6.5, color="#555555", ha="left", va="bottom")

    ax = fig.add_subplot(gs[1, 1])
    bio1_c = pd.to_numeric(env["chelsa_bio01"], errors="coerce") * 0.1 - 273.15
    bio12 = pd.to_numeric(env["chelsa_bio12"], errors="coerce")
    ok = np.isfinite(bio1_c) & np.isfinite(bio12)
    hb = ax.hexbin(bio1_c[ok], bio12[ok], gridsize=48, bins="log", mincnt=1, cmap="viridis")
    cb = fig.colorbar(hb, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label("Observations per hexbin (log scale)", fontsize=7)
    cb.ax.tick_params(labelsize=6)
    ax.set_xlabel("CHELSA BIO1 annual mean temperature (°C)", fontsize=8)
    ax.set_ylabel("CHELSA BIO12 annual precipitation (mm)", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.set_title("Environmental domain of the primary cohort\n46,276 spatially thinned observations", fontsize=9)
    panel_label(ax, "D")

    fig.suptitle("Geographic sampling and analytical domain", fontsize=13, y=0.995)
    save(fig, out, "Figure_2_geographic_sampling_and_analysis_domain")


def plot_s6(world, densities, out: Path, mappable_counts: dict[str, int]) -> None:
    stages = [
        ("all_detected", "A  Detector-positive", f"{mappable_counts['all_detected']:,} mapped / 406,582 total"),
        ("coordinate_usable", "B  Coordinate usable", "392,989 observations"),
        ("strict_10km", "C  ≤10 km accuracy", "297,293 observations"),
        ("primary", "D  Primary thinned", "46,276 observations"),
    ]
    vmax = float(np.log10(densities["n"].max() + 1))
    norm = Normalize(0, vmax)
    fig, axes = plt.subplots(2, 2, figsize=(12, 7.5))
    for ax, (stage, title, subtitle) in zip(axes.flat, stages):
        base_map(ax, world)
        x = densities.loc[densities["stage"] == stage]
        sc = ax.scatter(x["lon_cell"], x["lat_cell"], c=np.log10(x["n"] + 1), s=7, marker="s",
                        linewidths=0, cmap="viridis", norm=norm, zorder=2)
        ax.set_title(f"{title}\n{subtitle}", fontsize=9)
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
    frames = {
        "all_detected": load(args.all_detected),
        "coordinate_usable": load(args.coordinate_usable),
        "strict_10km": load(args.strict_10km),
        "primary": load(args.primary),
    }
    for name, frame in frames.items():
        validate(name, frame)

    env = load(args.environment)
    if len(env) != EXPECTED["primary"][0] or env["taxon_name"].nunique() != EXPECTED["primary"][1]:
        raise ValueError("Enriched environment does not match frozen primary cohort")
    if set(env["obs_id"].astype(str)) != set(frames["primary"]["obs_id"].astype(str)):
        raise ValueError("Environment and primary cohort obs_id sets differ")

    mappable_counts = {name: int(len(valid_coords(frame))) for name, frame in frames.items()}
    densities = pd.concat([density_grid(frame, name) for name, frame in frames.items()], ignore_index=True)
    assess = assessability_grid(frames["coordinate_usable"])
    densities.to_csv(out / "Figure2_S6_sampling_density_1deg.csv", index=False)
    assess.to_csv(out / "FigureS7_trait_assessability_2deg.csv", index=False)

    bio1_c = pd.to_numeric(env["chelsa_bio01"], errors="coerce") * 0.1 - 273.15
    bio12 = pd.to_numeric(env["chelsa_bio12"], errors="coerce")
    pd.DataFrame({"obs_id": env["obs_id"].astype(str), "bio1_c": bio1_c, "bio12_mm": bio12}).to_csv(
        out / "Figure2D_primary_environment_points.csv", index=False
    )

    world = world_layer(args.countries)
    plot_main(world, densities, frames["primary"], env, out, mappable_counts)
    plot_s6(world, densities, out, mappable_counts)
    plot_s7(world, assess, out)

    summary = {
        "frozen_counts": {k: {"observations": v[0], "taxa": v[1]} for k, v in EXPECTED.items()},
        "mappable_observations": mappable_counts,
        "sampling_density_cells_1deg": {s: int((densities["stage"] == s).sum()) for s in densities["stage"].unique()},
        "assessability_cells_2deg_n_ge_20": int(len(assess)),
        "assessability_cell_medians": {c: float(assess[c].median()) for c in ("orientation_usable", "colour_usable", "outline_usable", "mean_endpoint_usable")},
        "bio1_c_range": [float(bio1_c.min()), float(bio1_c.max())],
        "bio12_mm_range": [float(bio12.min()), float(bio12.max())],
        "interpretation": "Visualization only; no new inferential family or biological claim.",
    }
    (out / "figure2_geography_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    if len(list(out.iterdir())) != 13:
        raise RuntimeError(f"Unexpected release file count: {len(list(out.iterdir()))} != 13")


if __name__ == "__main__":
    main()
