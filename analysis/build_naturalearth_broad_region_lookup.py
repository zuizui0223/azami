#!/usr/bin/env python3
"""Assign reproducible broad geographic regions from Natural Earth polygons.

The output is a diagnostic lookup, not a taxonomic or ecological classification.
Polygon-contained assignments are deterministic. Points that require a nearest-
polygon fallback are flagged for review rather than silently treated as exact.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--observations", required=True, type=Path)
    p.add_argument("--countries", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--report", required=True, type=Path)
    p.add_argument("--observation-id", default="obs_id")
    p.add_argument("--latitude", default="latitude")
    p.add_argument("--longitude", default="longitude")
    p.add_argument("--nearest-max-km", type=float, default=100.0)
    return p.parse_args()


def pick_column(columns: list[str], candidates: tuple[str, ...], label: str) -> str:
    lookup = {c.upper(): c for c in columns}
    for candidate in candidates:
        if candidate.upper() in lookup:
            return lookup[candidate.upper()]
    raise ValueError(f"Natural Earth file has no recognized {label} column: {columns}")


def main() -> int:
    args = parse_args()
    observations = pd.read_csv(args.observations, low_memory=False)
    required = [args.observation_id, args.latitude, args.longitude]
    missing = [c for c in required if c not in observations.columns]
    if missing:
        raise ValueError(f"Observation table lacks columns: {missing}")
    if observations[args.observation_id].duplicated().any():
        raise ValueError("Observation table must be unique by observation id")

    work = observations[required].copy()
    work[args.latitude] = pd.to_numeric(work[args.latitude], errors="coerce")
    work[args.longitude] = pd.to_numeric(work[args.longitude], errors="coerce")
    if work[[args.latitude, args.longitude]].isna().any().any():
        raise ValueError("Region lookup requires finite coordinates for every row")

    countries = gpd.read_file(args.countries).to_crs(4326)
    continent_col = pick_column(
        list(countries.columns),
        ("CONTINENT", "continent"),
        "continent",
    )
    country_col = pick_column(
        list(countries.columns),
        ("ADMIN", "SOVEREIGNT", "NAME", "NAME_EN"),
        "country name",
    )
    polygons = countries[[continent_col, country_col, "geometry"]].copy()
    polygons = polygons.rename(
        columns={continent_col: "broad_region", country_col: "country_name"}
    )
    polygons = polygons.loc[
        polygons.geometry.notna()
        & polygons["broad_region"].notna()
        & ~polygons["broad_region"].astype(str).str.contains("Seven seas", case=False, na=False)
    ].copy()

    points = gpd.GeoDataFrame(
        work,
        geometry=gpd.points_from_xy(work[args.longitude], work[args.latitude]),
        crs=4326,
    )
    joined = gpd.sjoin(
        points,
        polygons,
        how="left",
        predicate="within",
    )
    if joined.index.duplicated().any():
        duplicated = joined.index[joined.index.duplicated(keep=False)].unique().tolist()
        raise ValueError(f"Boundary join produced multiple polygon matches for rows: {duplicated[:10]}")
    joined = joined.sort_index()
    joined["assignment_method"] = np.where(
        joined["broad_region"].notna(), "polygon_within", "unmapped"
    )
    joined["nearest_distance_km"] = 0.0

    # Natural Earth's medium-scale polygons omit some tiny islands. For those
    # only, use a projected nearest-country fallback with an explicit distance
    # threshold and flag the result for review.
    missing_mask = joined["broad_region"].isna()
    if missing_mask.any():
        projected_crs = "EPSG:8857"  # Equal Earth, metres
        missing_points = joined.loc[missing_mask, ["geometry"]].to_crs(projected_crs)
        projected_polygons = polygons.to_crs(projected_crs)
        nearest = gpd.sjoin_nearest(
            missing_points,
            projected_polygons,
            how="left",
            distance_col="distance_m",
        )
        if nearest.index.duplicated().any():
            nearest = nearest.loc[~nearest.index.duplicated(keep="first")]
        for idx, row in nearest.iterrows():
            distance_km = float(row.get("distance_m", np.nan)) / 1000.0
            if np.isfinite(distance_km) and distance_km <= args.nearest_max_km:
                joined.at[idx, "broad_region"] = row["broad_region"]
                joined.at[idx, "country_name"] = row["country_name"]
                joined.at[idx, "assignment_method"] = "nearest_country_fallback"
                joined.at[idx, "nearest_distance_km"] = distance_km

    joined["broad_region"] = joined["broad_region"].fillna("UNMAPPED")
    joined["country_name"] = joined["country_name"].fillna("")
    joined["review_status"] = np.where(
        joined["assignment_method"].eq("polygon_within"),
        "routine",
        "review_required",
    )

    output = pd.DataFrame(
        {
            args.observation_id: joined[args.observation_id].astype(str),
            "broad_region": joined["broad_region"].astype(str),
            "country_name": joined["country_name"].astype(str),
            "assignment_method": joined["assignment_method"].astype(str),
            "nearest_distance_km": pd.to_numeric(joined["nearest_distance_km"], errors="coerce"),
            "review_status": joined["review_status"].astype(str),
        }
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, index=False, encoding="utf-8-sig")

    report = {
        "n_observations": int(len(output)),
        "n_regions": int(output.loc[output.broad_region.ne("UNMAPPED"), "broad_region"].nunique()),
        "region_counts": output["broad_region"].value_counts().to_dict(),
        "assignment_methods": output["assignment_method"].value_counts().to_dict(),
        "n_review_required": int(output["review_status"].eq("review_required").sum()),
        "n_unmapped": int(output["broad_region"].eq("UNMAPPED").sum()),
        "nearest_max_km": args.nearest_max_km,
        "region_definition": "Natural Earth admin-0 CONTINENT field; nearest-country fallback only for small-island polygon omissions",
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
