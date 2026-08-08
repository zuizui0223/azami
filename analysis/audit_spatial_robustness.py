#!/usr/bin/env python3
"""Audit residual spatial autocorrelation and broad-region robustness.

This is a diagnostic layer over frozen Chapter 1 outputs. It does not refit or
replace the accepted SPDE-INLA models. Input rows must contain coordinates,
model residuals and an explicit broad-region label supplied by the frozen
artifact or a reviewed mapping step.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

from azami_ch1.provenance import write_json
from azami_ch1.tabular import require_columns, require_complete_text


def lonlat_to_unit_xyz(lonlat: np.ndarray) -> np.ndarray:
    """Map longitude/latitude degrees to the 3-D unit sphere for global kNN.

    Moran's I below uses binary k-nearest-neighbour weights, so Euclidean chord
    distance on the unit sphere gives a stable global neighbour ordering without
    treating degrees of longitude as equal distances at all latitudes or breaking
    neighbourhoods across the antimeridian.
    """
    lonlat = np.asarray(lonlat, dtype=float)
    if lonlat.ndim != 2 or lonlat.shape[1] != 2:
        raise ValueError("lonlat must be an n x 2 array ordered longitude, latitude")
    lon = np.deg2rad(lonlat[:, 0])
    lat = np.deg2rad(lonlat[:, 1])
    cos_lat = np.cos(lat)
    return np.column_stack(
        [
            cos_lat * np.cos(lon),
            cos_lat * np.sin(lon),
            np.sin(lat),
        ]
    )


def knn_neighbors(xy: np.ndarray, k: int) -> np.ndarray:
    xy = np.asarray(xy, dtype=float)
    if len(xy) < k + 2:
        return np.empty((0, k), dtype=int)
    tree = cKDTree(xy)
    _, indices = tree.query(xy, k=k + 1)
    return np.asarray(indices[:, 1:], dtype=int)


def morans_i_from_neighbors(values: np.ndarray, neighbours: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    n = len(values)
    if neighbours.size == 0 or neighbours.shape[0] != n:
        return float("nan")
    k = neighbours.shape[1]
    z = values - np.nanmean(values)
    numerator = np.sum(z[:, None] * z[neighbours])
    denominator = np.sum(z**2)
    return float((n / (n * k)) * numerator / denominator) if denominator else float("nan")


def morans_i(values: np.ndarray, xy: np.ndarray, k: int) -> float:
    return morans_i_from_neighbors(values, knn_neighbors(xy, k))


def permutation_p(
    values: np.ndarray,
    xy: np.ndarray,
    k: int,
    permutations: int,
    seed: int,
) -> tuple[float, float]:
    neighbours = knn_neighbors(xy, k)
    observed = morans_i_from_neighbors(values, neighbours)
    if not np.isfinite(observed):
        return observed, float("nan")
    rng = np.random.default_rng(seed)
    null = np.array(
        [
            morans_i_from_neighbors(rng.permutation(values), neighbours)
            for _ in range(permutations)
        ]
    )
    p = (1 + np.sum(np.abs(null) >= abs(observed))) / (permutations + 1)
    return observed, float(p)


def region_summary(df: pd.DataFrame, region: str, taxon: str) -> pd.DataFrame:
    counts = df.groupby(region, dropna=False).size().rename("n_rows").reset_index()
    taxa = df.groupby(region, dropna=False)[taxon].nunique().rename("n_taxa").reset_index()
    out = counts.merge(taxa, on=region)
    out["row_fraction"] = out["n_rows"] / len(df)
    return out.sort_values("n_rows", ascending=False)


def leave_one_region_out_rank(
    df: pd.DataFrame,
    region: str,
    taxon: str,
    value: str,
) -> pd.DataFrame:
    full = df.groupby(taxon)[value].mean()
    rows: list[dict[str, object]] = []
    for held_out in sorted(df[region].dropna().astype(str).unique()):
        reduced = df[df[region].astype(str) != held_out].groupby(taxon)[value].mean()
        shared = full.index.intersection(reduced.index)
        rho = (
            spearmanr(full.loc[shared], reduced.loc[shared]).statistic
            if len(shared) >= 3
            else np.nan
        )
        rows.append(
            {
                "held_out_region": held_out,
                "n_shared_taxa": len(shared),
                "spearman_rho": rho,
            }
        )
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--latitude", default="latitude")
    parser.add_argument("--longitude", default="longitude")
    parser.add_argument("--residual", default="residual")
    parser.add_argument("--region", default="broad_region")
    parser.add_argument("--taxon", default="taxon_name")
    parser.add_argument("--endpoint", default="endpoint")
    parser.add_argument("--value", default="residual")
    parser.add_argument("--k", type=int, default=8)
    parser.add_argument("--permutations", type=int, default=999)
    parser.add_argument("--seed", type=int, default=20260717)
    parser.add_argument("--max-rows-per-endpoint", type=int, default=10000)
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    require_columns(
        df,
        [
            args.latitude,
            args.longitude,
            args.residual,
            args.region,
            args.taxon,
            args.endpoint,
            args.value,
        ],
        "spatial diagnostic table",
    )
    require_complete_text(df, args.region, "broad-region labels")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    region_summary(df, args.region, args.taxon).to_csv(
        args.output_dir / "broad_region_coverage.csv",
        index=False,
    )

    moran_rows = []
    rank_tables = []
    for index, (endpoint, part) in enumerate(df.groupby(args.endpoint, dropna=False)):
        clean = part.dropna(
            subset=[args.latitude, args.longitude, args.residual]
        ).copy()
        if len(clean) > args.max_rows_per_endpoint:
            clean = clean.sample(
                args.max_rows_per_endpoint,
                random_state=args.seed + index,
            )
        lonlat = clean[[args.longitude, args.latitude]].to_numpy(float)
        xy = lonlat_to_unit_xyz(lonlat)
        values = clean[args.residual].to_numpy(float)
        statistic, p_value = permutation_p(
            values,
            xy,
            args.k,
            args.permutations,
            args.seed + index,
        )
        moran_rows.append(
            {
                "endpoint": endpoint,
                "n": len(clean),
                "k": args.k,
                "morans_i": statistic,
                "permutation_p": p_value,
            }
        )
        ranks = leave_one_region_out_rank(
            part.dropna(subset=[args.value]),
            args.region,
            args.taxon,
            args.value,
        )
        ranks.insert(0, "endpoint", endpoint)
        rank_tables.append(ranks)

    pd.DataFrame(moran_rows).to_csv(
        args.output_dir / "residual_morans_i.csv",
        index=False,
    )
    pd.concat(rank_tables, ignore_index=True).to_csv(
        args.output_dir / "leave_one_region_out_rank_stability.csv",
        index=False,
    )

    write_json(
        args.output_dir / "spatial_robustness_summary.json",
        {
            "input": str(args.input),
            "n_rows": int(len(df)),
            "n_taxa": int(df[args.taxon].nunique()),
            "n_regions": int(df[args.region].nunique()),
            "n_endpoints": int(df[args.endpoint].nunique()),
            "k": args.k,
            "permutations": args.permutations,
            "seed": args.seed,
            "neighbor_geometry": "kNN on 3-D unit-sphere chord distance from longitude/latitude",
            "permutation_rule": "neighbor graph fixed once per endpoint; residual values permuted over fixed graph",
            "scope": "diagnostic only; frozen models and claims unchanged",
        },
        include_generated_utc=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
