#!/usr/bin/env python3
"""Audit residual spatial autocorrelation and broad-region robustness.

This is a diagnostic layer over frozen Chapter 1 outputs. It does not refit or
replace the accepted SPDE-INLA models. Input rows must contain coordinates,
model residuals and an explicit broad-region label supplied by the frozen
artifact or a reviewed mapping step.
"""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import spearmanr


def morans_i(values: np.ndarray, xy: np.ndarray, k: int) -> float:
    n = len(values)
    if n < k + 2:
        return float("nan")
    z = values - np.nanmean(values)
    tree = cKDTree(xy)
    _, indices = tree.query(xy, k=k + 1)
    neighbours = indices[:, 1:]
    numerator = np.sum(z[:, None] * z[neighbours])
    denominator = np.sum(z**2)
    return float((n / (n * k)) * numerator / denominator) if denominator else float("nan")


def permutation_p(values: np.ndarray, xy: np.ndarray, k: int, permutations: int, seed: int) -> tuple[float, float]:
    observed = morans_i(values, xy, k)
    if not np.isfinite(observed):
        return observed, float("nan")
    rng = np.random.default_rng(seed)
    null = np.array([morans_i(rng.permutation(values), xy, k) for _ in range(permutations)])
    p = (1 + np.sum(np.abs(null) >= abs(observed))) / (permutations + 1)
    return observed, float(p)


def region_summary(df: pd.DataFrame, region: str, taxon: str) -> pd.DataFrame:
    counts = df.groupby(region, dropna=False).size().rename("n_rows").reset_index()
    taxa = df.groupby(region, dropna=False)[taxon].nunique().rename("n_taxa").reset_index()
    out = counts.merge(taxa, on=region)
    out["row_fraction"] = out["n_rows"] / len(df)
    return out.sort_values("n_rows", ascending=False)


def leave_one_region_out_rank(df: pd.DataFrame, region: str, taxon: str, value: str) -> pd.DataFrame:
    full = df.groupby(taxon)[value].mean()
    rows: list[dict[str, object]] = []
    for held_out in sorted(df[region].dropna().astype(str).unique()):
        reduced = df[df[region].astype(str) != held_out].groupby(taxon)[value].mean()
        shared = full.index.intersection(reduced.index)
        rho = spearmanr(full.loc[shared], reduced.loc[shared]).statistic if len(shared) >= 3 else np.nan
        rows.append({"held_out_region": held_out, "n_shared_taxa": len(shared), "spearman_rho": rho})
    return pd.DataFrame(rows)


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True, type=Path)
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--latitude", default="latitude")
    p.add_argument("--longitude", default="longitude")
    p.add_argument("--residual", default="residual")
    p.add_argument("--region", default="broad_region")
    p.add_argument("--taxon", default="taxon_name")
    p.add_argument("--endpoint", default="endpoint")
    p.add_argument("--value", default="residual")
    p.add_argument("--k", type=int, default=8)
    p.add_argument("--permutations", type=int, default=999)
    p.add_argument("--seed", type=int, default=20260717)
    p.add_argument("--max-rows-per-endpoint", type=int, default=10000)
    args = p.parse_args()

    df = pd.read_csv(args.input)
    required = {args.latitude, args.longitude, args.residual, args.region, args.taxon, args.endpoint, args.value}
    missing = sorted(required - set(df.columns))
    if missing:
        raise SystemExit(f"Missing required columns: {missing}")
    if df[args.region].isna().any() or (df[args.region].astype(str).str.strip() == "").any():
        raise SystemExit("Broad-region labels must be explicit and complete; automatic geographic guessing is not allowed")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    region_summary(df, args.region, args.taxon).to_csv(args.output_dir / "broad_region_coverage.csv", index=False)

    moran_rows = []
    rank_tables = []
    for i, (endpoint, part) in enumerate(df.groupby(args.endpoint, dropna=False)):
        clean = part.dropna(subset=[args.latitude, args.longitude, args.residual]).copy()
        if len(clean) > args.max_rows_per_endpoint:
            clean = clean.sample(args.max_rows_per_endpoint, random_state=args.seed + i)
        xy = clean[[args.longitude, args.latitude]].to_numpy(float)
        values = clean[args.residual].to_numpy(float)
        statistic, p_value = permutation_p(values, xy, args.k, args.permutations, args.seed + i)
        moran_rows.append({"endpoint": endpoint, "n": len(clean), "k": args.k, "morans_i": statistic, "permutation_p": p_value})
        ranks = leave_one_region_out_rank(part.dropna(subset=[args.value]), args.region, args.taxon, args.value)
        ranks.insert(0, "endpoint", endpoint)
        rank_tables.append(ranks)

    pd.DataFrame(moran_rows).to_csv(args.output_dir / "residual_morans_i.csv", index=False)
    pd.concat(rank_tables, ignore_index=True).to_csv(args.output_dir / "leave_one_region_out_rank_stability.csv", index=False)

    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "input": str(args.input),
        "n_rows": int(len(df)),
        "n_taxa": int(df[args.taxon].nunique()),
        "n_regions": int(df[args.region].nunique()),
        "n_endpoints": int(df[args.endpoint].nunique()),
        "k": args.k,
        "permutations": args.permutations,
        "seed": args.seed,
        "scope": "diagnostic only; frozen models and claims unchanged",
    }
    (args.output_dir / "spatial_robustness_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
