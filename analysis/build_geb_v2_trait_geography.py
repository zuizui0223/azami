#!/usr/bin/env python3
"""Describe the global continuous-trait geography for GEB v2.

This is pattern-first descriptive analysis. It summarizes coverage, spatial
extent and observation-level within-taxon variance for every primary/candidate
continuous endpoint, then builds taxon-level PCA morphospaces when enough taxa
are complete. Candidate observation-level variance is not a head/photo nested
variance decomposition and is labelled separately.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--observation-long", required=True)
    p.add_argument("--species-long", required=True)
    p.add_argument("--contract", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--cell-degrees", type=float, default=5.0)
    p.add_argument("--minimum-pca-taxa", type=int, default=20)
    return p.parse_args()


def as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def within_taxon_fraction(frame: pd.DataFrame) -> float:
    if len(frame) < 2:
        return float("nan")
    values = frame["value"].to_numpy(float)
    total = float(np.sum((values - np.mean(values)) ** 2))
    if total <= 0:
        return float("nan")
    means = frame.groupby("taxon_name")["value"].transform("mean").to_numpy(float)
    within = float(np.sum((values - means) ** 2))
    return within / total


def spatial_cells(frame: pd.DataFrame, width: float) -> int:
    if "latitude" not in frame or "longitude" not in frame:
        return 0
    lat = pd.to_numeric(frame["latitude"], errors="coerce")
    lon = pd.to_numeric(frame["longitude"], errors="coerce")
    valid = lat.notna() & lon.notna()
    if not valid.any():
        return 0
    ilat = np.floor((lat[valid].to_numpy(float) + 90.0) / width).astype(int)
    ilon = np.floor((lon[valid].to_numpy(float) + 180.0) / width).astype(int)
    return len(set(zip(ilat.tolist(), ilon.tolist())))


def pca_tables(wide: pd.DataFrame, minimum_taxa: int, scope: str) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    complete = wide.dropna().copy()
    if len(complete) < minimum_taxa or complete.shape[1] < 2:
        return pd.DataFrame(), pd.DataFrame(), {
            "scope": scope, "status": "insufficient_complete_taxa",
            "n_complete_taxa": int(len(complete)), "n_endpoints": int(complete.shape[1]),
        }
    x = complete.to_numpy(float)
    means = x.mean(axis=0)
    sds = x.std(axis=0, ddof=0)
    keep = np.isfinite(sds) & (sds > 0)
    if keep.sum() < 2:
        return pd.DataFrame(), pd.DataFrame(), {
            "scope": scope, "status": "insufficient_variable_endpoints",
            "n_complete_taxa": int(len(complete)), "n_endpoints": int(keep.sum()),
        }
    columns = np.asarray(complete.columns)[keep]
    z = (x[:, keep] - means[keep]) / sds[keep]
    u, s, vt = np.linalg.svd(z, full_matrices=False)
    variance = s ** 2
    explained = variance / variance.sum()
    n_pc = min(3, vt.shape[0])
    scores = pd.DataFrame(
        u[:, :n_pc] * s[:n_pc], index=complete.index,
        columns=[f"PC{i+1}" for i in range(n_pc)],
    ).reset_index().rename(columns={complete.index.name or "index": "taxon_name"})
    loadings = pd.DataFrame({"endpoint_id": columns})
    for i in range(n_pc):
        loadings[f"PC{i+1}_loading"] = vt[i, :]
    report = {
        "scope": scope, "status": "ok", "n_complete_taxa": int(len(complete)),
        "n_endpoints": int(len(columns)),
        "explained_variance": {f"PC{i+1}": float(explained[i]) for i in range(n_pc)},
    }
    return scores, loadings, report


def main() -> None:
    args = parse_args()
    obs = pd.read_csv(args.observation_long, low_memory=False)
    species = pd.read_csv(args.species_long, low_memory=False)
    contract = pd.read_csv(args.contract, dtype=str, keep_default_na=False)
    inferential_ids = set(contract.loc[contract["analysis_tier"].isin(["primary", "candidate"]), "endpoint_id"])
    obs = obs[obs["endpoint_id"].isin(inferential_ids)].copy()
    species = species[species["endpoint_id"].isin(inferential_ids)].copy()
    obs["measurement_available_bool"] = as_bool(obs["measurement_available"])
    species["measurement_available_bool"] = as_bool(species["measurement_available"])
    obs["value"] = pd.to_numeric(obs["value"], errors="coerce")
    species["value"] = pd.to_numeric(species["value"], errors="coerce")

    summaries: list[dict[str, object]] = []
    for endpoint_id, group in obs.groupby("endpoint_id", sort=True):
        usable = group[group["measurement_available_bool"] & group["value"].notna()].copy()
        meta = contract[contract["endpoint_id"].eq(endpoint_id)].iloc[0]
        values = usable["value"].to_numpy(float)
        summaries.append({
            "endpoint_id": endpoint_id,
            "module": meta["module"],
            "analysis_tier": meta["analysis_tier"],
            "unit": meta["unit"],
            "n_observations": int(len(usable)),
            "n_taxa": int(usable["taxon_name"].nunique()),
            "n_spatial_cells_5deg": spatial_cells(usable, args.cell_degrees),
            "median": float(np.median(values)) if len(values) else np.nan,
            "q05": float(np.quantile(values, 0.05)) if len(values) else np.nan,
            "q95": float(np.quantile(values, 0.95)) if len(values) else np.nan,
            "observation_level_within_taxon_fraction": within_taxon_fraction(usable[["taxon_name", "value"]].dropna()),
            "variance_scope": "observation-level taxon decomposition; not the nested taxon-photo-head variance analysis",
        })
    summary = pd.DataFrame(summaries)

    usable_species = species[species["measurement_available_bool"] & species["value"].notna()].copy()
    wide_all = usable_species.pivot(index="taxon_name", columns="endpoint_id", values="value")
    pca_reports: list[dict[str, object]] = []
    all_scores, all_loadings, all_report = pca_tables(wide_all, args.minimum_pca_taxa, "all_inferential_endpoints")
    pca_reports.append(all_report)
    score_frames: list[pd.DataFrame] = []
    loading_frames: list[pd.DataFrame] = []
    if not all_scores.empty:
        all_scores.insert(0, "scope", "all_inferential_endpoints")
        all_loadings.insert(0, "scope", "all_inferential_endpoints")
        score_frames.append(all_scores)
        loading_frames.append(all_loadings)

    for module, endpoints in contract[contract["analysis_tier"].isin(["primary", "candidate"])].groupby("module"):
        ids = [value for value in endpoints["endpoint_id"] if value in wide_all.columns]
        if len(ids) < 2:
            continue
        scores, loadings, report = pca_tables(wide_all[ids], args.minimum_pca_taxa, f"module:{module}")
        pca_reports.append(report)
        if not scores.empty:
            scores.insert(0, "scope", f"module:{module}")
            loadings.insert(0, "scope", f"module:{module}")
            score_frames.append(scores)
            loading_frames.append(loadings)

    scores_out = pd.concat(score_frames, ignore_index=True, sort=False) if score_frames else pd.DataFrame()
    loadings_out = pd.concat(loading_frames, ignore_index=True, sort=False) if loading_frames else pd.DataFrame()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    summary.to_csv(out / "geb_v2_endpoint_geography_summary.csv", index=False)
    scores_out.to_csv(out / "geb_v2_taxon_morphospace_scores.csv", index=False)
    loadings_out.to_csv(out / "geb_v2_taxon_morphospace_loadings.csv", index=False)
    report = {
        "n_registered_inferential_endpoints": int(len(inferential_ids)),
        "n_endpoints_with_observation_measurements": int((summary["n_observations"] > 0).sum()) if not summary.empty else 0,
        "n_primary_endpoints_with_measurements": int(((summary["analysis_tier"] == "primary") & (summary["n_observations"] > 0)).sum()) if not summary.empty else 0,
        "n_candidate_endpoints_with_measurements": int(((summary["analysis_tier"] == "candidate") & (summary["n_observations"] > 0)).sum()) if not summary.empty else 0,
        "pca_scopes": pca_reports,
        "interpretation": (
            "Pattern-first GEB v2 phenotype geography. Observation-level within-taxon fractions for candidate traits "
            "are descriptive and do not replace the frozen nested head-photo-taxon decomposition for primary traits."
        ),
    }
    (out / "geb_v2_trait_geography_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
