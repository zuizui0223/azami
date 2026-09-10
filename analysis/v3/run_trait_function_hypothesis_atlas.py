#!/usr/bin/env python3
"""Hypothesis-driven species-level Chapter 1 reanalysis.

Uses the frozen 46,276 strict-spatial observation cohort and all 27 endpoint
measurements, but only tests trait-environment pairs listed in the explicit
hypothesis registry.  Deferred endpoints stay visible in the coverage ledger.

The registry was written after prior Chapter 1 outcomes were known, so this is
not a preregistered confirmatory analysis.  It is a biologically structured
reanalysis intended to replace the undifferentiated 27 x environment atlas as
a candidate main-analysis framing.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from analysis.run_geb_v2_full27_environment_atlas import bh_adjust, inferential_units, stable_rng
from analysis.run_geb_v2_full27_historical_sensitivity import load_trees, tree_covariance
from analysis.run_geb_v2_full27_spatial_sensitivity import moran_test
from analysis.v3.run_full46276_multivariable_vif10 import (
    fit_phylo_multivariable,
    freedman_lane_coefficient,
    load_and_materialize_traits,
    standardized_design,
    standardized_response,
    unit_taxon_data,
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--endpoint-contract", required=True, type=Path)
    p.add_argument("--traits-long", required=True, type=Path)
    p.add_argument("--environment", required=True, type=Path)
    p.add_argument("--recovered-display", required=True, type=Path)
    p.add_argument("--registry", required=True, type=Path)
    p.add_argument("--tree-dir", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--minimum-trait-observations-per-taxon", type=int, default=5)
    p.add_argument("--minimum-taxa", type=int, default=30)
    p.add_argument("--permutations", type=int, default=9999)
    p.add_argument("--moran-permutations", type=int, default=999)
    p.add_argument("--seed", type=int, default=20260910)
    return p.parse_args()


def load_registry(path: Path, units: list[dict[str, Any]]) -> dict[str, Any]:
    registry = json.loads(path.read_text(encoding="utf-8"))
    tests = registry["tests"]
    if len({row["hypothesis_id"] for row in tests}) != len(tests):
        raise ValueError("hypothesis_id must be unique")
    unit_ids = {unit["unit_id"] for unit in units}
    tested = {row["unit_id"] for row in tests}
    deferred = {row["unit_id"] for row in registry["deferred_units"]}
    unknown = (tested | deferred) - unit_ids
    if unknown:
        raise ValueError(f"registry contains unknown inferential units: {sorted(unknown)}")
    uncovered = unit_ids - (tested | deferred)
    if uncovered:
        raise ValueError(f"registry does not account for all 26 inferential units: {sorted(uncovered)}")
    overlap = tested & deferred
    if overlap:
        raise ValueError(f"units cannot be both tested and deferred: {sorted(overlap)}")
    return registry


def load_environment(path: Path, predictors: list[str]) -> pd.DataFrame:
    cols = ["obs_id", "taxon_name", "latitude", "longitude", *predictors]
    env = pd.read_csv(path, usecols=cols, low_memory=False)
    env["obs_id"] = env["obs_id"].astype(str)
    env["taxon_name"] = env["taxon_name"].astype(str)
    if len(env) != 46276 or env["obs_id"].nunique() != 46276:
        raise ValueError(f"Expected 46,276 unique frozen environment observations; found {len(env)}")
    for col in ["latitude", "longitude", *predictors]:
        env[col] = pd.to_numeric(env[col], errors="coerce")
    for predictor in predictors:
        if env[predictor].notna().mean() < 0.98:
            raise ValueError(f"environment coverage below 0.98 for {predictor}")
    return env


def taxon_environment(env: pd.DataFrame, predictors: list[str]) -> pd.DataFrame:
    cols = [*predictors, "latitude", "longitude"]
    return env.groupby("taxon_name")[cols].median(numeric_only=True)


def direction_matches(beta: float, expected: str) -> bool:
    if expected == "positive":
        return bool(beta > 0)
    if expected == "negative":
        return bool(beta < 0)
    return True


def fit_base(
    registry: dict[str, Any],
    units_by_id: dict[str, dict[str, Any]],
    traits: pd.DataFrame,
    env_taxon: pd.DataFrame,
    min_trait_obs: int,
    min_taxa: int,
    permutations: int,
    seed: int,
) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    rows: list[dict[str, Any]] = []
    cache: dict[str, pd.DataFrame] = {}
    for test in registry["tests"]:
        unit = units_by_id[test["unit_id"]]
        predictor = test["predictor"]
        data = unit_taxon_data(unit, traits, env_taxon, min_trait_obs).dropna(
            subset=[predictor, "latitude", "longitude"]
        )
        cache[test["hypothesis_id"]] = data
        base = {
            **test,
            "module": unit["module"],
            "analysis_tier": unit["analysis_tier"],
            "validation_status": unit["validation_status"],
            "inferential_unit": unit["inferential_unit"],
            "member_endpoint_ids": "|".join(unit["member_endpoint_ids"]),
            "n_taxa": int(len(data)),
        }
        if len(data) < min_taxa:
            rows.append({**base, "status": "insufficient_taxa"})
            continue
        design, names = standardized_design(data, [predictor], include_space=False)
        response = standardized_response(data, unit["member_endpoint_ids"])
        beta, p_value, _ = freedman_lane_coefficient(
            response,
            design,
            names.index(predictor),
            permutations,
            stable_rng(seed, "hypothesis_base", test["hypothesis_id"]),
        )
        out = {**base, "status": "ok", "p_value": p_value}
        if unit["inferential_unit"] == "linear_endpoint":
            out["beta_std"] = float(beta[0])
            out["matches_expected_direction"] = direction_matches(float(beta[0]), test["expected_direction"])
        else:
            out["beta_sine_std"] = float(beta[0])
            out["beta_cosine_std"] = float(beta[1])
            out["effect_magnitude"] = float(np.linalg.norm(beta))
            out["effect_direction_degrees"] = float(math.degrees(math.atan2(beta[0], beta[1])) % 360.0)
            out["matches_expected_direction"] = True
        rows.append(out)

    result = pd.DataFrame(rows)
    result["q_fdr_bh_within_hypothesis_family"] = np.nan
    for family, part in result[result["status"].eq("ok")].groupby("family", sort=True):
        result.loc[part.index, "q_fdr_bh_within_hypothesis_family"] = bh_adjust(part["p_value"])
    return result, cache


def fit_spatial(
    base: pd.DataFrame,
    units_by_id: dict[str, dict[str, Any]],
    cache: dict[str, pd.DataFrame],
    permutations: int,
    moran_permutations: int,
    seed: int,
) -> pd.DataFrame:
    selected = base[
        base["status"].eq("ok")
        & base["q_fdr_bh_within_hypothesis_family"].lt(0.05)
        & base["matches_expected_direction"].fillna(True).astype(bool)
    ]
    rows: list[dict[str, Any]] = []
    for _, row in selected.iterrows():
        unit = units_by_id[row["unit_id"]]
        predictor = row["predictor"]
        data = cache[row["hypothesis_id"]]
        design, names = standardized_design(data, [predictor], include_space=True)
        response = standardized_response(data, unit["member_endpoint_ids"])
        beta, p_value, residual = freedman_lane_coefficient(
            response,
            design,
            names.index(predictor),
            permutations,
            stable_rng(seed, "hypothesis_spatial", row["hypothesis_id"]),
        )
        residual_for_moran = residual[:, 0] if residual.shape[1] == 1 else np.linalg.norm(residual, axis=1)
        moran_i, moran_p, moran_n = moran_test(
            residual_for_moran,
            data["latitude"].to_numpy(float),
            data["longitude"].to_numpy(float),
            moran_permutations,
            5000,
            stable_rng(seed, "hypothesis_moran", row["hypothesis_id"]),
        )
        out = {
            "hypothesis_id": row["hypothesis_id"],
            "family": row["family"],
            "unit_id": row["unit_id"],
            "predictor": predictor,
            "role": row["role"],
            "module": row["module"],
            "analysis_tier": row["analysis_tier"],
            "inferential_unit": row["inferential_unit"],
            "member_endpoint_ids": row["member_endpoint_ids"],
            "expected_direction": row["expected_direction"],
            "n_taxa": int(len(data)),
            "base_p_value": float(row["p_value"]),
            "base_q_fdr": float(row["q_fdr_bh_within_hypothesis_family"]),
            "spatial_p_value": float(p_value),
            "residual_morans_i": float(moran_i),
            "residual_morans_p_value": float(moran_p),
            "residual_morans_n": int(moran_n),
        }
        same_direction = True
        expected_ok = True
        if unit["inferential_unit"] == "linear_endpoint":
            out["base_beta_std"] = float(row["beta_std"])
            out["spatial_beta_std"] = float(beta[0])
            same_direction = bool(np.sign(beta[0]) == np.sign(float(row["beta_std"])))
            expected_ok = direction_matches(float(beta[0]), row["expected_direction"])
            out["same_direction_as_base"] = same_direction
        else:
            out["spatial_beta_sine_std"] = float(beta[0])
            out["spatial_beta_cosine_std"] = float(beta[1])
            out["spatial_effect_magnitude"] = float(np.linalg.norm(beta))
            out["same_direction_as_base"] = np.nan
        out["spatial_pass"] = bool(
            p_value < 0.05
            and np.isfinite(moran_p)
            and moran_p >= 0.05
            and same_direction
            and expected_ok
        )
        rows.append(out)
    return pd.DataFrame(rows)


def fit_phylogeny(
    spatial: pd.DataFrame,
    units_by_id: dict[str, dict[str, Any]],
    cache: dict[str, pd.DataFrame],
    tree_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if spatial.empty:
        return pd.DataFrame(), pd.DataFrame()
    selected = spatial[spatial["spatial_pass"].fillna(False).astype(bool)]
    trees = load_trees(tree_dir)
    resources = []
    for scenario, replicate, tree in trees:
        names, covariance = tree_covariance(tree)
        resources.append((scenario, replicate, names, covariance))
    if len(resources) != 52:
        raise ValueError(f"Expected 52 frozen placement trees, found {len(resources)}")

    model_rows: list[dict[str, Any]] = []
    for _, row in selected.iterrows():
        unit = units_by_id[row["unit_id"]]
        predictor = row["predictor"]
        data = cache[row["hypothesis_id"]].copy()
        data["taxon_key"] = data["taxon_name"].astype(str).str.strip().str.replace(" ", "_", regex=False)
        data = data.set_index("taxon_key", drop=False)
        for scenario, replicate, names, covariance in resources:
            lookup = {name: i for i, name in enumerate(names)}
            common = [name for name in names if name in data.index]
            base = {
                "hypothesis_id": row["hypothesis_id"],
                "family": row["family"],
                "unit_id": row["unit_id"],
                "predictor": predictor,
                "role": row["role"],
                "module": row["module"],
                "analysis_tier": row["analysis_tier"],
                "inferential_unit": row["inferential_unit"],
                "scenario": scenario,
                "replicate": replicate,
                "n_taxa": len(common),
            }
            if len(common) < 30:
                model_rows.append({**base, "status": "insufficient_taxa"})
                continue
            ordered = data.loc[common].copy()
            design, names_design = standardized_design(ordered, [predictor], include_space=False)
            response = standardized_response(ordered, unit["member_endpoint_ids"])
            positions = [lookup[name] for name in common]
            cov = covariance[np.ix_(positions, positions)]
            try:
                fit = fit_phylo_multivariable(response, design, cov, names_design.index(predictor))
                beta = fit.pop("beta")
                out = {**base, "status": "ok", **fit}
                if unit["inferential_unit"] == "linear_endpoint":
                    out["pgls_beta_std"] = float(beta[0])
                    out["same_direction_as_spatial"] = bool(
                        np.sign(beta[0]) == np.sign(float(row["spatial_beta_std"]))
                    )
                    out["matches_expected_direction"] = direction_matches(float(beta[0]), row["expected_direction"])
                else:
                    out["pgls_beta_sine_std"] = float(beta[0])
                    out["pgls_beta_cosine_std"] = float(beta[1])
                    out["pgls_effect_magnitude"] = float(np.linalg.norm(beta))
                    out["same_direction_as_spatial"] = True
                    out["matches_expected_direction"] = True
                model_rows.append(out)
            except Exception as exc:
                model_rows.append({**base, "status": f"failed:{type(exc).__name__}:{exc}"})

    models = pd.DataFrame(model_rows)
    summaries: list[dict[str, Any]] = []
    if not models.empty:
        for hypothesis_id, part in models.groupby("hypothesis_id", sort=True):
            ok = part[part["status"].eq("ok")]
            source = selected[selected["hypothesis_id"].eq(hypothesis_id)].iloc[0]
            direction_ok = bool(ok["same_direction_as_spatial"].fillna(False).all()) if not ok.empty else False
            expected_ok = bool(ok["matches_expected_direction"].fillna(False).all()) if not ok.empty else False
            passed = bool(
                len(ok) == 52
                and ok["p_value"].lt(0.05).all()
                and direction_ok
                and expected_ok
            )
            summaries.append({
                "hypothesis_id": hypothesis_id,
                "family": source["family"],
                "unit_id": source["unit_id"],
                "predictor": source["predictor"],
                "role": source["role"],
                "module": source["module"],
                "analysis_tier": source["analysis_tier"],
                "n_successful_trees": int(len(ok)),
                "n_trees_p_lt_0_05": int(ok["p_value"].lt(0.05).sum()) if not ok.empty else 0,
                "maximum_p_value": float(ok["p_value"].max()) if not ok.empty else np.nan,
                "lambda_min": float(ok["lambda"].min()) if not ok.empty else np.nan,
                "lambda_max": float(ok["lambda"].max()) if not ok.empty else np.nan,
                "phylogenetic_placement_pass": passed,
            })
    return models, pd.DataFrame(summaries)


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    contract = pd.read_csv(args.endpoint_contract, dtype=str, keep_default_na=False)
    if len(contract) != 27:
        raise ValueError(f"Expected 27 endpoint rows, found {len(contract)}")
    units = inferential_units(contract)
    if len(units) != 26:
        raise ValueError(f"Expected 26 inferential units, found {len(units)}")
    units_by_id = {u["unit_id"]: u for u in units}
    registry = load_registry(args.registry, units)
    predictors = sorted({row["predictor"] for row in registry["tests"]})
    environment = load_environment(args.environment, predictors)
    traits, recovered_counts = load_and_materialize_traits(
        args.traits_long, args.recovered_display, contract, set(environment["obs_id"])
    )
    taxon_lookup = environment.set_index("obs_id")["taxon_name"]
    missing = traits["taxon_name"].isna()
    if missing.any():
        traits.loc[missing, "taxon_name"] = traits.loc[missing, "obs_id"].map(taxon_lookup)
    env_taxon = taxon_environment(environment, predictors)

    base, cache = fit_base(
        registry, units_by_id, traits, env_taxon,
        args.minimum_trait_observations_per_taxon, args.minimum_taxa,
        args.permutations, args.seed,
    )
    spatial = fit_spatial(
        base, units_by_id, cache,
        args.permutations, args.moran_permutations, args.seed,
    )
    phylo_models, phylo_summary = fit_phylogeny(
        spatial, units_by_id, cache, args.tree_dir
    )

    base.to_csv(args.out_dir / "hypothesis_base.csv", index=False, encoding="utf-8-sig")
    spatial.to_csv(args.out_dir / "hypothesis_spatial.csv", index=False, encoding="utf-8-sig")
    phylo_models.to_csv(args.out_dir / "hypothesis_phylo_models.csv", index=False, encoding="utf-8-sig")
    phylo_summary.to_csv(args.out_dir / "hypothesis_phylo_summary.csv", index=False, encoding="utf-8-sig")

    coverage_rows = []
    tested_units = {row["unit_id"] for row in registry["tests"]}
    deferred_reason = {row["unit_id"]: row["reason"] for row in registry["deferred_units"]}
    for unit in units:
        coverage_rows.append({
            "unit_id": unit["unit_id"],
            "member_endpoint_ids": "|".join(unit["member_endpoint_ids"]),
            "module": unit["module"],
            "analysis_tier": unit["analysis_tier"],
            "ecological_test_status": "tested" if unit["unit_id"] in tested_units else "deferred",
            "deferred_reason": deferred_reason.get(unit["unit_id"], ""),
        })
    pd.DataFrame(coverage_rows).to_csv(args.out_dir / "hypothesis_unit_coverage.csv", index=False, encoding="utf-8-sig")

    family_summary = {}
    for family, part in base.groupby("family", sort=True):
        qpass = part[part["q_fdr_bh_within_hypothesis_family"].lt(0.05)]
        spatial_pass_ids = set(spatial.loc[spatial["spatial_pass"].fillna(False).astype(bool), "hypothesis_id"]) if not spatial.empty else set()
        phylo_pass_ids = set(phylo_summary.loc[phylo_summary["phylogenetic_placement_pass"].fillna(False).astype(bool), "hypothesis_id"]) if not phylo_summary.empty else set()
        family_summary[family] = {
            "registered_tests": int(len(part)),
            "base_fdr_pass": int(len(qpass)),
            "spatial_pass": int(sum(h in spatial_pass_ids for h in part["hypothesis_id"])),
            "phylogenetic_placement_pass": int(sum(h in phylo_pass_ids for h in part["hypothesis_id"])),
        }

    final_rows = []
    if not phylo_summary.empty:
        final_rows = phylo_summary[phylo_summary["phylogenetic_placement_pass"].fillna(False).astype(bool)].to_dict("records")
    summary = {
        "analysis_id": registry["analysis_id"],
        "registry_status": registry["status"],
        "cohort_observations": 46276,
        "cohort_taxa": 259,
        "registered_endpoints": 27,
        "inferential_units": 26,
        "registered_hypothesis_tests": len(registry["tests"]),
        "deferred_inferential_units": len(registry["deferred_units"]),
        "environment_predictors_used": predictors,
        "recovered_endpoint_counts": recovered_counts,
        "family_summary": family_summary,
        "final_phylogenetic_pass_rows": final_rows,
        "claim_boundary": registry["honesty_note"],
    }
    (args.out_dir / "hypothesis_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, default=str) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
