#!/usr/bin/env python3
"""Align GEB v2 within-taxon and among-taxon climate structure.

The among-taxon layer uses the same four CHELSA predictors and the same
registry-defined inferential units as the GEB v2 within-taxon screen. Linear
traits are summarized by taxon medians and tested with standardized taxon-level
slopes. Circular hue is kept as one joint sine/cosine inferential unit and tested
by the magnitude of its two standardized slopes. P values are Monte Carlo
permutation tests that shuffle the taxon-level predictor across the same taxa.

This analysis is a cross-scale pattern comparison. Taxon exchangeability is a
permutation null, not a phylogenetic correction, and neither scale establishes
plasticity, adaptation, selection, or mechanism.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

DEFAULT_PREDICTORS = ["chelsa_bio01", "chelsa_bio04", "chelsa_bio12", "chelsa_bio15"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--traits-long", required=True, type=Path)
    p.add_argument("--environment", required=True, type=Path)
    p.add_argument("--within-coefficients", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--predictors", nargs="+", default=DEFAULT_PREDICTORS)
    p.add_argument("--primary-min-trait-observations", type=int, default=5)
    p.add_argument("--sensitivity-min-trait-observations", type=int, default=2)
    p.add_argument("--minimum-taxa", type=int, default=20)
    p.add_argument("--permutations", type=int, default=10000)
    p.add_argument("--seed", type=int, default=20260827)
    return p.parse_args()


def as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def bh_adjust(values: pd.Series) -> pd.Series:
    p = np.asarray(values, dtype=float)
    n = len(p)
    if n == 0:
        return pd.Series(dtype=float, index=values.index)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    out = np.empty(n, dtype=float)
    out[order] = adjusted
    return pd.Series(out, index=values.index)


def stable_rng(seed: int, *parts: str) -> np.random.Generator:
    payload = "|".join([str(seed), *map(str, parts)]).encode("utf-8")
    digest = hashlib.sha256(payload).digest()
    derived = int.from_bytes(digest[:8], "little", signed=False)
    return np.random.default_rng(derived)


def standardize(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    sd = float(np.std(values, ddof=0))
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError("No finite between-taxon variation")
    return (values - float(np.mean(values))) / sd


def slope(y: np.ndarray, x: np.ndarray) -> float:
    denom = float(np.dot(x, x))
    if denom <= 0:
        raise ValueError("Predictor has zero standardized sum of squares")
    return float(np.dot(x, y) / denom)


def permutation_p_linear(
    y: np.ndarray,
    x: np.ndarray,
    permutations: int,
    rng: np.random.Generator,
) -> tuple[float, float]:
    yz = standardize(y)
    xz = standardize(x)
    observed = slope(yz, xz)
    exceed = 0
    for _ in range(permutations):
        xp = rng.permutation(xz)
        if abs(slope(yz, xp)) >= abs(observed) - 1e-15:
            exceed += 1
    return observed, float((exceed + 1) / (permutations + 1))


def permutation_p_circular(
    sine: np.ndarray,
    cosine: np.ndarray,
    x: np.ndarray,
    permutations: int,
    rng: np.random.Generator,
) -> tuple[float, float, float, float, float]:
    ys = standardize(sine)
    yc = standardize(cosine)
    xz = standardize(x)
    beta_sine = slope(ys, xz)
    beta_cosine = slope(yc, xz)
    magnitude = float(math.hypot(beta_sine, beta_cosine))
    exceed = 0
    for _ in range(permutations):
        xp = rng.permutation(xz)
        if math.hypot(slope(ys, xp), slope(yc, xp)) >= magnitude - 1e-15:
            exceed += 1
    direction = float(math.degrees(math.atan2(beta_sine, beta_cosine)) % 360.0)
    return beta_sine, beta_cosine, magnitude, direction, float((exceed + 1) / (permutations + 1))


def load_inputs(
    traits_path: Path,
    environment_path: Path,
    predictors: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    traits = pd.read_csv(
        traits_path,
        usecols=[
            "obs_id", "taxon_name", "endpoint_id", "module", "analysis_tier",
            "analysis_eligible", "circular_group", "value",
        ],
        low_memory=False,
    )
    environment = pd.read_csv(
        environment_path,
        usecols=["obs_id", "taxon_name", *predictors],
        low_memory=False,
    )
    if traits.duplicated(["obs_id", "endpoint_id"]).any():
        raise ValueError("Trait table must be unique by obs_id/endpoint_id")
    if environment["obs_id"].astype(str).duplicated().any():
        raise ValueError("Environment table must be unique by obs_id")
    traits["analysis_eligible_bool"] = as_bool(traits["analysis_eligible"])
    traits["circular_group"] = traits["circular_group"].fillna("").astype(str).str.strip()
    traits["value"] = pd.to_numeric(traits["value"], errors="coerce")
    for predictor in predictors:
        environment[predictor] = pd.to_numeric(environment[predictor], errors="coerce")
    return traits, environment


def taxon_environment(environment: pd.DataFrame, predictors: list[str]) -> pd.DataFrame:
    counts = environment.groupby("taxon_name").size().rename("n_strict_observations")
    medians = environment.groupby("taxon_name")[predictors].median(numeric_only=True)
    return medians.join(counts).reset_index()


def inferential_units(traits: pd.DataFrame) -> list[dict[str, Any]]:
    eligible = traits[
        traits["analysis_eligible_bool"]
        & traits["analysis_tier"].isin(["primary", "candidate"])
        & traits["value"].notna()
    ].copy()
    meta = eligible.drop_duplicates("endpoint_id")[["endpoint_id", "module", "analysis_tier", "circular_group"]]
    units: list[dict[str, Any]] = []
    for _, row in meta[meta["circular_group"].eq("")].iterrows():
        units.append({
            "endpoint_id": row["endpoint_id"],
            "module": row["module"],
            "analysis_tier": row["analysis_tier"],
            "inferential_unit": "linear_endpoint",
            "members": [row["endpoint_id"]],
        })
    for group, part in meta[meta["circular_group"].ne("")].groupby("circular_group", sort=True):
        members = sorted(part["endpoint_id"].tolist())
        if len(members) != 2 or not any("sin" in x for x in members) or not any("cos" in x for x in members):
            raise ValueError(f"Circular group {group!r} must contain one sine and one cosine endpoint")
        units.append({
            "endpoint_id": group,
            "module": part.iloc[0]["module"],
            "analysis_tier": part.iloc[0]["analysis_tier"],
            "inferential_unit": "circular_joint",
            "members": members,
        })
    return sorted(units, key=lambda x: (x["analysis_tier"], x["module"], x["endpoint_id"]))


def linear_taxon_table(
    traits: pd.DataFrame,
    endpoint_id: str,
    taxon_env: pd.DataFrame,
    min_trait_observations: int,
) -> pd.DataFrame:
    part = traits[
        traits["endpoint_id"].eq(endpoint_id)
        & traits["analysis_eligible_bool"]
        & traits["value"].notna()
    ][["taxon_name", "value"]].copy()
    counts = part.groupby("taxon_name").size().rename("n_trait_observations")
    medians = part.groupby("taxon_name")["value"].median().rename("trait_value")
    summary = pd.concat([medians, counts], axis=1).reset_index()
    summary = summary[summary["n_trait_observations"].ge(min_trait_observations)]
    return summary.merge(taxon_env, on="taxon_name", how="inner", validate="one_to_one")


def circular_taxon_table(
    traits: pd.DataFrame,
    members: list[str],
    taxon_env: pd.DataFrame,
    min_trait_observations: int,
) -> pd.DataFrame:
    sine = next(x for x in members if "sin" in x)
    cosine = next(x for x in members if "cos" in x)
    part = traits[
        traits["endpoint_id"].isin([sine, cosine])
        & traits["analysis_eligible_bool"]
        & traits["value"].notna()
    ][["obs_id", "taxon_name", "endpoint_id", "value"]].copy()
    wide = part.pivot(index=["obs_id", "taxon_name"], columns="endpoint_id", values="value")
    wide = wide.dropna(subset=[sine, cosine]).reset_index()
    counts = wide.groupby("taxon_name").size().rename("n_trait_observations")
    medians = wide.groupby("taxon_name")[[sine, cosine]].median()
    summary = medians.join(counts).reset_index()
    summary = summary[summary["n_trait_observations"].ge(min_trait_observations)]
    return summary.merge(taxon_env, on="taxon_name", how="inner", validate="one_to_one")


def run_scope(
    traits: pd.DataFrame,
    environment: pd.DataFrame,
    predictors: list[str],
    min_trait_observations: int,
    minimum_taxa: int,
    permutations: int,
    seed: int,
    scope: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    taxon_env = taxon_environment(environment, predictors)
    rows: list[dict[str, Any]] = []
    taxon_rows: list[pd.DataFrame] = []
    for unit in inferential_units(traits):
        if unit["inferential_unit"] == "linear_endpoint":
            table = linear_taxon_table(traits, unit["endpoint_id"], taxon_env, min_trait_observations)
        else:
            table = circular_taxon_table(traits, unit["members"], taxon_env, min_trait_observations)
        table = table.copy()
        table.insert(0, "scope", scope)
        table.insert(1, "endpoint_id", unit["endpoint_id"])
        taxon_rows.append(table)
        for predictor in predictors:
            needed = [predictor, "trait_value"] if unit["inferential_unit"] == "linear_endpoint" else [predictor, *unit["members"]]
            complete = table.dropna(subset=needed).copy()
            base = {
                "scope": scope,
                "min_trait_observations_per_taxon": min_trait_observations,
                "endpoint_id": unit["endpoint_id"],
                "module": unit["module"],
                "analysis_tier": unit["analysis_tier"],
                "inferential_unit": unit["inferential_unit"],
                "predictor": predictor,
                "n_taxa": int(len(complete)),
            }
            if len(complete) < minimum_taxa:
                rows.append({**base, "status": "insufficient_taxa"})
                continue
            rng = stable_rng(seed, scope, unit["endpoint_id"], predictor)
            try:
                if unit["inferential_unit"] == "linear_endpoint":
                    beta, p_value = permutation_p_linear(
                        complete["trait_value"].to_numpy(float),
                        complete[predictor].to_numpy(float),
                        permutations,
                        rng,
                    )
                    rows.append({
                        **base,
                        "status": "ok",
                        "p_value_method": f"taxon_label_permutation_two_sided_{permutations}",
                        "beta_std_among": beta,
                        "p_value": p_value,
                    })
                else:
                    sine = next(x for x in unit["members"] if "sin" in x)
                    cosine = next(x for x in unit["members"] if "cos" in x)
                    bs, bc, magnitude, direction, p_value = permutation_p_circular(
                        complete[sine].to_numpy(float),
                        complete[cosine].to_numpy(float),
                        complete[predictor].to_numpy(float),
                        permutations,
                        rng,
                    )
                    rows.append({
                        **base,
                        "status": "ok",
                        "p_value_method": f"taxon_label_permutation_joint_circular_{permutations}",
                        "beta_sine_std_among": bs,
                        "beta_cosine_std_among": bc,
                        "effect_magnitude_among": magnitude,
                        "effect_direction_degrees_among": direction,
                        "p_value": p_value,
                    })
            except Exception as error:
                rows.append({**base, "status": f"failed:{type(error).__name__}:{error}"})
    result = pd.DataFrame(rows)
    result["q_fdr_bh_within_tier_scope"] = np.nan
    ok = result["status"].eq("ok") & pd.to_numeric(result["p_value"], errors="coerce").notna()
    for (_, _tier), index in result[ok].groupby(["scope", "analysis_tier"]).groups.items():
        result.loc[index, "q_fdr_bh_within_tier_scope"] = bh_adjust(result.loc[index, "p_value"].astype(float))
    result["fdr_significant_0_05"] = result["q_fdr_bh_within_tier_scope"].lt(0.05)
    taxon_summary = pd.concat(taxon_rows, ignore_index=True, sort=False) if taxon_rows else pd.DataFrame()
    return result, taxon_summary


def build_alignment(among: pd.DataFrame, within_path: Path, primary_scope: str) -> pd.DataFrame:
    within = pd.read_csv(within_path, low_memory=False)
    within = within[within["status"].eq("ok")].copy()
    keep = [
        "endpoint_id", "module", "analysis_tier", "inferential_unit", "predictor",
        "beta_std_within", "beta_sine_std_within", "beta_cosine_std_within",
        "effect_magnitude", "effect_direction_degrees", "p_value",
        "q_fdr_bh_within_tier", "fdr_significant_0_05",
    ]
    for column in keep:
        if column not in within:
            within[column] = np.nan
    within = within[keep].rename(columns={
        "p_value": "p_value_within",
        "q_fdr_bh_within_tier": "q_fdr_bh_within",
        "fdr_significant_0_05": "within_fdr_significant_0_05",
    })
    selected = among[among["scope"].eq(primary_scope) & among["status"].eq("ok")].copy().rename(columns={
        "p_value": "p_value_among",
        "q_fdr_bh_within_tier_scope": "q_fdr_bh_among",
        "fdr_significant_0_05": "among_fdr_significant_0_05",
    })
    merged = within.merge(
        selected,
        on=["endpoint_id", "module", "analysis_tier", "inferential_unit", "predictor"],
        how="outer",
        validate="one_to_one",
        suffixes=("_within", "_among"),
    )
    merged["within_fdr_significant_0_05"] = merged["within_fdr_significant_0_05"].fillna(False).astype(bool)
    merged["among_fdr_significant_0_05"] = merged["among_fdr_significant_0_05"].fillna(False).astype(bool)
    merged["cross_scale_class"] = np.select(
        [
            merged["within_fdr_significant_0_05"] & merged["among_fdr_significant_0_05"],
            merged["within_fdr_significant_0_05"] & ~merged["among_fdr_significant_0_05"],
            ~merged["within_fdr_significant_0_05"] & merged["among_fdr_significant_0_05"],
        ],
        ["both_scales", "within_only", "among_only"],
        default="neither",
    )
    merged["linear_sign_concordant"] = pd.Series(pd.NA, index=merged.index, dtype="boolean")
    linear = merged["inferential_unit"].eq("linear_endpoint")
    bw = pd.to_numeric(merged["beta_std_within"], errors="coerce")
    ba = pd.to_numeric(merged["beta_std_among"], errors="coerce")
    comparable = linear & bw.notna() & ba.notna() & (bw != 0) & (ba != 0)
    concordant = np.sign(bw[comparable]).eq(np.sign(ba[comparable]))
    merged.loc[comparable, "linear_sign_concordant"] = concordant.to_numpy(dtype=bool)
    return merged


def main() -> int:
    args = parse_args()
    if min(args.primary_min_trait_observations, args.sensitivity_min_trait_observations) < 1:
        raise ValueError("Trait-observation thresholds must be positive")
    if args.permutations < 99:
        raise ValueError("At least 99 permutations are required")
    traits, environment = load_inputs(args.traits_long, args.environment, args.predictors)
    primary_scope = f"among_taxon_min{args.primary_min_trait_observations}"
    sensitivity_scope = f"among_taxon_min{args.sensitivity_min_trait_observations}"
    primary, primary_taxa = run_scope(
        traits, environment, args.predictors, args.primary_min_trait_observations,
        args.minimum_taxa, args.permutations, args.seed, primary_scope,
    )
    frames = [primary]
    taxa_frames = [primary_taxa]
    if sensitivity_scope != primary_scope:
        sensitivity, sensitivity_taxa = run_scope(
            traits, environment, args.predictors, args.sensitivity_min_trait_observations,
            args.minimum_taxa, args.permutations, args.seed, sensitivity_scope,
        )
        frames.append(sensitivity)
        taxa_frames.append(sensitivity_taxa)
    among = pd.concat(frames, ignore_index=True, sort=False)
    taxon_values = pd.concat(taxa_frames, ignore_index=True, sort=False)
    alignment = build_alignment(among, args.within_coefficients, primary_scope)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    among.to_csv(args.out_dir / "geb_v2_among_taxon_chelsa_coefficients.csv", index=False)
    alignment.to_csv(args.out_dir / "geb_v2_within_among_alignment.csv", index=False)
    taxon_values.to_csv(args.out_dir / "geb_v2_among_taxon_trait_summaries.csv", index=False)

    primary_ok = among[among["scope"].eq(primary_scope) & among["status"].eq("ok")]
    report = {
        "primary_scope": primary_scope,
        "sensitivity_scope": sensitivity_scope if sensitivity_scope != primary_scope else None,
        "predictors": args.predictors,
        "permutations": args.permutations,
        "seed": args.seed,
        "minimum_taxa": args.minimum_taxa,
        "n_inferential_units": int(primary_ok["endpoint_id"].nunique()),
        "n_primary_tier_units": int(primary_ok.loc[primary_ok["analysis_tier"].eq("primary"), "endpoint_id"].nunique()),
        "n_candidate_tier_units": int(primary_ok.loc[primary_ok["analysis_tier"].eq("candidate"), "endpoint_id"].nunique()),
        "n_primary_scope_models": int(len(primary_ok)),
        "n_primary_scope_fdr_signals": int(primary_ok["fdr_significant_0_05"].sum()),
        "n_cross_scale_both": int(alignment["cross_scale_class"].eq("both_scales").sum()),
        "n_cross_scale_within_only": int(alignment["cross_scale_class"].eq("within_only").sum()),
        "n_cross_scale_among_only": int(alignment["cross_scale_class"].eq("among_only").sum()),
        "multiplicity": "BH separately within primary and candidate tiers for each among-taxon replication scope; circular hue is one joint test per predictor",
        "taxon_environment_summary": "median CHELSA across the frozen strict spatial observations for each source-assigned taxon",
        "taxon_trait_summary": "median endpoint value among analysis-eligible observations; primary scope requires the frozen minimum usable observations per taxon",
        "interpretation_boundary": "Cross-scale observational pattern comparison. Among-taxon taxon-label permutations do not remove phylogenetic non-independence; within-taxon and among-taxon support do not establish plasticity, adaptation, selection, or mechanism."
    }
    (args.out_dir / "geb_v2_within_among_alignment_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())