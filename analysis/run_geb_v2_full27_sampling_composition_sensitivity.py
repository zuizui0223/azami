#!/usr/bin/env python3
"""Audit sampling-composition stability for the selected full-27 atlas rows.

This is a retrospective, selected-row sensitivity analysis.  It does not
create a new discovery family or replace the complete atlas.  Linear rows are
checked for sign preservation; circular rows are checked for positive vector
alignment (an angular difference below 90 degrees).  No post-selection
p-values are produced.
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-contract", required=True, type=Path)
    parser.add_argument("--sensitivity-contract", required=True, type=Path)
    parser.add_argument("--traits-long", required=True, type=Path)
    parser.add_argument("--environment", required=True, type=Path)
    parser.add_argument("--atlas-within", required=True, type=Path)
    parser.add_argument("--atlas-among", required=True, type=Path)
    parser.add_argument("--regions", required=True, type=Path)
    parser.add_argument("--native-status", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    return parser.parse_args()


def truth(values: pd.Series) -> pd.Series:
    return values.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stable_top_taxa(environment: pd.DataFrame, n: int) -> list[str]:
    """Return the most observed taxa with an explicit alphabetical tie-break."""
    counts = (
        environment.groupby("taxon_name", as_index=False)
        .size()
        .rename(columns={"size": "n_observations"})
        .sort_values(["n_observations", "taxon_name"], ascending=[False, True])
    )
    return counts.head(n)["taxon_name"].astype(str).tolist()


def equal_taxon_weights(taxa: pd.Series) -> np.ndarray:
    counts = taxa.groupby(taxa).transform("size").to_numpy(float)
    if np.any(counts <= 0):
        raise ValueError("Taxon counts must be positive")
    return 1.0 / counts


def standardize(values: np.ndarray, weights: np.ndarray | None = None) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if weights is None:
        mean = float(np.mean(values))
        variance = float(np.mean((values - mean) ** 2))
    else:
        weights = np.asarray(weights, dtype=float)
        mean = float(np.average(values, weights=weights))
        variance = float(np.average((values - mean) ** 2, weights=weights))
    if not np.isfinite(variance) or variance <= 0:
        raise ValueError("No finite variation")
    return (values - mean) / math.sqrt(variance)


def demean_standardize(
    frame: pd.DataFrame,
    column: str,
    weights: np.ndarray | None = None,
) -> np.ndarray:
    values = pd.to_numeric(frame[column], errors="coerce")
    centred = values - values.groupby(frame["taxon_name"]).transform("mean")
    return standardize(centred.to_numpy(float), weights)


def weighted_slope(y: np.ndarray, x: np.ndarray, weights: np.ndarray | None) -> float:
    if weights is None:
        weights = np.ones(len(x), dtype=float)
    denominator = float(np.dot(weights * x, x))
    if not np.isfinite(denominator) or denominator <= 0:
        raise ValueError("Predictor has no finite weighted variation")
    return float(np.dot(weights * x, y) / denominator)


def circular_alignment(
    baseline_sine: float,
    baseline_cosine: float,
    estimate_sine: float,
    estimate_cosine: float,
) -> float:
    baseline = np.asarray([baseline_sine, baseline_cosine], dtype=float)
    estimate = np.asarray([estimate_sine, estimate_cosine], dtype=float)
    denominator = float(np.linalg.norm(baseline) * np.linalg.norm(estimate))
    if not np.isfinite(denominator) or denominator <= 0:
        return float("nan")
    return float(np.dot(baseline, estimate) / denominator)


def load_inputs(args: argparse.Namespace) -> tuple[
    dict[str, Any], pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame
]:
    analysis_contract = json.loads(args.analysis_contract.read_text(encoding="utf-8"))
    predictors = list(analysis_contract["environment"]["predictors"])
    trait_columns = [
        "obs_id",
        "taxon_name",
        "endpoint_id",
        "measurement_available",
        "value",
    ]
    traits = pd.read_csv(args.traits_long, usecols=trait_columns, low_memory=False)
    environment = pd.read_csv(
        args.environment,
        usecols=["obs_id", "taxon_name", "latitude", "longitude", *predictors],
        low_memory=False,
    )
    regions = pd.read_csv(args.regions, usecols=["obs_id", "broad_region"], low_memory=False)
    native = pd.read_csv(
        args.native_status,
        usecols=["obs_id", "taxon_name", "native_range_status"],
        low_memory=False,
    )
    for frame in (traits, environment, regions, native):
        frame["obs_id"] = frame["obs_id"].astype(str)
    if environment["obs_id"].duplicated().any():
        raise ValueError("Environment must be unique by obs_id")
    if regions["obs_id"].duplicated().any():
        raise ValueError("Region lookup must be unique by obs_id")
    if native["obs_id"].duplicated().any():
        raise ValueError("Native-status lookup must be unique by obs_id")
    expected_ids = set(environment["obs_id"])
    if set(regions["obs_id"]) != expected_ids:
        raise ValueError("Region lookup must cover the exact environment observation universe")
    if set(native["obs_id"]) != expected_ids:
        raise ValueError("Native-status lookup must cover the exact environment observation universe")
    if traits.duplicated(["obs_id", "endpoint_id"]).any():
        raise ValueError("Trait table must be unique by obs_id and endpoint_id")
    available = truth(traits["measurement_available"])
    traits["value"] = pd.to_numeric(traits["value"], errors="coerce")
    traits = traits[available & traits["value"].notna()].copy()
    metadata = environment.merge(regions, on="obs_id", validate="one_to_one")
    metadata = metadata.merge(
        native.rename(columns={"taxon_name": "native_taxon_name"}),
        on="obs_id",
        validate="one_to_one",
    )
    if not metadata["taxon_name"].astype(str).eq(metadata["native_taxon_name"].astype(str)).all():
        raise ValueError("Native-status taxon names disagree with the environment table")
    metadata = metadata.drop(columns="native_taxon_name")
    for predictor in predictors:
        metadata[predictor] = pd.to_numeric(metadata[predictor], errors="coerce")
    return analysis_contract, traits, environment, metadata, native


def unit_observations(
    traits: pd.DataFrame,
    metadata: pd.DataFrame,
    member_endpoint_ids: list[str],
) -> pd.DataFrame:
    part = traits[traits["endpoint_id"].isin(member_endpoint_ids)][
        ["obs_id", "taxon_name", "endpoint_id", "value"]
    ]
    if part.empty or set(part["endpoint_id"]) != set(member_endpoint_ids):
        return pd.DataFrame()
    wide = part.pivot(
        index=["obs_id", "taxon_name"], columns="endpoint_id", values="value"
    ).reset_index()
    wide.columns.name = None
    return wide.merge(
        metadata.drop(columns="taxon_name"),
        on="obs_id",
        how="inner",
        validate="one_to_one",
    )


def fit_within(
    data: pd.DataFrame,
    members: list[str],
    predictor: str,
    inferential_unit: str,
    minimum_observations: int,
    minimum_taxa: int,
    equal_weight: bool = False,
) -> dict[str, Any]:
    complete = data.dropna(subset=[*members, predictor]).copy()
    counts = complete.groupby("taxon_name").size()
    complete = complete[complete["taxon_name"].isin(counts[counts >= 2].index)].copy()
    n_observations = int(len(complete))
    n_taxa = int(complete["taxon_name"].nunique())
    base = {"n_observations": n_observations, "n_taxa": n_taxa}
    if n_observations < minimum_observations or n_taxa < minimum_taxa:
        return {**base, "status": "insufficient_support"}
    weights = equal_taxon_weights(complete["taxon_name"]) if equal_weight else None
    try:
        x = demean_standardize(complete, predictor, weights)
        if inferential_unit == "linear_endpoint":
            y = demean_standardize(complete, members[0], weights)
            beta = weighted_slope(y, x, weights)
            return {
                **base,
                "status": "ok",
                "beta_std_sensitivity": beta,
                "effect_magnitude_sensitivity": abs(beta),
            }
        sine = next(member for member in members if "sin" in member)
        cosine = next(member for member in members if "cos" in member)
        ys = demean_standardize(complete, sine, weights)
        yc = demean_standardize(complete, cosine, weights)
        beta_sine = weighted_slope(ys, x, weights)
        beta_cosine = weighted_slope(yc, x, weights)
        return {
            **base,
            "status": "ok",
            "beta_sine_sensitivity": beta_sine,
            "beta_cosine_sensitivity": beta_cosine,
            "effect_magnitude_sensitivity": float(math.hypot(beta_sine, beta_cosine)),
            "effect_direction_degrees_sensitivity": float(
                math.degrees(math.atan2(beta_sine, beta_cosine)) % 360.0
            ),
        }
    except ValueError as error:
        return {**base, "status": f"failed:{type(error).__name__}:{error}"}


def fit_among(
    trait_data: pd.DataFrame,
    metadata: pd.DataFrame,
    members: list[str],
    predictor: str,
    inferential_unit: str,
    minimum_trait_observations: int,
    minimum_taxa: int,
) -> dict[str, Any]:
    complete_traits = trait_data.dropna(subset=members).copy()
    counts = complete_traits.groupby("taxon_name").size().rename("n_trait_observations")
    medians = complete_traits.groupby("taxon_name")[members].median().join(counts)
    medians = medians[medians["n_trait_observations"].ge(minimum_trait_observations)]
    env_medians = metadata.groupby("taxon_name")[[predictor]].median(numeric_only=True)
    data = medians.join(env_medians, how="inner").dropna(subset=[*members, predictor])
    n_taxa = int(len(data))
    base = {
        "n_observations": int(
            complete_traits[complete_traits["taxon_name"].isin(data.index)].shape[0]
        ),
        "n_taxa": n_taxa,
    }
    if n_taxa < minimum_taxa:
        return {**base, "status": "insufficient_support"}
    try:
        x = standardize(data[predictor].to_numpy(float))
        if inferential_unit == "linear_endpoint":
            y = standardize(data[members[0]].to_numpy(float))
            beta = weighted_slope(y, x, None)
            return {
                **base,
                "status": "ok",
                "beta_std_sensitivity": beta,
                "effect_magnitude_sensitivity": abs(beta),
            }
        sine = next(member for member in members if "sin" in member)
        cosine = next(member for member in members if "cos" in member)
        beta_sine = weighted_slope(standardize(data[sine].to_numpy(float)), x, None)
        beta_cosine = weighted_slope(standardize(data[cosine].to_numpy(float)), x, None)
        return {
            **base,
            "status": "ok",
            "beta_sine_sensitivity": beta_sine,
            "beta_cosine_sensitivity": beta_cosine,
            "effect_magnitude_sensitivity": float(math.hypot(beta_sine, beta_cosine)),
            "effect_direction_degrees_sensitivity": float(
                math.degrees(math.atan2(beta_sine, beta_cosine)) % 360.0
            ),
        }
    except ValueError as error:
        return {**base, "status": f"failed:{type(error).__name__}:{error}"}


def baseline_effect(row: pd.Series) -> dict[str, float]:
    if str(row["inferential_unit"]) == "linear_endpoint":
        beta = float(row["beta_std"])
        return {
            "baseline_beta_std": beta,
            "baseline_effect_magnitude": abs(beta),
        }
    sine_column = "beta_sine_std" if "beta_sine_std" in row.index else "beta_sine_std_among"
    cosine_column = "beta_cosine_std" if "beta_cosine_std" in row.index else "beta_cosine_std_among"
    beta_sine = float(row[sine_column])
    beta_cosine = float(row[cosine_column])
    return {
        "baseline_beta_sine": beta_sine,
        "baseline_beta_cosine": beta_cosine,
        "baseline_effect_magnitude": float(math.hypot(beta_sine, beta_cosine)),
        "baseline_effect_direction_degrees": float(
            math.degrees(math.atan2(beta_sine, beta_cosine)) % 360.0
        ),
    }


def add_stability(row: dict[str, Any], inferential_unit: str) -> dict[str, Any]:
    result = dict(row)
    if result.get("status") != "ok":
        result["direction_stable"] = False
        return result
    baseline_magnitude = float(result["baseline_effect_magnitude"])
    sensitivity_magnitude = float(result["effect_magnitude_sensitivity"])
    result["effect_magnitude_ratio_to_baseline"] = (
        sensitivity_magnitude / baseline_magnitude if baseline_magnitude > 0 else np.nan
    )
    if inferential_unit == "linear_endpoint":
        baseline = float(result["baseline_beta_std"])
        estimate = float(result["beta_std_sensitivity"])
        result["direction_stable"] = bool(
            baseline != 0 and estimate != 0 and np.sign(baseline) == np.sign(estimate)
        )
        return result
    alignment = circular_alignment(
        float(result["baseline_beta_sine"]),
        float(result["baseline_beta_cosine"]),
        float(result["beta_sine_sensitivity"]),
        float(result["beta_cosine_sensitivity"]),
    )
    result["circular_vector_alignment_cosine"] = alignment
    result["direction_stable"] = bool(np.isfinite(alignment) and alignment > 0)
    return result


def scenario_definitions(top_taxa: list[str], regions: list[str], scale: str) -> list[dict[str, Any]]:
    scenarios: list[dict[str, Any]] = []
    if scale == "within_taxon":
        scenarios.append(
            {
                "sensitivity_family": "equal_taxon_weight",
                "scenario": "equal_total_weight_per_taxon",
                "excluded_taxa": [],
                "excluded_region": "",
                "native_only": False,
                "equal_taxon_weight": True,
            }
        )
    for taxon in top_taxa:
        scenarios.append(
            {
                "sensitivity_family": "dominant_taxon_omission",
                "scenario": f"omit:{taxon}",
                "excluded_taxa": [taxon],
                "excluded_region": "",
                "native_only": False,
                "equal_taxon_weight": False,
            }
        )
    scenarios.append(
        {
            "sensitivity_family": "dominant_taxon_omission",
            "scenario": "omit_top2_joint",
            "excluded_taxa": top_taxa[:2],
            "excluded_region": "",
            "native_only": False,
            "equal_taxon_weight": False,
        }
    )
    for region in regions:
        scenarios.append(
            {
                "sensitivity_family": "leave_one_broad_region_out",
                "scenario": f"omit_region:{region}",
                "excluded_taxa": [],
                "excluded_region": region,
                "native_only": False,
                "equal_taxon_weight": False,
            }
        )
    scenarios.append(
        {
            "sensitivity_family": "native_only",
            "scenario": "native_only",
            "excluded_taxa": [],
            "excluded_region": "",
            "native_only": True,
            "equal_taxon_weight": False,
        }
    )
    return scenarios


def apply_scenario(
    unit_data: pd.DataFrame,
    metadata: pd.DataFrame,
    scenario: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    trait_scope = unit_data.copy()
    metadata_scope = metadata.copy()
    excluded_taxa = set(scenario["excluded_taxa"])
    if excluded_taxa:
        trait_scope = trait_scope[~trait_scope["taxon_name"].isin(excluded_taxa)]
        metadata_scope = metadata_scope[~metadata_scope["taxon_name"].isin(excluded_taxa)]
    excluded_region = str(scenario["excluded_region"])
    if excluded_region:
        trait_scope = trait_scope[~trait_scope["broad_region"].eq(excluded_region)]
        metadata_scope = metadata_scope[~metadata_scope["broad_region"].eq(excluded_region)]
    if bool(scenario["native_only"]):
        trait_scope = trait_scope[trait_scope["native_range_status"].eq("native")]
        metadata_scope = metadata_scope[metadata_scope["native_range_status"].eq("native")]
    return trait_scope.copy(), metadata_scope.copy()


def summarize_scenarios(scenarios: pd.DataFrame) -> pd.DataFrame:
    keys = ["scale", "unit_id", "predictor"]
    rows: list[dict[str, Any]] = []
    for key, part in scenarios.groupby(keys, sort=True):
        row: dict[str, Any] = dict(zip(keys, key))
        row["inferential_unit"] = str(part.iloc[0]["inferential_unit"])
        row["baseline_effect_magnitude"] = float(part.iloc[0]["baseline_effect_magnitude"])
        all_evaluable = True
        all_stable = True
        any_unstable = False
        for family, family_part in part.groupby("sensitivity_family", sort=True):
            ok = family_part["status"].eq("ok")
            stable = family_part.loc[ok, "direction_stable"].astype(bool)
            row[f"{family}_n_scenarios"] = int(len(family_part))
            row[f"{family}_n_evaluable"] = int(ok.sum())
            row[f"{family}_all_evaluable"] = bool(ok.all())
            row[f"{family}_all_directions_stable_where_evaluable"] = bool(
                len(stable) > 0 and stable.all()
            )
            row[f"{family}_minimum_effect_magnitude_ratio"] = (
                float(
                    pd.to_numeric(
                        family_part.loc[ok, "effect_magnitude_ratio_to_baseline"],
                        errors="coerce",
                    ).min()
                )
                if ok.any()
                else np.nan
            )
            all_evaluable = all_evaluable and bool(ok.all())
            all_stable = all_stable and bool(len(stable) > 0 and stable.all())
            any_unstable = any_unstable or bool((~stable).any())
        row["all_declared_scenarios_evaluable"] = all_evaluable
        row["all_directions_stable_where_evaluable"] = all_stable
        if any_unstable:
            row["sampling_composition_stability_class"] = "direction_unstable"
        elif all_evaluable and all_stable:
            row["sampling_composition_stability_class"] = "stable_all_declared_scenarios"
        else:
            row["sampling_composition_stability_class"] = "stable_where_evaluable_incomplete"
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> int:
    args = parse_args()
    analysis_contract, traits, environment, metadata, _ = load_inputs(args)
    within = pd.read_csv(args.atlas_within, low_memory=False)
    among = pd.read_csv(args.atlas_among, low_memory=False)
    selected_within = within[
        within["status"].eq("ok")
        & pd.to_numeric(within["q_fdr_bh_global_family"], errors="coerce").lt(0.05)
    ].copy()
    selected_among = among[
        among["scope"].eq("among_taxon_min5")
        & among["status"].eq("ok")
        & pd.to_numeric(among["q_fdr_bh_global_family"], errors="coerce").lt(0.05)
    ].copy()
    rules = json.loads(args.sensitivity_contract.read_text(encoding="utf-8"))
    if rules["status"] != "retrospective_selected_row_sensitivity":
        raise ValueError("Sampling-composition sensitivity must remain explicitly retrospective")
    top_taxa = stable_top_taxa(environment, int(rules["dominant_taxa_top_n"]))
    regions = sorted(
        region
        for region in metadata["broad_region"].dropna().astype(str).unique()
        if region != "UNMAPPED"
    )
    within_rules = analysis_contract["within_taxon"]
    among_rules = analysis_contract["among_taxon"]
    output_rows: list[dict[str, Any]] = []

    for scale, selected in (("within_taxon", selected_within), ("among_taxon", selected_among)):
        for _, atlas_row in selected.iterrows():
            members = str(atlas_row["member_endpoint_ids"]).split("|")
            unit_data = unit_observations(traits, metadata, members)
            if unit_data.empty:
                raise ValueError(f"No measurements for selected unit {atlas_row['unit_id']}")
            baseline = baseline_effect(atlas_row)
            for scenario in scenario_definitions(top_taxa, regions, scale):
                trait_scope, metadata_scope = apply_scenario(unit_data, metadata, scenario)
                if scale == "within_taxon":
                    fit = fit_within(
                        trait_scope,
                        members,
                        str(atlas_row["predictor"]),
                        str(atlas_row["inferential_unit"]),
                        int(within_rules["minimum_observations"]),
                        int(within_rules["minimum_taxa"]),
                        bool(scenario["equal_taxon_weight"]),
                    )
                else:
                    fit = fit_among(
                        trait_scope,
                        metadata_scope,
                        members,
                        str(atlas_row["predictor"]),
                        str(atlas_row["inferential_unit"]),
                        int(among_rules["primary_minimum_trait_observations_per_taxon"]),
                        int(among_rules["minimum_taxa"]),
                    )
                result = {
                    "scale": scale,
                    "unit_id": str(atlas_row["unit_id"]),
                    "member_endpoint_ids": str(atlas_row["member_endpoint_ids"]),
                    "module": str(atlas_row["module"]),
                    "inferential_unit": str(atlas_row["inferential_unit"]),
                    "predictor": str(atlas_row["predictor"]),
                    "environment_block": str(atlas_row["environment_block"]),
                    "sensitivity_family": scenario["sensitivity_family"],
                    "scenario": scenario["scenario"],
                    "excluded_taxa": "|".join(scenario["excluded_taxa"]),
                    "excluded_region": scenario["excluded_region"],
                    "native_only": scenario["native_only"],
                    "equal_taxon_weight": scenario["equal_taxon_weight"],
                    **baseline,
                    **fit,
                }
                output_rows.append(add_stability(result, str(atlas_row["inferential_unit"])))

    scenarios = pd.DataFrame(output_rows)
    summary = summarize_scenarios(scenarios)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    scenarios.to_csv(
        args.out_dir / "v2_full27_sampling_composition_scenarios.csv",
        index=False,
        encoding="utf-8-sig",
    )
    summary.to_csv(
        args.out_dir / "v2_full27_sampling_composition_summary.csv",
        index=False,
        encoding="utf-8-sig",
    )
    status_counts = summary["sampling_composition_stability_class"].value_counts().to_dict()
    report = {
        "analysis_id": "geb_v2_full27_sampling_composition_sensitivity_v1",
        "status": rules["status"],
        "entry_rule": rules["entry_rule"],
        "n_selected_within_pairs": int(len(selected_within)),
        "n_selected_among_pairs": int(len(selected_among)),
        "n_scenario_rows": int(len(scenarios)),
        "n_broad_regions_omitted_separately": int(len(regions)),
        "broad_regions": regions,
        "top_taxa": top_taxa,
        "stability_class_counts": status_counts,
        "n_direction_unstable_pairs": int(
            summary["sampling_composition_stability_class"].eq("direction_unstable").sum()
        ),
        "method": rules,
        "claim_boundary": (
            "Retrospective directional stability under declared sampling-composition "
            "perturbations; not a correction to a probability sample, not a new discovery "
            "family and not evidence of causal environmental response."
        ),
        "input_sha256": {
            "analysis_contract": sha256_file(args.analysis_contract),
            "sensitivity_contract": sha256_file(args.sensitivity_contract),
            "traits_long": sha256_file(args.traits_long),
            "environment": sha256_file(args.environment),
            "atlas_within": sha256_file(args.atlas_within),
            "atlas_among": sha256_file(args.atlas_among),
            "regions": sha256_file(args.regions),
            "native_status": sha256_file(args.native_status),
        },
    }
    (args.out_dir / "v2_full27_sampling_composition_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
