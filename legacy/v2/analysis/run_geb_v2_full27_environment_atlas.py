#!/usr/bin/env python3
"""Run the all-registered-endpoint by all-environment GEB v2 atlas.

This lane starts from every row in the frozen 27-endpoint GEB v2 contract.  It
does not use the old nine-endpoint result set, analysis tier, significance, or
measurement completeness to select endpoints.  Unexecuted endpoints remain in
the output as explicit no-measurement rows.  All measured endpoints, including
validation-only image proxies, are modelled but retain their validation labels.

The analysis is retrospective and exploratory.  It maps scale-specific
associations; it does not remove spatial or phylogenetic non-independence and it
does not identify causal environmental mechanisms.
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
import statsmodels.api as sm


EXPECTED_PREDICTORS = [
    "chelsa_bio01",
    "chelsa_bio04",
    "chelsa_bio12",
    "chelsa_bio15",
    "chelsa_rsds_mean",
    "chelsa_vpd_mean",
    "chelsa_sfcwind_mean",
    "chelsa_gsp",
    "chelsa_npp",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contract", required=True, type=Path)
    parser.add_argument("--analysis-contract", required=True, type=Path)
    parser.add_argument("--traits-long", required=True, type=Path)
    parser.add_argument("--environment", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    return parser.parse_args()


def as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stable_rng(seed: int, *parts: str) -> np.random.Generator:
    payload = "|".join([str(seed), *map(str, parts)]).encode("utf-8")
    digest = hashlib.sha256(payload).digest()
    return np.random.default_rng(int.from_bytes(digest[:8], "little", signed=False))


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


def standardize(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    sd = float(np.std(values, ddof=0))
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError("No finite variation")
    return (values - float(np.mean(values))) / sd


def demean_standardize(frame: pd.DataFrame, column: str) -> np.ndarray:
    values = pd.to_numeric(frame[column], errors="coerce")
    centred = values - values.groupby(frame["taxon_name"]).transform("mean")
    sd = float(centred.std(ddof=0))
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError(f"No within-taxon variation in {column}")
    return (centred / sd).to_numpy(float)


def slope(y: np.ndarray, x: np.ndarray) -> float:
    denominator = float(np.dot(x, x))
    if denominator <= 0:
        raise ValueError("Predictor has zero sum of squares")
    return float(np.dot(x, y) / denominator)


def within_taxon_permute(values: np.ndarray, taxa: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    result = values.copy()
    for taxon in np.unique(taxa):
        index = np.flatnonzero(taxa == taxon)
        result[index] = rng.permutation(result[index])
    return result


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
    remaining = permutations
    batch_size = 256
    while remaining:
        size = min(batch_size, remaining)
        matrix = np.broadcast_to(xz, (size, len(xz))).copy()
        permuted = rng.permuted(matrix, axis=1)
        simulated = permuted @ yz / float(np.dot(xz, xz))
        exceed += int(np.sum(np.abs(simulated) >= abs(observed) - 1e-15))
        remaining -= size
    return observed, float((exceed + 1) / (permutations + 1))


def permutation_p_circular_among(
    sine: np.ndarray,
    cosine: np.ndarray,
    x: np.ndarray,
    permutations: int,
    rng: np.random.Generator,
) -> dict[str, float]:
    ys = standardize(sine)
    yc = standardize(cosine)
    xz = standardize(x)
    beta_sine = slope(ys, xz)
    beta_cosine = slope(yc, xz)
    magnitude = float(math.hypot(beta_sine, beta_cosine))
    exceed = 0
    remaining = permutations
    batch_size = 256
    denominator = float(np.dot(xz, xz))
    while remaining:
        size = min(batch_size, remaining)
        matrix = np.broadcast_to(xz, (size, len(xz))).copy()
        permuted = rng.permuted(matrix, axis=1)
        bs = permuted @ ys / denominator
        bc = permuted @ yc / denominator
        exceed += int(np.sum(np.hypot(bs, bc) >= magnitude - 1e-15))
        remaining -= size
    return {
        "beta_sine_std_among": beta_sine,
        "beta_cosine_std_among": beta_cosine,
        "effect_magnitude": magnitude,
        "effect_direction_degrees": float(math.degrees(math.atan2(beta_sine, beta_cosine)) % 360.0),
        "p_value": float((exceed + 1) / (permutations + 1)),
    }


def permutation_p_circular_within(
    sine: np.ndarray,
    cosine: np.ndarray,
    x: np.ndarray,
    taxa: np.ndarray,
    permutations: int,
    rng: np.random.Generator,
) -> dict[str, float]:
    """Joint circular test with predictor permutations restricted within taxa.

    The permutation matrix is accumulated taxon by taxon in batches.  This is
    the same randomisation scheme as repeatedly shuffling ``x`` within every
    taxon, but avoids millions of small pandas/NumPy indexing operations for
    the full observation universe.
    """
    beta_sine = slope(sine, x)
    beta_cosine = slope(cosine, x)
    magnitude = float(math.hypot(beta_sine, beta_cosine))
    denominator = float(np.dot(x, x))
    groups = [np.flatnonzero(taxa == taxon) for taxon in np.unique(taxa)]
    exceed = 0
    remaining = permutations
    batch_size = 256
    while remaining:
        size = min(batch_size, remaining)
        numerator_sine = np.zeros(size, dtype=float)
        numerator_cosine = np.zeros(size, dtype=float)
        for index in groups:
            matrix = np.broadcast_to(x[index], (size, len(index))).copy()
            permuted = rng.permuted(matrix, axis=1)
            numerator_sine += permuted @ sine[index]
            numerator_cosine += permuted @ cosine[index]
        simulated = np.hypot(
            numerator_sine / denominator,
            numerator_cosine / denominator,
        )
        exceed += int(np.sum(simulated >= magnitude - 1e-15))
        remaining -= size
    return {
        "beta_sine_std": beta_sine,
        "beta_cosine_std": beta_cosine,
        "effect_magnitude": magnitude,
        "effect_direction_degrees": float(
            math.degrees(math.atan2(beta_sine, beta_cosine)) % 360.0
        ),
        "p_value": float((exceed + 1) / (permutations + 1)),
        "p_value_method": f"within_taxon_predictor_permutation_{permutations}",
    }


def load_contracts(contract_path: Path, analysis_contract_path: Path) -> tuple[pd.DataFrame, dict[str, Any]]:
    contract = pd.read_csv(contract_path, dtype=str, keep_default_na=False)
    analysis_contract = json.loads(analysis_contract_path.read_text(encoding="utf-8"))
    if contract["endpoint_id"].duplicated().any():
        raise ValueError("Endpoint contract must be unique by endpoint_id")
    expected_n = int(analysis_contract["trait_universe"]["expected_registered_endpoints"])
    if len(contract) != expected_n:
        raise ValueError(f"Expected {expected_n} registered endpoints, found {len(contract)}")
    predictors = analysis_contract["environment"]["predictors"]
    if predictors != EXPECTED_PREDICTORS:
        raise ValueError(f"Analysis contract predictors differ from the frozen nine: {predictors}")
    if analysis_contract["trait_universe"].get("old_v1_nine_endpoint_subset_used") is not False:
        raise ValueError("The v1 endpoint subset must not be used")
    return contract, analysis_contract


def load_data(
    traits_path: Path,
    environment_path: Path,
    contract: pd.DataFrame,
    analysis_contract: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    trait_columns = [
        "obs_id",
        "taxon_name",
        "endpoint_id",
        "module",
        "analysis_tier",
        "measurement_available",
        "value",
    ]
    traits = pd.read_csv(traits_path, usecols=trait_columns, low_memory=False)
    environment_columns = ["obs_id", "taxon_name", "latitude", "longitude", *EXPECTED_PREDICTORS]
    environment = pd.read_csv(environment_path, usecols=environment_columns, low_memory=False)
    traits["obs_id"] = traits["obs_id"].astype(str)
    environment["obs_id"] = environment["obs_id"].astype(str)
    if traits.duplicated(["obs_id", "endpoint_id"]).any():
        raise ValueError("Trait table must be unique by obs_id/endpoint_id")
    if environment["obs_id"].duplicated().any():
        raise ValueError("Environment table must be unique by obs_id")
    unknown = sorted(set(traits["endpoint_id"]) - set(contract["endpoint_id"]))
    if unknown:
        raise ValueError(f"Trait table contains endpoints outside the 27-row contract: {unknown}")
    expected_environment_rows = int(analysis_contract["environment"]["expected_observations"])
    if len(environment) != expected_environment_rows:
        raise ValueError(
            f"Expected {expected_environment_rows} full-cohort environment rows, found {len(environment)}"
        )
    traits["measurement_available_bool"] = as_bool(traits["measurement_available"])
    traits["value"] = pd.to_numeric(traits["value"], errors="coerce")
    traits = traits[traits["measurement_available_bool"] & traits["value"].notna()].copy()
    for predictor in EXPECTED_PREDICTORS:
        environment[predictor] = pd.to_numeric(environment[predictor], errors="coerce")
        coverage = float(environment[predictor].notna().mean())
        if coverage < 0.98:
            raise ValueError(f"{predictor} coverage {coverage:.4f} is below 0.98")
    return traits, environment


def environment_block_map(analysis_contract: dict[str, Any]) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for block in analysis_contract["environment"]["blocks"]:
        for predictor in block["predictors"]:
            if predictor in mapping:
                raise ValueError(f"Predictor appears in multiple blocks: {predictor}")
            mapping[predictor] = block["block_id"]
    if list(mapping) != EXPECTED_PREDICTORS:
        raise ValueError("Environment blocks do not cover the exact frozen predictor order")
    return mapping


def inferential_units(contract: pd.DataFrame) -> list[dict[str, Any]]:
    work = contract.copy()
    work["circular_group"] = work["circular_group"].fillna("").astype(str).str.strip()
    units: list[dict[str, Any]] = []
    for row in work[work["circular_group"].eq("")].to_dict("records"):
        units.append(
            {
                "unit_id": row["endpoint_id"],
                "member_endpoint_ids": [row["endpoint_id"]],
                "module": row["module"],
                "analysis_tier": row["analysis_tier"],
                "validation_status": row["validation_status"],
                "inferential_unit": "linear_endpoint",
            }
        )
    for group, part in work[work["circular_group"].ne("")].groupby("circular_group", sort=True):
        members = part["endpoint_id"].tolist()
        if len(members) != 2 or not any("sin" in x for x in members) or not any("cos" in x for x in members):
            raise ValueError(f"Circular group {group!r} must contain one sine and one cosine endpoint")
        units.append(
            {
                "unit_id": group,
                "member_endpoint_ids": members,
                "module": part.iloc[0]["module"],
                "analysis_tier": part.iloc[0]["analysis_tier"],
                "validation_status": part.iloc[0]["validation_status"],
                "inferential_unit": "circular_joint",
            }
        )
    return sorted(units, key=lambda row: (row["module"], row["unit_id"]))


def build_inventory(contract: pd.DataFrame, traits: pd.DataFrame, environment: pd.DataFrame) -> pd.DataFrame:
    measured = (
        traits.groupby("endpoint_id", as_index=False)
        .agg(
            n_measurements_available=("value", "size"),
            n_taxa_measured=("taxon_name", "nunique"),
            n_observations_measured=("obs_id", "nunique"),
        )
    )
    env_ids = set(environment["obs_id"])
    matched = traits[traits["obs_id"].isin(env_ids)]
    matched_counts = (
        matched.groupby("endpoint_id", as_index=False)
        .agg(
            n_measurements_joined_full_environment=("value", "size"),
            n_taxa_joined_full_environment=("taxon_name", "nunique"),
        )
    )
    inventory = contract[
        [
            "endpoint_id",
            "module",
            "biological_construct",
            "unit",
            "analysis_tier",
            "validation_status",
            "circular_group",
            "interpretation",
            "not_interpretable_as",
        ]
    ].copy()
    inventory = inventory.merge(measured, on="endpoint_id", how="left", validate="one_to_one")
    inventory = inventory.merge(matched_counts, on="endpoint_id", how="left", validate="one_to_one")
    count_columns = [column for column in inventory if column.startswith("n_")]
    inventory[count_columns] = inventory[count_columns].fillna(0).astype(int)
    inventory["registered_in_v2_starting_27"] = True
    inventory["measurement_status"] = np.where(
        inventory["n_measurements_available"].gt(0), "measured", "unexecuted_no_measurement"
    )
    inventory["model_selection_rule"] = "measurement_available_only_tier_not_used"
    return inventory


def variance_fraction(values: np.ndarray, taxa: np.ndarray) -> tuple[float, float, float]:
    frame = pd.DataFrame({"value": np.asarray(values, dtype=float), "taxon_name": taxa})
    grand_mean = float(frame["value"].mean())
    taxon_means = frame.groupby("taxon_name")["value"].transform("mean")
    total_ss = float(np.square(frame["value"] - grand_mean).sum())
    within_ss = float(np.square(frame["value"] - taxon_means).sum())
    among_ss = max(0.0, total_ss - within_ss)
    fraction = within_ss / total_ss if total_ss > 0 else float("nan")
    return total_ss, within_ss, fraction


def build_variance_decomposition(
    contract: pd.DataFrame,
    traits: pd.DataFrame,
    seed: int,
    bootstrap_repeats: int = 500,
) -> pd.DataFrame:
    """Quantify observation-level visible variance hidden by taxon means.

    The equal-replication sensitivity samples two observations from every taxon
    having at least two measurements.  This is a visible-image variance audit,
    not a genetic or plastic variance decomposition.
    """
    rows: list[dict[str, Any]] = []
    for endpoint in contract.to_dict("records"):
        endpoint_id = endpoint["endpoint_id"]
        part = traits[traits["endpoint_id"].eq(endpoint_id)][["taxon_name", "value"]].dropna().copy()
        base = {
            "endpoint_id": endpoint_id,
            "module": endpoint["module"],
            "analysis_tier_metadata_only": endpoint["analysis_tier"],
            "validation_status": endpoint["validation_status"],
            "variance_grain": "observation_below_taxon_mean",
            "interpretation_boundary": "visible_image_variance_not_genetic_variance_or_plasticity",
        }
        if part.empty:
            rows.append(
                {
                    **base,
                    "status": "unexecuted_no_measurement",
                    "n_observations": 0,
                    "n_taxa": 0,
                }
            )
            continue
        counts = part.groupby("taxon_name").size()
        total_ss, within_ss, fraction = variance_fraction(
            part["value"].to_numpy(float), part["taxon_name"].astype(str).to_numpy()
        )
        repeated_taxa = counts[counts >= 2].index
        sensitivity = part[part["taxon_name"].isin(repeated_taxa)].copy()
        bootstrap_values: list[float] = []
        if len(repeated_taxa) >= 2:
            rng = stable_rng(seed, "equal_replication_variance", endpoint_id)
            grouped = {taxon: group.index.to_numpy() for taxon, group in sensitivity.groupby("taxon_name")}
            for _ in range(bootstrap_repeats):
                sampled_index = np.concatenate(
                    [rng.choice(index, size=2, replace=False) for index in grouped.values()]
                )
                sampled = sensitivity.loc[sampled_index]
                _, _, sampled_fraction = variance_fraction(
                    sampled["value"].to_numpy(float),
                    sampled["taxon_name"].astype(str).to_numpy(),
                )
                bootstrap_values.append(sampled_fraction)
        rows.append(
            {
                **base,
                "status": "ok",
                "n_observations": int(len(part)),
                "n_taxa": int(part["taxon_name"].nunique()),
                "n_taxa_with_at_least_two": int(len(repeated_taxa)),
                "median_observations_per_taxon": float(counts.median()),
                "total_sum_of_squares": total_ss,
                "within_taxon_sum_of_squares": within_ss,
                "fraction_below_taxon_means": fraction,
                "fraction_retained_between_taxa": 1.0 - fraction if np.isfinite(fraction) else np.nan,
                "equal_replication_observations_per_taxon": 2,
                "equal_replication_bootstrap_repeats": bootstrap_repeats,
                "equal_replication_fraction_median": (
                    float(np.median(bootstrap_values)) if bootstrap_values else np.nan
                ),
                "equal_replication_fraction_low_95": (
                    float(np.quantile(bootstrap_values, 0.025)) if bootstrap_values else np.nan
                ),
                "equal_replication_fraction_high_95": (
                    float(np.quantile(bootstrap_values, 0.975)) if bootstrap_values else np.nan
                ),
            }
        )
    return pd.DataFrame(rows)


def build_trait_geography(
    contract: pd.DataFrame,
    traits: pd.DataFrame,
    environment: pd.DataFrame,
    cell_degrees: float = 5.0,
) -> pd.DataFrame:
    coordinates = environment[["obs_id", "latitude", "longitude"]]
    joined = traits.merge(coordinates, on="obs_id", how="inner", validate="many_to_one")
    rows: list[dict[str, Any]] = []
    for endpoint in contract.to_dict("records"):
        endpoint_id = endpoint["endpoint_id"]
        part = joined[joined["endpoint_id"].eq(endpoint_id)].dropna(
            subset=["value", "latitude", "longitude"]
        )
        base = {
            "endpoint_id": endpoint_id,
            "module": endpoint["module"],
            "analysis_tier_metadata_only": endpoint["analysis_tier"],
            "validation_status": endpoint["validation_status"],
            "spatial_cell_width_degrees": cell_degrees,
        }
        if part.empty:
            rows.append({**base, "status": "unexecuted_no_measurement", "n_observations": 0, "n_taxa": 0})
            continue
        ilat = np.floor((part["latitude"].to_numpy(float) + 90.0) / cell_degrees).astype(int)
        ilon = np.floor((part["longitude"].to_numpy(float) + 180.0) / cell_degrees).astype(int)
        values = part["value"].to_numpy(float)
        rows.append(
            {
                **base,
                "status": "ok",
                "n_observations": int(len(part)),
                "n_taxa": int(part["taxon_name"].nunique()),
                "n_spatial_cells": int(len(set(zip(ilat.tolist(), ilon.tolist())))),
                "latitude_min": float(part["latitude"].min()),
                "latitude_max": float(part["latitude"].max()),
                "longitude_min": float(part["longitude"].min()),
                "longitude_max": float(part["longitude"].max()),
                "value_q05": float(np.quantile(values, 0.05)),
                "value_median": float(np.median(values)),
                "value_q95": float(np.quantile(values, 0.95)),
            }
        )
    return pd.DataFrame(rows)


def unit_observations(unit: dict[str, Any], traits: pd.DataFrame, environment: pd.DataFrame) -> pd.DataFrame:
    members = unit["member_endpoint_ids"]
    part = traits[traits["endpoint_id"].isin(members)][
        ["obs_id", "taxon_name", "endpoint_id", "value"]
    ].copy()
    if part.empty or set(part["endpoint_id"]) != set(members):
        return pd.DataFrame()
    wide = part.pivot(index=["obs_id", "taxon_name"], columns="endpoint_id", values="value").reset_index()
    wide.columns.name = None
    env = environment.drop(columns="taxon_name")
    joined = wide.merge(env, on="obs_id", how="inner", validate="one_to_one")
    return joined


def base_result_row(
    unit: dict[str, Any], predictor: str, block_map: dict[str, str], scale: str, scope: str
) -> dict[str, Any]:
    return {
        "scale": scale,
        "scope": scope,
        "unit_id": unit["unit_id"],
        "member_endpoint_ids": "|".join(unit["member_endpoint_ids"]),
        "n_registered_endpoints_in_unit": len(unit["member_endpoint_ids"]),
        "module": unit["module"],
        "analysis_tier_metadata_only": unit["analysis_tier"],
        "validation_status": unit["validation_status"],
        "inferential_unit": unit["inferential_unit"],
        "predictor": predictor,
        "environment_block": block_map[predictor],
    }


def fit_linear_within(data: pd.DataFrame, response: str, predictor: str) -> dict[str, Any]:
    y = demean_standardize(data, response)
    x = demean_standardize(data, predictor)
    result = sm.OLS(y, x[:, None]).fit(
        cov_type="cluster", cov_kwds={"groups": data["taxon_name"].to_numpy()}
    )
    return {
        "beta_std": float(result.params[0]),
        "standard_error": float(result.bse[0]),
        "confidence_low_95": float(result.conf_int()[0, 0]),
        "confidence_high_95": float(result.conf_int()[0, 1]),
        "p_value": float(result.pvalues[0]),
        "p_value_method": "taxon_clustered_ols",
    }


def fit_circular_within(
    data: pd.DataFrame,
    members: list[str],
    predictor: str,
    permutations: int,
    rng: np.random.Generator,
) -> dict[str, Any]:
    sine = next(value for value in members if "sin" in value)
    cosine = next(value for value in members if "cos" in value)
    ys = demean_standardize(data, sine)
    yc = demean_standardize(data, cosine)
    x = demean_standardize(data, predictor)
    taxa = data["taxon_name"].astype(str).to_numpy()
    return permutation_p_circular_within(ys, yc, x, taxa, permutations, rng)


def run_within(
    units: list[dict[str, Any]],
    traits: pd.DataFrame,
    environment: pd.DataFrame,
    analysis_contract: dict[str, Any],
    block_map: dict[str, str],
) -> pd.DataFrame:
    rules = analysis_contract["within_taxon"]
    minimum_observations = int(rules["minimum_observations"])
    minimum_taxa = int(rules["minimum_taxa"])
    permutations = int(rules["permutations"])
    seed = int(analysis_contract["seed"])
    rows: list[dict[str, Any]] = []
    for unit in units:
        data = unit_observations(unit, traits, environment)
        for predictor in EXPECTED_PREDICTORS:
            base = base_result_row(unit, predictor, block_map, "within_taxon", "within_taxon")
            if data.empty:
                rows.append({**base, "status": "unexecuted_no_measurement", "n_observations": 0, "n_taxa": 0})
                continue
            complete = data.dropna(subset=[*unit["member_endpoint_ids"], predictor]).copy()
            counts = complete.groupby("taxon_name").size()
            complete = complete[complete["taxon_name"].isin(counts[counts >= 2].index)].copy()
            n_observations = int(len(complete))
            n_taxa = int(complete["taxon_name"].nunique())
            base.update({"n_observations": n_observations, "n_taxa": n_taxa})
            if n_observations < minimum_observations or n_taxa < minimum_taxa:
                rows.append({**base, "status": "insufficient_support"})
                continue
            try:
                if unit["inferential_unit"] == "linear_endpoint":
                    fit = fit_linear_within(complete, unit["member_endpoint_ids"][0], predictor)
                else:
                    rng = stable_rng(seed, "within", unit["unit_id"], predictor)
                    fit = fit_circular_within(
                        complete, unit["member_endpoint_ids"], predictor, permutations, rng
                    )
                rows.append({**base, "status": "ok", **fit})
            except Exception as error:
                rows.append({**base, "status": f"failed:{type(error).__name__}:{error}"})
    return pd.DataFrame(rows)


def taxon_environment(environment: pd.DataFrame) -> pd.DataFrame:
    counts = environment.groupby("taxon_name").size().rename("n_environment_observations")
    medians = environment.groupby("taxon_name")[EXPECTED_PREDICTORS].median(numeric_only=True)
    return medians.join(counts).reset_index()


def unit_taxon_traits(unit: dict[str, Any], traits: pd.DataFrame, minimum_observations: int) -> pd.DataFrame:
    members = unit["member_endpoint_ids"]
    part = traits[traits["endpoint_id"].isin(members)][
        ["obs_id", "taxon_name", "endpoint_id", "value"]
    ].copy()
    if part.empty or set(part["endpoint_id"]) != set(members):
        return pd.DataFrame()
    wide = part.pivot(index=["obs_id", "taxon_name"], columns="endpoint_id", values="value").reset_index()
    wide.columns.name = None
    wide = wide.dropna(subset=members)
    counts = wide.groupby("taxon_name").size().rename("n_trait_observations")
    medians = wide.groupby("taxon_name")[members].median()
    summary = medians.join(counts).reset_index()
    return summary[summary["n_trait_observations"].ge(minimum_observations)].copy()


def run_among_scope(
    units: list[dict[str, Any]],
    traits: pd.DataFrame,
    environment: pd.DataFrame,
    analysis_contract: dict[str, Any],
    block_map: dict[str, str],
    minimum_trait_observations: int,
    scope: str,
) -> pd.DataFrame:
    rules = analysis_contract["among_taxon"]
    minimum_taxa = int(rules["minimum_taxa"])
    permutations = int(rules["permutations"])
    seed = int(analysis_contract["seed"])
    env_taxon = taxon_environment(environment)
    rows: list[dict[str, Any]] = []
    for unit in units:
        taxon_traits = unit_taxon_traits(unit, traits, minimum_trait_observations)
        data = (
            taxon_traits.merge(env_taxon, on="taxon_name", how="inner", validate="one_to_one")
            if not taxon_traits.empty
            else pd.DataFrame()
        )
        for predictor in EXPECTED_PREDICTORS:
            base = base_result_row(unit, predictor, block_map, "among_taxon", scope)
            base["minimum_trait_observations_per_taxon"] = minimum_trait_observations
            if data.empty:
                rows.append({**base, "status": "unexecuted_no_measurement", "n_observations": 0, "n_taxa": 0})
                continue
            complete = data.dropna(subset=[*unit["member_endpoint_ids"], predictor]).copy()
            n_taxa = int(len(complete))
            base.update({"n_observations": n_taxa, "n_taxa": n_taxa})
            if n_taxa < minimum_taxa:
                rows.append({**base, "status": "insufficient_support"})
                continue
            rng = stable_rng(seed, scope, unit["unit_id"], predictor)
            try:
                if unit["inferential_unit"] == "linear_endpoint":
                    response = unit["member_endpoint_ids"][0]
                    beta, p_value = permutation_p_linear(
                        complete[response].to_numpy(float),
                        complete[predictor].to_numpy(float),
                        permutations,
                        rng,
                    )
                    fit = {
                        "beta_std": beta,
                        "p_value": p_value,
                        "p_value_method": f"taxon_label_permutation_two_sided_{permutations}",
                    }
                else:
                    sine = next(value for value in unit["member_endpoint_ids"] if "sin" in value)
                    cosine = next(value for value in unit["member_endpoint_ids"] if "cos" in value)
                    fit = permutation_p_circular_among(
                        complete[sine].to_numpy(float),
                        complete[cosine].to_numpy(float),
                        complete[predictor].to_numpy(float),
                        permutations,
                        rng,
                    )
                    fit["p_value_method"] = f"taxon_label_permutation_joint_circular_{permutations}"
                rows.append({**base, "status": "ok", **fit})
            except Exception as error:
                rows.append({**base, "status": f"failed:{type(error).__name__}:{error}"})
    return pd.DataFrame(rows)


def add_global_fdr(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    result["q_fdr_bh_global_family"] = np.nan
    ok = result["status"].eq("ok") & pd.to_numeric(result.get("p_value"), errors="coerce").notna()
    if ok.any():
        result.loc[ok, "q_fdr_bh_global_family"] = bh_adjust(result.loc[ok, "p_value"].astype(float))
    result["fdr_significant_0_05"] = result["q_fdr_bh_global_family"].lt(0.05)
    result["result_retained_regardless_of_q"] = True
    return result


def build_cross_scale(within: pd.DataFrame, among: pd.DataFrame, primary_scope: str) -> pd.DataFrame:
    keys = [
        "unit_id",
        "member_endpoint_ids",
        "module",
        "inferential_unit",
        "predictor",
        "environment_block",
    ]
    left = within.copy().rename(
        columns={
            "status": "status_within",
            "n_observations": "n_observations_within",
            "n_taxa": "n_taxa_within",
            "beta_std": "beta_std_within",
            "p_value": "p_value_within",
            "q_fdr_bh_global_family": "q_fdr_bh_within",
            "fdr_significant_0_05": "fdr_significant_within",
        }
    )
    right = among[among["scope"].eq(primary_scope)].copy().rename(
        columns={
            "status": "status_among",
            "n_observations": "n_observations_among",
            "n_taxa": "n_taxa_among",
            "beta_std": "beta_std_among",
            "p_value": "p_value_among",
            "q_fdr_bh_global_family": "q_fdr_bh_among",
            "fdr_significant_0_05": "fdr_significant_among",
        }
    )
    keep_left = keys + [
        "analysis_tier_metadata_only",
        "validation_status",
        "status_within",
        "n_observations_within",
        "n_taxa_within",
        "beta_std_within",
        "p_value_within",
        "q_fdr_bh_within",
        "fdr_significant_within",
    ]
    keep_right = keys + [
        "status_among",
        "n_observations_among",
        "n_taxa_among",
        "beta_std_among",
        "p_value_among",
        "q_fdr_bh_among",
        "fdr_significant_among",
    ]
    for column in keep_left:
        if column not in left:
            left[column] = np.nan
    for column in keep_right:
        if column not in right:
            right[column] = np.nan
    merged = left[keep_left].merge(right[keep_right], on=keys, how="outer", validate="one_to_one")
    comparable = merged["status_within"].eq("ok") & merged["status_among"].eq("ok")
    within_sig = merged["fdr_significant_within"].fillna(False).astype(bool)
    among_sig = merged["fdr_significant_among"].fillna(False).astype(bool)
    merged["cross_scale_class"] = np.select(
        [
            comparable & within_sig & among_sig,
            comparable & within_sig & ~among_sig,
            comparable & ~within_sig & among_sig,
            comparable,
        ],
        ["both_scales", "within_only", "among_only", "neither"],
        default="not_comparable",
    )
    merged["linear_sign_concordant"] = pd.Series(pd.NA, index=merged.index, dtype="boolean")
    linear = comparable & merged["inferential_unit"].eq("linear_endpoint")
    bw = pd.to_numeric(merged["beta_std_within"], errors="coerce")
    ba = pd.to_numeric(merged["beta_std_among"], errors="coerce")
    valid = linear & bw.notna() & ba.notna() & bw.ne(0) & ba.ne(0)
    merged.loc[valid, "linear_sign_concordant"] = (
        np.sign(bw[valid]).eq(np.sign(ba[valid])).to_numpy(dtype=bool)
    )
    return merged


def main() -> int:
    args = parse_args()
    contract, analysis_contract = load_contracts(args.contract, args.analysis_contract)
    block_map = environment_block_map(analysis_contract)
    traits, environment = load_data(args.traits_long, args.environment, contract, analysis_contract)
    units = inferential_units(contract)
    if len(units) != len(contract) - 1:
        raise ValueError(f"Expected 26 inferential units from 27 endpoints, found {len(units)}")

    inventory = build_inventory(contract, traits, environment)
    variance = build_variance_decomposition(contract, traits, int(analysis_contract["seed"]))
    geography = build_trait_geography(contract, traits, environment)
    within = add_global_fdr(run_within(units, traits, environment, analysis_contract, block_map))
    among_rules = analysis_contract["among_taxon"]
    primary_min = int(among_rules["primary_minimum_trait_observations_per_taxon"])
    sensitivity_min = int(among_rules["sensitivity_minimum_trait_observations_per_taxon"])
    primary_scope = f"among_taxon_min{primary_min}"
    sensitivity_scope = f"among_taxon_min{sensitivity_min}"
    among_frames = [
        add_global_fdr(
            run_among_scope(
                units, traits, environment, analysis_contract, block_map, primary_min, primary_scope
            )
        )
    ]
    if sensitivity_min != primary_min:
        among_frames.append(
            add_global_fdr(
                run_among_scope(
                    units,
                    traits,
                    environment,
                    analysis_contract,
                    block_map,
                    sensitivity_min,
                    sensitivity_scope,
                )
            )
        )
    among = pd.concat(among_frames, ignore_index=True, sort=False)
    cross_scale = build_cross_scale(within, among, primary_scope)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    inventory.to_csv(args.out_dir / "v2_full27_endpoint_inventory.csv", index=False, encoding="utf-8-sig")
    variance.to_csv(args.out_dir / "v2_full27_variance_decomposition.csv", index=False, encoding="utf-8-sig")
    geography.to_csv(args.out_dir / "v2_full27_trait_geography.csv", index=False, encoding="utf-8-sig")
    within.to_csv(args.out_dir / "v2_full27_environment_within.csv", index=False, encoding="utf-8-sig")
    among.to_csv(args.out_dir / "v2_full27_environment_among.csv", index=False, encoding="utf-8-sig")
    cross_scale.to_csv(args.out_dir / "v2_full27_environment_cross_scale.csv", index=False, encoding="utf-8-sig")

    unexecuted = inventory.loc[
        inventory["measurement_status"].eq("unexecuted_no_measurement"), "endpoint_id"
    ].tolist()
    report = {
        "analysis_id": analysis_contract["analysis_id"],
        "status": analysis_contract["status"],
        "n_registered_endpoints_start": int(len(contract)),
        "n_measured_endpoints_model_entry": int(inventory["measurement_status"].eq("measured").sum()),
        "n_unexecuted_endpoints_retained": int(len(unexecuted)),
        "unexecuted_endpoints": unexecuted,
        "n_inferential_units": int(len(units)),
        "n_measured_variance_decompositions": int(variance["status"].eq("ok").sum()),
        "fraction_below_taxon_means_range": [
            float(variance.loc[variance["status"].eq("ok"), "fraction_below_taxon_means"].min()),
            float(variance.loc[variance["status"].eq("ok"), "fraction_below_taxon_means"].max()),
        ],
        "n_environment_predictors": int(len(EXPECTED_PREDICTORS)),
        "environment_predictors": EXPECTED_PREDICTORS,
        "n_environment_observations": int(len(environment)),
        "n_environment_taxa": int(environment["taxon_name"].nunique()),
        "environment_coverage": {
            predictor: float(environment[predictor].notna().mean()) for predictor in EXPECTED_PREDICTORS
        },
        "within_status_counts": within["status"].value_counts(dropna=False).to_dict(),
        "among_status_counts_by_scope": {
            scope: part["status"].value_counts(dropna=False).to_dict()
            for scope, part in among.groupby("scope", sort=True)
        },
        "within_fdr_signals": int(within["fdr_significant_0_05"].sum()),
        "among_fdr_signals_by_scope": {
            scope: int(part["fdr_significant_0_05"].sum())
            for scope, part in among.groupby("scope", sort=True)
        },
        "cross_scale_class_counts": cross_scale["cross_scale_class"].value_counts().to_dict(),
        "multiplicity": analysis_contract["multiplicity"],
        "trait_selection": "all 27 registered; model entry by finite measurement_available only; analysis_tier is metadata only",
        "old_v1_nine_endpoint_subset_used": False,
        "claim_boundary": analysis_contract["claim_boundary"],
        "input_sha256": {
            "endpoint_contract": sha256_file(args.contract),
            "analysis_contract": sha256_file(args.analysis_contract),
            "traits_long": sha256_file(args.traits_long),
            "environment": sha256_file(args.environment),
        },
    }
    (args.out_dir / "v2_full27_environment_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
