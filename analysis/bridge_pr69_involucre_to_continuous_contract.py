#!/usr/bin/env python3
"""Bridge frozen PR #69 involucre outputs into the continuous-trait contract.

This script does not remeasure images. It performs three bounded operations:

1. map the three already executed PR #69 contour formulas to their sole
   canonical endpoint identifiers in the category-free contract;
2. rebuild observation- and species-grain extended tables without converting
   missingness to a category or biological zero; and
3. rerun the locked PR #69 resolution/sharpness models from their frozen
   observation table and, when supplied, compare them with the committed
   coefficient table.

The former ``spine_relative_length_max_proxy`` name is treated only as a
provenance alias for ``bract_projection_maximum``. It is not counted as a
second biological trait.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


ENDPOINT_ALIASES = {
    "involucre_projection_roughness": "bract_projection_roughness",
    "involucre_spread_fraction": "bract_spread_fraction",
    "spine_relative_length_max_proxy": "bract_projection_maximum",
}

EXPECTED_FORMULAS = {
    "bract_projection_roughness": "radial_residual_sd_over_equivalent_radius",
    "bract_spread_fraction": "radial_bins_above_four_percent_radius",
    "bract_projection_maximum": "positive_radial_residual_max_over_radius",
}

ENVIRONMENT_COLUMNS = (
    "env_chelsa_bio01_native",
    "env_chelsa_bio04_native",
    "env_chelsa_bio12_native",
    "env_chelsa_bio15_native",
)

COMPARISON_COLUMNS = (
    "beta_standardized",
    "se_cluster_taxon",
    "ci_low",
    "ci_high",
    "p_value",
    "q_fdr_bh_climate_12",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contract", required=True)
    parser.add_argument("--involucre-observations", required=True)
    parser.add_argument("--primary-observation", required=True)
    parser.add_argument("--expected-models")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--comparison-tolerance", type=float, default=1e-9)
    return parser.parse_args()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_resolution_module(path: Path):
    spec = importlib.util.spec_from_file_location("locked_involucre_resolution_audit", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to import locked audit module from {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def validate_inputs(
    contract: pd.DataFrame,
    observations: pd.DataFrame,
    primary: pd.DataFrame,
) -> dict[str, Any]:
    required_observation = {
        "obs_id",
        "taxon_name",
        "n_usable_heads",
        "coordinate_precision_tier",
        "log_min_dimension",
        "log1p_sharpness",
        "resolution_stratum",
        *ENDPOINT_ALIASES,
        *ENVIRONMENT_COLUMNS,
    }
    missing = sorted(required_observation.difference(observations.columns))
    if missing:
        raise ValueError(f"Frozen involucre observations are missing columns: {missing}")
    required_primary = {"obs_id", "taxon_name"}
    missing_primary = sorted(required_primary.difference(primary.columns))
    if missing_primary:
        raise ValueError(f"Primary observation table is missing columns: {missing_primary}")
    if observations["obs_id"].astype(str).duplicated().any():
        raise ValueError("Frozen involucre observations must be unique by obs_id")
    if primary["obs_id"].astype(str).duplicated().any():
        raise ValueError("Primary observation table must be unique by obs_id")

    contract_rows = contract.set_index("endpoint_id")
    for endpoint, expected_formula in EXPECTED_FORMULAS.items():
        if endpoint not in contract_rows.index:
            raise ValueError(f"Contract is missing canonical endpoint {endpoint!r}")
        actual_formula = str(contract_rows.loc[endpoint, "formula_id"])
        if actual_formula != expected_formula:
            raise ValueError(
                f"Formula mismatch for {endpoint}: expected={expected_formula!r}, actual={actual_formula!r}"
            )
        if str(contract_rows.loc[endpoint, "analysis_tier"]) != "candidate":
            raise ValueError(f"Bridged endpoint {endpoint!r} must remain candidate tier")

    observation_keys = observations[["obs_id", "taxon_name"]].astype(str)
    primary_keys = primary[["obs_id", "taxon_name"]].astype(str)
    joined = observation_keys.merge(
        primary_keys,
        on="obs_id",
        how="left",
        suffixes=("_involucre", "_primary"),
        validate="one_to_one",
        indicator=True,
    )
    outside = joined["_merge"].ne("both")
    if outside.any():
        examples = joined.loc[outside, "obs_id"].head().tolist()
        raise ValueError(f"Involucre observations fall outside the primary atlas: {examples}")
    mismatch = joined["taxon_name_involucre"].ne(joined["taxon_name_primary"])
    if mismatch.any():
        examples = joined.loc[mismatch, "obs_id"].head().tolist()
        raise ValueError(f"Taxon names disagree for overlapping obs_id values: {examples}")

    numeric = observations[["n_usable_heads", *ENDPOINT_ALIASES, *ENVIRONMENT_COLUMNS]].apply(
        pd.to_numeric, errors="coerce"
    )
    if numeric["n_usable_heads"].isna().any() or numeric["n_usable_heads"].le(0).any():
        raise ValueError("Every bridged observation must have at least one usable source head")
    return {
        "n_primary_observations": int(primary["obs_id"].astype(str).nunique()),
        "n_involucre_observations": int(observations["obs_id"].astype(str).nunique()),
        "n_involucre_taxa": int(observations["taxon_name"].astype(str).nunique()),
        "n_involucre_outside_primary": 0,
        "n_taxon_name_mismatches": 0,
        "endpoint_nonmissing_counts": {
            ENDPOINT_ALIASES[source]: int(numeric[source].notna().sum())
            for source in ENDPOINT_ALIASES
        },
    }


def build_extended_tables(
    contract: pd.DataFrame,
    observations: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    contract_rows = contract.set_index("endpoint_id")
    observation = observations[["obs_id", "taxon_name"]].copy()
    observation["obs_id"] = observation["obs_id"].astype(str)
    observation["taxon_name"] = observation["taxon_name"].astype(str)
    n_usable = pd.to_numeric(observations["n_usable_heads"], errors="raise").astype(int)

    for source, endpoint in ENDPOINT_ALIASES.items():
        row = contract_rows.loc[endpoint]
        value_column = str(row["observation_variable"])
        status_column = str(row["source_status_column"])
        observation[value_column] = pd.to_numeric(observations[source], errors="coerce")
        observation[status_column] = np.where(observation[value_column].notna(), n_usable, 0)

    species_rows: list[dict[str, Any]] = []
    for taxon_name, group in observation.groupby("taxon_name", sort=True, dropna=False):
        payload: dict[str, Any] = {"taxon_name": taxon_name}
        for endpoint in ENDPOINT_ALIASES.values():
            row = contract_rows.loc[endpoint]
            observation_value = str(row["observation_variable"])
            species_value = str(row["species_variable"])
            observation_status = str(row["source_status_column"])
            species_status = observation_status.replace(
                "n_usable_heads_", "n_observations_usable_", 1
            )
            usable = group.loc[
                group[observation_status].gt(0), observation_value
            ].dropna()
            payload[species_value] = float(usable.median()) if len(usable) else np.nan
            payload[species_status] = int(len(usable))
        species_rows.append(payload)
    species = pd.DataFrame(species_rows)

    environment = observations[["obs_id", "coordinate_precision_tier", *ENVIRONMENT_COLUMNS]].copy()
    environment["obs_id"] = environment["obs_id"].astype(str)
    environment_le10km = environment.loc[
        environment["coordinate_precision_tier"].isin(["high_le_1km", "moderate_1_to_10km"])
    ].copy()
    return observation, species, environment, environment_le10km


def compare_models(
    actual: pd.DataFrame,
    expected: pd.DataFrame,
    tolerance: float,
) -> dict[str, Any]:
    keys = ["endpoint", "cohort", "predictor"]
    expected = expected.copy()
    actual = actual.copy()
    for frame in (actual, expected):
        for column in keys:
            frame[column] = frame[column].astype(str)
    joined = actual.merge(expected, on=keys, how="outer", suffixes=("_actual", "_expected"), indicator=True)
    if joined["_merge"].ne("both").any():
        missing = joined.loc[joined["_merge"].ne("both"), keys + ["_merge"]].head().to_dict("records")
        raise ValueError(f"Rerun and expected model rows differ: {missing}")
    maxima: dict[str, float] = {}
    for column in COMPARISON_COLUMNS:
        actual_values = pd.to_numeric(joined[f"{column}_actual"], errors="coerce")
        expected_values = pd.to_numeric(joined[f"{column}_expected"], errors="coerce")
        both_missing = actual_values.isna() & expected_values.isna()
        mismatch_missing = actual_values.isna() ^ expected_values.isna()
        if mismatch_missing.any():
            raise ValueError(f"Missingness differs in model column {column}")
        difference = (actual_values - expected_values).abs().mask(both_missing)
        maximum = float(difference.max()) if difference.notna().any() else 0.0
        maxima[column] = maximum
        if maximum > tolerance:
            raise ValueError(
                f"Rerun differs from committed models in {column}: max_abs_difference={maximum}"
            )
    return {
        "status": "exact_within_tolerance",
        "comparison_tolerance": tolerance,
        "n_compared_rows": int(len(joined)),
        "max_abs_difference_by_column": maxima,
    }


def main() -> None:
    args = parse_args()
    contract_path = Path(args.contract)
    observation_path = Path(args.involucre_observations)
    primary_path = Path(args.primary_observation)
    contract = pd.read_csv(contract_path, dtype=str, keep_default_na=False)
    observations = pd.read_csv(observation_path, dtype={"obs_id": str}, low_memory=False)
    primary = pd.read_csv(primary_path, dtype={"obs_id": str}, low_memory=False)
    input_audit = validate_inputs(contract, observations, primary)
    extended_observation, extended_species, environment, environment_le10km = build_extended_tables(
        contract, observations
    )

    locked_module = _load_resolution_module(Path(__file__).with_name("run_involucre_resolution_audit.py"))
    coefficients, statuses = locked_module.run_models(observations, minimum=100)
    headline_audit = locked_module.build_audit(coefficients)
    comparison: dict[str, Any] = {"status": "not_requested"}
    if args.expected_models:
        expected_path = Path(args.expected_models)
        expected = pd.read_csv(expected_path, low_memory=False)
        comparison = compare_models(coefficients, expected, args.comparison_tolerance)
    else:
        expected_path = None

    coefficients = coefficients.copy()
    statuses = statuses.copy()
    headline_audit = headline_audit.copy()
    coefficients["source_endpoint_alias"] = coefficients["endpoint"]
    statuses["source_endpoint_alias"] = statuses["endpoint"]
    headline_audit["source_endpoint_alias"] = headline_audit["endpoint"]
    coefficients["endpoint"] = coefficients["endpoint"].map(ENDPOINT_ALIASES)
    statuses["endpoint"] = statuses["endpoint"].map(ENDPOINT_ALIASES)
    headline_audit["endpoint"] = headline_audit["endpoint"].map(ENDPOINT_ALIASES)
    headline_audit["claim_status"] = np.where(
        headline_audit["strict_resolution_control_retained"],
        "passes_pr69_resolution_gate_only",
        "withdrawn_by_pr69_resolution_gate",
    )
    headline_audit["submission_eligible"] = False
    headline_audit["submission_blocker"] = (
        "candidate endpoints require botanical reference validation; PR69 resolution gate also applies"
    )

    output = Path(args.out_dir)
    output.mkdir(parents=True, exist_ok=True)
    extended_observation.to_csv(
        output / "extended_continuous_observation_level.csv", index=False, encoding="utf-8-sig"
    )
    extended_species.to_csv(
        output / "extended_continuous_species_level.csv", index=False, encoding="utf-8-sig"
    )
    environment.to_csv(output / "continuous_environment_subset.csv", index=False, encoding="utf-8-sig")
    environment_le10km.to_csv(
        output / "continuous_environment_subset_le10km.csv", index=False, encoding="utf-8-sig"
    )
    coefficients.to_csv(output / "canonical_involucre_adjusted_models.csv", index=False, encoding="utf-8-sig")
    statuses.to_csv(output / "canonical_involucre_model_status.csv", index=False, encoding="utf-8-sig")
    headline_audit.to_csv(output / "canonical_involucre_headline_audit.csv", index=False, encoding="utf-8-sig")

    report = {
        "status": "complete_frozen_output_bridge",
        "scope": (
            "Reuses executed PR69 observation-level measurements; does not remeasure source images "
            "and does not fill unexecuted endpoints."
        ),
        "input_sha256": {
            "contract": _sha256(contract_path),
            "involucre_observations": _sha256(observation_path),
            "primary_observation": _sha256(primary_path),
            **({"expected_models": _sha256(expected_path)} if expected_path else {}),
        },
        **input_audit,
        "n_bridged_canonical_endpoints": len(ENDPOINT_ALIASES),
        "n_environment_observations_all_precision": int(len(environment)),
        "n_environment_observations_le10km": int(len(environment_le10km)),
        "endpoint_aliases": ENDPOINT_ALIASES,
        "model_reproduction": comparison,
        "n_endpoints_passing_pr69_resolution_gate": int(
            headline_audit["strict_resolution_control_retained"].sum()
        ),
        "submission_interpretation": (
            "No bridged candidate endpoint is promoted to a manuscript claim without both the frozen "
            "PR69 bias-control rule and independent continuous botanical reference validation."
        ),
    }
    (output / "bridge_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
