#!/usr/bin/env python3
"""Build one long, category-free trait table from primary and extended outputs."""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


FORBIDDEN_CATEGORY_COLUMNS = {
    "ai_candidate_state",
    "analysis_state_ai_all",
    "analysis_state_ai_conservative",
    "observation_ai_all_state",
    "observation_ai_conservative_state",
    "head_orientation_binary",
    "trait_state",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contract", required=True)
    parser.add_argument("--primary-observation", required=True)
    parser.add_argument("--primary-species", required=True)
    parser.add_argument("--extended-observation")
    parser.add_argument("--extended-species")
    parser.add_argument("--out-dir", required=True)
    return parser.parse_args()


def _load_validator(path: Path):
    spec = importlib.util.spec_from_file_location("continuous_trait_contract_validator", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to import contract validator from {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.validate_contract


def _assert_unique(frame: pd.DataFrame, key: str, label: str) -> None:
    if key not in frame:
        raise ValueError(f"{label} is missing key {key!r}")
    if frame[key].astype(str).duplicated().any():
        example = frame.loc[frame[key].astype(str).duplicated(keep=False), key].astype(str).head().tolist()
        raise ValueError(f"{label} must be unique by {key}; examples={example}")


def _merge_optional(left: pd.DataFrame, right_path: str | None, key: str, label: str) -> pd.DataFrame:
    if not right_path:
        return left
    right = pd.read_csv(right_path, low_memory=False)
    _assert_unique(right, key, label)
    overlap = sorted((set(left.columns) & set(right.columns)) - {key, "taxon_name"})
    if overlap:
        raise ValueError(f"{label} duplicates measurement columns already present: {overlap}")
    if key == "obs_id" and "taxon_name" in right:
        right = right.drop(columns="taxon_name")
    merged = left.merge(right, on=key, how="left", validate="one_to_one")
    if len(merged) != len(left):
        raise RuntimeError(f"{label} merge changed row count")
    return merged


def _build_long(frame: pd.DataFrame, contract: pd.DataFrame, grain: str) -> pd.DataFrame:
    key = "obs_id" if grain == "observation" else "taxon_name"
    variable_column = "observation_variable" if grain == "observation" else "species_variable"
    metadata = [column for column in (key, "taxon_name", "latitude", "longitude", "positional_accuracy", "coordinate_precision_tier") if column in frame]
    rows: list[pd.DataFrame] = []
    for endpoint in contract.to_dict("records"):
        variable = str(endpoint[variable_column])
        status = str(endpoint["source_status_column"])
        if grain == "species" and status not in frame.columns:
            observation_count = status.replace("n_usable_heads_", "n_observations_usable_")
            if observation_count in frame.columns:
                status = observation_count
        required = {key, variable, status}
        missing = sorted(required.difference(frame.columns))
        if missing:
            # Extended endpoints may be absent before the expanded image run. Preserve
            # the contract and make the missing execution explicit in the audit table.
            continue
        part = frame[metadata + [variable, status]].copy()
        part = part.rename(columns={variable: "value", status: "n_usable_source_heads"})
        part["value"] = pd.to_numeric(part["value"], errors="coerce")
        part["n_usable_source_heads"] = pd.to_numeric(part["n_usable_source_heads"], errors="coerce").fillna(0).astype(int)
        for column, value in endpoint.items():
            part[column] = value
        part["grain"] = grain
        part["measurement_available"] = np.isfinite(part["value"]) & part["n_usable_source_heads"].gt(0)
        part["analysis_eligible"] = part["measurement_available"] & part["analysis_tier"].isin(["primary", "candidate"])
        part["exclusion_reason"] = np.where(
            part["measurement_available"],
            np.where(part["analysis_eligible"], "", "tier_not_in_inferential_screen"),
            "no_usable_continuous_measurement",
        )
        rows.append(part)
    if not rows:
        raise ValueError(f"No contract endpoints were available at {grain} grain")
    output = pd.concat(rows, ignore_index=True, sort=False)
    if output.duplicated([key, "endpoint_id"]).any():
        raise ValueError(f"Continuous long table has duplicate {key}/endpoint_id rows")
    return output


def build_universe(
    contract: pd.DataFrame,
    primary_observation: pd.DataFrame,
    primary_species: pd.DataFrame,
    extended_observation_path: str | None = None,
    extended_species_path: str | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    for label, frame, key in (
        ("primary observation", primary_observation, "obs_id"),
        ("primary species", primary_species, "taxon_name"),
    ):
        _assert_unique(frame, key, label)
        leaked = sorted(FORBIDDEN_CATEGORY_COLUMNS.intersection(frame.columns))
        if leaked:
            raise ValueError(f"Legacy categorical columns are forbidden in {label}: {leaked}")

    observation = _merge_optional(primary_observation, extended_observation_path, "obs_id", "extended observation")
    species = _merge_optional(primary_species, extended_species_path, "taxon_name", "extended species")
    observation_long = _build_long(observation, contract, "observation")
    species_long = _build_long(species, contract, "species")

    available = set(observation_long["endpoint_id"]).intersection(species_long["endpoint_id"])
    expected = set(contract["endpoint_id"])
    missing = sorted(expected - available)
    measured_observation = set(
        observation_long.loc[observation_long["measurement_available"], "endpoint_id"]
    )
    measured_species = set(
        species_long.loc[species_long["measurement_available"], "endpoint_id"]
    )
    inferential = set(contract.loc[contract["analysis_tier"].isin(["primary", "candidate"]), "endpoint_id"])
    inferential_without_measurements = sorted(inferential - measured_observation)
    report = {
        "status": "ready" if not missing and not inferential_without_measurements else "partial",
        "n_contract_endpoints": int(len(contract)),
        "n_available_endpoints": int(len(available)),
        "n_missing_unexecuted_endpoints": int(len(missing)),
        "missing_unexecuted_endpoints": missing,
        "n_endpoints_with_observation_measurements": int(len(measured_observation)),
        "n_endpoints_with_species_measurements": int(len(measured_species)),
        "n_inferential_endpoints_without_measurements": int(len(inferential_without_measurements)),
        "inferential_endpoints_without_measurements": inferential_without_measurements,
        "n_observations": int(primary_observation["obs_id"].astype(str).nunique()),
        # Observation and species-summary grains can legitimately contain
        # different taxon counts when the upstream species table enforces a
        # minimum number of observations. Keep the historical ``n_taxa`` field
        # for compatibility but expose both grains explicitly.
        "n_taxa": int(primary_species["taxon_name"].astype(str).nunique()),
        "n_observation_taxa": int(primary_observation["taxon_name"].astype(str).nunique()),
        "n_species_summary_taxa": int(primary_species["taxon_name"].astype(str).nunique()),
        "n_observation_measurements_available": int(observation_long["measurement_available"].sum()),
        "n_observation_measurements_analysis_eligible": int(observation_long["analysis_eligible"].sum()),
        "category_columns_in_output": sorted(FORBIDDEN_CATEGORY_COLUMNS.intersection(observation_long.columns)),
        "semantic_status": (
            "One row per grain and continuous endpoint. Missingness records failed or unexecuted measurement; "
            "it is never converted to a category or biological zero."
        ),
    }
    return observation_long, species_long, report


def main() -> None:
    args = parse_args()
    contract_path = Path(args.contract)
    contract = pd.read_csv(contract_path, dtype=str, keep_default_na=False)
    validator = _load_validator(Path(__file__).with_name("87_validate_continuous_trait_contract.py"))
    validator(contract)
    primary_observation = pd.read_csv(args.primary_observation, low_memory=False)
    primary_species = pd.read_csv(args.primary_species, low_memory=False)
    observation_long, species_long, report = build_universe(
        contract,
        primary_observation,
        primary_species,
        args.extended_observation,
        args.extended_species,
    )
    output = Path(args.out_dir)
    output.mkdir(parents=True, exist_ok=True)
    observation_long.to_csv(output / "continuous_trait_universe_observation_long.csv", index=False, encoding="utf-8-sig")
    species_long.to_csv(output / "continuous_trait_universe_species_long.csv", index=False, encoding="utf-8-sig")
    contract.to_csv(output / "continuous_trait_contract_frozen.csv", index=False, encoding="utf-8-sig")
    (output / "continuous_trait_universe_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
