#!/usr/bin/env python3
"""Run the locked WCVP/TDWG level-3 native-range sensitivity."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from urllib.parse import quote_plus

import numpy as np
import pandas as pd
import requests
import statsmodels.api as sm
from shapely.geometry import Point, shape
from shapely.strtree import STRtree
from statsmodels.stats.multitest import multipletests


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))
import run_bias_control_reanalysis as primary  # noqa: E402


GBIF_API = "https://api.gbif.org/v1"
USER_AGENT = "azami-ch1-native-range-sensitivity/1.0"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--observation", required=True)
    parser.add_argument("--claims", required=True)
    parser.add_argument("--contract", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--timeout-seconds", type=float, default=30.0)
    parser.add_argument("--max-retries", type=int, default=4)
    parser.add_argument("--sleep-seconds", type=float, default=0.04)
    parser.add_argument("--min-species-model-n", type=int, default=5)
    return parser.parse_args()


def normalized(value: Any) -> str:
    return " ".join(str(value or "").strip().split())


def canonical(result: dict[str, Any]) -> str:
    return normalized(result.get("canonicalName") or result.get("scientificName"))


def accepted_key(result: dict[str, Any]) -> int | None:
    status = normalized(result.get("taxonomicStatus")).upper()
    if status == "ACCEPTED" and not bool(result.get("synonym", False)):
        value = result.get("key")
    else:
        value = result.get("acceptedKey")
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def resolve_exact_name(name: str, payload: dict[str, Any]) -> dict[str, Any]:
    exact = [
        result for result in list(payload.get("results") or [])
        if canonical(result).casefold() == name.casefold()
    ]
    targets = sorted({key for result in exact if (key := accepted_key(result)) is not None})
    if len(targets) != 1:
        return {
            "input_name": name,
            "resolution_status": "unresolved",
            "accepted_key": "",
            "accepted_name": "",
            "n_exact_records": len(exact),
            "n_distinct_accepted_keys": len(targets),
        }
    key = targets[0]
    names = {
        normalized(result.get("accepted") or result.get("acceptedScientificName") or result.get("scientificName"))
        for result in exact if accepted_key(result) == key
    }
    names.discard("")
    return {
        "input_name": name,
        "resolution_status": "resolved_unique_accepted_key",
        "accepted_key": key,
        "accepted_name": sorted(names)[0] if names else "",
        "n_exact_records": len(exact),
        "n_distinct_accepted_keys": 1,
    }


def request_json(url: str, timeout: float, retries: int) -> dict[str, Any]:
    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            response = requests.get(
                url, timeout=timeout,
                headers={"Accept": "application/json", "User-Agent": USER_AGENT},
            )
            response.raise_for_status()
            return response.json()
        except Exception as error:  # network failures are never converted to absence
            last_error = error
            time.sleep(min(2 ** attempt, 8))
    raise RuntimeError(f"Request failed after {retries} attempts: {url}: {last_error}")


def fetch_name_resolution(
    names: list[str], dataset_key: str, timeout: float, retries: int, sleep: float,
) -> pd.DataFrame:
    rows = []
    for index, name in enumerate(names):
        url = (
            f"{GBIF_API}/species/search?dataset_key={quote_plus(dataset_key)}"
            f"&q={quote_plus(name)}&limit=100"
        )
        rows.append(resolve_exact_name(name, request_json(url, timeout, retries)))
        if index + 1 < len(names):
            time.sleep(sleep)
    return pd.DataFrame(rows)


def fetch_distributions(
    resolution: pd.DataFrame, timeout: float, retries: int, sleep: float,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    keys = sorted({int(value) for value in resolution["accepted_key"] if normalized(value)})
    for key_index, key in enumerate(keys):
        offset = 0
        while True:
            payload = request_json(
                f"{GBIF_API}/species/{key}/distributions?limit=1000&offset={offset}",
                timeout, retries,
            )
            results = list(payload.get("results") or [])
            for result in results:
                location_id = normalized(result.get("locationId"))
                code = location_id.split(":", 1)[-1].upper() if location_id.startswith("TDWG:") else ""
                means = normalized(result.get("establishmentMeans")).upper()
                occurrence = normalized(result.get("occurrenceStatus")).upper()
                if occurrence in {"DOUBTFUL", "ABSENT"}:
                    status = "uncertain"
                elif means in {"", "NATIVE"}:
                    status = "native"
                elif means == "INTRODUCED":
                    status = "introduced"
                else:
                    status = "uncertain"
                rows.append({
                    "accepted_key": key,
                    "tdwg_level3_code": code,
                    "locality": normalized(result.get("locality")),
                    "establishment_means": means,
                    "occurrence_status": occurrence,
                    "classified_status": status if code else "uncertain",
                })
            if bool(payload.get("endOfRecords", True)):
                break
            offset += len(results)
            if not results:
                raise RuntimeError(f"Distribution pagination stalled for accepted key {key}")
        if key_index + 1 < len(keys):
            time.sleep(sleep)
    return pd.DataFrame(rows)


def collapse_distribution_status(distributions: pd.DataFrame) -> dict[int, dict[str, str]]:
    lookup: dict[int, dict[str, str]] = {}
    if distributions.empty:
        return lookup
    for (key, code), group in distributions.groupby(["accepted_key", "tdwg_level3_code"]):
        statuses = set(group["classified_status"])
        status = next(iter(statuses)) if len(statuses) == 1 else "ambiguous"
        lookup.setdefault(int(key), {})[str(code)] = status
    return lookup


def load_pinned_level3(contract: dict[str, Any], timeout: float, retries: int) -> tuple[dict[str, Any], str]:
    source = contract["tdwg_level3"]
    url = (
        "https://raw.githubusercontent.com/tdwg/wgsrpd/"
        f"{source['commit']}/{source['path']}"
    )
    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            response = requests.get(url, timeout=timeout, headers={"User-Agent": USER_AGENT})
            response.raise_for_status()
            content = response.content
            return json.loads(content.decode("utf-8")), hashlib.sha256(content).hexdigest()
        except Exception as error:
            last_error = error
            time.sleep(min(2 ** attempt, 8))
    raise RuntimeError(f"Pinned TDWG geometry request failed: {last_error}")


def assign_level3(frame: pd.DataFrame, geojson: dict[str, Any]) -> pd.Series:
    geometries = [shape(feature["geometry"]) for feature in geojson["features"]]
    codes = [normalized(feature["properties"].get("LEVEL3_COD")).upper() for feature in geojson["features"]]
    tree = STRtree(geometries)
    points = np.array([
        Point(float(longitude), float(latitude))
        for longitude, latitude in zip(frame["longitude"], frame["latitude"])
    ], dtype=object)
    pairs = tree.query(points, predicate="within")
    matches: dict[int, set[str]] = {}
    if pairs.size:
        for point_index, geometry_index in zip(pairs[0], pairs[1]):
            matches.setdefault(int(point_index), set()).add(codes[int(geometry_index)])
    assigned = []
    for index in range(len(frame)):
        values = matches.get(index, set())
        assigned.append(next(iter(values)) if len(values) == 1 else "")
    return pd.Series(assigned, index=frame.index, dtype="string")


def classify_observations(
    frame: pd.DataFrame,
    resolution: pd.DataFrame,
    distributions: pd.DataFrame,
    geojson: dict[str, Any],
) -> pd.DataFrame:
    result = frame.copy()
    result["latitude"] = pd.to_numeric(result["latitude"], errors="coerce")
    result["longitude"] = pd.to_numeric(result["longitude"], errors="coerce")
    if result[["latitude", "longitude"]].isna().any().any():
        raise ValueError("Every native-range row must have finite coordinates")
    result["tdwg_level3_code"] = assign_level3(result, geojson)
    mapping = resolution.set_index("input_name")[["resolution_status", "accepted_key", "accepted_name"]]
    result = result.join(mapping, on="taxon_name", validate="many_to_one")
    lookup = collapse_distribution_status(distributions)

    statuses = []
    for row in result[["resolution_status", "accepted_key", "tdwg_level3_code"]].to_dict("records"):
        if row["resolution_status"] != "resolved_unique_accepted_key":
            statuses.append("unresolved_taxon")
        elif not normalized(row["tdwg_level3_code"]):
            statuses.append("unmapped_tdwg")
        else:
            statuses.append(
                lookup.get(int(row["accepted_key"]), {}).get(str(row["tdwg_level3_code"]), "unlisted")
            )
    result["native_range_status"] = statuses
    return result


def fit_native_models(frame: pd.DataFrame, minimum: int) -> pd.DataFrame:
    native = primary.add_cyclic_date(frame.loc[frame["native_range_status"].eq("native")].copy())
    rows = []
    for trait in primary.TRAITS:
        for predictor in primary.PREDICTORS:
            rows.append({
                "trait": trait,
                "predictor": predictor,
                **primary.fit_model(native, trait, predictor, minimum, "taxon_mean_only"),
            })
    table = pd.DataFrame(rows)
    table["q_fdr_bh"] = np.nan
    table["fdr_significant_0_05"] = False
    successful = table["status"].eq("ok")
    reject, q_values, _, _ = multipletests(
        table.loc[successful, "p_value"].to_numpy(float), alpha=0.05, method="fdr_bh"
    )
    table.loc[successful, "q_fdr_bh"] = q_values
    table.loc[successful, "fdr_significant_0_05"] = reject
    return table


def headline_audit(models: pd.DataFrame, claims: dict[str, Any]) -> pd.DataFrame:
    expected = primary.expected_supported_rows(claims)
    rows = []
    for claim in expected.to_dict("records"):
        match = models.loc[
            models["trait"].eq(claim["trait"])
            & models["predictor"].eq(claim["predictor"])
        ].iloc[0]
        expected_sign = float(np.sign(claim["expected_beta"]))
        retained = bool(
            match["status"] == "ok"
            and np.sign(float(match["beta_std_within"])) == expected_sign
            and bool(match["fdr_significant_0_05"])
        )
        rows.append({
            **claim,
            "native_only_beta": float(match["beta_std_within"]),
            "native_only_ci_low": float(match["ci_low"]),
            "native_only_ci_high": float(match["ci_high"]),
            "native_only_q": float(match["q_fdr_bh"]),
            "native_only_n_observations": int(match["n_observations"]),
            "native_only_n_taxa": int(match["n_taxa"]),
            "native_range_robust": retained,
        })
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    contract_path = Path(args.contract)
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    if contract.get("status") != "locked_before_native_range_outcome_execution":
        raise ValueError("Native-range contract is not locked")
    metadata = request_json(
        f"{GBIF_API}/dataset/{contract['wcvp']['dataset_key']}",
        args.timeout_seconds, args.max_retries,
    )
    if normalized(metadata.get("doi")).casefold() != contract["wcvp"]["dataset_doi"].casefold():
        raise ValueError("WCVP dataset DOI differs from locked contract")
    if normalized(metadata.get("modified")) != contract["wcvp"]["dataset_modified"]:
        raise ValueError("WCVP dataset modification timestamp differs from locked contract")

    frame = pd.read_csv(args.observation, low_memory=False)
    required = {"obs_id", "taxon_name", "latitude", "longitude", "observed_on", *primary.TRAITS, *primary.PREDICTORS}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Observation table is missing columns: {sorted(missing)}")
    if frame["obs_id"].astype(str).duplicated().any():
        raise ValueError("Observation table must be unique by obs_id")

    output = Path(args.out_dir)
    output.mkdir(parents=True, exist_ok=True)
    names = sorted(frame["taxon_name"].dropna().map(normalized).unique())
    resolution = fetch_name_resolution(
        names, contract["wcvp"]["dataset_key"],
        args.timeout_seconds, args.max_retries, args.sleep_seconds,
    )
    distributions = fetch_distributions(
        resolution, args.timeout_seconds, args.max_retries, args.sleep_seconds,
    )
    geojson, geometry_sha256 = load_pinned_level3(
        contract, args.timeout_seconds, args.max_retries,
    )
    classified = classify_observations(frame, resolution, distributions, geojson)
    models = fit_native_models(classified, args.min_species_model_n)
    claims = json.loads(Path(args.claims).read_text(encoding="utf-8"))
    audit = headline_audit(models, claims)

    status_summary = (
        classified.groupby("native_range_status", as_index=False)
        .agg(n_observations=("obs_id", "size"), n_taxa=("taxon_name", "nunique"))
        .sort_values("n_observations", ascending=False)
    )
    taxon_summary = (
        classified.groupby(["taxon_name", "native_range_status"], as_index=False)
        .size().rename(columns={"size": "n_observations"})
        .sort_values(["taxon_name", "n_observations"], ascending=[True, False])
    )

    resolution.to_csv(output / "wcvp_name_resolution.csv", index=False)
    distributions.to_csv(output / "wcvp_distribution_records.csv", index=False)
    classified[[
        "obs_id", "taxon_name", "latitude", "longitude", "tdwg_level3_code",
        "resolution_status", "accepted_key", "accepted_name", "native_range_status",
    ]].to_csv(output / "observation_native_status.csv", index=False)
    status_summary.to_csv(output / "native_status_summary.csv", index=False)
    taxon_summary.to_csv(output / "native_status_by_taxon.csv", index=False)
    models.to_csv(output / "native_range_models.csv", index=False)
    audit.to_csv(output / "native_range_headline_audit.csv", index=False)
    source_manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "wcvp_dataset_key": contract["wcvp"]["dataset_key"],
        "wcvp_dataset_doi": metadata.get("doi"),
        "wcvp_dataset_modified": metadata.get("modified"),
        "wcvp_license": metadata.get("license"),
        "tdwg_commit": contract["tdwg_level3"]["commit"],
        "tdwg_git_blob_sha": contract["tdwg_level3"]["git_blob_sha"],
        "tdwg_geojson_sha256": geometry_sha256,
    }
    (output / "native_range_source_manifest.json").write_text(
        json.dumps(source_manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    report = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "complete_native_range_sensitivity",
        "contract": str(contract_path),
        "n_observations": int(len(classified)),
        "n_taxa": int(classified["taxon_name"].nunique()),
        "n_resolved_taxa": int(resolution["resolution_status"].eq("resolved_unique_accepted_key").sum()),
        "native_status_counts": dict(zip(status_summary["native_range_status"], status_summary["n_observations"].astype(int))),
        "n_headline_rows_native_range_robust": int(audit["native_range_robust"].sum()),
        "n_headline_rows": int(len(audit)),
        "claim_boundary": "Native-only sensitivity does not establish niche response, adaptation or introduction history.",
    }
    (output / "native_range_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
