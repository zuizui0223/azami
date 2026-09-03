#!/usr/bin/env python3
"""Rebuild the frozen WCVP/TDWG native-status table without manuscript inputs.

This is the source-classification subset of the locked 2026-08-26 native-range
sensitivity.  It deliberately omits the old claim/model layer and reconstructs
only the observation-level table needed by the frozen Chapter 1 v2 sampling-
composition audit.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path
from typing import Any
from urllib.parse import quote_plus

import numpy as np
import pandas as pd
import requests
from shapely.geometry import Point, shape
from shapely.strtree import STRtree

GBIF_API = "https://api.gbif.org/v1"
USER_AGENT = "azami-ch1-frozen-native-status/1.0"
EXPECTED_OBSERVATION_SHA256 = "2172e3570f684770d0f919ecd81265c8460574e287bc4fb057db4f719cab7bb0"
EXPECTED_OUTPUT_SHA256 = "c01eeb9ff245d7f73da1a12fa4eede904dd9770467655f20e3d85de2ac8dd84a"
EXPECTED_TDWG_GEOJSON_SHA256 = "c172bcf6aba20e19477adc60aebf0023068f0175dca8480f0760b090dcf64840"
EXPECTED_ROWS = 46276
EXPECTED_TAXA = 259
EXPECTED_RESOLVED_TAXA = 245
EXPECTED_STATUS_COUNTS = {
    "native": 27066,
    "introduced": 10554,
    "unresolved_taxon": 5491,
    "unmapped_tdwg": 2100,
    "unlisted": 1065,
}
OUTPUT_COLUMNS = [
    "obs_id",
    "taxon_name",
    "latitude",
    "longitude",
    "tdwg_level3_code",
    "resolution_status",
    "accepted_key",
    "accepted_name",
    "native_range_status",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--observation", required=True, type=Path)
    parser.add_argument(
        "--contract",
        type=Path,
        default=Path("analysis/ch1/native_range_sensitivity_contract.json"),
    )
    parser.add_argument("--out-csv", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--timeout-seconds", type=float, default=30.0)
    parser.add_argument("--max-retries", type=int, default=4)
    parser.add_argument("--sleep-seconds", type=float, default=0.04)
    parser.add_argument(
        "--expected-observation-sha256",
        default=EXPECTED_OBSERVATION_SHA256,
    )
    parser.add_argument(
        "--expected-output-sha256",
        default=EXPECTED_OUTPUT_SHA256,
    )
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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
        result
        for result in list(payload.get("results") or [])
        if canonical(result).casefold() == name.casefold()
    ]
    targets = sorted(
        {key for result in exact if (key := accepted_key(result)) is not None}
    )
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
        normalized(
            result.get("accepted")
            or result.get("acceptedScientificName")
            or result.get("scientificName")
        )
        for result in exact
        if accepted_key(result) == key
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
                url,
                timeout=timeout,
                headers={"Accept": "application/json", "User-Agent": USER_AGENT},
            )
            response.raise_for_status()
            return response.json()
        except Exception as error:
            last_error = error
            time.sleep(min(2**attempt, 8))
    raise RuntimeError(f"Request failed after {retries} attempts: {url}: {last_error}")


def fetch_name_resolution(
    names: list[str],
    dataset_key: str,
    timeout: float,
    retries: int,
    sleep: float,
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
    resolution: pd.DataFrame,
    timeout: float,
    retries: int,
    sleep: float,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    keys = sorted(
        {int(value) for value in resolution["accepted_key"] if normalized(value)}
    )
    for key_index, key in enumerate(keys):
        offset = 0
        while True:
            payload = request_json(
                f"{GBIF_API}/species/{key}/distributions?limit=1000&offset={offset}",
                timeout,
                retries,
            )
            results = list(payload.get("results") or [])
            for result in results:
                location_id = normalized(result.get("locationId"))
                code = (
                    location_id.split(":", 1)[-1].upper()
                    if location_id.startswith("TDWG:")
                    else ""
                )
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
                rows.append(
                    {
                        "accepted_key": key,
                        "tdwg_level3_code": code,
                        "locality": normalized(result.get("locality")),
                        "establishment_means": means,
                        "occurrence_status": occurrence,
                        "classified_status": status if code else "uncertain",
                    }
                )
            if bool(payload.get("endOfRecords", True)):
                break
            offset += len(results)
            if not results:
                raise RuntimeError(
                    f"Distribution pagination stalled for accepted key {key}"
                )
        if key_index + 1 < len(keys):
            time.sleep(sleep)
    return pd.DataFrame(rows)


def collapse_distribution_status(
    distributions: pd.DataFrame,
) -> dict[int, dict[str, str]]:
    lookup: dict[int, dict[str, str]] = {}
    if distributions.empty:
        return lookup
    for (key, code), group in distributions.groupby(
        ["accepted_key", "tdwg_level3_code"]
    ):
        statuses = set(group["classified_status"])
        status = next(iter(statuses)) if len(statuses) == 1 else "ambiguous"
        lookup.setdefault(int(key), {})[str(code)] = status
    return lookup


def load_pinned_level3(
    contract: dict[str, Any], timeout: float, retries: int
) -> tuple[dict[str, Any], str]:
    source = contract["tdwg_level3"]
    url = (
        "https://raw.githubusercontent.com/tdwg/wgsrpd/"
        f"{source['commit']}/{source['path']}"
    )
    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            response = requests.get(
                url, timeout=timeout, headers={"User-Agent": USER_AGENT}
            )
            response.raise_for_status()
            content = response.content
            digest = hashlib.sha256(content).hexdigest()
            return json.loads(content.decode("utf-8")), digest
        except Exception as error:
            last_error = error
            time.sleep(min(2**attempt, 8))
    raise RuntimeError(f"Pinned TDWG geometry request failed: {last_error}")


def assign_level3(frame: pd.DataFrame, geojson: dict[str, Any]) -> pd.Series:
    geometries = [shape(feature["geometry"]) for feature in geojson["features"]]
    codes = [
        normalized(feature["properties"].get("LEVEL3_COD")).upper()
        for feature in geojson["features"]
    ]
    tree = STRtree(geometries)
    points = np.array(
        [
            Point(float(longitude), float(latitude))
            for longitude, latitude in zip(frame["longitude"], frame["latitude"])
        ],
        dtype=object,
    )
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
    mapping = resolution.set_index("input_name")[
        ["resolution_status", "accepted_key", "accepted_name"]
    ]
    result = result.join(mapping, on="taxon_name", validate="many_to_one")
    lookup = collapse_distribution_status(distributions)

    statuses = []
    for row in result[
        ["resolution_status", "accepted_key", "tdwg_level3_code"]
    ].to_dict("records"):
        if row["resolution_status"] != "resolved_unique_accepted_key":
            statuses.append("unresolved_taxon")
        elif not normalized(row["tdwg_level3_code"]):
            statuses.append("unmapped_tdwg")
        else:
            statuses.append(
                lookup.get(int(row["accepted_key"]), {}).get(
                    str(row["tdwg_level3_code"]), "unlisted"
                )
            )
    result["native_range_status"] = statuses
    return result


def main() -> int:
    args = parse_args()
    observation_sha = sha256_file(args.observation)
    if observation_sha != args.expected_observation_sha256:
        raise SystemExit(
            "strict-spatial observation SHA-256 mismatch: "
            f"{observation_sha}; expected {args.expected_observation_sha256}"
        )

    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    if contract.get("status") != "locked_before_native_range_outcome_execution":
        raise SystemExit("Native-range contract is not the frozen locked contract")

    metadata = request_json(
        f"{GBIF_API}/dataset/{contract['wcvp']['dataset_key']}",
        args.timeout_seconds,
        args.max_retries,
    )
    if normalized(metadata.get("doi")).casefold() != contract["wcvp"][
        "dataset_doi"
    ].casefold():
        raise SystemExit("WCVP dataset DOI differs from the frozen contract")
    if normalized(metadata.get("modified")) != contract["wcvp"]["dataset_modified"]:
        raise SystemExit(
            "WCVP dataset modification timestamp differs from the frozen contract"
        )

    source = pd.read_csv(
        args.observation,
        usecols=["obs_id", "taxon_name", "latitude", "longitude"],
        low_memory=False,
    )
    if len(source) != EXPECTED_ROWS or source["obs_id"].nunique() != EXPECTED_ROWS:
        raise SystemExit("Strict-spatial source row identity/count is not frozen")
    if source["taxon_name"].nunique() != EXPECTED_TAXA:
        raise SystemExit("Strict-spatial source taxon count is not frozen")

    names = sorted(source["taxon_name"].dropna().map(normalized).unique())
    resolution = fetch_name_resolution(
        names,
        contract["wcvp"]["dataset_key"],
        args.timeout_seconds,
        args.max_retries,
        args.sleep_seconds,
    )
    distributions = fetch_distributions(
        resolution,
        args.timeout_seconds,
        args.max_retries,
        args.sleep_seconds,
    )
    geojson, tdwg_sha = load_pinned_level3(
        contract, args.timeout_seconds, args.max_retries
    )
    if tdwg_sha != EXPECTED_TDWG_GEOJSON_SHA256:
        raise SystemExit(
            f"Pinned TDWG GeoJSON SHA-256 mismatch: {tdwg_sha}; "
            f"expected {EXPECTED_TDWG_GEOJSON_SHA256}"
        )

    classified = classify_observations(source, resolution, distributions, geojson)
    status_counts = {
        str(key): int(value)
        for key, value in classified["native_range_status"].value_counts().items()
    }
    resolved_taxa = int(
        resolution["resolution_status"].eq("resolved_unique_accepted_key").sum()
    )

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    classified[OUTPUT_COLUMNS].to_csv(args.out_csv, index=False)
    output_sha = sha256_file(args.out_csv)
    checks = {
        "observation_sha256_matches": observation_sha
        == args.expected_observation_sha256,
        "tdwg_geojson_sha256_matches": tdwg_sha == EXPECTED_TDWG_GEOJSON_SHA256,
        "n_rows_matches": int(len(classified)) == EXPECTED_ROWS,
        "n_taxa_matches": int(classified["taxon_name"].nunique()) == EXPECTED_TAXA,
        "n_resolved_taxa_matches": resolved_taxa == EXPECTED_RESOLVED_TAXA,
        "native_status_counts_match": status_counts == EXPECTED_STATUS_COUNTS,
        "output_sha256_matches": output_sha == args.expected_output_sha256,
    }
    payload = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "source_observation": str(args.observation),
        "source_observation_sha256": observation_sha,
        "contract": str(args.contract),
        "wcvp_dataset_key": contract["wcvp"]["dataset_key"],
        "wcvp_dataset_doi": metadata.get("doi"),
        "wcvp_dataset_modified": metadata.get("modified"),
        "tdwg_commit": contract["tdwg_level3"]["commit"],
        "tdwg_geojson_sha256": tdwg_sha,
        "n_rows": int(len(classified)),
        "n_taxa": int(classified["taxon_name"].nunique()),
        "n_resolved_taxa": resolved_taxa,
        "native_status_counts": status_counts,
        "output_sha256": output_sha,
        "expected_output_sha256": args.expected_output_sha256,
        "checks": checks,
        "claim_boundary": "This reconstructs frozen native-status classification only; it does not infer adaptation, niche response, or introduction history.",
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    if payload["status"] != "PASS":
        raise SystemExit("Frozen native-status reconstruction failed closed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
