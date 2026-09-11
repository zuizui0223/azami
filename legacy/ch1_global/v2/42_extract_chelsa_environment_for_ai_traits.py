#!/usr/bin/env python3
"""Attach predeclared CHELSA V2.1 bioclim predictors to AI trait observation units.

The input is an observation-level table generated from fully automated trait
measurements. This script does not change trait states, select traits, infer
causality, or collapse multistate outcomes. It only joins four predeclared
macroclimate predictors, coordinate-quality fields, and spatial blocking fields.

Raster values remain in source-native numeric units. Raster metadata (scale,
offset, nodata, CRS, transform, tags) are exported so later model code can
standardize predictors without assuming undocumented unit conversions.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import rasterio
from pyproj import Transformer


REQUIRED_LONG = {
    "obs_id", "trait_id", "taxon_name", "latitude", "longitude",
    "coordinate_usable_for_environment", "positional_accuracy",
    "observation_ai_all_state", "observation_ai_conservative_state",
}
REQUIRED_WIDE = {
    "obs_id", "taxon_name", "latitude", "longitude",
    "coordinate_usable_for_environment", "positional_accuracy",
}
GLOBAL_SAMPLING_SCOPE = (
    "Research-grade global Cirsium image-inference sample: species-rank, non-captive, "
    "public-coordinate iNaturalist observations; one primary photo per observation; "
    "species and 5-degree spatial-block caps defined upstream. Trait rows include only "
    "observations with at least one bootstrap-detector visible-capitulum candidate. "
    "This is not a random global occurrence sample."
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract CHELSA climate predictors for automated global AI trait tables.")
    parser.add_argument("--observation-long", required=True)
    parser.add_argument("--observation-wide", required=True)
    parser.add_argument("--registry", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--min-complete-fraction", type=float, default=0.98)
    parser.add_argument("--sample-batch-size", type=int, default=512)
    return parser.parse_args()


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def as_bool(value: Any) -> bool:
    return text(value).lower() in {"true", "1", "yes"}


def require_columns(frame: pd.DataFrame, required: set[str], label: str) -> None:
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"{label} missing required columns: {missing}")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def json_safe(value: Any) -> Any:
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    return str(value)


def coordinate_precision_tier(value: Any) -> str:
    try:
        numeric = float(text(value))
    except ValueError:
        return "unknown"
    if not math.isfinite(numeric) or numeric < 0:
        return "unknown"
    if numeric <= 1_000:
        return "high_le_1km"
    if numeric <= 10_000:
        return "moderate_1_to_10km"
    return "low_gt_10km"


def spatial_block(latitude: Any, longitude: Any, degrees: int) -> str:
    try:
        lat, lon = float(latitude), float(longitude)
    except (TypeError, ValueError):
        return ""
    if not (math.isfinite(lat) and math.isfinite(lon) and -90 <= lat <= 90 and -180 <= lon <= 180):
        return ""
    return f"lat{math.floor((lat + 90) / degrees):03d}_lon{math.floor((lon + 180) / degrees):03d}"


def open_remote_raster(url: str):
    """Open source GeoTIFF with GDAL range requests and conservative retry settings."""
    remote_path = url if url.startswith("/vsicurl/") else f"/vsicurl/{url}"
    environment = rasterio.Env(
        GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
        CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif,.tiff",
        GDAL_HTTP_MULTIRANGE="YES",
        GDAL_HTTP_TIMEOUT="120",
        GDAL_HTTP_MAX_RETRY="4",
        GDAL_HTTP_RETRY_DELAY="3",
        VSI_CACHE="TRUE",
        VSI_CACHE_SIZE="50000000",
    )
    environment.__enter__()
    try:
        dataset = rasterio.open(remote_path)
    except Exception:
        environment.__exit__(None, None, None)
        raise
    return environment, dataset


def sample_predictor(url: str, lons: np.ndarray, lats: np.ndarray, batch_size: int) -> tuple[np.ndarray, dict[str, Any]]:
    environment, src = open_remote_raster(url)
    try:
        if src.count != 1:
            raise ValueError(f"Expected one-band raster, found {src.count}: {url}")
        transformer = Transformer.from_crs("EPSG:4326", src.crs, always_xy=True)
        xs, ys = transformer.transform(lons.tolist(), lats.tolist())
        coords = list(zip(xs, ys))
        values = np.full(len(coords), np.nan, dtype=float)
        for start in range(0, len(coords), batch_size):
            stop = min(start + batch_size, len(coords))
            samples = src.sample(coords[start:stop], indexes=1, masked=True)
            for offset, sample in enumerate(samples):
                item = sample[0]
                if np.ma.is_masked(item):
                    continue
                try:
                    numeric = float(item)
                except (TypeError, ValueError):
                    continue
                if math.isfinite(numeric):
                    values[start + offset] = numeric
        metadata = {
            "url": url,
            "driver": src.driver,
            "crs": str(src.crs),
            "width": int(src.width),
            "height": int(src.height),
            "dtype": str(src.dtypes[0]),
            "nodata": json_safe(src.nodata),
            "scales": json_safe(src.scales),
            "offsets": json_safe(src.offsets),
            "transform": json_safe(tuple(src.transform)),
            "dataset_tags": json_safe(src.tags()),
            "band_1_tags": json_safe(src.tags(1)),
        }
        return values, metadata
    finally:
        src.close()
        environment.__exit__(None, None, None)


def numeric_summary(values: pd.Series) -> dict[str, Any]:
    numeric = pd.to_numeric(values, errors="coerce")
    finite = numeric[np.isfinite(numeric)]
    if finite.empty:
        return {"n_nonmissing": 0, "min": None, "p01": None, "median": None, "p99": None, "max": None}
    return {
        "n_nonmissing": int(len(finite)),
        "min": float(finite.min()),
        "p01": float(finite.quantile(0.01)),
        "median": float(finite.median()),
        "p99": float(finite.quantile(0.99)),
        "max": float(finite.max()),
    }


def build_feasibility(long: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for mode, state_column in [
        ("all_ensemble_winner", "observation_ai_all_state"),
        ("conservative_strict_consensus", "observation_ai_conservative_state"),
    ]:
        work = long.loc[long["environment_primary_complete"].eq(True)].copy()
        work["trait_state"] = work[state_column].map(text)
        work = work.loc[work["trait_state"].ne("")].copy()
        for (trait_id, state), group in work.groupby(["trait_id", "trait_state"], sort=True):
            rows.append({
                "measurement_mode": mode,
                "trait_id": trait_id,
                "trait_state": state,
                "n_observations_environment_complete": int(group["obs_id"].nunique()),
                "n_species_environment_complete": int(group["taxon_name"].nunique()),
                "n_spatial_blocks_5deg": int(group["spatial_block_5deg"].nunique()),
                "n_spatial_blocks_10deg": int(group["spatial_block_10deg"].nunique()),
                "n_high_precision_le_1km": int(group.loc[group["coordinate_precision_tier"].eq("high_le_1km"), "obs_id"].nunique()),
                "n_moderate_precision_le_10km": int(group.loc[group["coordinate_precision_tier"].isin(["high_le_1km", "moderate_1_to_10km"]), "obs_id"].nunique()),
                "interpretation": "Descriptive feasibility only. No state is automatically retained, collapsed, or excluded by this table.",
            })
    return pd.DataFrame(rows).sort_values(["measurement_mode", "trait_id", "n_observations_environment_complete"], ascending=[True, True, False]).reset_index(drop=True)


def main() -> None:
    args = parse_args()
    if not 0 < args.min_complete_fraction <= 1:
        raise ValueError("--min-complete-fraction must be in (0, 1]")
    if args.sample_batch_size < 1:
        raise ValueError("--sample-batch-size must be positive")

    long_path = Path(args.observation_long)
    wide_path = Path(args.observation_wide)
    registry_path = Path(args.registry)
    for path in (long_path, wide_path, registry_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    long = pd.read_csv(long_path, dtype=str, keep_default_na=False)
    wide = pd.read_csv(wide_path, dtype=str, keep_default_na=False)
    require_columns(long, REQUIRED_LONG, "observation-long table")
    require_columns(wide, REQUIRED_WIDE, "observation-wide table")
    if long.duplicated(["obs_id", "trait_id"]).any():
        raise ValueError("observation-long table has duplicate obs_id/trait_id rows")
    if wide["obs_id"].duplicated().any():
        raise ValueError("observation-wide table has duplicate obs_id rows")
    if set(long["obs_id"].map(text)) != set(wide["obs_id"].map(text)):
        raise ValueError("observation-long and observation-wide tables do not cover the same obs_id set")

    registry = json.loads(registry_path.read_text(encoding="utf-8"))
    predictors = registry.get("predictors", [])
    if not predictors or not all(item.get("id") and item.get("url") for item in predictors):
        raise ValueError("Predictor registry must provide id and url for every predictor")
    ids = [str(item["id"]) for item in predictors]
    if len(ids) != len(set(ids)):
        raise ValueError("Predictor registry contains duplicate IDs")

    environment = wide[["obs_id", "taxon_name", "latitude", "longitude", "coordinate_usable_for_environment", "positional_accuracy"]].copy()
    environment["latitude_numeric"] = pd.to_numeric(environment["latitude"], errors="coerce")
    environment["longitude_numeric"] = pd.to_numeric(environment["longitude"], errors="coerce")
    environment["coordinate_usable_bool"] = environment["coordinate_usable_for_environment"].map(as_bool)
    coordinate_valid = (
        environment["coordinate_usable_bool"]
        & environment["latitude_numeric"].between(-90, 90, inclusive="both")
        & environment["longitude_numeric"].between(-180, 180, inclusive="both")
    )
    if not coordinate_valid.any():
        raise ValueError("No valid coordinates are available for CHELSA extraction")
    environment["environment_coordinate_eligible"] = coordinate_valid
    environment["coordinate_precision_tier"] = environment["positional_accuracy"].map(coordinate_precision_tier)
    environment["spatial_block_5deg"] = [spatial_block(lat, lon, 5) for lat, lon in zip(environment["latitude_numeric"], environment["longitude_numeric"])]
    environment["spatial_block_10deg"] = [spatial_block(lat, lon, 10) for lat, lon in zip(environment["latitude_numeric"], environment["longitude_numeric"])]

    eligible_indices = environment.index[coordinate_valid].to_numpy()
    lons = environment.loc[eligible_indices, "longitude_numeric"].to_numpy(dtype=float)
    lats = environment.loc[eligible_indices, "latitude_numeric"].to_numpy(dtype=float)
    raster_metadata: dict[str, Any] = {}
    qc_rows: list[dict[str, Any]] = []
    for predictor in predictors:
        predictor_id = str(predictor["id"])
        values, metadata = sample_predictor(str(predictor["url"]), lons, lats, args.sample_batch_size)
        column = f"env_{predictor_id}_native"
        environment[column] = np.nan
        environment.loc[eligible_indices, column] = values
        raster_metadata[predictor_id] = metadata
        summary = numeric_summary(environment.loc[coordinate_valid, column])
        complete_fraction = summary["n_nonmissing"] / int(coordinate_valid.sum())
        qc_rows.append({
            "predictor_id": predictor_id,
            "ecological_axis": predictor.get("ecological_axis", ""),
            "source_url": predictor["url"],
            "n_coordinate_eligible": int(coordinate_valid.sum()),
            "complete_fraction_among_coordinate_eligible": complete_fraction,
            **summary,
        })
        if complete_fraction < args.min_complete_fraction:
            raise RuntimeError(
                f"CHELSA extraction coverage for {predictor_id} is {complete_fraction:.3f}, below --min-complete-fraction={args.min_complete_fraction:.3f}"
            )

    predictor_columns = [f"env_{item['id']}_native" for item in predictors]
    environment["environment_primary_complete"] = environment[predictor_columns].notna().all(axis=1)
    environment["environment_dataset"] = registry.get("dataset", {}).get("name", "CHELSA")
    environment["environment_baseline_period"] = registry.get("dataset", {}).get("baseline_period", "")
    environment["sampling_scope"] = GLOBAL_SAMPLING_SCOPE
    environment = environment.drop(columns=["coordinate_usable_bool"])

    long_enriched = long.drop(columns=["sampling_scope"], errors="ignore").merge(
        environment.drop(columns=["taxon_name"], errors="ignore"), on="obs_id", how="left", validate="many_to_one"
    )
    long_enriched["sampling_scope"] = GLOBAL_SAMPLING_SCOPE
    wide_enriched = wide.drop(columns=["sampling_scope"], errors="ignore").merge(
        environment.drop(columns=["taxon_name", "latitude", "longitude", "coordinate_usable_for_environment", "positional_accuracy"], errors="ignore"),
        on="obs_id", how="left", validate="one_to_one"
    )
    wide_enriched["sampling_scope"] = GLOBAL_SAMPLING_SCOPE

    if long_enriched["environment_primary_complete"].isna().any() or wide_enriched["environment_primary_complete"].isna().any():
        raise ValueError("Environment join did not cover every observation")
    long_enriched["environment_primary_complete"] = long_enriched["environment_primary_complete"].map(as_bool)
    wide_enriched["environment_primary_complete"] = wide_enriched["environment_primary_complete"].map(as_bool)

    feasibility = build_feasibility(long_enriched)
    correlation = environment.loc[environment["environment_primary_complete"], predictor_columns].corr(method="spearman", min_periods=20)
    correlation.index.name = "predictor"
    correlation = correlation.reset_index()

    output = Path(args.out_dir)
    output.mkdir(parents=True, exist_ok=True)
    environment.to_csv(output / "global_ai_observation_environment_metadata.csv", index=False, encoding="utf-8-sig")
    long_enriched.to_csv(output / "global_ai_trait_observation_environment_long.csv", index=False, encoding="utf-8-sig")
    wide_enriched.to_csv(output / "global_ai_trait_observation_environment_wide.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(qc_rows).to_csv(output / "chelsa_environment_extraction_qc.csv", index=False, encoding="utf-8-sig")
    correlation.to_csv(output / "chelsa_primary_predictor_spearman_correlation.csv", index=False, encoding="utf-8-sig")
    feasibility.to_csv(output / "trait_outcome_feasibility_after_environment_join.csv", index=False, encoding="utf-8-sig")

    provenance = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "semantic_status": "Environmental covariates joined to fully automated model-derived trait measurements. This artifact supports association analysis only and does not establish adaptive causality.",
        "sampling_scope": GLOBAL_SAMPLING_SCOPE,
        "n_observations_input": int(len(wide)),
        "n_coordinate_eligible": int(coordinate_valid.sum()),
        "n_environment_primary_complete": int(environment["environment_primary_complete"].sum()),
        "predictor_registry_version": registry.get("registry_version", ""),
        "predictor_ids": ids,
        "source_native_values": True,
        "raster_metadata": raster_metadata,
        "input_sha256": {
            "observation_long": sha256_file(long_path),
            "observation_wide": sha256_file(wide_path),
            "registry": sha256_file(registry_path),
        },
        "spatial_policy": registry.get("spatial_policy", {}),
        "coordinate_policy": registry.get("coordinate_policy", {}),
        "outcome_policy": "Trait states are retained as provided by all-ensemble and conservative measurement fields. No binary collapse or outcome selection occurs in environmental extraction.",
    }
    (output / "environment_analysis_provenance.json").write_text(json.dumps(provenance, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({
        "n_observations_input": provenance["n_observations_input"],
        "n_coordinate_eligible": provenance["n_coordinate_eligible"],
        "n_environment_primary_complete": provenance["n_environment_primary_complete"],
        "predictors": ids,
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
