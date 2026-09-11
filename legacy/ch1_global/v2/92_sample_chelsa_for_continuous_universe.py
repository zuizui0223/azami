#!/usr/bin/env python3
"""Attach predeclared CHELSA predictors without running a trait model.

Environmental extraction is deliberately separated from inference so every
registered continuous endpoint reaches the same observation-level predictor
table.  The input cohort is never filtered on a measured phenotype.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import rasterio
from pyproj import Transformer


PREDICTORS = {
    "chelsa_bio01": "https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/bioclim/bio01/1981-2010/CHELSA_bio01_1981-2010_V.2.1.tif",
    "chelsa_bio04": "https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/bioclim/bio04/1981-2010/CHELSA_bio04_1981-2010_V.2.1.tif",
    "chelsa_bio12": "https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/bioclim/bio12/1981-2010/CHELSA_bio12_1981-2010_V.2.1.tif",
    "chelsa_bio15": "https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/bioclim/bio15/1981-2010/CHELSA_bio15_1981-2010_V.2.1.tif",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--observation", required=True)
    parser.add_argument("--out-csv", required=True)
    parser.add_argument("--report", required=True)
    parser.add_argument("--sample-batch-size", type=int, default=1024)
    parser.add_argument("--minimum-coverage", type=float, default=0.98)
    return parser.parse_args()


def sample_raster(
    url: str,
    longitudes: np.ndarray,
    latitudes: np.ndarray,
    batch_size: int,
) -> tuple[np.ndarray, dict[str, Any]]:
    with rasterio.Env(
        GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
        CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif,.tiff",
        GDAL_HTTP_MULTIRANGE="YES",
        GDAL_HTTP_TIMEOUT="180",
        GDAL_HTTP_MAX_RETRY="5",
        GDAL_HTTP_RETRY_DELAY="3",
        VSI_CACHE="TRUE",
        VSI_CACHE_SIZE="100000000",
    ):
        with rasterio.open(f"/vsicurl/{url}") as source:
            transformer = Transformer.from_crs("EPSG:4326", source.crs, always_xy=True)
            xs, ys = transformer.transform(longitudes.tolist(), latitudes.tolist())
            coordinates = list(zip(xs, ys))
            values = np.full(len(coordinates), np.nan, dtype=float)
            for start in range(0, len(coordinates), batch_size):
                stop = min(start + batch_size, len(coordinates))
                for offset, sample in enumerate(
                    source.sample(coordinates[start:stop], indexes=1, masked=True)
                ):
                    value = sample[0]
                    if np.ma.is_masked(value):
                        continue
                    value = float(value)
                    if math.isfinite(value):
                        values[start + offset] = value
            metadata = {
                "url": url,
                "crs": str(source.crs),
                "dtype": str(source.dtypes[0]),
                "nodata": None if source.nodata is None else float(source.nodata),
                "scales": [float(value) for value in source.scales],
                "offsets": [float(value) for value in source.offsets],
                "transform": [float(value) for value in list(source.transform)[:6]],
            }
    return values, metadata


def main() -> None:
    args = parse_args()
    if args.sample_batch_size < 1:
        raise ValueError("--sample-batch-size must be positive")
    if not 0 < args.minimum_coverage <= 1:
        raise ValueError("--minimum-coverage must be in (0,1]")
    observation = pd.read_csv(args.observation, low_memory=False)
    required = {"obs_id", "taxon_name", "latitude", "longitude"}
    missing = sorted(required.difference(observation.columns))
    if missing:
        raise ValueError(f"Observation table is missing columns: {missing}")
    if observation["obs_id"].astype(str).duplicated().any():
        raise ValueError("Observation table must be unique by obs_id")
    latitude = pd.to_numeric(observation["latitude"], errors="coerce").to_numpy(float)
    longitude = pd.to_numeric(observation["longitude"], errors="coerce").to_numpy(float)
    if not (np.isfinite(latitude) & np.isfinite(longitude)).all():
        raise ValueError("Observation table contains invalid coordinates")

    environment = observation[["obs_id", "taxon_name", "latitude", "longitude"]].copy()
    raster_metadata: dict[str, Any] = {}
    coverage: dict[str, float] = {}
    for predictor, url in PREDICTORS.items():
        values, metadata = sample_raster(
            url, longitude, latitude, args.sample_batch_size
        )
        environment[predictor] = values
        raster_metadata[predictor] = metadata
        coverage[predictor] = float(np.isfinite(values).mean())
        if coverage[predictor] < args.minimum_coverage:
            raise RuntimeError(
                f"{predictor} extraction coverage {coverage[predictor]:.4f} "
                f"is below {args.minimum_coverage:.4f}"
            )

    output = Path(args.out_csv)
    output.parent.mkdir(parents=True, exist_ok=True)
    environment.to_csv(output, index=False, encoding="utf-8-sig")
    report = {
        "n_observations": int(len(environment)),
        "n_taxa": int(environment["taxon_name"].nunique()),
        "predictor_coverage": coverage,
        "minimum_required_coverage": float(args.minimum_coverage),
        "raster_metadata": raster_metadata,
        "selection_rule": "No phenotype-based filtering; predictors attached to the frozen strict spatial cohort.",
    }
    report_path = Path(args.report)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
