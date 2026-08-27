#!/usr/bin/env python3
"""Attach frozen CHELSA process-extension variables to an existing cohort.

Source identities and aggregation rules are read from a frozen JSON contract.
Every source must meet the frozen coverage threshold; otherwise extraction fails
rather than silently substituting a post-hoc predictor. A source may be a single
COG or a predeclared list of monthly climatology COGs aggregated by arithmetic
mean at each frozen observation coordinate.
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


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--environment", required=True, type=Path)
    p.add_argument("--source-contract", required=True, type=Path)
    p.add_argument("--out-csv", required=True, type=Path)
    p.add_argument("--report", required=True, type=Path)
    p.add_argument("--sample-batch-size", type=int, default=1024)
    return p.parse_args()


def sample_raster(url: str, lon: np.ndarray, lat: np.ndarray, batch_size: int) -> tuple[np.ndarray, dict[str, Any]]:
    print(f"Sampling CHELSA source: {url}", flush=True)
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
        with rasterio.open(f"/vsicurl/{url}") as src:
            transformer = Transformer.from_crs("EPSG:4326", src.crs, always_xy=True)
            xs, ys = transformer.transform(lon.tolist(), lat.tolist())
            coords = list(zip(xs, ys))
            vals = np.full(len(coords), np.nan, dtype=float)
            for start in range(0, len(coords), batch_size):
                stop = min(start + batch_size, len(coords))
                for offset, sample in enumerate(src.sample(coords[start:stop], indexes=1, masked=True)):
                    v = sample[0]
                    if np.ma.is_masked(v):
                        continue
                    v = float(v)
                    if math.isfinite(v):
                        vals[start + offset] = v
            meta = {
                "url": url,
                "crs": str(src.crs),
                "dtype": str(src.dtypes[0]),
                "nodata": None if src.nodata is None else float(src.nodata),
                "scales": [float(x) for x in src.scales],
                "offsets": [float(x) for x in src.offsets],
            }
    return vals, meta


def sample_source(
    source: dict[str, Any],
    lon: np.ndarray,
    lat: np.ndarray,
    batch_size: int,
    minimum_coverage: float,
) -> tuple[np.ndarray, dict[str, Any]]:
    if "url" in source and "urls" in source:
        raise ValueError(f"Source {source['column']} cannot define both url and urls")
    if "url" in source:
        values, meta = sample_raster(source["url"], lon, lat, batch_size)
        cov = float(np.isfinite(values).mean())
        if cov < minimum_coverage:
            raise RuntimeError(
                f"{source['column']} source coverage {cov:.4f} below frozen minimum {minimum_coverage:.4f}"
            )
        return values, {"aggregation": "single_cog", "component_coverage": [cov], "components": [meta]}

    urls = source.get("urls")
    if not isinstance(urls, list) or not urls:
        raise ValueError(f"Source {source['column']} must define url or non-empty urls")
    aggregation = source.get("aggregation")
    if aggregation != "mean":
        raise ValueError(f"Source {source['column']} multi-COG aggregation must be frozen as mean")
    matrices: list[np.ndarray] = []
    component_meta: list[dict[str, Any]] = []
    component_coverage: list[float] = []
    for url in urls:
        values, meta = sample_raster(str(url), lon, lat, batch_size)
        cov = float(np.isfinite(values).mean())
        if cov < minimum_coverage:
            raise RuntimeError(
                f"{source['column']} component coverage {cov:.4f} below frozen minimum {minimum_coverage:.4f}: {url}"
            )
        matrices.append(values)
        component_meta.append(meta)
        component_coverage.append(cov)
    matrix = np.vstack(matrices)
    with np.errstate(invalid="ignore"):
        aggregated = np.nanmean(matrix, axis=0)
    return aggregated, {
        "aggregation": "arithmetic_mean_across_predeclared_monthly_climatologies",
        "n_components": len(urls),
        "component_coverage": component_coverage,
        "components": component_meta,
    }


def main() -> int:
    args = parse_args()
    env = pd.read_csv(args.environment, low_memory=False)
    contract = json.loads(args.source_contract.read_text())
    required = {"obs_id", "taxon_name", "latitude", "longitude"}
    missing = required.difference(env.columns)
    if missing:
        raise ValueError(f"Environment table missing columns: {sorted(missing)}")
    if env["obs_id"].astype(str).duplicated().any():
        raise ValueError("Environment table must be unique by obs_id")
    lon = pd.to_numeric(env["longitude"], errors="coerce").to_numpy(float)
    lat = pd.to_numeric(env["latitude"], errors="coerce").to_numpy(float)
    if not (np.isfinite(lon) & np.isfinite(lat)).all():
        raise ValueError("Invalid coordinates")
    min_cov = float(contract.get("minimum_coverage", 0.98))
    metadata, coverage = {}, {}
    out = env.copy()
    for source in contract["sources"]:
        column = source["column"]
        values, meta = sample_source(source, lon, lat, args.sample_batch_size, min_cov)
        out[column] = values
        coverage[column] = float(np.isfinite(values).mean())
        metadata[column] = meta
        if coverage[column] < min_cov:
            raise RuntimeError(f"{column} coverage {coverage[column]:.4f} below frozen minimum {min_cov:.4f}")
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out_csv, index=False)
    report = {
        "n_observations": int(len(out)),
        "n_taxa": int(out["taxon_name"].nunique()),
        "contract": contract,
        "coverage": coverage,
        "raster_metadata": metadata,
        "selection_rule": "phenotype-blind extraction on the frozen GEB-v2 strict spatial cohort",
        "aggregation_rule": "only aggregation declared before outcome inspection in the source contract is permitted",
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
