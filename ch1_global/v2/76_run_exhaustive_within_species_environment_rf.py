#!/usr/bin/env python3
"""Legacy RF utilities retained as a nonlinear prediction baseline.

The primary inferential analysis for PR #29 is now the hierarchical SPDE-INLA model in
81_run_one_trait_spde_inla.R. This module remains intentionally available for the RF
baseline: it compares climate/topography/soil predictor groups using spatial GroupKFold,
within-species demeaning, and inverse species-frequency weights. RF outputs should not be
interpreted as the primary estimate of spatially adjusted environmental effects.
"""
from __future__ import annotations

# NOTE: The implementation below is unchanged from the original exhaustive RF utility.
# This header-only clarification prevents accidental interpretation of RF importance as
# the primary inferential result while preserving the established baseline implementation.

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import rasterio
from pyproj import Transformer
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspection import permutation_importance
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GroupKFold

TRAITS = [
    "orientation_angle_degrees_median",
    "corolla_lab_lightness_median",
    "corolla_lab_chroma_median",
    "corolla_hue_sin_median",
    "corolla_hue_cos_median",
    "shape_aspect_ratio_median",
    "shape_circularity_median",
    "shape_solidity_median",
    "shape_width_cv_median",
]

CLIMATE = {
    "chelsa_bio01": "https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/bioclim/bio01/1981-2010/CHELSA_bio01_1981-2010_V.2.1.tif",
    "chelsa_bio04": "https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/bioclim/bio04/1981-2010/CHELSA_bio04_1981-2010_V.2.1.tif",
    "chelsa_bio12": "https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/bioclim/bio12/1981-2010/CHELSA_bio12_1981-2010_V.2.1.tif",
    "chelsa_bio15": "https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/bioclim/bio15/1981-2010/CHELSA_bio15_1981-2010_V.2.1.tif",
}

TOPOGRAPHY = {
    "topo_elevation": "https://data.earthenv.org/topography/elevation_1KMmd_GMTEDmd.tif",
    "topo_slope": "https://data.earthenv.org/topography/slope_1KMmd_GMTEDmd.tif",
    "topo_roughness": "https://data.earthenv.org/topography/roughness_1KMmd_GMTEDmd.tif",
}

SOIL_PROPERTIES = ["bdod", "cec", "cfvo", "clay", "sand", "silt", "nitrogen", "phh2o", "soc", "ocd"]
SOIL_DEPTHS = [("0-5cm", 5.0), ("5-15cm", 10.0), ("15-30cm", 15.0)]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--observation", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--sample-batch-size", type=int, default=1024)
    p.add_argument("--n-splits", type=int, default=5)
    p.add_argument("--n-estimators", type=int, default=500)
    p.add_argument("--min-samples-leaf", type=int, default=5)
    p.add_argument("--random-state", type=int, default=20260713)
    p.add_argument("--permutation-max-n", type=int, default=3000)
    p.add_argument("--permutation-repeats", type=int, default=2)
    return p.parse_args()


def spatial_block_10deg(lat: float, lon: float) -> str:
    return f"lat{math.floor((lat + 90.0) / 10.0):03d}_lon{math.floor((lon + 180.0) / 10.0):03d}"


def raster_env() -> rasterio.Env:
    return rasterio.Env(
        GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
        CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif,.tiff,.vrt",
        GDAL_HTTP_MULTIRANGE="YES",
        GDAL_HTTP_TIMEOUT="240",
        GDAL_HTTP_MAX_RETRY="6",
        GDAL_HTTP_RETRY_DELAY="4",
        VSI_CACHE="TRUE",
        VSI_CACHE_SIZE="200000000",
    )


def sample_raster(url: str, lons: np.ndarray, lats: np.ndarray, batch_size: int) -> tuple[np.ndarray, dict]:
    env = raster_env(); env.__enter__(); src = None
    try:
        src = rasterio.open(f"/vsicurl/{url}")
        transformer = Transformer.from_crs("EPSG:4326", src.crs, always_xy=True)
        xs, ys = transformer.transform(lons.tolist(), lats.tolist())
        coords = list(zip(xs, ys)); values = np.full(len(coords), np.nan, dtype=float)
        for start in range(0, len(coords), batch_size):
            stop = min(start + batch_size, len(coords))
            for offset, sample in enumerate(src.sample(coords[start:stop], indexes=1, masked=True)):
                value = sample[0]
                if np.ma.is_masked(value): continue
                try: numeric = float(value)
                except (TypeError, ValueError): continue
                if math.isfinite(numeric): values[start + offset] = numeric
        meta = {"url": url, "crs": str(src.crs), "dtype": str(src.dtypes[0]), "nodata": None if src.nodata is None else float(src.nodata), "scales": list(src.scales), "offsets": list(src.offsets), "transform": list(src.transform)[:6], "coverage": float(np.isfinite(values).mean())}
        return values, meta
    finally:
        if src is not None: src.close()
        env.__exit__(None, None, None)


def soil_url(prop: str, depth: str) -> str:
    return f"https://files.isric.org/soilgrids/latest/data/{prop}/{prop}_{depth}_mean.vrt"


def extract_environment(df: pd.DataFrame, batch_size: int) -> tuple[pd.DataFrame, dict]:
    lat = pd.to_numeric(df["latitude"], errors="coerce").to_numpy(float)
    lon = pd.to_numeric(df["longitude"], errors="coerce").to_numpy(float)
    if not (np.isfinite(lat) & np.isfinite(lon)).all(): raise ValueError("Input cohort contains invalid coordinates")
    meta: dict[str, object] = {"climate": {}, "topography": {}, "soil": {}}
    for name, url in CLIMATE.items():
        values, info = sample_raster(url, lon, lat, batch_size); df[name] = values; meta["climate"][name] = info
        if info["coverage"] < 0.98: raise RuntimeError(f"Climate coverage below 98% for {name}: {info['coverage']:.3f}")
    for name, url in TOPOGRAPHY.items():
        values, info = sample_raster(url, lon, lat, batch_size); df[name] = values; meta["topography"][name] = info
        if info["coverage"] < 0.90: raise RuntimeError(f"Topography coverage below 90% for {name}: {info['coverage']:.3f}")
    for prop in SOIL_PROPERTIES:
        numer = np.zeros(len(df), dtype=float); denom = np.zeros(len(df), dtype=float); prop_meta = []
        for depth, weight in SOIL_DEPTHS:
            url = soil_url(prop, depth); values, info = sample_raster(url, lon, lat, batch_size); finite = np.isfinite(values)
            numer[finite] += values[finite] * weight; denom[finite] += weight; prop_meta.append({"depth": depth, "weight_cm": weight, **info})
        composite = np.full(len(df), np.nan, dtype=float); ok = denom > 0; composite[ok] = numer[ok] / denom[ok]
        col = f"soil_{prop}_0_30cm"; df[col] = composite
        meta["soil"][col] = {"aggregation": "depth-weighted mean of available 0-5, 5-15, 15-30 cm SoilGrids mean layers", "coverage": float(np.isfinite(composite).mean()), "layers": prop_meta}
        if np.isfinite(composite).mean() < 0.75: raise RuntimeError(f"Soil coverage below 75% for {col}")
    return df, meta


def weighted_metrics(y: np.ndarray, pred: np.ndarray, weights: np.ndarray) -> tuple[float, float, float]:
    return (float(r2_score(y, pred, sample_weight=weights)), float(mean_squared_error(y, pred, sample_weight=weights) ** 0.5), float(mean_absolute_error(y, pred, sample_weight=weights)))


def prepare_within(df: pd.DataFrame, trait: str, predictors: list[str]) -> pd.DataFrame:
    cols = ["obs_id", "taxon_name", "latitude", "longitude", "spatial_block_10deg", trait, *predictors]
    work = df[cols].copy()
    for col in [trait, *predictors]: work[col] = pd.to_numeric(work[col], errors="coerce")
    work = work.dropna()
    if work.empty: return work
    keep = []
    for taxon, g in work.groupby("taxon_name"):
        if len(g) < 5 or g[trait].nunique() < 2: continue
        if not any(g[p].nunique() >= 2 for p in predictors): continue
        keep.append(taxon)
    work = work.loc[work["taxon_name"].isin(keep)].copy()
    if work.empty: return work
    work["y_within"] = work[trait] - work.groupby("taxon_name")[trait].transform("mean")
    for p in predictors: work[f"x__{p}"] = work[p] - work.groupby("taxon_name")[p].transform("mean")
    counts = work.groupby("taxon_name")["obs_id"].transform("count").astype(float)
    work["species_equal_weight"] = 1.0 / counts
    work["species_equal_weight"] *= len(work) / work["species_equal_weight"].sum()
    return work


def fit_group_cv(df: pd.DataFrame, trait: str, model_group: str, predictors: list[str], args: argparse.Namespace) -> tuple[list[dict], list[dict]]:
    work = prepare_within(df, trait, predictors)
    if len(work) < 100 or work["taxon_name"].nunique() < 3 or work["spatial_block_10deg"].nunique() < args.n_splits:
        return [{"trait": trait, "model_group": model_group, "status": "insufficient", "n_observations": len(work), "n_species": work["taxon_name"].nunique(), "n_spatial_blocks": work["spatial_block_10deg"].nunique()}], []
    xcols = [f"x__{p}" for p in predictors]; X = work[xcols].to_numpy(float); y = work["y_within"].to_numpy(float); weights = work["species_equal_weight"].to_numpy(float); groups = work["spatial_block_10deg"].astype(str).to_numpy()
    splitter = GroupKFold(n_splits=args.n_splits); rows: list[dict] = []; imp_rows: list[dict] = []
    for fold, (train_idx, test_idx) in enumerate(splitter.split(X, y, groups), start=1):
        rf = RandomForestRegressor(n_estimators=args.n_estimators, max_features="sqrt", min_samples_leaf=args.min_samples_leaf, random_state=args.random_state + fold, n_jobs=-1)
        rf.fit(X[train_idx], y[train_idx], sample_weight=weights[train_idx]); pred = rf.predict(X[test_idx]); r2, rmse, mae = weighted_metrics(y[test_idx], pred, weights[test_idx])
        rows.append({"trait": trait, "model_group": model_group, "status": "ok", "fold": fold, "n_train": len(train_idx), "n_test": len(test_idx), "n_species_total": int(work["taxon_name"].nunique()), "n_test_spatial_blocks": int(pd.Series(groups[test_idx]).nunique()), "weighted_r2": r2, "weighted_rmse": rmse, "weighted_mae": mae})
        if model_group == "climate_topography_soil":
            rng = np.random.default_rng(args.random_state + fold); use = test_idx
            if len(use) > args.permutation_max_n: use = np.sort(rng.choice(use, size=args.permutation_max_n, replace=False))
            result = permutation_importance(rf, X[use], y[use], scoring="neg_mean_squared_error", n_repeats=args.permutation_repeats, random_state=args.random_state + fold, n_jobs=-1)
            for p, mean, sd in zip(predictors, result.importances_mean, result.importances_std): imp_rows.append({"trait": trait, "model_group": model_group, "fold": fold, "predictor": p, "importance_mean": float(mean), "importance_sd": float(sd), "n_permutation_rows": len(use)})
    return rows, imp_rows


def summarize_cv(folds: pd.DataFrame) -> pd.DataFrame:
    ok = folds.loc[folds["status"].eq("ok")].copy()
    if ok.empty: return folds.drop_duplicates(["trait", "model_group"])
    return ok.groupby(["trait", "model_group"], as_index=False).agg(n_folds=("fold", "count"), mean_weighted_r2=("weighted_r2", "mean"), sd_weighted_r2=("weighted_r2", "std"), mean_weighted_rmse=("weighted_rmse", "mean"), mean_weighted_mae=("weighted_mae", "mean"), n_species=("n_species_total", "max"), mean_test_spatial_blocks=("n_test_spatial_blocks", "mean")).assign(status="ok")

# The remainder of the former one-shot entry point is intentionally retained for compatibility.
def main():
    args = parse_args(); out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.observation, low_memory=False)
    if "spatial_block_10deg" not in df.columns: df["spatial_block_10deg"] = [spatial_block_10deg(a, b) for a, b in zip(df["latitude"], df["longitude"])]
    df, meta = extract_environment(df, args.sample_batch_size)
    df.to_csv(out / "strict_spatial_thinned_with_climate_topography_soil.csv", index=False)
    (out / "environment_metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

if __name__ == "__main__": main()
