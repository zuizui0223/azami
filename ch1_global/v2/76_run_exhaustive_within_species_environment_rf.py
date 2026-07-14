#!/usr/bin/env python3
"""RF utilities used by the nonlinear prediction baseline in PR #29.

The primary inferential analysis is the hierarchical SPDE-INLA model in
81_run_one_trait_spde_inla.R. RF is retained for nonlinear predictive comparison only.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from pyproj import Transformer
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspection import permutation_importance
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GroupKFold

TRAITS = [
    "orientation_angle_degrees_median", "corolla_lab_lightness_median",
    "corolla_lab_chroma_median", "corolla_hue_sin_median", "corolla_hue_cos_median",
    "shape_aspect_ratio_median", "shape_circularity_median", "shape_solidity_median",
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


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--observation", required=True); p.add_argument("--out-dir", required=True)
    p.add_argument("--sample-batch-size", type=int, default=1024); p.add_argument("--n-splits", type=int, default=5)
    p.add_argument("--n-estimators", type=int, default=500); p.add_argument("--min-samples-leaf", type=int, default=5)
    p.add_argument("--random-state", type=int, default=20260713); p.add_argument("--permutation-max-n", type=int, default=3000)
    p.add_argument("--permutation-repeats", type=int, default=2); return p.parse_args()


def spatial_block_10deg(lat, lon):
    return f"lat{math.floor((lat + 90.0) / 10.0):03d}_lon{math.floor((lon + 180.0) / 10.0):03d}"


def raster_env():
    return rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR", CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif,.tiff,.vrt", GDAL_HTTP_MULTIRANGE="YES", GDAL_HTTP_TIMEOUT="240", GDAL_HTTP_MAX_RETRY="6", GDAL_HTTP_RETRY_DELAY="4", VSI_CACHE="TRUE", VSI_CACHE_SIZE="200000000")


def sample_raster(url, lons, lats, batch_size):
    env = raster_env(); env.__enter__(); src = None
    try:
        src = rasterio.open(f"/vsicurl/{url}"); transformer = Transformer.from_crs("EPSG:4326", src.crs, always_xy=True)
        xs, ys = transformer.transform(lons.tolist(), lats.tolist()); coords = list(zip(xs, ys)); values = np.full(len(coords), np.nan)
        for start in range(0, len(coords), batch_size):
            stop = min(start + batch_size, len(coords))
            for offset, sample in enumerate(src.sample(coords[start:stop], indexes=1, masked=True)):
                value = sample[0]
                if np.ma.is_masked(value): continue
                try: numeric = float(value)
                except (TypeError, ValueError): continue
                if math.isfinite(numeric): values[start + offset] = numeric
        return values, {"url": url, "crs": str(src.crs), "coverage": float(np.isfinite(values).mean())}
    finally:
        if src is not None: src.close()
        env.__exit__(None, None, None)


def soil_url(prop, depth):
    return f"https://files.isric.org/soilgrids/latest/data/{prop}/{prop}_{depth}_mean.vrt"


def extract_environment(df, batch_size):
    lat = pd.to_numeric(df["latitude"], errors="coerce").to_numpy(float); lon = pd.to_numeric(df["longitude"], errors="coerce").to_numpy(float)
    if not (np.isfinite(lat) & np.isfinite(lon)).all(): raise ValueError("Input cohort contains invalid coordinates")
    meta = {"climate": {}, "topography": {}, "soil": {}}
    for name, url in CLIMATE.items():
        values, info = sample_raster(url, lon, lat, batch_size); df[name] = values; meta["climate"][name] = info
    for name, url in TOPOGRAPHY.items():
        values, info = sample_raster(url, lon, lat, batch_size); df[name] = values; meta["topography"][name] = info
    for prop in SOIL_PROPERTIES:
        numer = np.zeros(len(df)); denom = np.zeros(len(df)); prop_meta = []
        for depth, weight in SOIL_DEPTHS:
            values, info = sample_raster(soil_url(prop, depth), lon, lat, batch_size); finite = np.isfinite(values)
            numer[finite] += values[finite] * weight; denom[finite] += weight; prop_meta.append({"depth": depth, "weight_cm": weight, **info})
        composite = np.full(len(df), np.nan); ok = denom > 0; composite[ok] = numer[ok] / denom[ok]
        col = f"soil_{prop}_0_30cm"; df[col] = composite; meta["soil"][col] = {"coverage": float(np.isfinite(composite).mean()), "layers": prop_meta}
    return df, meta


def weighted_metrics(y, pred, weights):
    return float(r2_score(y, pred, sample_weight=weights)), float(mean_squared_error(y, pred, sample_weight=weights) ** 0.5), float(mean_absolute_error(y, pred, sample_weight=weights))


def prepare_within(df, trait, predictors):
    cols = ["obs_id", "taxon_name", "latitude", "longitude", "spatial_block_10deg", trait, *predictors]; work = df[cols].copy()
    for col in [trait, *predictors]: work[col] = pd.to_numeric(work[col], errors="coerce")
    work = work.dropna()
    keep = [taxon for taxon, g in work.groupby("taxon_name") if len(g) >= 5 and g[trait].nunique() >= 2 and any(g[p].nunique() >= 2 for p in predictors)]
    work = work.loc[work["taxon_name"].isin(keep)].copy()
    if work.empty: return work
    work["y_within"] = work[trait] - work.groupby("taxon_name")[trait].transform("mean")
    for p in predictors: work[f"x__{p}"] = work[p] - work.groupby("taxon_name")[p].transform("mean")
    counts = work.groupby("taxon_name")["obs_id"].transform("count").astype(float); work["species_equal_weight"] = 1.0 / counts; work["species_equal_weight"] *= len(work) / work["species_equal_weight"].sum()
    return work


def fit_group_cv(df, trait, model_group, predictors, args):
    work = prepare_within(df, trait, predictors)
    if len(work) < 100 or work["taxon_name"].nunique() < 3 or work["spatial_block_10deg"].nunique() < args.n_splits:
        return [{"trait": trait, "model_group": model_group, "status": "insufficient", "n_observations": len(work), "n_species": work["taxon_name"].nunique(), "n_spatial_blocks": work["spatial_block_10deg"].nunique()}], []
    xcols = [f"x__{p}" for p in predictors]; X = work[xcols].to_numpy(float); y = work["y_within"].to_numpy(float); weights = work["species_equal_weight"].to_numpy(float); groups = work["spatial_block_10deg"].astype(str).to_numpy()
    rows = []; imp_rows = []; splitter = GroupKFold(n_splits=args.n_splits)
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


def summarize_cv(folds):
    ok = folds.loc[folds["status"].eq("ok")].copy()
    if ok.empty: return folds.drop_duplicates(["trait", "model_group"])
    return ok.groupby(["trait", "model_group"], as_index=False).agg(n_folds=("fold", "count"), mean_weighted_r2=("weighted_r2", "mean"), sd_weighted_r2=("weighted_r2", "std"), mean_weighted_rmse=("weighted_rmse", "mean"), mean_weighted_mae=("weighted_mae", "mean"), n_species=("n_species_total", "max"), mean_test_spatial_blocks=("n_test_spatial_blocks", "mean")).assign(status="ok")


def main():
    args = parse_args(); out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True); df = pd.read_csv(args.observation, low_memory=False)
    if "spatial_block_10deg" not in df.columns: df["spatial_block_10deg"] = [spatial_block_10deg(a, b) for a, b in zip(df["latitude"], df["longitude"])]
    df, meta = extract_environment(df, args.sample_batch_size); df.to_csv(out / "strict_spatial_thinned_with_climate_topography_soil.csv", index=False); (out / "environment_metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

if __name__ == "__main__": main()
