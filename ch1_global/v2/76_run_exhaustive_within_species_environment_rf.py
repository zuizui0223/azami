#!/usr/bin/env python3
"""Run exhaustive within-species RF models across climate, topography and SoilGrids.

The input is the strict spatially thinned all-photo observation cohort. Environmental
predictors are sampled only after trait inference. RF models compare four predeclared
predictor groups using 10-degree spatial GroupKFold. Outcomes and predictors are
species-demeaned before fitting so the analysis targets within-species structure.
Species contribute inverse-frequency sample weights so hyper-abundant species do not
numerically dominate model fitting or scoring.
"""
from __future__ import annotations

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

# EarthEnv 1-km global topography, GMTED2010-derived layers.
TOPOGRAPHY = {
    "topo_elevation": "https://data.earthenv.org/topography/elevation_1KMmd_GMTEDmd.tif",
    "topo_slope": "https://data.earthenv.org/topography/slope_1KMmd_GMTEDmd.tif",
    "topo_roughness": "https://data.earthenv.org/topography/roughness_1KMmd_GMTEDmd.tif",
}

# SoilGrids 2.0 mean predictions. A depth-weighted 0-30 cm composite is built from
# 0-5, 5-15 and 15-30 cm layers with thickness weights 5, 10 and 15 cm.
SOIL_PROPERTIES = [
    "bdod", "cec", "cfvo", "clay", "sand", "silt", "nitrogen", "phh2o", "soc", "ocd"
]
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
    env = raster_env()
    env.__enter__()
    src = None
    try:
        src = rasterio.open(f"/vsicurl/{url}")
        transformer = Transformer.from_crs("EPSG:4326", src.crs, always_xy=True)
        xs, ys = transformer.transform(lons.tolist(), lats.tolist())
        coords = list(zip(xs, ys))
        values = np.full(len(coords), np.nan, dtype=float)
        for start in range(0, len(coords), batch_size):
            stop = min(start + batch_size, len(coords))
            for offset, sample in enumerate(src.sample(coords[start:stop], indexes=1, masked=True)):
                value = sample[0]
                if np.ma.is_masked(value):
                    continue
                try:
                    numeric = float(value)
                except (TypeError, ValueError):
                    continue
                if math.isfinite(numeric):
                    values[start + offset] = numeric
        meta = {
            "url": url,
            "crs": str(src.crs),
            "dtype": str(src.dtypes[0]),
            "nodata": None if src.nodata is None else float(src.nodata),
            "scales": list(src.scales),
            "offsets": list(src.offsets),
            "transform": list(src.transform)[:6],
            "coverage": float(np.isfinite(values).mean()),
        }
        return values, meta
    finally:
        if src is not None:
            src.close()
        env.__exit__(None, None, None)


def soil_url(prop: str, depth: str) -> str:
    return f"https://files.isric.org/soilgrids/latest/data/{prop}/{prop}_{depth}_mean.vrt"


def extract_environment(df: pd.DataFrame, batch_size: int) -> tuple[pd.DataFrame, dict]:
    lat = pd.to_numeric(df["latitude"], errors="coerce").to_numpy(float)
    lon = pd.to_numeric(df["longitude"], errors="coerce").to_numpy(float)
    if not (np.isfinite(lat) & np.isfinite(lon)).all():
        raise ValueError("Input cohort contains invalid coordinates")
    meta: dict[str, object] = {"climate": {}, "topography": {}, "soil": {}}

    for name, url in CLIMATE.items():
        values, info = sample_raster(url, lon, lat, batch_size)
        df[name] = values
        meta["climate"][name] = info
        if info["coverage"] < 0.98:
            raise RuntimeError(f"Climate coverage below 98% for {name}: {info['coverage']:.3f}")

    for name, url in TOPOGRAPHY.items():
        values, info = sample_raster(url, lon, lat, batch_size)
        df[name] = values
        meta["topography"][name] = info
        if info["coverage"] < 0.90:
            raise RuntimeError(f"Topography coverage below 90% for {name}: {info['coverage']:.3f}")

    for prop in SOIL_PROPERTIES:
        numer = np.zeros(len(df), dtype=float)
        denom = np.zeros(len(df), dtype=float)
        prop_meta = []
        for depth, weight in SOIL_DEPTHS:
            url = soil_url(prop, depth)
            values, info = sample_raster(url, lon, lat, batch_size)
            finite = np.isfinite(values)
            numer[finite] += values[finite] * weight
            denom[finite] += weight
            prop_meta.append({"depth": depth, "weight_cm": weight, **info})
        composite = np.full(len(df), np.nan, dtype=float)
        ok = denom > 0
        composite[ok] = numer[ok] / denom[ok]
        col = f"soil_{prop}_0_30cm"
        df[col] = composite
        meta["soil"][col] = {
            "aggregation": "depth-weighted mean of available 0-5, 5-15, 15-30 cm SoilGrids mean layers",
            "coverage": float(np.isfinite(composite).mean()),
            "layers": prop_meta,
        }
        if np.isfinite(composite).mean() < 0.75:
            raise RuntimeError(f"Soil coverage below 75% for {col}")
    return df, meta


def weighted_metrics(y: np.ndarray, pred: np.ndarray, weights: np.ndarray) -> tuple[float, float, float]:
    return (
        float(r2_score(y, pred, sample_weight=weights)),
        float(mean_squared_error(y, pred, sample_weight=weights) ** 0.5),
        float(mean_absolute_error(y, pred, sample_weight=weights)),
    )


def prepare_within(df: pd.DataFrame, trait: str, predictors: list[str]) -> pd.DataFrame:
    cols = ["obs_id", "taxon_name", "latitude", "longitude", "spatial_block_10deg", trait, *predictors]
    work = df[cols].copy()
    for col in [trait, *predictors]:
        work[col] = pd.to_numeric(work[col], errors="coerce")
    work = work.dropna()
    if work.empty:
        return work
    # Require within-species variation in the response and at least one predictor.
    keep = []
    for taxon, g in work.groupby("taxon_name"):
        if len(g) < 5 or g[trait].nunique() < 2:
            continue
        if not any(g[p].nunique() >= 2 for p in predictors):
            continue
        keep.append(taxon)
    work = work.loc[work["taxon_name"].isin(keep)].copy()
    if work.empty:
        return work
    work["y_within"] = work[trait] - work.groupby("taxon_name")[trait].transform("mean")
    for p in predictors:
        work[f"x__{p}"] = work[p] - work.groupby("taxon_name")[p].transform("mean")
    counts = work.groupby("taxon_name")["obs_id"].transform("count").astype(float)
    work["species_equal_weight"] = 1.0 / counts
    # Normalize weights to mean 1 for numerical stability.
    work["species_equal_weight"] *= len(work) / work["species_equal_weight"].sum()
    return work


def fit_group_cv(
    df: pd.DataFrame,
    trait: str,
    model_group: str,
    predictors: list[str],
    args: argparse.Namespace,
) -> tuple[list[dict], list[dict]]:
    work = prepare_within(df, trait, predictors)
    if len(work) < 100 or work["taxon_name"].nunique() < 3 or work["spatial_block_10deg"].nunique() < args.n_splits:
        return [{
            "trait": trait, "model_group": model_group, "status": "insufficient",
            "n_observations": len(work), "n_species": work["taxon_name"].nunique(),
            "n_spatial_blocks": work["spatial_block_10deg"].nunique(),
        }], []

    xcols = [f"x__{p}" for p in predictors]
    X = work[xcols].to_numpy(float)
    y = work["y_within"].to_numpy(float)
    weights = work["species_equal_weight"].to_numpy(float)
    groups = work["spatial_block_10deg"].astype(str).to_numpy()
    splitter = GroupKFold(n_splits=args.n_splits)
    rows: list[dict] = []
    imp_rows: list[dict] = []

    for fold, (train_idx, test_idx) in enumerate(splitter.split(X, y, groups), start=1):
        rf = RandomForestRegressor(
            n_estimators=args.n_estimators,
            max_features="sqrt",
            min_samples_leaf=args.min_samples_leaf,
            random_state=args.random_state + fold,
            n_jobs=-1,
        )
        rf.fit(X[train_idx], y[train_idx], sample_weight=weights[train_idx])
        pred = rf.predict(X[test_idx])
        r2, rmse, mae = weighted_metrics(y[test_idx], pred, weights[test_idx])
        rows.append({
            "trait": trait,
            "model_group": model_group,
            "status": "ok",
            "fold": fold,
            "n_train": len(train_idx),
            "n_test": len(test_idx),
            "n_species_total": int(work["taxon_name"].nunique()),
            "n_test_spatial_blocks": int(pd.Series(groups[test_idx]).nunique()),
            "weighted_r2": r2,
            "weighted_rmse": rmse,
            "weighted_mae": mae,
        })

        # Permutation importance only for the full environment model. Subsample held-out
        # observations deterministically to control runtime.
        if model_group == "climate_topography_soil":
            rng = np.random.default_rng(args.random_state + fold)
            use = test_idx
            if len(use) > args.permutation_max_n:
                use = np.sort(rng.choice(use, size=args.permutation_max_n, replace=False))
            result = permutation_importance(
                rf,
                X[use],
                y[use],
                scoring="neg_mean_squared_error",
                n_repeats=args.permutation_repeats,
                random_state=args.random_state + fold,
                n_jobs=-1,
            )
            for p, mean_imp, sd_imp in zip(predictors, result.importances_mean, result.importances_std):
                imp_rows.append({
                    "trait": trait,
                    "fold": fold,
                    "predictor": p,
                    "permutation_importance_mse_mean": float(mean_imp),
                    "permutation_importance_mse_sd": float(sd_imp),
                    "n_permutation_test": int(len(use)),
                })
    return rows, imp_rows


def summarize_cv(folds: pd.DataFrame) -> pd.DataFrame:
    ok = folds.loc[folds["status"].eq("ok")].copy()
    if ok.empty:
        return pd.DataFrame()
    return (
        ok.groupby(["trait", "model_group"], as_index=False)
        .agg(
            folds=("fold", "nunique"),
            n_species=("n_species_total", "max"),
            mean_weighted_r2=("weighted_r2", "mean"),
            sd_weighted_r2=("weighted_r2", "std"),
            mean_weighted_rmse=("weighted_rmse", "mean"),
            mean_weighted_mae=("weighted_mae", "mean"),
        )
    )


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.observation, low_memory=False)
    required = {"obs_id", "taxon_name", "latitude", "longitude", *TRAITS}
    missing = sorted(required.difference(df.columns))
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    if df["obs_id"].astype(str).duplicated().any():
        raise ValueError("Observation table must have unique obs_id")

    lat = pd.to_numeric(df["latitude"], errors="coerce")
    lon = pd.to_numeric(df["longitude"], errors="coerce")
    if not (lat.notna() & lon.notna()).all():
        raise ValueError("Input contains invalid coordinates")
    df["spatial_block_10deg"] = [spatial_block_10deg(a, b) for a, b in zip(lat, lon)]

    df, env_meta = extract_environment(df, args.sample_batch_size)
    env_path = out_dir / "strict_spatial_thinned_with_climate_topography_soil.csv"
    df.to_csv(env_path, index=False, encoding="utf-8-sig")

    climate = list(CLIMATE)
    topo = list(TOPOGRAPHY)
    soil = [f"soil_{p}_0_30cm" for p in SOIL_PROPERTIES]
    groups = {
        "climate": climate,
        "climate_topography": climate + topo,
        "climate_soil": climate + soil,
        "climate_topography_soil": climate + topo + soil,
    }

    fold_rows: list[dict] = []
    importance_rows: list[dict] = []
    for trait in TRAITS:
        for group_name, predictors in groups.items():
            rows, imps = fit_group_cv(df, trait, group_name, predictors, args)
            fold_rows.extend(rows)
            importance_rows.extend(imps)

    folds = pd.DataFrame(fold_rows)
    folds.to_csv(out_dir / "rf_spatial_cv_folds.csv", index=False, encoding="utf-8-sig")
    summary = summarize_cv(folds)
    summary.to_csv(out_dir / "rf_spatial_cv_summary.csv", index=False, encoding="utf-8-sig")
    importance = pd.DataFrame(importance_rows)
    importance.to_csv(out_dir / "rf_full_model_permutation_importance.csv", index=False, encoding="utf-8-sig")

    if not summary.empty:
        wide = summary.pivot(index="trait", columns="model_group", values="mean_weighted_r2")
        ablation_rows = []
        for trait, row in wide.iterrows():
            base = row.get("climate", np.nan)
            ablation_rows.append({
                "trait": trait,
                "r2_climate": base,
                "r2_climate_topography": row.get("climate_topography", np.nan),
                "r2_climate_soil": row.get("climate_soil", np.nan),
                "r2_full": row.get("climate_topography_soil", np.nan),
                "delta_topography_over_climate": row.get("climate_topography", np.nan) - base,
                "delta_soil_over_climate": row.get("climate_soil", np.nan) - base,
                "delta_full_over_climate": row.get("climate_topography_soil", np.nan) - base,
            })
        ablation = pd.DataFrame(ablation_rows)
    else:
        ablation = pd.DataFrame()
    ablation.to_csv(out_dir / "rf_environment_group_ablation.csv", index=False, encoding="utf-8-sig")

    coverage = []
    for name in climate + topo + soil:
        coverage.append({
            "predictor": name,
            "n_nonmissing": int(pd.to_numeric(df[name], errors="coerce").notna().sum()),
            "complete_fraction": float(pd.to_numeric(df[name], errors="coerce").notna().mean()),
        })
    pd.DataFrame(coverage).to_csv(out_dir / "environment_predictor_coverage.csv", index=False, encoding="utf-8-sig")

    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "input_observations": int(len(df)),
        "input_species": int(df["taxon_name"].nunique()),
        "traits": TRAITS,
        "predictor_groups": groups,
        "cv": {
            "method": "GroupKFold by 10-degree spatial block",
            "n_splits": args.n_splits,
            "species_weighting": "inverse complete-case observations per species, normalized to mean 1",
            "within_transform": "species-demeaned outcome and predictors",
        },
        "rf": {
            "n_estimators": args.n_estimators,
            "max_features": "sqrt",
            "min_samples_leaf": args.min_samples_leaf,
            "random_state": args.random_state,
        },
        "environment_metadata": env_meta,
        "outputs": {
            "cv_fold_rows": int(len(folds)),
            "cv_summary_rows": int(len(summary)),
            "permutation_importance_rows": int(len(importance)),
            "ablation_rows": int(len(ablation)),
        },
        "interpretation_boundary": "RF is predictive association analysis. Variable importance and performance gains are not causal effects.",
    }
    (out_dir / "exhaustive_within_species_environment_rf_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps({
        "input_observations": report["input_observations"],
        "input_species": report["input_species"],
        "cv_summary_rows": report["outputs"]["cv_summary_rows"],
        "ablation": ablation.to_dict("records"),
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
