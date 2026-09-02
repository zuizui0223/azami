#!/usr/bin/env python3
"""Run exhaustive within-species climate diagnostics and fixed-effects models.

Input is the post-prediction strict spatially thinned observation table. Climate is
sampled from the four predeclared CHELSA V2.1 rasters. No trait-based sampling is
performed. Models use species-demeaned outcome and predictor values, then global
standardization of the demeaned values, with species-clustered standard errors.
"""
from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
import statsmodels.api as sm
from pyproj import Transformer
from statsmodels.stats.multitest import multipletests

PREDICTORS = {
    "chelsa_bio01": "https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/bioclim/bio01/1981-2010/CHELSA_bio01_1981-2010_V.2.1.tif",
    "chelsa_bio04": "https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/bioclim/bio04/1981-2010/CHELSA_bio04_1981-2010_V.2.1.tif",
    "chelsa_bio12": "https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/bioclim/bio12/1981-2010/CHELSA_bio12_1981-2010_V.2.1.tif",
    "chelsa_bio15": "https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/bioclim/bio15/1981-2010/CHELSA_bio15_1981-2010_V.2.1.tif",
}
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


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--observation", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--min-species-model-n", type=int, default=5)
    p.add_argument("--sample-batch-size", type=int, default=1024)
    return p.parse_args()


def sample_raster(url: str, lons: np.ndarray, lats: np.ndarray, batch_size: int) -> tuple[np.ndarray, dict]:
    env = rasterio.Env(
        GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
        CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif,.tiff",
        GDAL_HTTP_MULTIRANGE="YES",
        GDAL_HTTP_TIMEOUT="180",
        GDAL_HTTP_MAX_RETRY="5",
        GDAL_HTTP_RETRY_DELAY="3",
        VSI_CACHE="TRUE",
        VSI_CACHE_SIZE="100000000",
    )
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
                value = float(value)
                if math.isfinite(value):
                    values[start + offset] = value
        metadata = {
            "url": url,
            "crs": str(src.crs),
            "dtype": str(src.dtypes[0]),
            "nodata": None if src.nodata is None else float(src.nodata),
            "scales": list(src.scales),
            "offsets": list(src.offsets),
            "transform": list(src.transform)[:6],
        }
        return values, metadata
    finally:
        if src is not None:
            src.close()
        env.__exit__(None, None, None)


def finite_range(series: pd.Series) -> float:
    x = pd.to_numeric(series, errors="coerce").dropna()
    return float(x.max() - x.min()) if len(x) else np.nan


def build_species_diagnostics(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for taxon, g in df.groupby("taxon_name", sort=True):
        row = {
            "taxon_name": taxon,
            "n_observations": int(len(g)),
            "n_analysis_cells": int(g["analysis_cell"].nunique()) if "analysis_cell" in g else int(len(g)),
        }
        for pred in PREDICTORS:
            row[f"n_{pred}"] = int(pd.to_numeric(g[pred], errors="coerce").notna().sum())
            row[f"range_{pred}"] = finite_range(g[pred])
        for trait in TRAITS:
            row[f"n_{trait}"] = int(pd.to_numeric(g[trait], errors="coerce").notna().sum())
        rows.append(row)
    out = pd.DataFrame(rows)
    for threshold in (30, 50, 100, 200, 500):
        out[f"eligible_n_ge_{threshold}"] = out["n_observations"].ge(threshold)
    return out


def fit_within_model(df: pd.DataFrame, trait: str, predictor: str, min_n: int) -> dict:
    work = df[["taxon_name", trait, predictor]].copy()
    work[trait] = pd.to_numeric(work[trait], errors="coerce")
    work[predictor] = pd.to_numeric(work[predictor], errors="coerce")
    work = work.dropna()
    keep_species = []
    for taxon, g in work.groupby("taxon_name"):
        if len(g) < min_n:
            continue
        if g[trait].nunique() < 2 or g[predictor].nunique() < 2:
            continue
        keep_species.append(taxon)
    work = work.loc[work["taxon_name"].isin(keep_species)].copy()
    if len(work) < 20 or work["taxon_name"].nunique() < 2:
        return {"trait": trait, "predictor": predictor, "status": "insufficient"}

    work["y_dm"] = work[trait] - work.groupby("taxon_name")[trait].transform("mean")
    work["x_dm"] = work[predictor] - work.groupby("taxon_name")[predictor].transform("mean")
    y_sd = float(work["y_dm"].std(ddof=1))
    x_sd = float(work["x_dm"].std(ddof=1))
    if not (math.isfinite(y_sd) and math.isfinite(x_sd) and y_sd > 0 and x_sd > 0):
        return {"trait": trait, "predictor": predictor, "status": "zero_variance"}
    y = work["y_dm"] / y_sd
    x = work[["x_dm"]] / x_sd
    fit = sm.OLS(y, x).fit(cov_type="cluster", cov_kwds={"groups": work["taxon_name"]})
    return {
        "trait": trait,
        "predictor": predictor,
        "status": "ok",
        "beta_std_within": float(fit.params["x_dm"]),
        "se_cluster_species": float(fit.bse["x_dm"]),
        "z": float(fit.tvalues["x_dm"]),
        "p_value": float(fit.pvalues["x_dm"]),
        "ci_low": float(fit.conf_int().loc["x_dm", 0]),
        "ci_high": float(fit.conf_int().loc["x_dm", 1]),
        "n_observations": int(len(work)),
        "n_species": int(work["taxon_name"].nunique()),
        "min_species_complete_n": int(work.groupby("taxon_name").size().min()),
        "median_species_complete_n": float(work.groupby("taxon_name").size().median()),
    }


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

    lat = pd.to_numeric(df["latitude"], errors="coerce").to_numpy(float)
    lon = pd.to_numeric(df["longitude"], errors="coerce").to_numpy(float)
    valid = np.isfinite(lat) & np.isfinite(lon)
    if not valid.all():
        raise ValueError("Strict spatial cohort unexpectedly contains invalid coordinates")

    raster_meta = {}
    for pred, url in PREDICTORS.items():
        values, meta = sample_raster(url, lon, lat, args.sample_batch_size)
        df[pred] = values
        raster_meta[pred] = meta
        if np.isfinite(values).mean() < 0.98:
            raise RuntimeError(f"{pred} extraction coverage below 98%")

    climate_cols = list(PREDICTORS)
    df.to_csv(out_dir / "strict_spatial_thinned_with_chelsa.csv", index=False, encoding="utf-8-sig")

    diagnostics = build_species_diagnostics(df)
    diagnostics.to_csv(out_dir / "species_within_analysis_diagnostics.csv", index=False, encoding="utf-8-sig")

    model_rows = [fit_within_model(df, trait, pred, args.min_species_model_n) for trait in TRAITS for pred in PREDICTORS]
    models = pd.DataFrame(model_rows)
    ok = models["status"].eq("ok")
    models["q_fdr_bh"] = np.nan
    models["fdr_significant_0_05"] = False
    if ok.any():
        reject, qvals, _, _ = multipletests(models.loc[ok, "p_value"].to_numpy(float), alpha=0.05, method="fdr_bh")
        models.loc[ok, "q_fdr_bh"] = qvals
        models.loc[ok, "fdr_significant_0_05"] = reject
    models.to_csv(out_dir / "exhaustive_within_species_climate_models.csv", index=False, encoding="utf-8-sig")

    eligibility = {}
    for threshold in (30, 50, 100, 200, 500):
        eligibility[str(threshold)] = int(diagnostics["n_observations"].ge(threshold).sum())
    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "input_observations": int(len(df)),
        "input_species": int(df["taxon_name"].nunique()),
        "climate_complete_observations": int(df[climate_cols].notna().all(axis=1).sum()),
        "species_by_minimum_observation_threshold": eligibility,
        "models_total": int(len(models)),
        "models_ok": int(ok.sum()),
        "fdr_significant_models": int(models["fdr_significant_0_05"].sum()),
        "significant_models": models.loc[models["fdr_significant_0_05"]].to_dict("records"),
        "model_definition": {
            "source_cohort": "strict_spatial_thinned_observations: <=10 km positional accuracy, one observation per taxon x 0.25-degree cell",
            "sampling": "all-photo prediction first; no trait-based sampling; no 40-observation cap",
            "within_transform": "species-demean outcome and predictor",
            "standardization": "global SD after within-species demeaning",
            "uncertainty": "species-clustered robust standard errors",
            "minimum_complete_observations_per_species_per_model": args.min_species_model_n,
            "multiplicity": "Benjamini-Hochberg FDR across 36 primary trait-predictor models",
        },
        "raster_metadata": raster_meta,
    }
    (out_dir / "exhaustive_within_species_climate_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
