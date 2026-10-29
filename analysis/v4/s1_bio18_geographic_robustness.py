"""Source-authenticated exploratory BIO18 × image angle geographic robustness.

This is a NEW geographically controlled analysis of the original private
S1 observations; it does NOT re-estimate or replace the frozen 80-test Holm
family. Real BIO18, BIO05 and elevation are sampled at source pixels.
Outputs contain AGGREGATE coefficients and sector counts only.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t

from analysis.v4.s1_rain_native_geographic_controls import (
    geo_features, longitude_sector, GEO_COLS, SECTOR_NAMES
)
from analysis.v4.s1_s2_extreme_environment import (
    retrieve, open_raster, sample_grid
)

SOURCE_SHA = "342d78e9a8a960e8d9bdcf0bc2cdacd13bee48214a8b962419d5a71d5aef46ed"
ENDPOINT = "orientation_image_vertical_angle"
EXTRA = ("chelsa_bio18", "chelsa_bio05", "topo_elevation")
ENV = ("chelsa_bio18", "chelsa_bio05", "chelsa_rsds_mean", "topo_elevation")
OLD_ENV = ("chelsa_bio01", "chelsa_bio04", "chelsa_bio12", "chelsa_bio15",
           "chelsa_rsds_mean")
SCALE = np.array([100., 10., 1., 1000.], float)
CONTROLS = ("none", "longitude_sector_fractions", "spherical_second_order")
MODES = ("original_photo_weighted", "taxon_equal_5degree_cells")


def require(condition, msg):
    if not condition:
        raise ValueError(msg)


def original_table(path: Path) -> pd.DataFrame:
    require(hashlib.sha256(path.read_bytes()).hexdigest() == SOURCE_SHA,
            "ORIGINAL_ENCRYPTED_S1_DECODED_SHA_DRIFT")
    d = pd.read_csv(path, dtype={"obs_id": str}, low_memory=False)
    required = {"obs_id", "taxon_name", "accepted_key", "latitude", "longitude",
                ENDPOINT, ENDPOINT + "__eligible", *OLD_ENV}
    require(required <= set(d), "S1_ORIGINAL_IDENTITIES_PHENOTYPES_CLIMATE_MISSING")
    require(len(d) == 313350 and d.obs_id.is_unique,
            "ORIGINAL_S1_NATIVE_POPULATION_OR_ID_CHANGED")
    require(np.isfinite(d[["latitude", "longitude"]].to_numpy(float)).all(),
            "NONFINITE_SOURCE_COORDINATES")
    require(d.latitude.between(-90, 90).all() and
            d.longitude.between(-180, 180).all(),
            "SOURCE_COORDINATES_INVALID")
    eligible = d[ENDPOINT+"__eligible"].fillna(False).astype(bool).to_numpy()
    old = eligible & d[ENDPOINT].between(0, 180).to_numpy()
    old &= np.isfinite(d[list(OLD_ENV)].to_numpy(float)).all(axis=1)
    require(int(old.sum()) == 170066 and
            int(d.loc[old, "accepted_key"].nunique()) == 244,
            "OLD_FIVE_CLIMATE_S1_ANGLE_SUPPORT_CHANGED")
    return d


def add_source_pixels(d: pd.DataFrame, temp: Path) -> tuple[pd.DataFrame, list]:
    """Query only original raster pixels. Never save plaintext observation rows."""
    temp.mkdir(parents=True, exist_ok=True)
    receipts = []
    for variable in EXTRA:
        archive, source = retrieve(variable, temp)
        raster = open_raster(archive, variable, temp)
        values, audit = sample_grid(raster, variable, d)
        require(len(values) == len(d), "RASTER_SOURCE_ID_ALIGNMENT_CHANGED")
        d[variable] = values
        receipts.append({"variable": variable,
                         "raster_sha256": source["download_sha256"],
                         "raster_url": source["url"],
                         "n_finite": audit["n_finite_extracted"],
                         "decode_scale": audit["scale"],
                         "decode_offset": audit["offset"],
                         "source_pixel": audit["sampled_from_original_raster_pixels"]})
        if raster != archive:
            raster.unlink(missing_ok=True)
        archive.unlink(missing_ok=True)
        print(f"Sampled original pixels: {variable}; finite={audit['n_finite_extracted']}",
              flush=True)
    return d, receipts


def prepare(d: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    need = {ENDPOINT, ENDPOINT+"__eligible", "accepted_key", "taxon_name",
            "latitude", "longitude", *OLD_ENV, *ENV}
    require(need <= set(d), "SOURCE_PHOTO_OR_FOUR_CLIMATE_MISSING")
    old = d[ENDPOINT+"__eligible"].fillna(False).astype(bool).to_numpy()
    old &= d[ENDPOINT].between(0, 180).to_numpy()
    old &= np.isfinite(d[list(OLD_ENV)].to_numpy(float)).all(axis=1)
    require(int(old.sum()) == 170066, "ORIGINAL_S1_FIVE_CLIMATE_COHORT_DRIFT")
    subset = old & np.isfinite(d[list(ENV)].to_numpy(float)).all(axis=1)
    x = d.loc[subset, ["accepted_key", "taxon_name", "latitude", "longitude",
                       ENDPOINT, *ENV]].copy()
    require(x.accepted_key.notna().all() and
            x.latitude.between(-90, 90).all() and
            x.longitude.between(-180, 180).all(),
            "SOURCE_TAXON_OR_LOCATION_NOT_VALID")
    x["accepted_key"] = x.accepted_key.astype(str)
    x["cell_lat"] = np.floor((x.latitude+90)/5).astype(int)
    x["cell_lon"] = np.floor((x.longitude+180)/5).astype(int)
    sector = longitude_sector(x.longitude.to_numpy(float))
    for k, label in enumerate(SECTOR_NAMES):
        x["sector_"+label] = (sector == k).astype(float)
    harmonics = geo_features(x.latitude.to_numpy(float),
                             x.longitude.to_numpy(float))
    for j, col in enumerate(GEO_COLS):
        x[col] = harmonics[:, j]
    return x, {
        "original_native_source_observations": int(len(d)),
        "original_five_climate_angle_eligible": int(old.sum()),
        "new_four_predictor_geo_audit_angle_rows": int(len(x)),
        "new_four_predictor_geo_audit_taxa_before_taxonomy_and_min_n":
            int(x.accepted_key.nunique()),
        "not_identical_to_frozen_all_12_climate_quality_subset": True,
    }


def taxon_roster(x: pd.DataFrame, mode: str) -> pd.DataFrame:
    require(mode in MODES, "UNREGISTERED_PHOTO_VS_CELL_WEIGHTING")
    measures = [ENDPOINT, *ENV, *GEO_COLS,
                *["sector_"+name for name in SECTOR_NAMES]]
    n = x.groupby("accepted_key").size().rename("photo_n")
    concepts = x.groupby("accepted_key").taxon_name.nunique().rename("n_names")
    cells = x.groupby("accepted_key")[["cell_lat", "cell_lon"]].apply(
        lambda v: len(v.drop_duplicates())).rename("cell_n")
    if mode == "original_photo_weighted":
        means = x.groupby("accepted_key")[measures].mean()
    else:
        cell = x.groupby(["accepted_key", "cell_lat", "cell_lon"], sort=True)[
            measures].mean().reset_index()
        means = cell.groupby("accepted_key")[measures].mean()
    y = means.join([n, concepts, cells]).reset_index()
    y = y.loc[(y.n_names == 1) & (y.photo_n >= 5)].copy()
    props = y[["sector_"+name for name in SECTOR_NAMES]].to_numpy(float)
    y["dominant_sector"] = np.argmax(props, axis=1)
    require((np.abs(props.sum(axis=1)-1) < 1.e-9).all(),
            "OBSERVED_SOURCE_SECTOR_WEIGHT_NOT_ONE")
    require(y.accepted_key.is_unique and len(y) >= 80, "TOO_FEW_DISTINCT_NATIVE_TAXA")
    return y.sort_values("accepted_key")


def design(x: pd.DataFrame, geo: str) -> np.ndarray:
    require(geo in CONTROLS, "UNREGISTERED_GEO_CONTROL")
    climate = (x[list(ENV)].to_numpy(float) -
               x[list(ENV)].to_numpy(float).mean(axis=0))/SCALE
    parts = [np.ones((len(x), 1)), climate]
    if geo == "longitude_sector_fractions":
        extra = x[["sector_"+name for name in SECTOR_NAMES[:2]]].to_numpy(float)
    elif geo == "spherical_second_order":
        extra = x[list(GEO_COLS)].to_numpy(float)
    else:
        extra = np.zeros((len(x), 0))
    if extra.shape[1]:
        sd = extra.std(axis=0)
        require((sd > 1.e-10).all(), "SPATIAL_TERMS_CONSTANT")
        parts.append((extra-extra.mean(axis=0))/sd)
    m = np.column_stack(parts)
    require(len(x) > m.shape[1]+8 and np.linalg.matrix_rank(m) == m.shape[1],
            "GEOGRAPHIC_DESIGN_NOT_ESTIMABLE")
    require(np.linalg.cond(m) < 1.e7, "GEOGRAPHIC_DESIGN_COLLINEAR")
    return m


def fit(x: pd.DataFrame, geo: str) -> dict:
    try:
        m = design(x, geo)
        y = x[ENDPOINT].to_numpy(float)
        inv = np.linalg.inv(m.T@m)
        b = inv@(m.T@y)
        resid = y-m@b
        h = np.einsum("ij,jk,ik->i", m, inv, m)
        require(np.max(h) < .999, "EXCESSIVE_SINGLE_TAXON_LEVERAGE")
        score = m*(resid/(1-h))[:, None]
        cov = inv@(score.T@score)@inv
        se = np.sqrt(max(0., cov[1,1]))
        q = float(t.ppf(.975, len(y)-m.shape[1]))
        return {"status":"estimated","n_taxa":int(len(x)),
                "BIO18_degrees_per_100mm_month":float(b[1]),
                "BIO18_HC3_pointwise_CI95":[float(b[1]-q*se),float(b[1]+q*se)],
                "n_design_parameters":int(m.shape[1]),
                "max_taxon_leverage":float(h.max()),
                "not_frozen_S1_WLS_or_Holm":True,
                "pointwise_nonphylogenetic_spatial_sensitivity":True}
    except (ValueError, np.linalg.LinAlgError) as exc:
        return {"status":"not_identifiable", "reason":str(exc)[:160],
                "n_taxa":int(len(x))}


def evaluate(roster: pd.DataFrame) -> dict:
    n_by_sector = {name:int((roster.dominant_sector==j).sum())
                   for j,name in enumerate(SECTOR_NAMES)}
    full = {name:fit(roster,name) for name in CONTROLS}
    leave = {}
    within = {}
    for j,name in enumerate(SECTOR_NAMES):
        without = roster.loc[roster.dominant_sector!=j]
        inside = roster.loc[roster.dominant_sector==j]
        leave[name] = fit(without, "none") if len(without)>=50 else {
            "status":"not_identifiable", "reason":"TOO_FEW_TAXA_AFTER_SECTOR_HOLDOUT"}
        within[name] = fit(inside, "none") if len(inside)>=25 else {
            "status":"not_identifiable", "reason":"TOO_FEW_TAXA_IN_SECTOR"}
    return {"n_taxa":int(len(roster)), "dominant_sector_n":n_by_sector,
            "full_cohort_geographic_adjustments":full,
            "leave_dominant_sector_out":leave,
            "within_fixed_longitude_sector":within,
            "not_independent_evolutionary_replication":True}


def selftest():
    """Guard against a misleading global relationship caused by region sorting."""
    rng = np.random.default_rng(20261029)
    n = 240
    sector = np.repeat(np.arange(3),n//3)
    z = rng.normal(0, 1, n)
    fake = pd.DataFrame({
       ENDPOINT: 60.+15.*sector+rng.normal(0,3,n),
       "chelsa_bio18": 80.+60.*sector+10*z,
       "chelsa_bio05": 25.+rng.normal(0,2,n),
       "chelsa_rsds_mean": 15.+rng.normal(0,1,n),
       "topo_elevation": 500.+rng.normal(0,100,n),
    })
    for k,name in enumerate(SECTOR_NAMES):
        fake["sector_"+name]=(sector==k).astype(float)
    naive = fit(fake,"none")
    adjusted = fit(fake,"longitude_sector_fractions")
    require(naive["status"] == adjusted["status"] == "estimated",
            "SYNTHETIC_GEO_MODEL_NOT_ESTIMABLE")
    require(naive["BIO18_degrees_per_100mm_month"] > 8 and
            abs(adjusted["BIO18_degrees_per_100mm_month"]) < 8,
            "SYNTHETIC_REGION_CONFOUNDING_NOT_RECOVERED")
    print("SYNTHETIC_SECTOR_CONFOUNDING_DIAGNOSTIC_OK")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--native",type=Path)
    p.add_argument("--temporary",type=Path)
    p.add_argument("--output",type=Path)
    p.add_argument("--selftest",action="store_true")
    args = p.parse_args()
    if args.selftest:
        selftest()
        return
    require(args.native and args.temporary and args.output,
            "ORIGINAL_NATIVE_INPUT_AND_PRIVATE_TEMP_REQUIRED")
    native = original_table(args.native)
    native, sources = add_source_pixels(native,args.temporary)
    image, scope = prepare(native)
    panels = {mode:evaluate(taxon_roster(image,mode)) for mode in MODES}
    output = {
      "schema":"azami_v4_followup_S1_BIO18_geographic_robustness_20261029_v1",
      "status":"SOURCE_VERIFIED_EXPLORATORY_GEO_SENSITIVITY",
      "source_SHA256":SOURCE_SHA,
      "scope":scope,
      "predictors":list(ENV),
      "climate_raster_source_receipts":sources,
      "photo_weighting_and_taxon_cell_equalization":panels,
      "meaning":"Test if between-taxon BIO18-angle association survives geographic composition",
      "no_original_observation_coordinates_names_or_ids_exported":True,
      "fixed_80_test_S1_primary_result_NOT_replaced_or_recorrected":True,
      "multiple_geographic_diagnostics_exploratory_not_new_discovery_family":True,
      "longitude_sectors_are_not_independent_evolutionary_lineages":True,
      "head_angle_is_camera_frame_not_gravity":True
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(output,indent=2,allow_nan=False)+"\n",
                           encoding="utf-8")
    print(json.dumps({"status":output["status"],
        "scope":scope,"panels":{
            k:v["full_cohort_geographic_adjustments"] for k,v in panels.items()
        }},ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
