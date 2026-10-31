"""Optional private WorldClim 2.1 elevation-derived terrain RF pilot.

This is NOT a refit of frozen GEB hypotheses. Five DEM-derived quantities are
observational site descriptors, not causal/topographic adaptation evidence.
Only source SHA, aggregate counts and held-out scores may be exported.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.v4.native_integration import FEATURES
from analysis.v4.rf_archived_atmosphere import (
    ATMOSPHERE23, PREDICTOR_GROUPS, verify_original_1981_2010_airfields,
)
from analysis.v4.rf_extended_environment import (
    BIO9, TERRAIN5, verify_extended_source,
)
from analysis.v4.rf_multitrait_predictability import (
    aggregate_endpoint, score_models, verified_source, sha256,
)

MODEL_SET = {
    "chelsa_vpd_wind23": ATMOSPHERE23,
    "chelsa_vpd_wind23_plus_worldclim_terrain28": ATMOSPHERE23 + TERRAIN5,
}
GROUPS = {**PREDICTOR_GROUPS, "terrain_worldclim_dem": TERRAIN5}
MIN_TPI_NEIGHBORS = 60     # among the 120 OTHER positions in an 11x11 window
RADIUS = 5                 # explicitly fixed landscape neighbourhood (~5 km)
TILE = 256
SOURCE_URL = "https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_30s_elev.zip"

assert len(MODEL_SET["chelsa_vpd_wind23_plus_worldclim_terrain28"]) == 28
assert len(set(MODEL_SET["chelsa_vpd_wind23_plus_worldclim_terrain28"])) == 28


def _rect(integral: np.ndarray, r0, c0, r1, c1):
    return integral[r1, c1] - integral[r0, c1] - integral[r1, c0] + integral[r0, c0]


def _terrain_from_pixels(tile: np.ndarray, local_r: np.ndarray,
                         local_c: np.ndarray, row_lat: np.ndarray,
                         lon_step: float, lat_step: float):
    """Take precise 30-sec pixels and calculate physical slope/ruggedness/TPI.

    In a geographic DEM, distance between longitude grid columns decreases
    with cos(latitude). Heading is *down-slope*; flat cells receive 0/0.
    Slope, aspect and RMS ruggedness require all nine 3x3 DEM cells.
    TPI is focal height minus the average OTHER usable 11x11 cells; if fewer
    than 60 neighbors exist, TPI is missing (not filled from another habitat).
    """
    from scipy.ndimage import uniform_filter
    arr = np.asarray(tile, float)
    r = np.asarray(local_r, int); c = np.asarray(local_c, int)
    lat = np.asarray(row_lat, float)
    n = len(r)
    if (arr.ndim != 2 or len(c) != n or len(lat) != n or
            not (np.isfinite(lat).all()) or lon_step <= 0 or lat_step <= 0 or
            (r - RADIUS < 0).any() or (c - RADIUS < 0).any() or
            (r + RADIUS >= arr.shape[0]).any() or
            (c + RADIUS >= arr.shape[1]).any()):
        raise ValueError("INVALID_DEM_TILE_NEIGHBOURHOOD")
    window = np.stack([arr[r + dr, c + dc]
                       for dr in (-1, 0, 1) for dc in (-1, 0, 1)])
    complete3 = np.isfinite(window).all(axis=0) & (np.abs(lat) < 85)
    dx = 111_320 * np.cos(np.deg2rad(lat)) * lon_step
    dy = 111_132 * lat_step
    # Horn derivatives with rows increasing southward; 0: top-left, 8: bottom-right.
    east_grad = (window[2] + 2*window[5] + window[8]
                 - window[0] - 2*window[3] - window[6]) / (8*dx)
    north_grad = (window[0] + 2*window[1] + window[2]
                  - window[6] - 2*window[7] - window[8]) / (8*dy)
    norm = np.hypot(east_grad, north_grad)
    slope = np.degrees(np.arctan(norm))
    neighbors = window[[0, 1, 2, 3, 5, 6, 7, 8]]
    rugged = np.sqrt(np.mean((neighbors - window[4][None, :])**2, axis=0))
    # north/east components of DOWNHILL unit heading: flat -> 0/0.
    northness = np.divide(-north_grad, norm, out=np.zeros(n), where=norm > 1e-12)
    eastness = np.divide(-east_grad, norm, out=np.zeros(n), where=norm > 1e-12)
    # Integral-image-like 11x11 totals using SciPy's fast separable uniform filter.
    valid = np.isfinite(arr)
    totals = uniform_filter(np.where(valid, arr, 0.), size=11,
                            mode="constant", cval=0.) * 121
    counts = uniform_filter(valid.astype(float), size=11,
                            mode="constant", cval=0.) * 121
    support = np.rint(counts[r, c]).astype(int) - np.isfinite(arr[r, c]).astype(int)
    other_sum = totals[r, c] - np.where(np.isfinite(arr[r, c]), arr[r, c], 0.)
    tpi = np.full(n, np.nan)
    support_ok = (support >= MIN_TPI_NEIGHBORS) & np.isfinite(arr[r, c])
    tpi[support_ok] = arr[r[support_ok], c[support_ok]] - other_sum[support_ok] / support[support_ok]
    result = np.column_stack([slope, rugged, tpi, northness, eastness])
    result[~complete3, :] = np.nan
    # TPI must separately satisfy 11x11 neighbor count (other channels unaffected).
    result[~support_ok, 2] = np.nan
    result[~np.isfinite(result)] = np.nan
    return result, np.asarray(arr[r, c], float), support


def sample_dem_terrain(source_tif: Path, expected_sha256: str,
                       lat: np.ndarray, lon: np.ndarray, *,
                       source_zip_sha256: str | None = None,
                       source_zip_file: Path | None = None,
                       require_global: bool = True):
    """Sample 3x3/11x11 topography; group 313k coordinates by DEM tiles.

    No interpolation, no source ID or location emitted. In-memory per-photo
    arrays exist only inside the protected runner and are subsequently removed.
    """
    import rasterio
    from rasterio.windows import Window

    if sha256(source_tif) != expected_sha256:
        raise ValueError("WORLDCLIM_DEM_TIFF_SHA_MISMATCH")
    if source_zip_file is not None and sha256(source_zip_file) != source_zip_sha256:
        raise ValueError("WORLDCLIM_ELEVATION_ZIP_SHA_NOT_PRIOR_SOURCE")
    lat = np.asarray(lat, float); lon = np.asarray(lon, float)
    if (len(lat) != len(lon) or len(lat) == 0 or
            not np.isfinite(lat).all() or not np.isfinite(lon).all() or
            (np.abs(lat) > 90).any() or (np.abs(lon) > 180).any()):
        raise ValueError("DEM_POINT_COORDINATE_INPUT_INVALID")
    n = len(lat)
    result = np.full((n, len(TERRAIN5)), np.nan)
    center = np.full(n, np.nan)
    with rasterio.open(source_tif) as src:
        transform = src.transform
        if (src.count != 1 or src.crs is None or src.crs.to_epsg() != 4326 or
                abs(transform.b) > 1e-12 or abs(transform.d) > 1e-12 or
                transform.a <= 0 or transform.e >= 0 or
                (require_global and (src.width < 20000 or src.height < 10000))):
            raise ValueError("WORLDCLIM_30S_GEOGRAPHIC_SINGLE_BAND_REQUIRED")
        scale = float(src.scales[0]); offset = float(src.offsets[0])
        if not (np.isfinite(scale) and scale > 0 and np.isfinite(offset)):
            raise ValueError("DEM_RASTER_UNIT_SCALE_INVALID")
        row = np.floor((lat - transform.f) / transform.e).astype(int)
        col = np.floor((lon - transform.c) / transform.a).astype(int)
        good = ((row >= 0) & (row < src.height) & (col >= 0) &
                (col < src.width) & (np.abs(lat) < 85))
        included = np.flatnonzero(good)
        blocks = (row[included] // TILE)*((src.width+TILE-1)//TILE) + col[included]//TILE
        order = np.argsort(blocks, kind="stable")
        ordered = included[order]
        group = blocks[order]
        breaks = np.r_[0, np.flatnonzero(np.diff(group))+1, len(ordered)]
        for left, right in zip(breaks[:-1], breaks[1:]):
            ids = ordered[left:right]
            tile_r = int(row[ids[0]]//TILE)*TILE
            tile_c = int(col[ids[0]]//TILE)*TILE
            window = Window(tile_c-RADIUS, tile_r-RADIUS,
                            TILE+2*RADIUS, TILE+2*RADIUS)
            raw = src.read(1, window=window, boundless=True, masked=True)
            z = raw.astype(float).filled(np.nan) * scale + offset
            if src.nodata is not None:
                z[np.isclose(z, float(src.nodata)*scale+offset, atol=0, rtol=0)] = np.nan
            local_r = row[ids] - tile_r + RADIUS
            local_c = col[ids] - tile_c + RADIUS
            # Use actual row-center latitude for east-west length.
            lat0 = transform.f + (row[ids]+0.5)*transform.e
            vals, elevations, _ = _terrain_from_pixels(z, local_r, local_c,
                                     lat0, transform.a, abs(transform.e))
            result[ids] = vals
            center[ids] = elevations
    summary = {"dem_tif_sha256":expected_sha256,
               "source_zip_sha256":source_zip_sha256,
               "source_url":SOURCE_URL,
               "pixel_radius_for_TPI":RADIUS,
               "TPI_min_land_neighbors":MIN_TPI_NEIGHBORS,
               "slope_method":"Horn 3x3 geographic physical pixel spacing",
               "roughness_method":"RMS focal-minus-8-neighbour elevation difference",
               "aspect_method":"downslope east/north heading unit components; flat=0,0",
               "TPI_method":"center minus available-other-cell mean in 11x11, >=60 valid",
               "unit_by_variable":dict(zip(TERRAIN5,("degree","m_RMS","m","unitless","unitless"))),
               "no_nodata_imputation_or_coordinate_export":True,
               "finite_source_records_by_trait":{
                   name:int(np.isfinite(result[:,i]).sum()) for i,name in enumerate(TERRAIN5)},
               "finite_all_terrain_records":int(np.isfinite(result).all(axis=1).sum()),
               "n_total_source_points":int(n)}
    return pd.DataFrame(result, columns=TERRAIN5),center,summary


def run_terrain(native:Path, receipt:Path, archived_seven:Path, archived_seven_sha256:str,
                extra_nine:Path, extra_nine_sha256:str, nine_manifest:Path,
                original_bio21_receipt:Path, original_metadata:Path,
                output_dir:Path, temporary:Path, *,
                constructs:tuple[str,...]=("angle","chroma","aspect"),
                min_per_cell:int=3, degrees:int=5,
                splits:tuple[str,...]=("unseen_taxa_balanced","unseen_longitude_sectors"),
                trees:int=80, repeats:int=2):
    """Paired 23-vs-28 on identical 28-complete cases, protected originals only."""
    if (not constructs or len(set(constructs)) != len(constructs) or
            not set(constructs).issubset(FEATURES) or min_per_cell not in (1,3) or
            degrees not in (1,2,5) or not splits or
            not set(splits).issubset({"unseen_taxa_balanced", "unseen_longitude_sectors"}) or
            len(set(splits))!=len(splits) or trees<8 or repeats<1):
        raise ValueError("TERRAIN_SOURCE_PILOT_SCOPE_NOT_DECLARED")
    decoded = json.loads(receipt.read_text())
    report_prev = json.loads(original_bio21_receipt.read_text())
    if (report_prev.get("source_observations") != 313350 or
            decoded.get("prepared_table_sha256") != report_prev.get("native_s1_table_sha256") or
            not report_prev.get("source_variables",{}).get("topo_elevation",{}).get("source_sha256")):
        raise ValueError("TERRAIN_DEM_PRIOR_SOURCE_PROVENANCE_INVALID")
    frame,prov=verified_source(native,receipt,archived_seven,archived_seven_sha256)
    frame,bio9=verify_extended_source(frame,extra_nine,extra_nine_sha256,nine_manifest,
                                      required_columns=BIO9)
    air=verify_original_1981_2010_airfields(frame,receipt,original_metadata)
    from analysis.v4.s1_s2_extreme_environment import (
        retrieve as retrieve_original_dem, open_raster,
    )
    temp=Path(temporary);temp.mkdir(parents=True,exist_ok=True)
    original_zip, download=retrieve_original_dem("topo_elevation",temp)
    expected_zip_sha=report_prev["source_variables"]["topo_elevation"]["source_sha256"]
    if download["download_sha256"] != expected_zip_sha:
        raise ValueError("TERRAIN_DEM_ZIP_DIFFERS_FROM_ARCHIVED_TOPO_ELEVATION")
    original_tif=open_raster(original_zip,"topo_elevation",temp)
    try:
        geom,center,topo=sample_dem_terrain(original_tif,sha256(original_tif),
                                           frame.latitude.to_numpy(float),
                                           frame.longitude.to_numpy(float),
                                           source_zip_file=original_zip,
                                           source_zip_sha256=expected_zip_sha)
    finally:
        original_tif.unlink(missing_ok=True)
        original_zip.unlink(missing_ok=True)
    a=frame["topo_elevation"].to_numpy(float)
    same=np.isfinite(a)&np.isfinite(center)
    n_same=int(same.sum());n_disagree=int((~np.isclose(a[same],center[same],atol=1e-4,rtol=0)).sum())
    if n_same < 250000 or n_disagree:
        raise ValueError("DEM_ELEVATION_NOT_IDENTICAL_TO_ARCHIVED_7_AXIS_SOURCE")
    private=temp/"terrain_private"
    private.mkdir(parents=True,exist_ok=True)
    geom.insert(0,"obs_id",frame.obs_id.to_numpy())
    raster_table=private/"worldclim_terrain5.csv.gz"
    geom.to_csv(raster_table,index=False,compression={"method":"gzip","mtime":0})
    table_sha=sha256(raster_table)
    manifest={"schema":"azami_rf_extra32_source_manifest_v1", "table_sha256":table_sha,
              "layers":{name:{"dataset":"WorldClim-elevation-derived",
                              "temporal_support":"static-topography",
                              "source_bundle_sha256":topo["dem_tif_sha256"],
                              "source_url":SOURCE_URL,
                              "unit":topo["unit_by_variable"][name],
                              "derivation":"Original 30sec DEM local 3x3 or 11x11 fixed-neighborhood terrain metric",
                              "sampling_method":"rasterio 30sec native pixel, non-interpolated"}
                        for name in TERRAIN5},
              "private_raw_obs_ids_and_coordinates_not_exported":True}
    man=private/"terrain_manifest.json"
    man.write_text(json.dumps(manifest,indent=2)+"\n")
    frame,extension=verify_extended_source(frame,raster_table,table_sha,man,
                                           required_columns=TERRAIN5)
    predictions={}
    predictor_sets=MODEL_SET
    for construct in constructs:
        try:
            cells=aggregate_endpoint(frame,construct,predictors=MODEL_SET[list(MODEL_SET)[1]],
                                     min_per_cell=min_per_cell,degrees=degrees)
        except ValueError as e:
            if "INSUFFICIENT_QUALIFIED_TAXA_OR_THREE_REGIONS" not in str(e):raise
            predictions[construct]={"status":"insufficient_sustained_species_geographic_support"}
            continue
        entries=[]
        for split in splits:
            if split=="unseen_longitude_sectors" and cells.groupby("region").taxon.nunique().min()<5:
                entries.append({"status":"insufficient_five_species_per_sector","split":split})
                continue
            entries.append(score_models(cells,construct,split,trees=trees,repeats=repeats,
                                        predictor_sets=predictor_sets,
                                        predictor_groups=GROUPS))
        predictions[construct]=entries
    out={"schema":"azami_original_s1_chelsa_atmosphere23_vs_worldclim_terrain28_exploratory_v1",
         "status":"original_source_terrain_derived_and_evaluated",
         "source":prov,"bio9":bio9,"air_verified":air,
         "terrain_provenance":{**topo,"matching_original_topo_elevation_records":n_same,
                               "mismatch_n":n_disagree,
                               "terrain_private_table_sha256":table_sha,
                               "terrain_private_manifest_sha256":sha256(man)},
         "paired_models":{k:list(v) for k,v in MODEL_SET.items()},
         "n_estimators":trees,"permutation_repeats":repeats,
         "design":{"constructs":list(constructs),"min_per_cell":min_per_cell,
                   "degrees":degrees,"splits":list(splits),"posthoc_sensitivity":True},
         "interpretation":"out-of-group photographic phenotype predictability, not causal topo effect",
         "frozen_original_GEB_S1_S2_H3c_unchanged":True,
         "no_raw_observations_coordinates_taxon_ids_or_predictions_exported":True,
         "results":predictions}
    output_dir.mkdir(parents=True,exist_ok=True)
    outfile=output_dir/(f"rf_23_to_28_terrain_{degrees}degree_min{min_per_cell}_aggregate.json")
    outfile.write_text(json.dumps(out,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    return out


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ("native","receipt","archived-seven","archived-seven-sha256",
                "bio-nine","bio-nine-sha256","nine-manifest","bio21-source-receipt",
                "original-metadata","outdir","temporary"):
        p.add_argument("--"+arg,type=str if arg.endswith("sha256") else Path,required=True)
    p.add_argument("--constructs",default="angle,chroma,aspect")
    p.add_argument("--degrees",type=int,choices=[1,2,5],default=5)
    p.add_argument("--min-per-cell",type=int,choices=[1,3],default=3)
    p.add_argument("--splits",default="unseen_taxa_balanced,unseen_longitude_sectors")
    p.add_argument("--trees",type=int,default=80);p.add_argument("--repeats",type=int,default=2)
    a=p.parse_args()
    run_terrain(a.native,a.receipt,a.archived_seven,a.archived_seven_sha256,
                a.bio_nine,a.bio_nine_sha256,a.nine_manifest,
                a.bio21_source_receipt,a.original_metadata,a.outdir,a.temporary,
                constructs=tuple(x.strip() for x in a.constructs.split(",") if x.strip()),
                degrees=a.degrees,min_per_cell=a.min_per_cell,
                splits=tuple(x.strip() for x in a.splits.split(",") if x.strip()),
                trees=a.trees,repeats=a.repeats)

if __name__=="__main__":
    main()