"""Extend V4 S1+S2 at original coordinates with CHELSA extreme BIOCLIM and elevation.

S1 source = checksum-verified 313350 native observation table. S2 source =
source-authenticated 297622 wild native/introduced observations in the TWO
approved focal species. Same obs_id and coordinates are deduplicated before
raster sampling; no original observation coordinates or identifiers are uploaded.

New source files/physical units are audited independently of biological traits.
No coarse grid average, guessed elevation from latitude, or imputed BIO18.
"""
from __future__ import annotations
import argparse
import gzip
import hashlib
import json
import shutil
import time
import zipfile
from pathlib import Path
import numpy as np
import pandas as pd

EXTREMES = ("chelsa_bio05","chelsa_bio06","chelsa_bio13","chelsa_bio14",
            "chelsa_bio17","chelsa_bio18","topo_elevation")
CHELSA = set(EXTREMES) - {"topo_elevation"}
CHELSA_TEMPLATE = ("https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/"
                   "bioclim/{bio}/1981-2010/CHELSA_{bio}_1981-2010_V.2.1.tif")
CHELSA_LEGACY = ("https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/"
                 "climatologies/1981-2010/bio/CHELSA_{bio}_1981-2010_V.2.1.tif")
ELEVATION_URL = "https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_30s_elev.zip"
CHUNK = 1024 * 1024
MAX_BYTES = 3_500_000_000
MIN_PAYLOAD_BYTES = 10_000_000


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as source:
        for data in iter(lambda: source.read(CHUNK * 8), b""):
            h.update(data)
    return h.hexdigest()


def validate_coordinates(native, focal):
    columns = ("obs_id", "latitude", "longitude")
    for name, frame in (("S1", native), ("S2", focal)):
        if not set(columns) <= set(frame) or frame.obs_id.duplicated().any():
            raise ValueError(name + "_MISSING_OR_DUPLICATE_SOURCE_OBS_ID")
        if len(frame) != (313350 if name == "S1" else 297622):
            raise ValueError(name + "_WRONG_ORIGINAL_DENOMINATOR")
        coords = frame[["latitude", "longitude"]].to_numpy(float)
        if not np.isfinite(coords).all() or np.any(np.abs(coords[:, 0]) > 90) or np.any(np.abs(coords[:, 1]) > 180):
            raise ValueError(name + "_INVALID_SOURCE_COORDINATE")
    if not ({"range_state", "taxon_name"} <= set(focal)):
        raise ValueError("S2_WCVP_RANGE_OR_TAXON_MISSING")
    if set(focal.taxon_name.unique()) != {"Cirsium arvense", "Cirsium vulgare"}:
        raise ValueError("S2_NOT_TWO_DECLARED_SPECIES")
    if set(focal.range_state.unique()) != {"native", "introduced"}:
        raise ValueError("S2_NOT_NATIVE_AND_INTRODUCED")
    cat = pd.concat([native[list(columns)], focal[list(columns)]], ignore_index=True)
    dup = cat.loc[cat.duplicated("obs_id", keep=False)].copy()
    if len(dup):
        bounds = dup.groupby("obs_id")[["latitude", "longitude"]].agg(["min","max"])
        for coordinate in ("latitude","longitude"):
            if ((bounds[(coordinate,"max")] - bounds[(coordinate,"min")]).abs() > 1e-8).any():
                raise ValueError("S1_S2_SAME_OBS_ID_DIFFERENT_COORDINATES")
    distinct = cat.drop_duplicates("obs_id", keep="first").reset_index(drop=True)
    if not 313350 <= len(distinct) <= 313350 + 156750:
        raise ValueError("UNEXPECTED_OVERLAP_OF_S1_AND_S2_ORIGINAL_CENSUS")
    return distinct


def urls(variable):
    if variable == "topo_elevation":
        return [ELEVATION_URL]
    if variable not in CHELSA:
        raise ValueError("UNREGISTERED_ENVIRONMENT_VARIABLE")
    bio = variable.split("_")[-1]
    return [CHELSA_TEMPLATE.format(bio=bio), CHELSA_LEGACY.format(bio=bio)]


def retrieve(variable, temp_dir):
    import requests
    paths = urls(variable)
    failures = []
    target = Path(temp_dir) / (variable + (".zip" if variable == "topo_elevation" else ".tif"))
    for url in paths:
        start = time.monotonic()
        n = 0
        target.unlink(missing_ok=True)
        try:
            with requests.get(url, headers={"User-Agent":"azami-v4-s1-s2-extreme-bioclim/20261024"},
                              timeout=(25,100), stream=True) as response:
                response.raise_for_status()
                with target.open("wb") as file:
                    for payload in response.iter_content(chunk_size=CHUNK):
                        if not payload:
                            continue
                        n += len(payload)
                        if n > MAX_BYTES or time.monotonic() - start > 1000:
                            raise RuntimeError("DOWNLOAD_SIZE_OR_TIME_BUDGET")
                        file.write(payload)
            if n < MIN_PAYLOAD_BYTES:
                raise RuntimeError("SOURCE_TOO_SMALL")
            return target, {
                "variable":variable, "url":url, "download_bytes":n,
                "download_sha256":sha256(target), "elapsed_seconds":round(time.monotonic()-start,2),
                "fallbacks_attempted":failures,
            }
        except Exception as exc:
            failures.append({"url":url, "error":str(exc)[:250]})
            target.unlink(missing_ok=True)
    raise RuntimeError("ENVIRONMENT_ACQUISITION_FAILED " + variable + " " + json.dumps(failures))


def open_raster(downloaded, variable, temporary):
    if variable != "topo_elevation":
        return downloaded
    with zipfile.ZipFile(downloaded) as archive:
        candidates = [m for m in archive.namelist() if m.endswith(".tif") and
                      "elev" in m.lower()]
        if len(candidates) != 1 or ".." in Path(candidates[0]).parts:
            raise ValueError("WORLDCLIM_ELEVATION_ZIP_MEMBER_INVALID")
        target = Path(temporary) / "worldclim_elevation_30s.tif"
        with archive.open(candidates[0]) as stream, target.open("wb") as dest:
            shutil.copyfileobj(stream, dest, CHUNK * 4)
        return target


def sample_grid(raster_path, variable, frame):
    import rasterio
    from rasterio.windows import Window
    with rasterio.open(raster_path) as source:
        if (source.count != 1 or source.crs is None or
            source.crs.to_epsg() != 4326 or source.width < 20000 or source.height < 12000):
            raise ValueError("NON_GLOBAL_OR_UNEXPECTED_RASTER_GEOMETRY " + variable)
        if not all(np.isfinite(source.transform)[:6]):
            raise ValueError("INVALID_RASTER_AFFINE")
        scale, offset = float(source.scales[0]), float(source.offsets[0])
        if not np.isfinite(scale) or scale <= 0 or not np.isfinite(offset):
            raise ValueError("INVALID_RASTER_DECODE_SCALE_OFFSET")
        lat = frame.latitude.to_numpy(float); lon = frame.longitude.to_numpy(float)
        # Unique pixel indices avoid repeatedly decoding identical raster blocks.
        inv = ~source.transform
        col_float, row_float = inv * (lon, lat)
        row, col = np.floor(row_float).astype(int), np.floor(col_float).astype(int)
        valid = (row >= 0)&(row < source.height)&(col >= 0)&(col < source.width)
        values = np.full(len(frame), np.nan, dtype=float)
        indices = np.flatnonzero(valid)
        if len(indices):
            # Read sorted 256-by-256 source tiles; only source pixels are sampled.
            tile_rows = row[indices]//256
            tile_cols = col[indices]//256
            tile_ids = tile_rows*((source.width + 255)//256) + tile_cols
            ordering = np.argsort(tile_ids, kind="stable")
            sorted_ids = indices[ordering]
            sorted_tile = tile_ids[ordering]
            ends = np.r_[0, np.flatnonzero(np.diff(sorted_tile))+1, len(sorted_ids)]
            for l, r in zip(ends[:-1], ends[1:]):
                ix = sorted_ids[l:r]
                y0 = (row[ix[0]]//256)*256
                x0 = (col[ix[0]]//256)*256
                height = min(256, source.height-y0)
                width = min(256, source.width-x0)
                block = source.read(1, window=Window(x0,y0,width,height), masked=True)
                block_values = block.astype(float).filled(np.nan)
                sample = block_values[row[ix]-y0, col[ix]-x0]
                if source.nodata is not None:
                    sample[np.isclose(sample, float(source.nodata), atol=0, rtol=0)] = np.nan
                # Guard against unsigned overflow encodings left undeclared as nodata.
                if np.dtype(source.dtypes[0]).kind == "u":
                    sample[sample == np.iinfo(np.dtype(source.dtypes[0])).max] = np.nan
                values[ix] = sample
        decoded = values * scale + offset
        if (np.isfinite(decoded).sum() < len(frame) * .75):
            raise ValueError("NEW_ENVIRONMENT_LAYER_COVERAGE_BELOW_75_PERCENT " + variable)
        receipt = {
            "variable":variable,"n_original_observation_locations":len(frame),
            "n_source_pixels_inside_extent":int(valid.sum()),
            "n_finite_extracted":int(np.isfinite(decoded).sum()),
            "coverage_fraction":round(float(np.isfinite(decoded).mean()),8),
            "crs":str(source.crs),"width":source.width,"height":source.height,
            "source_nodata":None if source.nodata is None else float(source.nodata),
            "source_dtype":source.dtypes[0],
            "scale":scale,"offset":offset,
            "variable_unit":"m" if variable == "topo_elevation"
                         else ("degC" if variable in ("chelsa_bio05","chelsa_bio06")
                               else "kg m-2 month-1"),
            "sampled_from_original_raster_pixels":True,
            "interpolation_or_cell_averaging":False,
        }
        return decoded, receipt


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--native", type=Path, required=True)
    p.add_argument("--focal", type=Path, required=True)
    p.add_argument("--out-private", type=Path, required=True)
    p.add_argument("--out-public", type=Path, required=True)
    p.add_argument("--temporary", type=Path, required=True)
    p.add_argument("--variables", default=",".join(EXTREMES))
    args = p.parse_args()
    selected = [v.strip() for v in args.variables.split(",") if v.strip()]
    if not selected or len(selected) != len(set(selected)) or any(v not in EXTREMES for v in selected):
        raise ValueError("UNREGISTERED_OR_DUPLICATE_VARIABLE")
    for target in (args.out_private,args.out_public,args.temporary):
        target.mkdir(parents=True,exist_ok=True)
    native = pd.read_csv(args.native,usecols=["obs_id","latitude","longitude"],dtype={"obs_id":str})
    focal = pd.read_csv(args.focal,usecols=["obs_id","latitude","longitude","range_state","taxon_name"],
                        dtype={"obs_id":str})
    source = validate_coordinates(native,focal)
    data = source[["obs_id"]].copy()
    receipt = {"analysis_id":"azami_v4_S1_S2_extreme_bioclim_20261024",
               "population":{"S1_native_observations":len(native),
                             "S2_both_range_observations":len(focal),
                             "deduplicated_original_locations":len(source)},
               "variables_requested":selected,
               "sources":[],"no_point_level_public_output":True}
    for variable in selected:
        archive,meta = retrieve(variable,args.temporary)
        decoded = open_raster(archive,variable,args.temporary)
        values,diagnostic = sample_grid(decoded,variable,source)
        data[variable] = values
        receipt["sources"].append({**meta,**diagnostic})
        if decoded != archive:
            decoded.unlink(missing_ok=True)
        archive.unlink(missing_ok=True)
        # Save progress privately after each layer; no public "all completed" receipt until entire set passes.
        data.to_csv(args.out_private/"environment_ext_partial.csv.gz",index=False,compression="gzip")
        print(json.dumps({"variable":variable,"finite":diagnostic["n_finite_extracted"],
                          "coverage":diagnostic["coverage_fraction"]}),flush=True)
    output=args.out_private/"environment_ext.csv.gz"
    (args.out_private/"environment_ext_partial.csv.gz").replace(output)
    receipt["private_rows"]=len(data)
    receipt["status"]="ALL_REQUESTED_REAL_RASTERS_SAMPLED"
    (args.out_public/"environment_source_receipt.json").write_text(
        json.dumps(receipt,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({"status":receipt["status"],"rows":len(data),
                      "sampled_variables":selected}),flush=True)


if __name__ == "__main__":
    main()
