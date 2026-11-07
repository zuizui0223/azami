"""Secure original-source CHELSA BIOCLIM 12->21 stage; no raw rows published.

Runs ONLY inside a protected environment that can decrypt the frozen S1 source.
This source-extraction script intentionally has no GITHUB upload action.
Downloads full original CHELSA V2.1 1981-2010 rasters sequentially,
authenticates each raster's identity and raw byte digest, samples original
1-km pixels at the archived S1 coordinates without interpolation, and leaves
only PRIVATE observation-level CSVs plus a restricted source manifest. The
public receipt contains source URLs, hashes and aggregate coverage, but no
photo/observer/taxon/location identifiers or individual predictions.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd

from analysis.v4.rf_extended_environment import BIO9, DATASET_FAMILIES, EXPECTED_PERIODS
from analysis.v4.rf_multitrait_predictability import SOURCE_N, sha256
from analysis.v4.s1_s2_extreme_lego import EXTRA
from analysis.v4.s1_s2_extreme_environment import (
    retrieve as retrieve_original_extra,
    open_raster,
    sample_grid,
    CHELSA_TEMPLATE,
    CHELSA_LEGACY,
)

# These are the *physical* display units, not claims of new independently
# measured ecology. BIO03 is dimensionless isothermality (%) rather than degC.
PHYSICAL_UNIT = {
    **{f"chelsa_bio{i:02d}": "degree_C" for i in (2, 7, 8, 9, 10, 11)},
    "chelsa_bio03": "percent_isothermality",
    "chelsa_bio16": "kg_m-2_month-1",
    "chelsa_bio19": "kg_m-2_month-1",
}

EXTENSION = tuple(BIO9)
SOURCE_VARS = (*EXTRA, *EXTENSION)
MAX_LAYER_BYTES = 3_500_000_000
MIN_LAYER_BYTES = 10_000_000
MAX_FETCH_SECONDS = 900


def chelsa_candidates(variable: str) -> tuple[str, str]:
    if variable not in BIO9:
        raise ValueError("UNREGISTERED_BIOCLIM_LAYER")
    # Existing repository CHELSA source URLs retain two independently
    # documented original-version hosts, never substituted WorldClim data.
    return (
        CHELSA_TEMPLATE.format(bio=variable.removeprefix("chelsa_")),
        CHELSA_LEGACY.format(bio=variable.removeprefix("chelsa_")),
    )


def retrieve_bio9(variable: str, temporary: Path):
    """Download one CHELSA original GeoTIFF; keep checksum and failed routes."""
    import requests
    temporary.mkdir(parents=True, exist_ok=True)
    dest = temporary / (variable + ".tif")
    errors = []
    for url in chelsa_candidates(variable):
        t0 = time.monotonic()
        dest.unlink(missing_ok=True)
        count = 0
        try:
            with requests.get(url, stream=True, timeout=(25, 90),
                              headers={"User-Agent": "azami-v4-rf-bioclim-original-v2.1/20261031"}) as response:
                response.raise_for_status()
                with dest.open("wb") as stream:
                    for part in response.iter_content(chunk_size=1024 * 1024):
                        if not part:
                            continue
                        count += len(part)
                        if count > MAX_LAYER_BYTES or time.monotonic() - t0 > MAX_FETCH_SECONDS:
                            raise RuntimeError("CHELSA_DOWNLOAD_EXCEEDS_FIXED_RESOURCE_BUDGET")
                        stream.write(part)
            if count < MIN_LAYER_BYTES:
                raise ValueError("SOURCE_RASTER_TOO_SMALL")
            return dest, {"variable": variable, "url": url,
                          "download_bytes": count,
                          "download_sha256": sha256(dest),
                          "failed_routes": errors}
        except Exception as error:
            errors.append({"url": url, "error_type": type(error).__name__})
            dest.unlink(missing_ok=True)
    raise RuntimeError("AUTHENTIC_BIOCLIM_DOWNLOAD_FAILED: " + variable +
                       "; routes=" + ",".join(x["error_type"] for x in errors))


def check_raster_identity(path: Path, variable: str, *, is_new: bool) -> dict:
    """Reject mismatched year, wrong variable, projection or impossible decode."""
    import rasterio
    from analysis.v4.s1_s2_extreme_environment import EXTREMES

    if is_new and variable not in BIO9:
        raise ValueError("BIO9_UNREGISTERED_SOURCE")
    if not is_new and variable not in EXTREMES:
        raise ValueError("ORIGINAL7_UNREGISTERED_SOURCE")
    with rasterio.open(path) as src:
        if (src.count != 1 or src.crs is None or src.crs.to_epsg() != 4326 or
                src.width < 20000 or src.height < 12000):
            raise ValueError("CHELSA_GLOBAL_RASTER_GEOMETRY_INCOMPATIBLE")
        tags = src.tags()
        variable_tag = tags.get("variable", "")
        if variable != "topo_elevation":
            expected = variable.removeprefix("chelsa_")
            if variable_tag and variable_tag.lstrip("0") != expected.lstrip("0"):
                raise ValueError("CHELSA_WRONG_BIOCLIM_VARIABLE_TAG")
            version = str(tags.get("version", ""))
            timestamp = str(tags.get("datetime", ""))
            if version and version != "2.1":
                raise ValueError("CHELSA_WRONG_VERSION")
            if timestamp and timestamp != "1981-2010":
                raise ValueError("CHELSA_WRONG_CLIMATOLOGY_PERIOD")
        scale, offset = float(src.scales[0]), float(src.offsets[0])
        if not (np.isfinite(scale) and scale > 0 and np.isfinite(offset)):
            raise ValueError("CHELSA_INVALID_RASTER_SCALE_OFFSET")
        return {"crs": str(src.crs), "width": int(src.width),
                "height": int(src.height), "scale": scale,
                "offset": offset, "source_variable_tag": variable_tag,
                "source_version_tag": tags.get("version"),
                "source_datetime_tag": tags.get("datetime")}


def _sanity(values, variable):
    """Catch obvious physical-unit errors; do NOT truncate extremes."""
    v = np.asarray(values, float)
    finite = v[np.isfinite(v)]
    if finite.size < 100 or finite.size < 0.75 * len(v):
        raise ValueError("BIOCLIM_SAMPLE_COVERAGE_BELOW_75_PERCENT_" + variable)
    p01, p99 = np.quantile(finite, [.01, .99])
    if variable in PHYSICAL_UNIT:
        if variable == "chelsa_bio03":
            if p01 < 0 or p99 > 105:
                raise ValueError("ISOTHERMALITY_PHYSICAL_UNITS_UNEXPECTED")
        elif variable in ("chelsa_bio16", "chelsa_bio19"):
            if p01 < 0 or p99 > 3000:
                raise ValueError("QUARTER_RAIN_PHYSICAL_UNITS_UNEXPECTED")
        elif p01 < -100 or p99 > 80:
            raise ValueError("QUARTER_TEMPERATURE_PHYSICAL_UNITS_UNEXPECTED")
    return {"finite": int(finite.size), "fraction_finite": float(finite.size / len(v)),
            "p01": float(p01), "p99": float(p99)}


def _write_private_csv(frame: pd.DataFrame, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False,
                 compression={"method": "gzip", "mtime": 0})
    return sha256(path)


def extract(native: Path, decoded_receipt: Path, private_dir: Path,
            public_dir: Path, temporary: Path, *, expected_n=SOURCE_N,
            original_retrieve=retrieve_original_extra,
            bio_retrieve=retrieve_bio9,
            raster_sample=sample_grid):
    """Authenticate S1 and sequentially extract 7 archived + 9 new BIO layers.

    For test isolation, retrieval/sample functions may be injected. The
    original release requires the full 313350 source and SHA contract.
    """
    receipt = json.loads(decoded_receipt.read_text(encoding="utf-8"))
    if (receipt.get("prepared_table_sha256") != sha256(native) or
            receipt.get("native_observations") != expected_n or
            receipt.get("environment_decoding_status") !=
                "physical_units_gsp_invalid_encoding_masked_v1"):
        raise ValueError("S1_ORIGINAL_DECODE_AND_SHA_GATE_FAILED")
    frame = pd.read_csv(native, dtype={"obs_id": str}, low_memory=False,
                        usecols=["obs_id", "latitude", "longitude"])
    if (len(frame) != expected_n or frame.obs_id.isna().any() or
            not frame.obs_id.is_unique or
            not np.isfinite(frame[["latitude", "longitude"]].to_numpy(float)).all() or
            not frame.latitude.between(-90, 90).all() or
            not frame.longitude.between(-180, 180).all()):
        raise ValueError("S1_SOURCE_GRAIN_OR_COORDINATES_INVALID")
    private_dir.mkdir(parents=True, exist_ok=True)
    public_dir.mkdir(parents=True, exist_ok=True)
    temporary.mkdir(parents=True, exist_ok=True)
    original = pd.DataFrame({"obs_id": frame.obs_id.copy()})
    added = pd.DataFrame({"obs_id": frame.obs_id.copy()})
    layers = {}
    receipts = {}
    for variable in SOURCE_VARS:
        if variable == "topo_elevation" or variable in EXTRA:
            target, provenance = original_retrieve(variable, temporary)
        else:
            target, provenance = bio_retrieve(variable, temporary)
        opened = open_raster(target, variable, temporary)
        try:
            identity = check_raster_identity(opened, variable, is_new=variable in BIO9)
            values, numeric_receipt = raster_sample(opened, variable, frame)
            support = _sanity(values, variable)
            if variable in EXTRA:
                original[variable] = values
            else:
                added[variable] = values
                layers[variable] = {
                    "dataset": DATASET_FAMILIES[variable],
                    "temporal_support": EXPECTED_PERIODS[variable],
                    "source_bundle_sha256": provenance["download_sha256"],
                    "source_url": provenance["url"],
                    "unit": PHYSICAL_UNIT[variable],
                    "derivation": ("Original CHELSA v2.1 un-interpolated original pixel; "
                                   "raw GeoTIFF multiplied by native scale and offset"),
                    "sampling_method": "rasterio-nearest-original-WGS84-coordinate",
                    "finite_sampled_observations": support["finite"],
                }
            receipts[variable] = {**identity,
                                   "source_sha256": provenance["download_sha256"],
                                   "source_url": provenance["url"],
                                   "download_bytes": provenance["download_bytes"],
                                   "finite_count": support["finite"],
                                   "coverage": support["fraction_finite"],
                                   "physical_unit": PHYSICAL_UNIT.get(variable,
                                                                      numeric_receipt.get("variable_unit"))}
            # Do not log confidential photo IDs, coordinates or trait values.
            print(json.dumps({"stage": "source_raster_sampled", "variable": variable,
                              "source_ok": True, "n_finite": support["finite"]}), flush=True)
        finally:
            opened.unlink(missing_ok=True)
            if opened != target:
                target.unlink(missing_ok=True)
    if set(original.columns) != {"obs_id", *EXTRA} or set(added.columns) != {"obs_id", *BIO9}:
        raise ValueError("BIO21_SOURCE_VARIABLE_UNIVERSE_DRIFT")
    orig_path = private_dir / "original_seven.csv.gz"
    nine_path = private_dir / "chelsa_extra_nine.csv.gz"
    original_sha = _write_private_csv(original, orig_path)
    extra_sha = _write_private_csv(added, nine_path)
    manifest = {
        "schema": "azami_rf_extra32_source_manifest_v1",
        "table_sha256": extra_sha,
        "layers": layers,
        "non_public_source_data_held_locally": True,
        "raster_bytes_checksum_verified_by_sampler": True,
        "source_photo_ids_and_original_coordinates_not_exportable": True,
    }
    man_path = private_dir / "chelsa_nine_manifest.json"
    man_path.write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n",
                        encoding="utf-8")
    public = {
        "schema": "azami_bio21_source_reconstruction_receipt_20261031_v1",
        "status": "verified_source_pixels_before_fitting",
        "native_s1_table_sha256": sha256(native),
        "original_seven_private_table_sha256": original_sha,
        "extra_nine_private_table_sha256": extra_sha,
        "nine_manifest_sha256": sha256(man_path),
        "source_observations": int(expected_n),
        "source_variables": receipts,
        "private_identifier_or_coordinate_rows_exported": False,
        "source_trait_predictions_exported": False,
    }
    (public_dir / "bio21_source_receipt.json").write_text(
        json.dumps(public, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return {"seven": orig_path, "seven_sha": original_sha,
            "nine": nine_path, "nine_sha": extra_sha,
            "manifest": man_path,
            "public": public}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("native", "decoded-receipt", "private-dir", "public-dir", "temporary"):
        p.add_argument("--" + name, type=Path, required=True)
    args = p.parse_args()
    extract(args.native, args.decoded_receipt, args.private_dir, args.public_dir,
            args.temporary)


if __name__ == "__main__":
    main()
