"""Official SoilGrids 0-5cm Q0.5 WebDAV source-availability audit.

NO photo, location, species, observation, or user data are used here. Six
public VRT master files are obtained, validated as *VRT* rather than a 200-OK
HTML error, and their exact bytes hashed. This is **NOT** a soil observation
extraction and NOT an RF result. A single public, fixed test point may be
sampled to check whether GDAL /vsicurl can access the underlying tiles.

The project GEB manuscript, S1/S2/H3c and existing 23/28 RF stay immutable.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from urllib.parse import urlparse
import xml.etree.ElementTree as ET

from analysis.v4.rf_extended_environment import SOIL6

ROOT = "https://files.isric.org/soilgrids/latest/data"
DEPTH = "0-5cm"
QUANTILE = "Q0.5"
MAX_VRT_BYTES = 16_000_000
# SoilGrids official integer-to-conventional-unit factors; do not infer from
# GDAL's Scale/Offset because the SoilGrids source defines mapped units.
CONVERSIONS = {
    "soil_phh2o_0_5cm": ("phh2o", 10, "pH"),
    "soil_soc_0_5cm": ("soc", 10, "g/kg"),
    "soil_clay_0_5cm": ("clay", 10, "percent_mass"),
    "soil_sand_0_5cm": ("sand", 10, "percent_mass"),
    "soil_bdod_0_5cm": ("bdod", 100, "kg/dm3"),
    "soil_cec_0_5cm": ("cec", 10, "cmol_c_per_kg"),
}
assert set(CONVERSIONS) == set(SOIL6)


def official_url(registered_name: str) -> str:
    if registered_name not in CONVERSIONS:
        raise ValueError("UNREGISTERED_SOILGRIDS_PROPERTY")
    property_id = CONVERSIONS[registered_name][0]
    return f"{ROOT}/{property_id}/{property_id}_{DEPTH}_{QUANTILE}.vrt"


def parsed_vrt_metadata(content: bytes, expected_property: str) -> dict:
    """Require a complete single-band GDAL VRT with raster sources and SRS."""
    if not content or len(content) > MAX_VRT_BYTES:
        raise ValueError("SOILGRIDS_VRT_MISSING_OR_OVER_BUDGET")
    # Defuse external entities, large node attributes, or bizarre XML files.
    if b"<!DOCTYPE" in content[:4000] or b"<!ENTITY" in content[:4000]:
        raise ValueError("SOILGRIDS_VRT_ENTITY_DEFINITION_REFUSED")
    try:
        root = ET.fromstring(content)
    except ET.ParseError as exc:
        raise ValueError("SOILGRIDS_RETURNED_NOT_VALID_VRT_XML") from exc
    if root.tag != "VRTDataset":
        raise ValueError("SOILGRIDS_RETURNED_NON_VRT_PAYLOAD")
    bands = root.findall("VRTRasterBand")
    if len(bands) != 1 or not root.get("rasterXSize") or not root.get("rasterYSize"):
        raise ValueError("SOILGRIDS_VRT_BAND_OR_GRID_UNEXPECTED")
    srs = root.findtext("SRS")
    transform = root.findtext("GeoTransform")
    if not srs or not transform or len(transform.split(",")) != 6:
        raise ValueError("SOILGRIDS_VRT_SRS_OR_GEOTRANSFORM_MISSING")
    paths = [x.text for x in root.iter("SourceFilename") if x.text]
    if not paths or not all(".tif" in x.lower() for x in paths):
        raise ValueError("SOILGRIDS_VRT_NO_TIFF_TILES")
    # The source VRT is allowed to use paths with relative filenames. Never
    # follow any of those paths in this public-only metadata verification.
    if not any(expected_property in p for p in paths):
        raise ValueError("SOILGRIDS_VRT_WRONG_PROPERTY_TILE_REFERENCE")
    return {"raster_width": int(root.get("rasterXSize")),
            "raster_height": int(root.get("rasterYSize")),
            "band_count": len(bands),
            "raster_tile_references": len(paths),
            "srs_present": True, "geotransform_present": True,
            "raw_source_tiles_not_downloaded": True}


def _download_vrt(url: str, *, max_bytes: int = MAX_VRT_BYTES) -> tuple[bytes, dict]:
    import requests
    parsed = urlparse(url)
    if parsed.scheme != "https" or parsed.netloc != "files.isric.org" or not url.endswith("_Q0.5.vrt"):
        raise ValueError("ONLY_OFFICIAL_SOILGRIDS_WEBDAV_URL_ALLOWED")
    with requests.get(url, stream=True, timeout=(12, 28), allow_redirects=False,
                      headers={"User-Agent": "azami-v4-source-audit/20261031"}) as response:
        response.raise_for_status()
        payload = bytearray()
        for part in response.iter_content(256 * 1024):
            payload.extend(part)
            if len(payload) > max_bytes:
                raise ValueError("SOILGRIDS_VRT_EXCEEDED_TRANSFER_BUDGET")
        return bytes(payload), {
            "http_status": response.status_code,
            "server_etag": response.headers.get("ETag"),
            "server_last_modified": response.headers.get("Last-Modified"),
            "content_type": response.headers.get("Content-Type")}


def fixed_public_sentinel(url: str, scale_divisor: int) -> dict:
    """Optional ONLY Amsterdam, NL public probe: no private S1 positions."""
    import rasterio
    import numpy as np
    from rasterio.warp import transform
    # /vsicurl? is the ISRIC documented way to access their own master VRT.
    virtual = "/vsicurl?max_retry=1&retry_delay=1&list_dir=no&url=" + url
    with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
                      GDAL_HTTP_MAX_RETRY="1", GDAL_HTTP_RETRY_DELAY="1",
                      GDAL_HTTP_TIMEOUT="12", GDAL_HTTP_CONNECTTIMEOUT="8"):
        with rasterio.open(virtual) as src:
            if src.count != 1 or src.crs is None:
                raise ValueError("SOILGRIDS_REMOTE_VRT_NOT_GEOREFERENCED")
            xs, ys = transform("EPSG:4326", src.crs, [4.8952], [52.3702])
            a = list(src.sample([(xs[0], ys[0])], indexes=1, masked=True))
            if len(a) != 1 or np.ma.is_masked(a[0][0]):
                return {"sentinel_sample_status": "nodata_at_public_test_point"}
            raw = float(a[0][0])
            if not np.isfinite(raw):
                return {"sentinel_sample_status": "invalid_at_public_test_point"}
            return {"sentinel_sample_status": "one_original_pixel_decoded",
                    "sentinel_physical_value": raw / scale_divisor,
                    "public_point_only": True}


def audit(*, fetcher=None, sampler=None, sample_property="soil_phh2o_0_5cm") -> dict:
    """Fail-closed per-source status and an overall readiness gate."""
    if fetcher is None:
        fetcher = _download_vrt
    if sample_property not in SOIL6:
        raise ValueError("INVALID_PUBLIC_SENTINEL_PROPERTY")
    results = {}
    for name in SOIL6:
        code, divisor, unit = CONVERSIONS[name]
        url = official_url(name)
        entry = {"source_url": url, "depth_interval": DEPTH,
                 "statistic": QUANTILE, "mapped_value_divisor": divisor,
                 "conventional_unit": unit, "data_family": "SoilGrids-2.0",
                 "source_is_modelled_soil_not_field_sampled": True}
        try:
            content, headers = fetcher(url)
            entry.update(headers)
            entry.update(parsed_vrt_metadata(content, code))
            entry["vrt_file_sha256"] = hashlib.sha256(content).hexdigest()
            entry["vrt_bytes"] = len(content)
            entry["status"] = "authentic_webdav_vrt_bytes_verified"
        except Exception as exc:
            entry["status"] = "unavailable_or_invalid_source"
            entry["failure_type"] = type(exc).__name__
            entry["failure_message_code"] = str(exc)[:180] if isinstance(exc, ValueError) else "network_failure"
        results[name] = entry
    all6 = all(r["status"] == "authentic_webdav_vrt_bytes_verified" for r in results.values())
    sentinel = {"status": "not_requested"}
    if sampler is not None and results[sample_property]["status"] == "authentic_webdav_vrt_bytes_verified":
        try:
            value = sampler(official_url(sample_property), CONVERSIONS[sample_property][1])
            sentinel = dict(value)
            sentinel["registered_property"] = sample_property
        except Exception as exc:
            sentinel = {"status": "sampling_failed", "failure_type": type(exc).__name__,
                        "no_source_observation_sampled": True}
    return {"schema": "azami_soilgrids_official_webdav_probe_v1",
            "status": "all_six_vrt_authenticated_no_original_s1_extracted" if all6 else "at_least_one_soil_vrt_unavailable",
            "requested_layer_count": 6,
            "authentic_vrt_layers": sum(r["status"] == "authentic_webdav_vrt_bytes_verified" for r in results.values()),
            "source": "ISRIC SoilGrids 2.0 public WebDAV VRT 0-5cm Q0.5",
            "no_private_source_access": True,
            "not_an_RF_model_or_soil_extraction": True,
            "public_sentinel": sentinel,
            "layers": results}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--public-sentinel", action="store_true")
    a = p.parse_args()
    out = audit(sampler=fixed_public_sentinel if a.public_sentinel else None)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(out, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
    print(json.dumps({"schema": out["schema"], "layers_verified": out["authentic_vrt_layers"],
                      "status": out["status"], "public_sentinel": out["public_sentinel"]}), flush=True)
    if out["authentic_vrt_layers"] != 6:
        raise SystemExit("SOILGRIDS_OFFICIAL_SOURCE_AUDIT_INCOMPLETE")


if __name__ == "__main__":
    main()