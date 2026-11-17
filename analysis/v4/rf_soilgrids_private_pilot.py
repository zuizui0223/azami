"""Protected SoilGrids 0-5cm Q0.5 pixel-sampling feasibility gate.

Uses checksum-authenticated original S1 inside GitHub Actions, selects a fixed
outcome-blind sample of source coordinates stratified into three broad regions,
reads original SoilGrids VRT pixels via /vsicurl, and emits only aggregate
validity, source SHA and physical-unit sanity ranges. It DOES NOT fit an RF and
DOES NOT upload source coordinates, IDs, taxon names, or pixel-level values.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.v4.rf_multitrait_predictability import sha256, SOURCE_N
from analysis.v4.rf_soilgrids_webdav_probe import (
    SOIL6, CONVERSIONS, official_url, parsed_vrt_metadata, _download_vrt,
)

# Hashes of the exact official 2026-10-31 six VRT source archives; a changed
# upstream 'latest' cannot silently alter the species/photo comparison.
VRT_PINNED_SHA256 = {
    "soil_phh2o_0_5cm": "f916ab12288eb895443f95830ce8dbdf72baccd3829819f5aae3cbf0a6f5c804",
    "soil_soc_0_5cm": "f5286f3f4cd502eacf0b0a6ec608f3b009580820017da6ba87300622824903a2",
    "soil_clay_0_5cm": "0faf370ff572b800e5d4ca090defba8882280fc6b9e81ef3bb29a2f53d040443",
    "soil_sand_0_5cm": "5e5279fc1c2903681712f386f05e6ab26f488cddd90709b87b4d35020f9abd2b",
    "soil_bdod_0_5cm": "9a13f440b5a8cfc961812ec74ddfdc4274c80468857658c5c1d4368b173dfe33",
    "soil_cec_0_5cm": "d26fca3d716511ff74cf8da638a780b038d88916ae0401b48b460d06309b7990",
}
PHYSICAL_BOUNDS = {
    "soil_phh2o_0_5cm": (0, 14),
    "soil_soc_0_5cm": (0, 600),
    "soil_clay_0_5cm": (0, 100),
    "soil_sand_0_5cm": (0, 100),
    "soil_bdod_0_5cm": (0, 3.5),
    "soil_cec_0_5cm": (0, 250),
}
SEED = "azami-original-s1-soilgrids-public-gate-20261031-v1"
assert set(VRT_PINNED_SHA256) == set(SOIL6) == set(PHYSICAL_BOUNDS)


def validated_roster(native: Path, receipt_path: Path, *, expected_n=SOURCE_N) -> pd.DataFrame:
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    if (sha256(native) != receipt.get("prepared_table_sha256") or
            receipt.get("native_observations") != expected_n or
            receipt.get("environment_decoding_status") !=
            "physical_units_gsp_invalid_encoding_masked_v1"):
        raise ValueError("S1_ORIGINAL_SOURCE_SHA_DECODE_DENOMINATOR_NOT_AUTHENTICATED")
    frame = pd.read_csv(native, dtype={"obs_id": str}, low_memory=False,
                        usecols=["obs_id", "latitude", "longitude"])
    if (len(frame) != expected_n or frame.obs_id.isna().any() or
            not frame.obs_id.is_unique):
        raise ValueError("ORIGINAL_S1_GRAIN_OR_IDENTITY_INVALID")
    coord = frame[["latitude", "longitude"]].to_numpy(float)
    if not (np.isfinite(coord).all() and
            (np.abs(coord[:, 0]) <= 90).all() and
            (np.abs(coord[:, 1]) <= 180).all()):
        raise ValueError("ORIGINAL_S1_COORDINATE_INVALID")
    return frame


def source_sentinels(frame: pd.DataFrame, *, per_sector=12) -> tuple[pd.DataFrame, dict]:
    """Outcome-blind, 3-sector deterministic unique-coordinate source sample."""
    if per_sector < 3 or per_sector > 32:
        raise ValueError("SOILGRIDS_PILOT_SECTOR_SAMPLE_BUDGET_INVALID")
    df = frame.loc[frame.latitude.abs() <= 80,
                   ["obs_id", "latitude", "longitude"]].copy()
    df["sector"] = np.select((df.longitude < -30, df.longitude < 60),
                              (0, 1), default=2).astype(int)
    # Stable cryptographic ranking of ids: not by exposure, outcome or location.
    df["rank"] = df.obs_id.map(lambda x:
        hashlib.blake2b((SEED+str(x)).encode("utf-8"), digest_size=8).hexdigest())
    chosen = []
    counts = {}
    for sector in (0, 1, 2):
        subset = df.loc[df.sector == sector].sort_values("rank", kind="stable")
        subset = subset.drop_duplicates(["latitude", "longitude"])
        if len(subset) < per_sector:
            raise ValueError("INSUFFICIENT_12_SOURCE_POINTS_IN_A_LONGITUDE_SECTOR")
        chosen.append(subset.head(per_sector))
        counts[str(sector)] = int(per_sector)
    out = pd.concat(chosen, ignore_index=True)
    if out[["latitude", "longitude"]].duplicated().any():
        raise ValueError("SOILGRIDS_SELECTED_COORDINATE_DUPLICATE")
    return out, {"sample_count": int(len(out)),
                 "source_sector_observations": counts,
                 "selection_uses_trait_or_environment_values": False,
                 "selected_ids_or_coordinates_exported": False}


def raster_points_reader(url: str, latitude, longitude) -> np.ndarray:
    """Exact native VRT pixels; no interpolation or point-by-point HTTP API."""
    import rasterio
    from rasterio.warp import transform
    virtual = "/vsicurl?max_retry=1&retry_delay=1&list_dir=no&url=" + url
    with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
                      GDAL_HTTP_MAX_RETRY="1", GDAL_HTTP_RETRY_DELAY="1",
                      GDAL_HTTP_TIMEOUT="15", GDAL_HTTP_CONNECTTIMEOUT="9"):
        with rasterio.open(virtual) as src:
            if src.count != 1 or src.crs is None:
                raise ValueError("SOILGRIDS_SOURCE_VRT_NOT_GEOREFERENCED")
            xs, ys = transform("EPSG:4326", src.crs,
                               np.asarray(longitude).tolist(),
                               np.asarray(latitude).tolist())
            vals = np.full(len(xs), np.nan)
            # Spatial sorting reduces repeated HTTP requests to the same TIFF.
            order = np.lexsort((np.asarray(ys), np.asarray(xs)))
            test_points = np.column_stack((xs, ys))
            for idx, a in zip(order, src.sample(test_points[order], indexes=1, masked=True)):
                raw = a[0]
                if not np.ma.is_masked(raw) and np.isfinite(float(raw)):
                    vals[idx] = float(raw)
            return vals


def audit_pilot(frame: pd.DataFrame, *, per_sector=12, fetcher=None, pixel_reader=None) -> dict:
    """All incoming rows private. Export only coarse summaries."""
    fetcher = _download_vrt if fetcher is None else fetcher
    pixel_reader = raster_points_reader if pixel_reader is None else pixel_reader
    sample, selection = source_sentinels(frame, per_sector=per_sector)
    values_by_property = {}
    outputs = {}
    for col in SOIL6:
        t0 = time.monotonic()
        code, factor, unit = CONVERSIONS[col]
        url = official_url(col)
        info = {"property": code, "unit": unit,
                "depth_interval": "0-5cm", "statistic": "Q0.5",
                "source_vrt_sha256_expected": VRT_PINNED_SHA256[col],
                "source_original_pixel_never_imputed": True}
        try:
            data, headers = fetcher(url)
            got_hash = hashlib.sha256(data).hexdigest()
            if got_hash != VRT_PINNED_SHA256[col]:
                raise ValueError("SOILGRIDS_SOURCE_VRT_SHA_CHANGED")
            info.update(parsed_vrt_metadata(data, code))
            raw = np.asarray(pixel_reader(url, sample.latitude.to_numpy(float),
                                              sample.longitude.to_numpy(float)), dtype=float)
            if raw.shape != (len(sample),):
                raise ValueError("SOILGRIDS_SAMPLE_VECTOR_LENGTH_DRIFT")
            if np.isinf(raw).any():
                raise ValueError("SOILGRIDS_INFINITE_RAW_PIXEL")
            v = raw / factor
            good = np.isfinite(v)
            lo, hi = PHYSICAL_BOUNDS[col]
            plausible = good & (v >= lo) & (v <= hi)
            info.update({"finite_original_pixel_n": int(good.sum()),
                         "finite_fraction": float(good.mean()),
                         "implausible_pixel_n": int((good & ~plausible).sum()),
                         "status": "sampled_authentic_original_soilgrids_pixels",
                         "elapsed_seconds": round(time.monotonic()-t0, 2)})
            if good.any():
                info["original_pixel_physical_p05_p50_p95"] = [
                    float(x) for x in np.quantile(v[good], [.05, .5, .95])]
            values_by_property[col] = v
        except Exception as error:
            info["status"] = "source_or_sample_rejected"
            info["failure_type"] = type(error).__name__
            info["failure_code"] = str(error)[:170] if isinstance(error, ValueError) else "network_or_raster_error"
        outputs[col] = info
        print(json.dumps({"layer": col, "status": info["status"],
                          "finite": info.get("finite_original_pixel_n", 0),
                          "elapsed": info.get("elapsed_seconds")}), flush=True)
    all_pass = all(x["status"] == "sampled_authentic_original_soilgrids_pixels" and
                   x.get("finite_fraction", 0) >= .5 and
                   x.get("implausible_pixel_n", 1) == 0 for x in outputs.values())
    joint = 0
    if len(values_by_property) == len(SOIL6):
        z = np.column_stack([values_by_property[k] for k in SOIL6])
        joint = int(np.isfinite(z).all(axis=1).sum())
    return {"schema": "azami_v4_protected_source_S1_soilgrids_six_pilot_v1",
            "status": "pilot_source_pixels_validated_not_full_cohort" if all_pass else "insufficient_or_inconsistent_soil_sample",
            "source_original_s1_n": int(len(frame)),
            "source_original_s1_checksum_protected": True,
            "selection": selection,
            "six_soil_properties_covered": int(sum(x["status"] ==
                "sampled_authentic_original_soilgrids_pixels" for x in outputs.values())),
            "complete_six_original_pixel_pilot_n": joint,
            "private_s1_coordinates_records_and_taxa_never_exported": True,
            "not_a_Random_Forest_or_confirmatory_ecological_result": True,
            "prior_GEB_S1_S2_H3c_and_RF_results_unmodified": True,
            "layers": outputs}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--native", type=Path, required=True)
    p.add_argument("--receipt", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--per-sector", type=int, default=12)
    a = p.parse_args()
    df = validated_roster(a.native, a.receipt)
    out = audit_pilot(df, per_sector=a.per_sector)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(out, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
    print(json.dumps({"status":out["status"], "joint_valid": out["complete_six_original_pixel_pilot_n"],
                      "n_layers":out["six_soil_properties_covered"]}),flush=True)
    if out["status"] != "pilot_source_pixels_validated_not_full_cohort":
        raise SystemExit("PROTECTED_S1_SOILGRIDS_SOURCE_PILOT_QUALITY_GATE_FAILED")

if __name__ == "__main__":
    main()