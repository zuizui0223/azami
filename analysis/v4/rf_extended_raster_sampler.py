"""Private-only extraction of 32 *preprepared* environmental rasters at S1 points.

The source recipe resolves each registered feature to ONE concrete raster and
pins its bytes. Monthly -> climatology, DEM -> terrain and 250m -> fraction
calculations must be prepared *before* this point sampler. No raster download or
third-party redistribution is attempted. Outputs contain exact source IDs and
coordinates only in the authorised working environment (never GitHub or SI).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.v4.rf_extended_environment import (
    DATASET_FAMILIES, EXPECTED_PERIODS, NEW32, SOIL6, HEX, added_columns,
)
from analysis.v4.rf_multitrait_predictability import SOURCE_N, sha256


SAMPLING = "rasterio-nearest-original-WGS84-coordinate"
SOURCE_LANDING_PAGES = {
    "CHELSA-bioclim-v2.1": "https://www.chelsa-climate.org/datasets/chelsa_bioclim",
    "TerraClimate": "https://www.climatologylab.org/terraclimate-variables.html",
    "WorldClim-elevation-derived": "https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_30s_elev.zip",
    "SoilGrids-v2.0": "https://isric.org/explore/soilgrids",
    "MOD13A3-v061": "https://ladsweb.modaps.eosdis.nasa.gov/missions-and-measurements/products/MOD13A3",
    "Copernicus-LC100-C3-2015": "https://land.copernicus.eu/en/products/global-dynamic-land-cover/copernicus-global-land-service-land-cover-100m-collection-3-epoch-2015-2019-globe",
}


def recipe_template(upper_tier: str):
    """Generate a NON-EXECUTABLE template with explicit missing-source gates."""
    required = added_columns(upper_tier)
    return {
        "schema": "azami_rf_extended32_raster_recipe_v1",
        "scope": upper_tier,
        "notice": "SOURCE FILES ARE NOT ACQUIRED. Complete, verify and pin all fields before sampling.",
        "layers": {
            name: {
                "dataset": DATASET_FAMILIES[name],
                "temporal_support": EXPECTED_PERIODS[name],
                "raster_file": "",
                "expected_sha256": "",
                "source_url": SOURCE_LANDING_PAGES[DATASET_FAMILIES[name]],
                "unit": "",
                "derivation": "",
                "scale": None,
                "offset": None,
                **({"depth_interval": "0-5cm", "statistic": "Q0.50"}
                   if name in SOIL6 else {}),
            }
            for name in required
        },
    }


def check_recipe(recipe, required_columns: tuple[str, ...] = NEW32):
    if (recipe.get("schema") != "azami_rf_extended32_raster_recipe_v1" or
            set(recipe.get("layers", {})) != set(required_columns)):
        raise ValueError("SOURCE_RASTER_RECIPE_MISSING_OR_UNEXPECTED_LAYERS")
    for name in required_columns:
        spec = recipe["layers"][name]
        if (spec.get("dataset") != DATASET_FAMILIES[name] or
                spec.get("temporal_support") != EXPECTED_PERIODS[name] or
                not isinstance(spec.get("raster_file"), str) or
                not HEX.fullmatch(str(spec.get("expected_sha256", ""))) or
                not all(isinstance(spec.get(v), str) and len(spec[v].strip()) > 2
                        for v in ("source_url", "unit", "derivation")) or
                not isinstance(spec.get("scale"), (float, int)) or
                not isinstance(spec.get("offset"), (float, int)) or
                not np.isfinite(float(spec["scale"])) or
                not np.isfinite(float(spec["offset"])) or
                spec["scale"] == 0):
            raise ValueError("INVALID_RASTER_SOURCE_CONTRACT: " + name)
        if name in SOIL6 and (spec.get("depth_interval") != "0-5cm" or
                              spec.get("statistic") != "Q0.50"):
            raise ValueError("SOIL_DEPTH_OR_STATISTIC_NOT_Q050_0_5CM: " + name)


def _validate_coordinates(source):
    if (source.obs_id.isna().any() or not source.obs_id.is_unique or
            not np.isfinite(source[["latitude", "longitude"]].to_numpy(float)).all() or
            not source.latitude.between(-90, 90).all() or
            not source.longitude.between(-180, 180).all()):
        raise ValueError("SOURCE_COORDINATE_OR_ID_INVALID")


def sample_layer(raster_file: Path, expected_sha256: str, lat, lon,
                 *, scale=1., offset=0.):
    """Return NaN for valid masked/no-data cells; never nearest-observation fill."""
    import rasterio
    from rasterio.warp import transform

    if sha256(raster_file) != expected_sha256:
        raise ValueError("RASTER_SOURCE_BYTES_CHANGED: " + Path(raster_file).name)
    with rasterio.open(raster_file) as dataset:
        if dataset.count != 1 or dataset.crs is None:
            raise ValueError("ONE_BAND_GEOREFERENCED_RASTER_REQUIRED")
        xs, ys = transform("EPSG:4326", dataset.crs, lon.tolist(), lat.tolist())
        values = np.full(len(lat), np.nan, dtype=float)
        for idx, pair in enumerate(dataset.sample(zip(xs, ys), indexes=1, masked=True)):
            raw = pair[0]
            if not np.ma.is_masked(raw):
                q = float(raw) * scale + offset
                if np.isfinite(q):
                    values[idx] = q
        return values


def extract(native_csv: Path, native_receipt: Path, recipe_path: Path,
            destination_table: Path, destination_manifest: Path,
            *, expected_n=SOURCE_N,
            required_columns: tuple[str, ...] = NEW32):
    """Produce a private 32-variable table and hash-locked companion manifest."""
    receipt = json.loads(native_receipt.read_text(encoding="utf-8"))
    if (sha256(native_csv) != receipt.get("prepared_table_sha256") or
            receipt.get("native_observations") != expected_n or
            receipt.get("environment_decoding_status") !=
            "physical_units_gsp_invalid_encoding_masked_v1"):
        raise ValueError("S1_NATIVE_ARCHIVE_AND_DECODE_RECEIPT_NOT_VERIFIED")
    source = pd.read_csv(native_csv, dtype={"obs_id": str}, low_memory=False,
                         usecols=["obs_id", "latitude", "longitude"])
    if len(source) != expected_n:
        raise ValueError("SOURCE_S1_OBSERVATION_COUNT_CHANGED")
    _validate_coordinates(source)
    recipe = json.loads(recipe_path.read_text(encoding="utf-8"))
    check_recipe(recipe, required_columns=required_columns)
    output = pd.DataFrame({"obs_id": source.obs_id.copy()})
    layers = {}
    lat = source.latitude.to_numpy(float)
    lon = source.longitude.to_numpy(float)
    for name in required_columns:
        spec = recipe["layers"][name]
        path = Path(spec["raster_file"])
        values = sample_layer(path, spec["expected_sha256"], lat, lon,
                              scale=float(spec["scale"]), offset=float(spec["offset"]))
        if len(values) != len(output):
            raise ValueError("ENVIRONMENT_RASTER_SAMPLE_COUNT_DRIFT")
        output[name] = values
        layers[name] = {
            "dataset": spec["dataset"],
            "temporal_support": spec["temporal_support"],
            "source_bundle_sha256": spec["expected_sha256"],
            "source_url": spec["source_url"], "unit": spec["unit"],
            "derivation": spec["derivation"], "sampling_method": SAMPLING,
            "finite_sampled_observations": int(np.isfinite(values).sum()),
            **({"depth_interval": "0-5cm", "statistic": "Q0.50"}
               if name in SOIL6 else {}),
        }
    destination_table.parent.mkdir(parents=True, exist_ok=True)
    destination_manifest.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(destination_table, index=False)
    manifest = {
        "schema": "azami_rf_extra32_source_manifest_v1",
        "table_sha256": sha256(destination_table),
        "layers": layers,
        "non_public_source_data_held_locally": True,
        "raster_bytes_checksum_verified_by_sampler": True,
        "source_photo_ids_and_original_coordinates_not_exportable": True,
    }
    destination_manifest.write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n",
                                    encoding="utf-8")
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--create-recipe-template", type=Path)
    parser.add_argument("--upper-tier", default="complete_bioclim21",
                        choices=["complete_bioclim21", "water_balance26",
                                 "wind28", "terrain33", "soil39", "vegetation44"])
    for name in ("native-table", "native-receipt", "recipe", "output-table", "output-manifest"):
        parser.add_argument("--" + name, type=Path)
    args = parser.parse_args()
    if args.create_recipe_template is not None:
        args.create_recipe_template.parent.mkdir(parents=True, exist_ok=True)
        args.create_recipe_template.write_text(
            json.dumps(recipe_template(args.upper_tier), indent=2) + "\n")
        return
    if any(getattr(args, arg.replace("-", "_")) is None for arg in
           ("native-table", "native-receipt", "recipe", "output-table", "output-manifest")):
        parser.error("source extraction needs native table/receipt, recipe, and output paths")
    extract(args.native_table, args.native_receipt, args.recipe,
            args.output_table, args.output_manifest,
            required_columns=added_columns(args.upper_tier))


if __name__ == "__main__":
    main()