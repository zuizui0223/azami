"""Synthetic rasters only: byte-pinned source extraction, masked pixels, and QC."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

rasterio = pytest.importorskip("rasterio")

from rasterio.transform import from_origin
from analysis.v4.rf_extended_environment import (
    NEW32, DATASET_FAMILIES, EXPECTED_PERIODS, SOIL6, verify_extended_source, added_columns,
)
from analysis.v4.rf_extended_raster_sampler import extract, sample_layer, check_recipe, recipe_template
from analysis.v4.rf_multitrait_predictability import sha256


def prepare(tmp_path):
    raster = tmp_path / "synthetic_environment.tif"
    with rasterio.open(raster, "w", driver="GTiff", width=2, height=2,
                       count=1, dtype="float32", crs="EPSG:4326",
                       transform=from_origin(-180, 90, 180, 90), nodata=-9999) as d:
        d.write(np.array([[1., 2.], [3., -9999.]], dtype="float32"), 1)
    source = tmp_path / "native.csv"
    pd.DataFrame({
        "obs_id": ["a", "b", "c", "d"],
        "latitude": [45, 45, -45, -45],
        "longitude": [-90, 90, -90, 90],
    }).to_csv(source, index=False)
    receipt_path = tmp_path / "receipt.json"
    receipt_path.write_text(json.dumps({
        "prepared_table_sha256": sha256(source), "native_observations": 4,
        "environment_decoding_status": "physical_units_gsp_invalid_encoding_masked_v1",
    }))
    fingerprint = sha256(raster)
    recipe = {
        "schema": "azami_rf_extended32_raster_recipe_v1",
        "layers": {name: {
            "dataset": DATASET_FAMILIES[name],
            "temporal_support": EXPECTED_PERIODS[name],
            "raster_file": str(raster),
            "expected_sha256": fingerprint,
            "source_url": "https://example.org/mock-source",
            "unit": "source physical unit",
            "derivation": "test only, synthetic single raster",
            "scale": 1.0, "offset": 0.0,
            **({"depth_interval": "0-5cm", "statistic": "Q0.50"}
               if name in SOIL6 else {}),
        } for name in NEW32},
    }
    recipe_path = tmp_path / "recipe.json"
    recipe_path.write_text(json.dumps(recipe))
    return raster, source, receipt_path, recipe_path


def test_pixel_source_sampler_nodata(tmp_path):
    raster, source, receipt, recipe = prepare(tmp_path)
    y = sample_layer(raster, sha256(raster), np.array([45, 45, -45, -45]),
                     np.array([-90, 90, -90, 90]))
    assert np.allclose(y[:3], [1, 2, 3])
    assert np.isnan(y[3])


def test_raster_source_extraction_hashes_and_private_row_join(tmp_path):
    raster, source, receipt, recipe = prepare(tmp_path)
    extra_csv = tmp_path / "private" / "extra.csv"
    extra_receipt = tmp_path / "private" / "manifest.json"
    manifest = extract(source, receipt, recipe, extra_csv, extra_receipt, expected_n=4)
    assert manifest["raster_bytes_checksum_verified_by_sampler"]
    assert manifest["layers"][NEW32[0]]["finite_sampled_observations"] == 3
    assert "obs_id" not in str(manifest["layers"][NEW32[0]])
    native = pd.read_csv(source, dtype={"obs_id": str})
    joined, provenance = verify_extended_source(native, extra_csv, sha256(extra_csv), extra_receipt)
    assert joined.obs_id.tolist() == native.obs_id.tolist()
    assert joined.loc[3, NEW32[0]] != joined.loc[3, NEW32[0]]  # NaN
    assert provenance["manifest_sha256"] == sha256(extra_receipt)


def test_raster_bytes_changed_fails_closed(tmp_path):
    raster, source, receipt, recipe = prepare(tmp_path)
    r = json.loads(recipe.read_text())
    r["layers"][NEW32[0]]["expected_sha256"] = "0"*64
    recipe.write_text(json.dumps(r))
    with pytest.raises(ValueError, match="RASTER_SOURCE_BYTES_CHANGED"):
        extract(source, receipt, recipe, tmp_path / "out.csv", tmp_path / "out.json", expected_n=4)


def test_recipe_missing_layer_rejected(tmp_path):
    raster, source, receipt, recipe = prepare(tmp_path)
    r = json.loads(recipe.read_text())
    del r["layers"][NEW32[0]]
    with pytest.raises(ValueError, match="MISSING_OR_UNEXPECTED"):
        check_recipe(r)


def test_s1_source_sha_guard(tmp_path):
    raster, source, receipt, recipe = prepare(tmp_path)
    p = json.loads(receipt.read_text())
    p["prepared_table_sha256"] = "0"*64
    receipt.write_text(json.dumps(p))
    with pytest.raises(ValueError, match="ARCHIVE_AND_DECODE"):
        extract(source, receipt, recipe, tmp_path / "out.csv", tmp_path / "out.json", expected_n=4)


def test_21_feature_template_is_explicitly_non_executable():
    blueprint = recipe_template("complete_bioclim21")
    assert len(blueprint["layers"]) == 9
    assert set(blueprint["layers"]) == {"chelsa_bio02", "chelsa_bio03",
                                             "chelsa_bio07", "chelsa_bio08",
                                             "chelsa_bio09", "chelsa_bio10",
                                             "chelsa_bio11", "chelsa_bio16",
                                             "chelsa_bio19"}
    with pytest.raises(ValueError, match="INVALID_RASTER_SOURCE_CONTRACT"):
        check_recipe(blueprint, required_columns=tuple(blueprint["layers"]))


def test_can_extract_bioclim21_without_other_23_future_layers(tmp_path):
    raster, source, receipt, recipe_path = prepare(tmp_path)
    recipe = json.loads(recipe_path.read_text())
    only = added_columns("complete_bioclim21")
    recipe["layers"] = {k: recipe["layers"][k] for k in only}
    recipe_path.write_text(json.dumps(recipe))
    table = tmp_path / "only9.csv"
    manifest_path = tmp_path / "only9.json"
    extract(source, receipt, recipe_path, table, manifest_path, expected_n=4,
            required_columns=only)
    assert set(pd.read_csv(table).columns) == {"obs_id", *only}
    assert len(json.loads(manifest_path.read_text())["layers"]) == 9