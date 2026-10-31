"""Synthetic-only source and provenance tests; NO original S1 or CHELSA download."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from analysis.v4.rf_bio21_source_real import (
    BIO9, EXTRA, SOURCE_VARS, PHYSICAL_UNIT,
    chelsa_candidates, extract,
)
from analysis.v4.rf_multitrait_predictability import sha256


def test_registered_variable_set_and_original_chelsa_routes():
    assert len(BIO9) == 9 and len(set(SOURCE_VARS)) == 16
    assert len(EXTRA) == 7
    assert PHYSICAL_UNIT["chelsa_bio03"] == "percent_isothermality"
    assert PHYSICAL_UNIT["chelsa_bio19"] == "kg_m-2_month-1"
    for k in BIO9:
        urls = chelsa_candidates(k)
        assert len(urls) == 2
        assert all("1981-2010" in u and "2.1.tif" in u for u in urls)
        assert all(k.removeprefix("chelsa_") in u for u in urls)
    with pytest.raises(ValueError, match="UNREGISTERED"):
        chelsa_candidates("chelsa_bio05")


def test_original_source_sha_and_identifier_guard(tmp_path):
    source = pd.DataFrame({"obs_id": ["a", "b", "c"],
                           "latitude": [12.1, 13.2, 14.3],
                           "longitude": [111., 112., 113.]})
    raw = tmp_path / "orig.csv.gz"
    source.to_csv(raw, index=False, compression="gzip")
    receipt = tmp_path / "receipt.json"
    receipt.write_text(json.dumps({"prepared_table_sha256": sha256(raw),
                                    "native_observations": 3,
                                    "environment_decoding_status":
                                    "physical_units_gsp_invalid_encoding_masked_v1"}))
    with pytest.raises(ValueError, match="S1_ORIGINAL_DECODE"):
        extract(raw, receipt, tmp_path / "private", tmp_path / "public", tmp_path / "temp",
                expected_n=4)
    receipt.write_text(json.dumps({"prepared_table_sha256": "0" * 64,
                                    "native_observations": 3,
                                    "environment_decoding_status":
                                    "physical_units_gsp_invalid_encoding_masked_v1"}))
    with pytest.raises(ValueError, match="S1_ORIGINAL_DECODE"):
        extract(raw, receipt, tmp_path / "private", tmp_path / "public", tmp_path / "temp",
                expected_n=3)


def test_private_7_plus_9_extract_and_aggregate_only_receipt(tmp_path, monkeypatch):
    import analysis.v4.rf_bio21_source_real as mod
    n = 120
    source = pd.DataFrame({"obs_id": [f"id_{i:04}" for i in range(n)],
                           "latitude": np.linspace(-25, 25, n),
                           "longitude": np.linspace(-65, 65, n)})
    raw = tmp_path / "source.csv.gz"
    source.to_csv(raw, index=False, compression="gzip")
    receipt = tmp_path / "receipt.json"
    receipt.write_text(json.dumps({"prepared_table_sha256": sha256(raw),
                                    "native_observations": n,
                                    "environment_decoding_status":
                                    "physical_units_gsp_invalid_encoding_masked_v1"}))
    counter = {name: i for i, name in enumerate(SOURCE_VARS)}

    def fake_fetch(var, folder):
        path = Path(folder) / f"{var}.tif"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(var.encode())
        return path, {"variable": var, "url": f"https://example.org/{var}.tif",
                      "download_bytes": len(var), "download_sha256": sha256(path)}

    def fake_sample(path, var, frame):
        return np.full(len(frame), counter[var] + 20., dtype=float), {"variable_unit": "test"}

    monkeypatch.setattr(mod, "open_raster", lambda f, var, folder: f)
    monkeypatch.setattr(mod, "check_raster_identity", lambda path, var, is_new:
                        {"crs": "EPSG:4326", "scale": 1, "offset": 0})
    out = extract(raw, receipt, tmp_path / "protected", tmp_path / "receipt_only",
                  tmp_path / "temp", expected_n=n,
                  original_retrieve=fake_fetch, bio_retrieve=fake_fetch,
                  raster_sample=fake_sample)
    seven = pd.read_csv(out["seven"], dtype={"obs_id": str})
    nine = pd.read_csv(out["nine"], dtype={"obs_id": str})
    assert seven.shape == (n, 8) and nine.shape == (n, 10)
    assert set(seven.columns) == {"obs_id", *EXTRA}
    assert set(nine.columns) == {"obs_id", *BIO9}
    assert list(nine.obs_id) == list(source.obs_id)
    assert out["seven_sha"] == sha256(out["seven"])
    assert out["nine_sha"] == sha256(out["nine"])
    manifest = json.loads(out["manifest"].read_text())
    assert set(manifest["layers"]) == set(BIO9)
    assert manifest["table_sha256"] == out["nine_sha"]
    public = json.loads((tmp_path / "receipt_only" / "bio21_source_receipt.json").read_text())
    assert public["source_observations"] == n
    assert len(public["source_variables"]) == 16
    assert "id_0000" not in json.dumps(public)
    assert not {"obs_id", "latitude", "longitude", "accepted_key", "user_id"} & set(public)
