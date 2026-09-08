import copy
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from analysis.v3.production_environment import (CONTRACT, SOURCE_COLUMNS, checkpoint,
    choose_process_variables, sample_blocks, source_frame, task_plan, validate_contract)
from analysis.v3.environment_matrix import load_contract, sha256_file


def specification():
    return json.loads(CONTRACT.read_text(encoding="utf-8"))


def environment():
    rng = np.random.default_rng(8945)
    contract = specification()
    data = pd.DataFrame(rng.normal(size=(1200, 15)), columns=contract["acquire"])
    data["accepted_key"] = np.repeat(["a", "b", "c"], 400)
    return data


def test_no_trait_fields_and_exact_candidate_accounting():
    contract = specification()
    validate_contract(contract)
    assert not any("trait" in name or "chroma" in name for name in SOURCE_COLUMNS)
    assert contract["source"]["observations"] == 319244
    assert contract["source"]["thinning_or_taxon_cap"] is False
    assert contract["model"]["primary_tests"] == 36
    assert set(contract["hypothesis_evidence"]) == set(contract["processes"])


def test_weighted_vif_removes_duplicate_proxy_not_protected_exposure():
    data = environment()
    data["cmi_month"] = data["pr_month"]
    result = choose_process_variables(data, specification())
    assert result["status"].startswith("FULL_NATIVE")
    assert "cmi_month" in result["removed"]
    assert "pr_month" in result["selected"]
    assert "tasmax_month" in result["selected"]
    assert set(result["processes"]) == set(specification()["processes"])


def test_protected_collinearity_blocks_rather_than_claiming_pass():
    data = environment()
    data["vpd_month"] = data["pr_month"]
    result = choose_process_variables(data, specification())
    assert result["status"] == "ENVIRONMENT_SELECTION_NOT_ESTIMABLE"
    assert result["selected"] == []


def test_fixed_complete_case_weights_resist_dominant_taxon_replication():
    data = environment()
    data.loc[0, "BIO12"] = np.nan
    first = choose_process_variables(data, specification())
    extra = data[data.accepted_key.eq("a")]
    repeated = pd.concat([data, extra, extra, extra], ignore_index=True)
    second = choose_process_variables(repeated, specification())
    assert first["complete_rows"] == len(data) - 1
    assert first["selected"] == second["selected"]
    for variable in first["selected"]:
        assert first["centers"][variable] == pytest.approx(second["centers"][variable])
        assert first["scales"][variable] == pytest.approx(second["scales"][variable])


@pytest.mark.parametrize("mode", ["constant", "all_missing"])
def test_unavailable_exposure_does_not_become_an_ecological_negative(mode):
    data = environment()
    data["pr_month"] = 1 if mode == "constant" else np.nan
    assert choose_process_variables(data, specification())["selected"] == []


def test_source_hash_membership_and_allowlisted_read(tmp_path):
    contract = copy.deepcopy(specification())
    data = pd.DataFrame({name: [0, 0] for name in SOURCE_COLUMNS})
    data["obs_id"], data["accepted_key"] = ["2", "1"], ["a", "b"]
    data["observation_month"], data["native_range_status"] = [1, 12], "native"
    data["secret_trait_chroma"] = [999, 1000]
    path = tmp_path / "source.csv"
    data.to_csv(path, index=False)
    contract["source"].update(enriched_csv_sha256=sha256_file(path), observations=2, taxa=2)
    loaded = source_frame(path, contract)
    assert "secret_trait_chroma" not in loaded
    assert loaded["obs_id"].tolist() == ["1", "2"]
    contract["source"]["enriched_csv_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="hash"):
        source_frame(path, contract)


def test_task_plan_preserves_all_observations_per_variable():
    frame = pd.DataFrame({"observation_month": [1, 1, 7, 12]})
    tasks = task_plan(frame, load_contract(), specification()["acquire"])
    assert len(tasks) == 8 * 3 + 7
    for variable in specification()["acquire"]:
        indices = np.concatenate([t["indices"] for t in tasks if t["variable"] == variable])
        assert sorted(indices) == list(range(len(frame)))


def raster_fixture(tmp_path):
    rasterio = pytest.importorskip("rasterio")
    path = tmp_path / "raster.tif"
    raw = np.arange(32 * 32, dtype="int16").reshape(32, 32)
    raw[4, 5] = -999
    with rasterio.open(path, "w", driver="GTiff", height=32, width=32, count=1,
                       dtype="int16", crs="EPSG:4326", nodata=-999, tiled=True,
                       blockxsize=16, blockysize=16,
                       transform=rasterio.transform.from_origin(0, 32, 1, 1)) as dataset:
        dataset.write(raw, 1)
        dataset.scales = [0.1]
        dataset.offsets = [-20]
    return path


def test_block_sampling_matches_reference_values_masks_scales_and_bounds(tmp_path):
    rasterio = pytest.importorskip("rasterio")
    path = raster_fixture(tmp_path)
    xy = np.array([[0.5, 31.5], [5.5, 27.5], [19.5, 12.5], [19.5, 12.5], [-1, 33]])
    with rasterio.open(path) as dataset:
        actual, report = sample_blocks(dataset, xy)
        expected = np.ma.concatenate(list(dataset.sample(xy, masked=True))).astype(float)
        expected = expected.filled(np.nan) * 0.1 - 20
    np.testing.assert_allclose(actual, expected, equal_nan=True)
    assert report["out_of_bounds"] == 1
    assert report["blocks_read"] == 2


def test_checkpoint_offline_resume_and_corruption(tmp_path):
    pytest.importorskip("rasterio")
    path = raster_fixture(tmp_path)
    directory = tmp_path / "checkpoint"
    directory.mkdir()
    frame = pd.DataFrame({"analysis_longitude": [0.5], "analysis_latitude": [31.5]})
    task = {"id": "x_01", "url": str(path), "indices": np.array([0])}
    identity = lambda url: {"url": url, "etag": "synthetic", "bytes": path.stat().st_size}
    first = checkpoint(task, frame, directory, "run", identity_reader=identity)
    def forbid_network(url):
        raise AssertionError("Committed checkpoints must replay offline")
    assert checkpoint(task, frame, directory, "run", identity_reader=forbid_network) == first
    (directory / "x_01.npz").write_bytes(b"corrupt")
    with pytest.raises(ValueError, match="hash"):
        checkpoint(task, frame, directory, "run", identity_reader=forbid_network)


def test_changed_remote_object_never_commits_checkpoint(tmp_path):
    pytest.importorskip("rasterio")
    path = raster_fixture(tmp_path)
    directory = tmp_path / "checkpoint"
    directory.mkdir()
    frame = pd.DataFrame({"analysis_longitude": [0.5], "analysis_latitude": [31.5]})
    task = {"id": "x_01", "url": str(path), "indices": np.array([0])}
    calls = iter(["before", "after"])
    with pytest.raises(ValueError, match="changed during"):
        checkpoint(task, frame, directory, "run", identity_reader=lambda url: {"etag": next(calls)})
    assert not list(directory.iterdir())
