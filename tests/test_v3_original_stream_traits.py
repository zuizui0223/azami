import csv
import hashlib
import io
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import pytest
from PIL import Image

from analysis.v3 import stream_original_traits as stream

from analysis.v3.stream_original_traits import (
    photo_schedule,
    qualified_endpoints,
    select_observation_ids,
)


ROOT = Path(__file__).resolve().parents[1]
DECISION = ROOT / "analysis" / "v3" / "measurement_qualification_decision_20260908.json"


def test_decision_exposes_exactly_14_original_stream_endpoints():
    decision = json.loads(DECISION.read_text(encoding="utf-8"))
    qualified, bbox = qualified_endpoints(decision)
    assert len(qualified) == 14
    assert "orientation_image_vertical_angle" in qualified
    assert "corolla_lab_chroma" in qualified
    assert "visible_floret_fraction" in qualified
    assert set(bbox) == {
        "orientation_image_vertical_angle",
        "capitulum_outline_aspect_ratio",
        "capitulum_outline_circularity",
        "capitulum_outline_solidity",
        "capitulum_width_profile_cv",
    }


def test_observation_selection_is_deterministic_and_sharded(tmp_path):
    cohort = tmp_path / "cohort.csv"
    with cohort.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["obs_id"])
        writer.writeheader()
        for obs in range(1, 101):
            writer.writerow({"obs_id": str(obs)})
    first = select_observation_ids(cohort, n=12, shard_index=0, shard_count=1)
    second = select_observation_ids(cohort, n=12, shard_index=0, shard_count=1)
    assert first == second and len(first) == 12
    shards = [set(select_observation_ids(cohort, n=0, shard_index=i, shard_count=4)) for i in range(4)]
    assert len(set.union(*shards)) == 100
    assert all(not (shards[i] & shards[j]) for i in range(4) for j in range(i+1, 4))


def test_photo_schedule_keeps_all_licensed_photos_for_selected_observations(tmp_path):
    metadata = tmp_path / "metadata.csv"
    fields = ["obs_id", "photo_id", "photo_license_code", "medium_image_url", "large_image_url"]
    rows = [
        {"obs_id":"1","photo_id":"11","photo_license_code":"cc-by","medium_image_url":"https://static.inaturalist.org/photos/11/medium.jpg","large_image_url":"https://static.inaturalist.org/photos/11/large.jpg"},
        {"obs_id":"1","photo_id":"12","photo_license_code":"cc0","medium_image_url":"https://static.inaturalist.org/photos/12/medium.jpg","large_image_url":"https://static.inaturalist.org/photos/12/large.jpg"},
        {"obs_id":"2","photo_id":"21","photo_license_code":"","medium_image_url":"https://static.inaturalist.org/photos/21/medium.jpg","large_image_url":"https://static.inaturalist.org/photos/21/large.jpg"},
        {"obs_id":"3","photo_id":"31","photo_license_code":"cc-by","medium_image_url":"https://static.inaturalist.org/photos/31/medium.jpg","large_image_url":"https://static.inaturalist.org/photos/31/large.jpg"},
    ]
    with metadata.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)
    schedule, report = photo_schedule(metadata, {"1", "2"})
    assert {row["photo_id"] for row in schedule} == {"11", "12"}
    assert report["selected_observations"] == 2
    assert report["scheduled_unique_photos"] == 2
    assert report["unavailable_states"]["license_unavailable"] == 1
    assert all("/original.jpg" in row["original_url"] for row in schedule)


def write_metadata(path, rows):
    fields = ["obs_id", "photo_id", "photo_license_code", "medium_image_url", "large_image_url"]
    fields.extend(sorted(set().union(*(row.keys() for row in rows)) - set(fields)))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def photo_row(obs="1", photo="11", license_code="cc-by", url=None):
    return {"obs_id": obs, "photo_id": photo, "photo_license_code": license_code,
            "medium_image_url": "", "large_image_url": url if url is not None else f"https://static.inaturalist.org/photos/{photo}/large.jpg"}


@pytest.mark.parametrize("conflict_license", ["", "all-rights-reserved", "cc0"])
def test_conflicting_license_blocks_photo_even_on_unselected_observation(tmp_path, conflict_license):
    metadata = tmp_path / "metadata.csv"
    write_metadata(metadata, [photo_row(), photo_row(obs="99", license_code=conflict_license)])
    links = []
    schedule, report = photo_schedule(metadata, {"1"}, links)
    assert schedule == []
    assert report["unavailable_states"]["metadata_conflict"] == 1
    assert links[0]["status"] == "metadata_conflict"
    assert json.loads(links[0]["known_photo_observation_ids_json"]) == ["1", "99"]
    assert len(json.loads(links[0]["known_photo_versions_json"])) == 2


def test_link_ledger_preserves_missing_license_url_photo_and_observation(tmp_path):
    metadata = tmp_path / "metadata.csv"
    write_metadata(metadata, [photo_row("1", "11", ""), photo_row("2", "22", url="invalid"),
                              photo_row("3", ""), photo_row("4", "44")])
    links = []
    schedule, report = photo_schedule(metadata, {"1", "2", "3", "4", "5"}, links)
    assert len(schedule) == 1
    assert len(links) == 5
    assert {row["status"] for row in links} == {"license_unavailable", "url_unavailable", "photo_id_missing", "scheduled", "metadata_link_missing"}
    assert report["retained_metadata_rows"] == 4
    assert report["retained_link_state_rows"] == 5


@pytest.fixture
def offline_run(tmp_path, monkeypatch):
    """Exercise real CSV/JSONL writers without network, YOLO or feature engines."""
    cohort, metadata, weights = (tmp_path / name for name in ("cohort.csv", "metadata.csv", "weights.pt"))
    cohort.write_text("obs_id\n1\n", encoding="utf-8")
    write_metadata(metadata, [{**photo_row(), "photo_attribution": "Offline Photographer (CC BY)"}])
    weights.write_bytes(b"offline-detector-fixture")
    from analysis.v3 import detect_cached_images
    monkeypatch.setattr(detect_cached_images, "MODEL_SHA", hashlib.sha256(weights.read_bytes()).hexdigest())
    session = SimpleNamespace(headers={}, close=lambda: None)
    monkeypatch.setitem(sys.modules, "requests", SimpleNamespace(Session=lambda: session))
    monkeypatch.setitem(sys.modules, "ultralytics", SimpleNamespace(YOLO=lambda _: object()))
    buffer = io.BytesIO()
    Image.new("RGB", (128, 128), (20, 110, 40)).save(buffer, format="PNG")
    payload = buffer.getvalue()
    monkeypatch.setattr(stream, "_download", lambda *_: payload)
    monkeypatch.setattr(stream, "_detector_records", lambda *_: [([10., 10., 40., 40.], .93, 0), ([60., 60., 100., 100.], .82, 0)])
    counter = {"calls": 0}

    def fake_measure(*_):
        counter["calls"] += 1
        endpoints = [{"endpoint_id": row["endpoint_id"], "unit": row["unit"], "value": float(counter["calls"]),
                      "original": float(counter["calls"]), "mirror": float(counter["calls"]),
                      "mirror_abs_difference": 0., "status": "usable"} for row in stream.features.registry()]
        return {"endpoints": endpoints, "raw": {"fixture": 17}, "legacy_qc": {"fixture": "usable"},
                "extended_combined": {}, "engine_errors": {}, "diagnostics": {
                    "head_min_dimension_px": 100, "head_laplacian_variance": 81.0,
                    "unavailable_diagnostic": float("nan"),
                    "paired_colour": {"non_head_context": {"support_status": "available", "lab_chroma": 37.41},
                                      "green_non_head_context": {"support_status": "insufficient_pixels", "lab_chroma": 22.0}}}}, {}
    monkeypatch.setattr(stream.features, "measure", fake_measure)
    def invoke(name="output"):
        out = tmp_path / name
        report = stream.run(metadata, cohort, weights, DECISION, out, pilot_observations=1)
        return out, report
    return SimpleNamespace(invoke=invoke, payload=payload, calls=counter, cohort=cohort, metadata=metadata, weights=weights)


def csv_rows(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def test_actual_stream_outputs_retain_all_raw_routes_context_shifts_and_identity(offline_run):
    out, report = offline_run.invoke()
    endpoint_rows = csv_rows(out / "endpoint_measurements_private.csv")
    assert len(endpoint_rows) == 54
    assert len({row["endpoint_id"] for row in endpoint_rows}) == 27
    assert sum(row["primary_measurement_eligible"] == "True" for row in endpoint_rows) == 28
    assert all(row["value"] for row in endpoint_rows if row["primary_measurement_eligible"] == "False")
    bbox = csv_rows(out / "bbox_measurements_private.csv")
    assert len(bbox) == 2 * 5 * 5
    assert {row["condition"] for row in bbox} == {"baseline", *(row["id"] for row in stream.BBOX_CONDITIONS)}
    assert len({row["value"] for row in bbox if row["head_index"] == "0"}) == 5
    diagnostics = [json.loads(line) for line in (out / "head_diagnostics_private.jsonl").read_text().splitlines()]
    assert diagnostics[0]["diagnostics"]["paired_colour"]["non_head_context"]["lab_chroma"] == 37.41
    assert diagnostics[0]["diagnostics"]["head_laplacian_variance"] == 81.0
    assert diagnostics[0]["diagnostics"]["unavailable_diagnostic"] is None
    assert "NaN" not in (out / "head_diagnostics_private.jsonl").read_text()
    links = csv_rows(out / "photo_observation_links_private.csv")
    assert json.loads(links[0]["source_metadata_json"])["photo_attribution"] == "Offline Photographer (CC BY)"
    provenance = diagnostics[0]["diagnostics"]["stream_measurement_provenance"]
    assert provenance["baseline_recipe"]["head_box"]
    assert provenance["raw"] == {"fixture": 17}
    assert len(provenance["bbox_conditions"]) == 5
    transfer = csv_rows(out / "transfer_private.csv")
    assert transfer[0]["source_byte_sha256"] == hashlib.sha256(offline_run.payload).hexdigest()
    assert len(transfer[0]["oriented_rgb_pixel_sha256"]) == 64
    detector = json.loads((out / "photo_detection_private.jsonl").read_text())
    assert detector["detections"][0]["box_xyxy"] == [10., 10., 40., 40.]
    assert detector["detections"][0]["confidence"] == .93
    assert report["transfer"]["success"] == 1 and report["transfer"]["error"] == 0
    assert len(report["all_endpoint_usable_heads"]) == 27
    assert len(report["endpoint_usable_heads"]) == 14
    assert report["production_execution_authorized"] is False
    for name, sha in report["numerical_file_sha256"].items():
        assert hashlib.sha256((out / name).read_bytes()).hexdigest() == sha
    assert not any(key in json.dumps(report) for key in ('"obs_id"', '"photo_id"', 'box_xyxy'))


def test_bbox_failure_preserves_baseline_but_blocks_primary_membership(offline_run, monkeypatch):
    original = stream.features.measure
    calls = {"n": 0}
    def fail_one_shift(*args):
        calls["n"] += 1
        if calls["n"] == 2:
            raise ValueError("offline-shift-failure")
        return original(*args)
    monkeypatch.setattr(stream.features, "measure", fail_one_shift)
    out, report = offline_run.invoke()
    rows = csv_rows(out / "endpoint_measurements_private.csv")
    orientation = next(row for row in rows if row["head_index"] == "0" and row["endpoint_id"] == "orientation_image_vertical_angle")
    assert orientation["value"] and orientation["status"] == "usable"
    assert orientation["primary_measurement_eligible"] == "False"
    bbox = csv_rows(out / "bbox_measurements_private.csv")
    assert len(bbox) == 50
    assert sum(row["status"] == "bbox_measurement_error" for row in bbox) == 5
    assert report["primary_measurement_eligible_heads"]["orientation_image_vertical_angle"] == 1


def test_head_error_retains_failure_slots_without_counting_photo_success(offline_run, monkeypatch):
    def fail(*_):
        raise ValueError("offline-head-failure")
    monkeypatch.setattr(stream.features, "measure", fail)
    out, report = offline_run.invoke()
    assert len(csv_rows(out / "endpoint_measurements_private.csv")) == 54
    assert {row["status"] for row in csv_rows(out / "endpoint_measurements_private.csv")} == {"head_measurement_error"}
    assert len(csv_rows(out / "bbox_measurements_private.csv")) == 50
    assert report["transfer"]["success"] == 0
    assert report["transfer"]["error"] == 1
    assert report["transfer"]["head_measurement_errors"] == 2
    assert report["transfer"]["bytes"] == len(offline_run.payload)


def test_decode_error_keeps_download_identity_and_bytes(offline_run, monkeypatch):
    def fail(*_):
        raise ValueError("offline-decode-failure")
    monkeypatch.setattr(stream, "_decode", fail)
    out, report = offline_run.invoke()
    rows = csv_rows(out / "transfer_private.csv")
    assert rows[0]["source_byte_sha256"] == hashlib.sha256(offline_run.payload).hexdigest()
    assert int(rows[0]["bytes"]) == len(offline_run.payload)
    assert report["transfer"]["success"] == 0 and report["transfer"]["error"] == 1
    assert len(csv_rows(out / "photo_observation_links_private.csv")) == 1


def test_download_error_keeps_failed_link_and_no_invented_trait(offline_run, monkeypatch):
    def fail(*_):
        raise ValueError("offline-download-failure")
    monkeypatch.setattr(stream, "_download", fail)
    out, report = offline_run.invoke()
    assert csv_rows(out / "endpoint_measurements_private.csv") == []
    assert len(csv_rows(out / "photo_observation_links_private.csv")) == 1
    assert report["transfer"]["success"] == 0 and report["transfer"]["error"] == 1


def test_full_production_is_blocked_before_weights_inputs_or_network(tmp_path, monkeypatch):
    monkeypatch.setattr(stream, "_download", lambda *_: pytest.fail("network must not be reached"))
    with pytest.raises(ValueError, match="Production stream is blocked"):
        stream.run(tmp_path / "absent", tmp_path / "absent", tmp_path / "absent", DECISION, tmp_path / "out", pilot_observations=0)
    assert not (tmp_path / "out").exists()


def test_private_outputs_cannot_be_written_under_tracked_source():
    with pytest.raises(ValueError, match="Private numerical outputs"):
        stream.run(Path("missing"), Path("missing"), Path("missing"), DECISION,
                   ROOT / "analysis" / "private-stream-test", pilot_observations=1)


def test_oriented_pixel_identity_encodes_dimensions_and_rgb():
    flat = Image.new("RGB", (2, 4), (10, 20, 30))
    wide = Image.new("RGB", (4, 2), (10, 20, 30))
    assert flat.tobytes() == wide.tobytes()
    assert stream.oriented_pixel_identity(flat) != stream.oriented_pixel_identity(wide)
    assert stream.oriented_pixel_identity(flat) == stream.oriented_pixel_identity(flat.copy())


def test_runtime_records_actual_metadata_not_requirement_pins(monkeypatch):
    versions = {name: "observed-test-version-9.7" for name in
                ("torch", "torchvision", "ultralytics", "numpy", "pandas", "pillow", "requests", "opencv-python")}
    def installed_version(name):
        if name not in versions:
            raise stream.importlib.metadata.PackageNotFoundError(name)
        return versions[name]
    monkeypatch.setattr(stream.importlib.metadata, "version", installed_version)
    monkeypatch.setattr(stream.platform, "python_version", lambda: "3.12.runtime-test")
    runtime = stream.software_runtime()
    assert runtime["python"] == "3.12.runtime-test"
    assert runtime["distribution_versions"]["torch"] == "observed-test-version-9.7"
    assert runtime["distribution_versions"]["opencv-python-headless"] is None
    assert runtime["missing_distribution_metadata"] == ["opencv-python-headless"]
    assert runtime["requirements_conformance_verified"] is False
    assert runtime["loaded_feature_module_versions"]["cv2"] == stream.features.cv2.__version__


def test_stream_execution_contract_persists_actual_runtime_snapshot(offline_run, monkeypatch):
    runtime = {"python": "actual-interpreter-test", "distribution_versions": {"torch": "actual-torch-test"},
               "requirements_conformance_verified": False}
    monkeypatch.setattr(stream, "software_runtime", lambda: runtime)
    out, report = offline_run.invoke()
    execution = json.loads((out / "execution_contract.json").read_text(encoding="utf-8"))
    assert execution["software_runtime"] == runtime
    assert report["execution_contract_sha256"] == hashlib.sha256((out / "execution_contract.json").read_bytes()).hexdigest()
