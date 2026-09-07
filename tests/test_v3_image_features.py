import importlib
import json
from pathlib import Path
import sqlite3
import sys

import cv2
import numpy as np
from PIL import Image
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
features = importlib.import_module("analysis.v3.image_features")
worker = importlib.import_module("analysis.v3.measure_cached_heads")


def scene():
    context = np.full((280, 280, 3), (20, 140, 20), dtype=np.uint8)
    cv2.ellipse(context, (140, 150), (40, 65), 0, 0, 360, (60, 90, 50), -1)
    cv2.ellipse(context, (140, 95), (35, 25), 0, 0, 360, (180, 50, 200), -1)
    return context[55:230, 85:195].copy(), context


def test_registry_all_slots_and_composition_visible():
    registry = features.registry()
    assert len(registry) == 27
    assert sum(bool(r["compositional_group"]) for r in registry) == 4
    assert any(r["endpoint_id"] == "visible_floret_fraction" for r in registry)


def test_context_excludes_all_detected_heads_not_only_focal():
    context = np.full((20, 30, 3), (20, 150, 20), dtype=np.uint8)
    available, green = features.context_masks(context, (10, 20, 40, 40), [(10, 20, 20, 30), (30, 30, 50, 50)])
    assert available.sum() == 400
    assert not available[:10, :10].any()
    assert not available[10:, 20:].any()
    assert np.array_equal(available, green)


def test_context_coordinate_mismatch_stops():
    with pytest.raises(ValueError, match="coordinates"):
        features.context_masks(np.zeros((5, 5, 3), np.uint8), (0, 0, 6, 6), [])


def test_empty_colour_is_missing_not_zero():
    result = features.colour_summary(np.zeros((5, 5, 3), np.uint8), np.zeros((5, 5), bool))
    assert result["n_pixels"] == 0
    assert result["lab_chroma"] is None


def test_same_colour_statistics_on_flower_and_background():
    image = np.full((20, 20, 3), (40, 70, 180), np.uint8)
    first, second = np.zeros((20, 20), bool), np.zeros((20, 20), bool)
    first[:10] = True
    second[10:] = True
    assert features.colour_summary(image, first) == features.colour_summary(image, second)


def test_all_endpoint_rows_and_raw_values_retained():
    head, context = scene()
    result, masks = features.measure(head, context, (0, 0, 280, 280), [(85, 55, 195, 230)])
    assert len(result["endpoints"]) == 27
    assert set(result["raw"]) == {"original", "mirror"}
    assert any(e["value"] is not None for e in result["endpoints"])
    assert any("low_surface_resolution" in e["status"] for e in result["endpoints"])
    assert len(masks) == 4
    json.dumps(result, allow_nan=False)


def test_repeated_measurement_identical():
    head, context = scene()
    one = features.measure(head, context, (0, 0, 280, 280), [(85, 55, 195, 230)])[0]
    two = features.measure(head, context, (0, 0, 280, 280), [(85, 55, 195, 230)])[0]
    assert one == two


def test_missing_and_bounded_values_not_clipped():
    row = {"lower_bound": "0", "upper_bound": "1"}
    assert "outside_registry_bounds" in features.bound_reasons([1.1, None], row)
    assert "measurement_missing" in features.bound_reasons([1.1, None], row)
    assert features.mean_pair(1.1, None) == 1.1
    assert features.mean_pair(None, None) is None


def fixture(tmp_path):
    root = tmp_path / "detector"
    root.mkdir()
    head, context = scene()
    crops = {}
    for name, image, box in (("head", head, [85, 55, 195, 230]), ("context", context, [0, 0, 280, 280])):
        path = root / f"{name}.png"
        sha = worker.save_png(Image.fromarray(image[:, :, ::-1]), path)
        crops[name] = {"path": path.name, "sha256": sha, "box": box}
    with sqlite3.connect(root / "detection.sqlite") as db:
        db.execute("CREATE TABLE detections(head_id TEXT,sha256 TEXT,det_index INTEGER,roi_status TEXT,crops_json TEXT)")
        db.executemany("INSERT INTO detections VALUES (?,?,?,?,?)", [("a:0", "a", 0, "roi_ready", json.dumps(crops)), ("a:1", "a", 1, "invalid_roi", "{}")])
    (root / "detection_report.json").write_text(json.dumps({"job_states": {"detected": 1}, "execution_contract": {"synthetic": True}, "detection_database_sha256": worker.digest(root / "detection.sqlite")}))
    return root


def test_worker_resumes_retains_invalid_and_is_not_accuracy(tmp_path):
    root = fixture(tmp_path)
    one = worker.run(root, tmp_path / "out", limit=1)
    assert one["status"] == "CACHED_27_BASELINE_PARTIAL"
    result = worker.run(root, tmp_path / "out")
    assert result["counts"]["endpoint_rows_retained"] == 54
    assert result["job_states"] == {"invalid_roi": 1, "measured": 1}
    assert result["independent_accuracy_estimated"] is False
    assert result["technical_perturbation_evaluation_completed"] is False
    repeat = worker.run(root, tmp_path / "out")
    assert repeat["counts"]["completed_this_invocation"] == 0
    assert repeat["measurement_database_sha256"] == result["measurement_database_sha256"]


def test_changed_crop_is_explicit_worker_failure(tmp_path):
    root = fixture(tmp_path)
    path = root / "head.png"
    path.write_bytes(b"changed synthetic input")
    result = worker.run(root, tmp_path / "out")
    assert result["job_states"] == {"error": 1, "invalid_roi": 1}
    assert result["counts"]["endpoint_rows_retained"] == 54
    assert result["status"] == "CACHED_27_BASELINE_COMPLETED_WITH_ERRORS"


def test_partial_detector_cannot_become_measurement_source(tmp_path):
    root = fixture(tmp_path)
    path = root / "detection_report.json"
    report = json.loads(path.read_text())
    report["job_states"]["pending"] = 1
    path.write_text(json.dumps(report))
    with pytest.raises(ValueError, match="Finish"):
        worker.run(root, tmp_path / "out")
