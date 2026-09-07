import hashlib
import importlib
import json
from pathlib import Path
import sqlite3
import sys

import pytest
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
detector = importlib.import_module("analysis.v3.detect_cached_images")


def database():
    db = sqlite3.connect(":memory:")
    db.execute("CREATE TABLE detections (head_id TEXT PRIMARY KEY, sha256 TEXT, det_index INTEGER, box_json TEXT, confidence REAL, class_id INTEGER, roi_status TEXT, reason TEXT, crops_json TEXT)")
    return db


def test_padding_rounding_and_clipping():
    assert detector.crop_box([10.2, 20.2, 30.4, 50.4], 100, 100, .12) == ((7, 16, 33, 55), False)
    assert detector.crop_box([0, 0, 10, 10], 10, 10, .8) == ((0, 0, 10, 10), True)


@pytest.mark.parametrize("box", [[1, 1, 0, 0], [1, 1, float("nan"), 3], [10, 10, 20, 20], [1, 2, 3]])
def test_invalid_boxes_remain_invalid(box):
    with pytest.raises(ValueError):
        detector.crop_box(box, 5, 5, .12)


def test_invalid_padding():
    with pytest.raises(ValueError):
        detector.crop_box([0, 0, 1, 1], 5, 5, -1)


def test_lossless_crops_and_all_proposals_retained(tmp_path):
    db = database()
    image = Image.new("RGB", (100, 100), "purple")
    sha = "a" * 64
    n = detector.write_detections(db, tmp_path, sha, image,
                                  [([20, 20, 40, 40], .8, 0), ([2, 2, 2, 5], .7, 0), ([20, 20, 40, 40], .6, 1)])
    assert n == 2
    rows = db.execute("SELECT * FROM detections ORDER BY det_index").fetchall()
    assert len(rows) == 3
    assert [r[6] for r in rows] == ["roi_ready", "invalid_roi", "invalid_roi"]
    for crop in json.loads(rows[0][8]).values():
        path = detector.readable(tmp_path / crop["path"])
        assert detector.digest(path) == crop["sha256"]
        with Image.open(path) as opened:
            assert opened.getpixel((0, 0)) == image.getpixel((0, 0))


def test_changed_crop_is_integrity_failure_not_bad_roi(tmp_path):
    db = database()
    sha = "a" * 64
    relative = tmp_path / "crops" / sha / "0000_head.png"
    detector.save_png(Image.new("RGB", (2, 2), "red"), relative)
    with pytest.raises(ValueError, match="do not overwrite"):
        detector.write_detections(db, tmp_path, sha, Image.new("RGB", (100, 100)), [([20, 20, 40, 40], .8, 0)])
    assert db.execute("SELECT COUNT(*) FROM detections").fetchone()[0] == 0


def test_repeat_png_bytes_identical(tmp_path):
    path = tmp_path / "crop.png"
    image = Image.new("RGB", (20, 30), "red")
    assert detector.save_png(image, path) == detector.save_png(image, path)


def image_row(tmp_path):
    image = Image.new("RGB", (20, 30), "red")
    path = tmp_path / "source.png"
    image.save(path)
    pixels = hashlib.sha256((20).to_bytes(8, "big") + (30).to_bytes(8, "big") + image.tobytes()).hexdigest()
    return (detector.digest(path), "source.png", 20, 30, pixels)


def test_decode_matches_byte_and_pixel_identity(tmp_path):
    row = image_row(tmp_path)
    assert detector.decode_object(tmp_path, row).size == (20, 30)


def test_changed_pixels_stop(tmp_path):
    row = image_row(tmp_path)
    with pytest.raises(ValueError, match="Decoded image differs"):
        detector.decode_object(tmp_path, (*row[:4], "0" * 64))


def test_changed_bytes_stop(tmp_path):
    row = image_row(tmp_path)
    with pytest.raises(ValueError, match="byte identity"):
        detector.decode_object(tmp_path, ("0" * 64, *row[1:]))


def test_escaped_input_path_stops(tmp_path):
    with pytest.raises(ValueError, match="escapes"):
        detector.decode_object(tmp_path, ("0" * 64, "../source.png", 20, 30, "0" * 64))
