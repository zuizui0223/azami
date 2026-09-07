import importlib
from pathlib import Path
import sys

import pytest
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
audit = importlib.import_module("analysis.v3.verify_cached_pipeline")


def test_exact_relationship_check():
    audit.require_equal({"a", "b"}, {"b", "a"}, "identity")
    with pytest.raises(ValueError, match="identity mismatch"):
        audit.require_equal({"a"}, {"b"}, "identity")


def mask(tmp_path, value=255):
    path = tmp_path / "mask.png"
    Image.new("L", (10, 20), value).save(path)
    return {"path": "mask.png", "sha256": audit.digest(path), "width": 10, "height": 20}


def test_valid_binary_mask(tmp_path):
    audit.verify_mask(tmp_path, mask(tmp_path))


def test_nonbinary_mask_fails(tmp_path):
    with pytest.raises(ValueError, match="binary"):
        audit.verify_mask(tmp_path, mask(tmp_path, value=12))


def test_changed_mask_fails(tmp_path):
    record = mask(tmp_path)
    record["sha256"] = "0" * 64
    with pytest.raises(ValueError, match="SHA-256"):
        audit.verify_mask(tmp_path, record)


def test_mask_dimensions_fail(tmp_path):
    record = mask(tmp_path)
    record["width"] = 20
    with pytest.raises(ValueError, match="dimensions"):
        audit.verify_mask(tmp_path, record)


def test_mask_escape_fails(tmp_path):
    record = mask(tmp_path)
    record["path"] = "../mask.png"
    with pytest.raises(ValueError, match="escapes"):
        audit.verify_mask(tmp_path, record)
