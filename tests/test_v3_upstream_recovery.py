"""Synthetic archive recovery checks; no network calls."""
import importlib
from pathlib import Path
import sys
import zipfile

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
recovery = importlib.import_module("analysis.v3.recover_upstream")


def archive_fixture(tmp_path, members=None):
    source = tmp_path / "input.zip"
    with zipfile.ZipFile(source, "w") as zipped:
        for name, value in (members or {"photo_metadata.csv": b"obs_id,photo_id\n1,2\n"}).items():
            zipped.writestr(name, value)
    return source, {"artifact_id": 1, "archive_bytes": source.stat().st_size,
                    "archive_sha256": recovery.digest(source)}


@pytest.mark.parametrize("name", ["../x.csv", "/x.csv", "a/../../x.csv", "C:/x.csv", "a\\x.csv"])
def test_archive_traversal_rejected(tmp_path, name):
    with pytest.raises(ValueError):
        recovery.member_target(tmp_path, name)


def test_exact_member_recovery_and_reuse(tmp_path):
    archive, spec = archive_fixture(tmp_path)
    out = tmp_path / "out"
    first = recovery.extract_verified(archive, out, spec)
    again = recovery.extract_verified(archive, out, spec)
    assert first == again
    assert first["members"][0]["sha256"] == recovery.digest(out / "photo_metadata.csv")


def test_changed_existing_member_is_not_overwritten(tmp_path):
    archive, spec = archive_fixture(tmp_path)
    out = tmp_path / "out"
    recovery.extract_verified(archive, out, spec)
    target = out / "photo_metadata.csv"
    target.write_bytes(b"different")
    with pytest.raises(ValueError, match="Existing extracted file differs"):
        recovery.extract_verified(archive, out, spec)
    assert target.read_bytes() == b"different"


def test_wrong_archive_hash_stops_before_extraction(tmp_path):
    archive, spec = archive_fixture(tmp_path)
    spec["archive_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="Archive identity"):
        recovery.extract_verified(archive, tmp_path / "out", spec)
    assert not (tmp_path / "out").exists()


def test_member_identity_is_independently_checked(tmp_path):
    archive, spec = archive_fixture(tmp_path)
    spec.update(required_member="photo_metadata.csv", required_member_sha256="0" * 64)
    with pytest.raises(ValueError, match="chunk identity"):
        recovery.extract_verified(archive, tmp_path / "out", spec)


def test_image_bytes_stay_in_original_archive(tmp_path):
    archive, spec = archive_fixture(tmp_path, {"photo.jpg": b"image", "run.json": b"{}"})
    report = recovery.extract_verified(archive, tmp_path / "out", spec)
    assert not (tmp_path / "out/photo.jpg").exists()
    assert archive.exists()
    assert report["members"][0]["extracted"] is False


def test_partial_extraction_is_not_silently_retried(tmp_path):
    archive, spec = archive_fixture(tmp_path)
    out = tmp_path / "out"
    out.mkdir()
    (out / "photo_metadata.csv.extracting").write_bytes(b"partial")
    with pytest.raises(ValueError, match="Partial extraction"):
        recovery.extract_verified(archive, out, spec)


def test_unsafe_member_does_not_escape(tmp_path):
    archive, spec = archive_fixture(tmp_path, {"../outside.csv": b"bad"})
    with pytest.raises(ValueError):
        recovery.extract_verified(archive, tmp_path / "out", spec)
    assert not (tmp_path / "outside.csv").exists()
