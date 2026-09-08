import json
from pathlib import Path

import pytest

from analysis.v3 import private_replay as replay
from analysis.v3.workflow import digest


def selection(tmp_path):
    source = tmp_path / "source.csv"
    source.write_bytes(b"id,value\n1,2.5\n")
    plan = tmp_path / "selection.json"
    plan.write_text(json.dumps({"schema_version": 1, "files": [
        {"name": "measurements/source.csv", "path": str(source.resolve()), "sha256": digest(source)},
        {"name": "identical/source.csv", "path": str(source.resolve()), "sha256": digest(source)}]}))
    return source, plan


def test_exact_private_snapshot_restore_and_content_deduplication(tmp_path):
    source, plan = selection(tmp_path)
    snapshot = tmp_path / "snapshot"
    report = replay.snapshot(plan, snapshot)
    assert report["files"] == 2 and report["unique_blobs"] == 1
    assert report["off_device_private_restore_verified"] is False
    assert report["production_image_execution_authorized"] is False
    restored = tmp_path / "restored"
    receipt = replay.restore(snapshot, restored, expected_manifest_sha256=report["snapshot_manifest_sha256"])
    assert receipt["files"] == 2 and receipt["source_files_deleted"] == 0
    assert (restored / "measurements/source.csv").read_bytes() == source.read_bytes()
    assert (restored / "identical/source.csv").read_bytes() == source.read_bytes()
    assert receipt["off_device_private_restore_verified"] is False


@pytest.mark.parametrize("name", ["../outside.csv", "/outside.csv", "C:/outside.csv", "a\\b.csv",
                                  "a//b.csv", "a/../b.csv", "a./b.csv", "NUL.csv", "a/photo.jpg",
                                  "private_restore_report.json", "incomplete_run.json", 123])
def test_unsafe_or_image_members_are_rejected(name):
    with pytest.raises(ValueError, match="Unsafe"):
        replay.checked_entries([{"name": name, "sha256": "a" * 64}])


def test_case_collisions_and_nested_file_collisions_rejected():
    for names in [("A.csv", "a.csv"), ("a.csv", "a.csv/b.csv")]:
        with pytest.raises(ValueError, match="colliding"):
            replay.checked_entries([{"name": n, "sha256": "a" * 64} for n in names])


def test_input_change_fails_without_publishing_a_manifest(tmp_path):
    source, plan = selection(tmp_path)
    source.write_bytes(b"changed")
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        replay.snapshot(plan, tmp_path / "bad")
    assert not (tmp_path / "bad/snapshot_manifest.json").exists()
    assert (tmp_path / "bad/incomplete_run.json").exists()


def test_changed_blob_or_manifest_cannot_restore_successfully(tmp_path):
    source, plan = selection(tmp_path)
    snapshot = tmp_path / "snapshot"
    report = replay.snapshot(plan, snapshot)
    with pytest.raises(ValueError, match="input identity"):
        replay.restore(snapshot, tmp_path / "bad-pin", expected_manifest_sha256="0" * 64)
    assert not (tmp_path / "bad-pin").exists()
    blob = next((snapshot / "blobs").iterdir())
    blob.write_bytes(b"corrupt")
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        replay.restore(snapshot, tmp_path / "bad-blob", expected_manifest_sha256=report["snapshot_manifest_sha256"])
    assert not (tmp_path / "bad-blob/private_restore_report.json").exists()
    assert (tmp_path / "bad-blob/incomplete_run.json").exists()


def test_existing_snapshot_and_restore_are_preserved(tmp_path):
    source, plan = selection(tmp_path)
    snapshot = tmp_path / "snapshot"
    report = replay.snapshot(plan, snapshot)
    with pytest.raises(ValueError, match="exists"):
        replay.snapshot(plan, snapshot)
    restored = tmp_path / "restored"
    replay.restore(snapshot, restored, expected_manifest_sha256=report["snapshot_manifest_sha256"])
    with pytest.raises(ValueError, match="exists"):
        replay.restore(snapshot, restored, expected_manifest_sha256=report["snapshot_manifest_sha256"])
    assert source.read_bytes() == (restored / "measurements/source.csv").read_bytes()


def test_numpy_checkpoint_bytes_can_be_preserved_without_loading_them(tmp_path):
    import numpy as np
    source = tmp_path / "checkpoint.npz"
    np.savez_compressed(source, indices=np.arange(3), values=np.array([1.0, np.nan, 3.0]))
    plan = tmp_path / "selection.json"
    plan.write_text(json.dumps({"schema_version": 1, "files": [
        {"name": "checkpoints/x.npz", "path": str(source.resolve()), "sha256": digest(source)}]}))
    report = replay.snapshot(plan, tmp_path / "packed")
    replay.restore(tmp_path / "packed", tmp_path / "restored", expected_manifest_sha256=report["snapshot_manifest_sha256"])
    assert (tmp_path / "restored/checkpoints/x.npz").read_bytes() == source.read_bytes()
