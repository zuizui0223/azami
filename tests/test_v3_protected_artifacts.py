import hashlib
import json
from types import SimpleNamespace
import zipfile

import pytest

from analysis.v3 import private_replay
from analysis.v3.protected_artifacts import DraftStore, pack, unpack, manifest_entries
from analysis.v3 import protected_artifacts as transport
from analysis.v3.workflow import digest


@pytest.fixture
def bundle(tmp_path):
    data = tmp_path / "data.csv"
    data.write_bytes(b"id,value\nprivate,2.75\n")
    selection = tmp_path / "selection.json"
    selection.write_text(json.dumps({"schema_version": 1, "files": [
        {"name": "measurement.csv", "path": str(data), "sha256": digest(data)},
        {"name": "alias.csv", "path": str(data), "sha256": digest(data)}]}))
    private_replay.snapshot(selection, tmp_path / "snapshot")
    archive = tmp_path / "numeric.zip"
    receipt = pack(tmp_path / "snapshot", archive)
    return SimpleNamespace(path=archive, receipt=receipt, data=data, snapshot=tmp_path / "snapshot")


def test_exact_numerical_roundtrip_and_duplicates(bundle, tmp_path):
    receipt = unpack(bundle.path, tmp_path / "unpacked", bundle.receipt)
    assert receipt["files"] == 2
    assert receipt["restored_bytes"] == 2 * bundle.data.stat().st_size
    assert (tmp_path / "unpacked/restored/measurement.csv").read_bytes() == bundle.data.read_bytes()
    assert receipt["production_image_execution_authorized"] is False


def test_changed_bundle_fails_before_extracting(bundle, tmp_path):
    with bundle.path.open("ab") as handle:
        handle.write(b"changed")
    with pytest.raises(ValueError, match="identity"):
        unpack(bundle.path, tmp_path / "bad", bundle.receipt)
    assert not (tmp_path / "bad").exists()


@pytest.mark.parametrize("name", ["../../out.csv", "photos/image.jpg", "unknown.txt", "blobs/" + "f" * 64])
def test_unlisted_members_are_rejected_even_if_zip_is_repinned(bundle, tmp_path, name):
    with zipfile.ZipFile(bundle.path, "a") as archive:
        archive.writestr(name, b"bad")
    receipt = {**bundle.receipt, "bundle_bytes": bundle.path.stat().st_size, "bundle_sha256": digest(bundle.path)}
    with pytest.raises(ValueError, match="Unexpected"):
        unpack(bundle.path, tmp_path / "bad", receipt)


def test_corrupt_blob_fails_even_after_outer_zip_rehash(bundle, tmp_path):
    corrupt = tmp_path / "changed.zip"
    with zipfile.ZipFile(bundle.path) as src, zipfile.ZipFile(corrupt, "x") as dst:
        for info in src.infolist():
            data = src.read(info.filename)
            dst.writestr(info.filename, data if info.filename == "snapshot_manifest.json" else b"z" * len(data))
    receipt = {**bundle.receipt, "bundle_bytes": corrupt.stat().st_size, "bundle_sha256": digest(corrupt)}
    with pytest.raises(ValueError, match="blob changed"):
        unpack(corrupt, tmp_path / "bad", receipt)
    assert not (tmp_path / "bad/restored/private_restore_report.json").exists()


def test_never_overwrites_pack_or_restore(bundle, tmp_path):
    with pytest.raises(ValueError, match="Preserve"):
        pack(bundle.snapshot, bundle.path)
    unpack(bundle.path, tmp_path / "out", bundle.receipt)
    with pytest.raises(ValueError, match="Preserve"):
        unpack(bundle.path, tmp_path / "out", bundle.receipt)


def test_original_images_are_not_allowed_in_numerical_manifest():
    row = {"name": "photo.jpg", "sha256": "0" * 64, "bytes": 4}
    with pytest.raises(ValueError, match="Unsafe"):
        manifest_entries(json.dumps({"schema_version": 1, "kind": "private_numerical_snapshot_v1", "raw_images_included": False, "files": [row]}).encode())


def store(public_code=404, draft=True, published=None, tag="private-v3-numerical-test", assets=None):
    release = {"id": 123, "draft": draft, "published_at": published, "tag_name": tag, "assets": assets or []}
    reply = SimpleNamespace(raise_for_status=lambda: None, json=lambda: release)
    def get(url, **kwargs):
        if url.endswith('/assets'):
            start = (kwargs['params']['page'] - 1) * 100
            return SimpleNamespace(raise_for_status=lambda: None, json=lambda: release['assets'][start:start + 100])
        return reply
    session = SimpleNamespace(get=get)
    anon = SimpleNamespace(get=lambda *a, **k: SimpleNamespace(status_code=public_code))
    return DraftStore({"repository": "zuizui0223/azami", "release_id": 123, "tag_name": "private-v3-numerical-test"}, session, anon)


def test_only_verified_hidden_draft_is_accepted():
    assert store().check()["draft"] is True


@pytest.mark.parametrize("kwargs", [{"public_code": 200}, {"public_code": 403}, {"draft": False}, {"published": "2026-01-01"}, {"tag": "public-v1"}])
def test_published_unknown_or_inaccessible_state_does_not_pass(kwargs):
    with pytest.raises(ValueError):
        store(**kwargs).check()


def test_asset_download_rejects_foreign_or_missing_id_before_transfer(tmp_path):
    with pytest.raises(ValueError, match="metadata"):
        store().download({"asset_id": 7}, tmp_path / "output.zip")
    assert not (tmp_path / "output.zip").exists()


def test_upload_does_not_replace_existing_assets(bundle):
    existing = store(assets=[{"id": 1, "name": "v3-replay-test.zip"}])
    with pytest.raises(ValueError, match="no overwrite"):
        existing.upload(bundle.path, "v3-replay-test.zip")


def test_all_pages_are_listed_before_asset_reuse_or_upload():
    assets = [{'id': i, 'name': f'v3-unit-{i}.zip'} for i in range(205)]
    assert store(assets=assets).check()['assets'] == assets


def test_embedded_release_list_does_not_hide_later_assets():
    target = store(assets=[{'id': 1}])
    original = target.session.get
    target.session.get = lambda url, **kwargs: (
        SimpleNamespace(raise_for_status=lambda: None, json=lambda: [{'id': 1}, {'id': 2}])
        if url.endswith('/assets') else original(url, **kwargs))
    assert [r['id'] for r in target.check()['assets']] == [1, 2]


def test_duplicate_paginated_identity_is_rejected():
    target = store(assets=[{'id': i} for i in range(101)] + [{'id': 0}])
    with pytest.raises(ValueError, match='Duplicate'):
        target.check()


def test_disappeared_embedded_asset_is_not_silently_accepted():
    target = store(assets=[{'id': 1}])
    original = target.session.get
    target.session.get = lambda url, **kwargs: (
        SimpleNamespace(raise_for_status=lambda: None, json=lambda: [])
        if url.endswith('/assets') else original(url, **kwargs))
    with pytest.raises(ValueError, match='incomplete'):
        target.check()


def test_complete_cloud_replay_keeps_private_packet_out_of_public_report(bundle, tmp_path, monkeypatch):
    from analysis.v3 import reconciled_stream_input
    def packet(path, expected, budget):
        assert path.read_bytes() == bundle.data.read_bytes()
        assert digest(path) == expected and budget == 128
        return {"selected": ["private-observation-id"], "report": {"selected_observations": 1, "request_candidates": 0}}
    monkeypatch.setattr(reconciled_stream_input, "pilot_input", packet)
    class OfflineStore:
        def __init__(self, contract):
            self.payload = None
        def download(self, asset, out):
            out.write_bytes(bundle.path.read_bytes() if asset["asset_id"] == 1 else self.payload)
        def upload(self, path, name):
            self.payload = path.read_bytes()
            return {"asset_id": 2, "asset_name": name}
        def close(self):
            pass
    monkeypatch.setattr(transport, "DraftStore", OfflineStore)
    monkeypatch.setenv("GITHUB_RUN_ID", "123")
    monkeypatch.setenv("GITHUB_RUN_ATTEMPT", "1")
    monkeypatch.setenv("GITHUB_SHA", "a" * 40)
    contract = tmp_path / "contract.json"
    contract.write_text(json.dumps({"status": "numerical_transport_only_no_image_or_ecological_authorization",
        "release_id": 9, "source_asset": {**bundle.receipt, "asset_id": 1},
        "schedule_member": "measurement.csv", "schedule_sha256": digest(bundle.data)}))
    result = transport.cloud_replay(contract, tmp_path / "run")
    assert result["source_files_restored"] == 2
    assert result["ecological_models_executed"] == result["image_requests_executed"] == 0
    assert "private-observation-id" not in (tmp_path / "run/public_report.json").read_text()
    assert "private-observation-id" in (tmp_path / "run/returned/restored/worker_packet_private.json").read_text()
