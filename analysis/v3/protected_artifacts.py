"""Protected numerical transport through an UNPUBLISHED GitHub draft release.

GitHub Actions remains the compute platform. Public workflow artifacts contain
only receipts. This module never creates/publishes/deletes a release, changes
visibility, retrieves source images, or authorizes ecological inference.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import stat
import zipfile

import requests

from .private_replay import checked_entries, restore, snapshot
from .workflow import digest, text_digest

REPOSITORY = "zuizui0223/azami"
MANIFEST = "snapshot_manifest.json"
MAX_ASSET_BYTES = 2_000_000_000


def require(condition, message):
    if not condition:
        raise ValueError(message)


def new_json(path: Path, value: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as handle:
        json.dump(value, handle, sort_keys=True, indent=2, allow_nan=False)
        handle.write("\n")


def manifest_entries(raw: bytes):
    require(len(raw) <= 10_000_000, "Numerical manifest too large")
    manifest = json.loads(raw)
    require(manifest.get("schema_version") == 1 and manifest.get("kind") == "private_numerical_snapshot_v1", "Unknown numerical manifest")
    require(manifest.get("raw_images_included") is False, "Source image archive is not permitted")
    entries = checked_entries(manifest["files"])
    sizes = {}
    for row in entries:
        require(isinstance(row.get("bytes"), int) and row["bytes"] >= 0, "Invalid numerical size")
        require(row["sha256"] not in sizes or sizes[row["sha256"]] == row["bytes"], "Inconsistent duplicate blob size")
        sizes[row["sha256"]] = row["bytes"]
    return entries, sizes


def pack(snapshot_dir: Path, out: Path) -> dict:
    require(not out.exists(), "Preserve previous numerical bundle")
    require(not (snapshot_dir / "incomplete_run.json").exists(), "Cannot pack incomplete snapshot")
    raw = (snapshot_dir / MANIFEST).read_bytes()
    entries, sizes = manifest_entries(raw)
    out.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(out, "x", compression=zipfile.ZIP_DEFLATED, compresslevel=1) as zipped:
        zipped.writestr(MANIFEST, raw)
        for sha, size in sorted(sizes.items()):
            path = snapshot_dir / "blobs" / sha
            require(path.is_file() and not path.is_symlink() and path.stat().st_size == size and digest(path) == sha, "Snapshot blob changed before packing")
            zipped.write(path, "blobs/" + sha)
    require(out.stat().st_size < MAX_ASSET_BYTES, "Bundle exceeds protected asset size limit; do not upload")
    return {"bundle_sha256": digest(out), "bundle_bytes": out.stat().st_size,
            "manifest_sha256": hashlib.sha256(raw).hexdigest(), "files": len(entries),
            "restored_bytes": sum(row["bytes"] for row in entries), "source_images_included": False}


def unpack(bundle: Path, out: Path, expected: dict) -> dict:
    require(not out.exists(), "Preserve previous numerical restoration")
    require(bundle.stat().st_size == expected["bundle_bytes"] and digest(bundle) == expected["bundle_sha256"], "Protected bundle identity differs")
    with zipfile.ZipFile(bundle) as zipped:
        infos = zipped.infolist()
        names = [info.filename for info in infos]
        require(names.count(MANIFEST) == 1 and len(names) == len(set(names)), "Duplicate or missing archive manifest/member")
        require(zipped.getinfo(MANIFEST).file_size <= 10_000_000, "Oversized numerical manifest")
        raw = zipped.read(MANIFEST)
        require(hashlib.sha256(raw).hexdigest() == expected["manifest_sha256"], "Protected manifest identity differs")
        entries, sizes = manifest_entries(raw)
        require(set(names) == {MANIFEST, *("blobs/" + sha for sha in sizes)}, "Unexpected numerical archive members")
        require(len(entries) == expected["files"] and sum(row["bytes"] for row in entries) == expected["restored_bytes"], "Numerical manifest counts differ")
        for info in infos:
            mode = info.external_attr >> 16
            require(not info.is_dir() and not stat.S_ISLNK(mode), "Archive directories/links not permitted")
            require(info.filename == MANIFEST or info.file_size == sizes[info.filename[6:]], "Numerical blob size differs")
        out.mkdir(parents=True)
        blobs = out / "snapshot"
        (blobs / "blobs").mkdir(parents=True)
        with (blobs / MANIFEST).open("xb") as handle:
            handle.write(raw)
        for sha in sorted(sizes):
            target = blobs / "blobs" / sha
            with zipped.open("blobs/" + sha) as src, target.open("xb") as dst:
                shutil.copyfileobj(src, dst, length=1024 * 1024)
            require(digest(target) == sha, "Downloaded numerical blob changed")
    receipt = restore(blobs, out / "restored", expected_manifest_sha256=expected["manifest_sha256"])
    require(receipt["files"] == expected["files"] and receipt["restored_bytes"] == expected["restored_bytes"], "Restoration count mismatch")
    return receipt


class DraftStore:
    """A fixed existing draft; refuses published or anonymously visible assets."""

    def __init__(self, contract: dict, session=None, anonymous=None):
        require(contract["repository"] == REPOSITORY, "Protected numerical repository differs")
        require(isinstance(contract["release_id"], int) and contract["release_id"] > 0, "Invalid draft release ID")
        require(contract["tag_name"].startswith("private-v3-numerical-"), "Not a private numerical draft")
        self.contract = contract
        self.base = f"https://api.github.com/repos/{REPOSITORY}"
        self.session = session or requests.Session()
        self.anonymous = anonymous or requests.Session()
        self.anonymous.trust_env = False
        if session is None:
            token = os.environ.get("GH_TOKEN")
            require(bool(token), "GitHub token required for protected artifact access")
            self.session.headers.update({"Authorization": "Bearer " + token, "X-GitHub-Api-Version": "2022-11-28"})

    def close(self):
        self.session.close()
        self.anonymous.close()

    def check(self) -> dict:
        url = self.base + "/releases/" + str(self.contract["release_id"])
        response = self.session.get(url, timeout=(15, 60))
        response.raise_for_status()
        release = response.json()
        require(release.get("draft") is True and release.get("published_at") is None, "Protected release is published or not draft; STOP")
        require(release.get("tag_name") == self.contract["tag_name"] and release.get("id") == self.contract["release_id"], "Draft identity differs")
        public = self.anonymous.get(url, timeout=(15, 60), allow_redirects=False)
        require(public.status_code == 404, "Draft is not verified hidden from anonymous access")
        # The embedded release list is not a complete inventory contract. Use
        # the documented paginated endpoint as production adds many chunks.
        assets, seen = [], set()
        for page in range(1, 1001):
            reply = self.session.get(url + "/assets", params={"per_page": 100, "page": page}, timeout=(15, 60))
            reply.raise_for_status()
            rows = reply.json()
            require(isinstance(rows, list) and len(rows) <= 100, "Invalid draft asset page")
            for row in rows:
                require(isinstance(row, dict) and type(row.get("id")) is int and row["id"] not in seen,
                        "Duplicate or invalid draft asset identity across pages; STOP")
                seen.add(row["id"])
                assets.append(row)
            if len(rows) < 100:
                break
        else:
            raise ValueError("Draft inventory exceeds bounded pagination; STOP")
        require(all(row["id"] in seen for row in release.get("assets", [])), "Draft inventory changed or is incomplete; STOP")
        release["assets"] = assets
        return release

    def download(self, asset: dict, out: Path):
        require(not out.exists(), "Preserve previous downloaded numerical asset")
        release = self.check()
        matches = [row for row in release["assets"] if row["id"] == asset["asset_id"]]
        require(len(matches) == 1 and matches[0]["name"] == asset["asset_name"] and matches[0]["state"] == "uploaded" and matches[0]["size"] == asset["bundle_bytes"], "Protected asset metadata differs")
        endpoint = self.base + "/releases/assets/" + str(asset["asset_id"])
        public = self.anonymous.get(endpoint, timeout=(15, 60), allow_redirects=False)
        require(public.status_code == 404, "Asset is not verified hidden from anonymous access")
        out.parent.mkdir(parents=True, exist_ok=True)
        with self.session.get(endpoint, headers={"Accept": "application/octet-stream"}, stream=True, timeout=(15, 180)) as response:
            response.raise_for_status()
            length = 0
            with out.open("xb") as handle:
                for block in response.iter_content(1024 * 1024):
                    length += len(block)
                    require(length <= asset["bundle_bytes"], "Asset exceeds pinned byte length")
                    handle.write(block)
        require(length == asset["bundle_bytes"] and digest(out) == asset["bundle_sha256"], "Downloaded asset identity differs")
        self.check()

    def upload(self, bundle: Path, name: str) -> dict:
        require(re.fullmatch(r"v3-[a-z0-9-]+\.zip", name) is not None, "Unsafe numerical asset name")
        release = self.check()
        require(not any(row["name"] == name for row in release["assets"]), "Preserve previous draft asset; no overwrite")
        size = bundle.stat().st_size
        require(size < MAX_ASSET_BYTES, "Numerical asset exceeds size limit")
        # Construct the host/path ourselves; never trust an upload_url in metadata.
        url = f"https://uploads.github.com/repos/{REPOSITORY}/releases/{self.contract['release_id']}/assets"
        with bundle.open("rb") as handle:
            response = self.session.post(url, params={"name": name}, headers={"Content-Type": "application/zip", "Content-Length": str(size)}, data=handle, timeout=(15, 180))
        response.raise_for_status()
        uploaded = response.json()
        release = self.check()
        require(any(row["id"] == uploaded["id"] and row["name"] == name and row["size"] == size for row in release["assets"]), "Draft upload not present after transfer")
        return {"asset_id": uploaded["id"], "asset_name": name}


def cloud_replay(contract_path: Path, out: Path) -> dict:
    require(not out.exists(), "Preserve previous replay output")
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    require(contract["status"] == "numerical_transport_only_no_image_or_ecological_authorization", "Unknown replay scope")
    out.mkdir(parents=True)
    store = DraftStore(contract)
    try:
        bundle = out / "input.zip"
        store.download(contract["source_asset"], bundle)
        restored = unpack(bundle, out / "input", contract["source_asset"])
        from .reconciled_stream_input import pilot_input
        root = out / "input" / "restored"
        packet = pilot_input(root / contract["schedule_member"], contract["schedule_sha256"], 128)
        private_output = out / "worker_packet_private.json"
        new_json(private_output, packet)
        selection = out / "output_selection.json"
        new_json(selection, {"schema_version": 1, "files": [
            {"name": "worker_packet_private.json", "path": str(private_output.resolve()), "sha256": digest(private_output)}]})
        snapshot(selection, out / "output_snapshot")
        output_bundle = out / "result.zip"
        output_asset = pack(out / "output_snapshot", output_bundle)
        run_id, attempt = os.environ.get("GITHUB_RUN_ID", "local"), os.environ.get("GITHUB_RUN_ATTEMPT", "1")
        require(run_id.isdigit() and attempt.isdigit(), "Cloud replay requires recorded GitHub run identity")
        output_asset.update(store.upload(output_bundle, f"v3-replay-{run_id}-{attempt}.zip"))
        # Verify a fresh download from GitHub, not merely the pre-upload file.
        store.download(output_asset, out / "returned.zip")
        returned = unpack(out / "returned.zip", out / "returned", output_asset)
        require(digest(out / "returned/restored/worker_packet_private.json") == digest(private_output), "Roundtrip worker packet differs")
        report = {
            "status": "GITHUB_DRAFT_NUMERICAL_RESTORE_AND_OUTPUT_ROUNDTRIP_VERIFIED",
            "repository": REPOSITORY, "release_id": contract["release_id"], "draft_verified": True,
            "anonymous_release_and_asset_requests": "404", "run_id": run_id, "run_attempt": attempt,
            "commit": os.environ.get("GITHUB_SHA"), "contract_sha256": digest(contract_path),
            "implementation_sha256_text_lf": text_digest(Path(__file__)),
            "source_files_restored": restored["files"], "source_bytes_restored": restored["restored_bytes"],
            "selected_observations": packet["report"]["selected_observations"],
            "request_candidates_not_executed": packet["report"]["request_candidates"],
            "returned_files_restored": returned["files"], "output_asset": output_asset,
            "source_images_persisted": 0, "image_requests_executed": 0, "ecological_models_executed": 0,
            "production_image_execution_authorized": False, "ecological_fitting_authorized": False,
            "limits": ["Draft access control is not encryption or permanent publication archiving; never publish this draft.",
                       "Transport/byte identity does not establish scientific validity or measurement coverage."]}
        new_json(out / "public_report.json", report)
        return report
    finally:
        store.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    packed = sub.add_parser("pack")
    packed.add_argument("--snapshot-dir", type=Path, required=True)
    packed.add_argument("--out", type=Path, required=True)
    replay = sub.add_parser("cloud-replay")
    replay.add_argument("--contract", type=Path, required=True)
    replay.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "pack":
        report = pack(args.snapshot_dir, args.out)
        new_json(args.out.with_suffix(".receipt.json"), report)
    else:
        report = cloud_replay(args.contract, args.out)
    print(json.dumps(report))


if __name__ == "__main__":
    main()
