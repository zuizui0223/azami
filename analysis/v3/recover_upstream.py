"""Recover pinned original chunks and historical all-photo processing locally.

No new source-photo retrieval and no raw-data publication. Exact existing files
are reusable; a changed existing file stops recovery instead of being overwritten.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path, PurePosixPath
import shutil
import zipfile

from .recover_revision_inputs import download
from .workflow import ROOT, digest


def member_target(root: Path, name: str) -> Path:
    parts = PurePosixPath(name)
    if not name or parts.is_absolute() or ".." in parts.parts or "\\" in name or ":" in name:
        raise ValueError("Unsafe archive member path")
    target = root.joinpath(*parts.parts).resolve()
    if not target.is_relative_to(root.resolve()):
        raise ValueError("Archive member outside recovery directory")
    return target


def extract_verified(archive: Path, out: Path, spec: dict) -> dict:
    if archive.stat().st_size != spec["archive_bytes"] or digest(archive) != spec["archive_sha256"]:
        raise ValueError("Archive identity mismatch: " + str(spec["artifact_id"]))
    items = []
    with zipfile.ZipFile(archive) as zipped:
        members = [item for item in zipped.infolist() if not item.is_dir()]
        if len({item.filename for item in members}) != len(members):
            raise ValueError("Duplicate archive member name")
        if sum(item.file_size for item in members) > 12_000_000_000:
            raise ValueError("Unexpected uncompressed recovery size")
        for item in members:
            if (item.external_attr >> 16) & 0o170000 == 0o120000:
                raise ValueError("Archive symlinks are not admitted")
            target = member_target(out, item.filename)
            # Preserve original archives; extract metadata/raw API/text records only.
            selected = target.suffix.lower() in {".csv", ".json", ".txt", ".ndjson", ".gz"}
            row = {"member": item.filename, "bytes": item.file_size, "extracted": selected}
            if selected:
                target.parent.mkdir(parents=True, exist_ok=True)
                temporary = target.with_name(target.name + ".extracting")
                if temporary.exists():
                    raise ValueError("Partial extraction exists; inspect before retrying: " + item.filename)
                if target.exists():
                    # Compare the existing bytes with the immutable archive, without
                    # trusting a prior receipt or replacing an existing local file.
                    import hashlib
                    expected = hashlib.sha256()
                    with zipped.open(item) as handle:
                        for block in iter(lambda: handle.read(1024 * 1024), b""):
                            expected.update(block)
                    if target.stat().st_size != item.file_size or digest(target) != expected.hexdigest():
                        raise ValueError("Existing extracted file differs: " + item.filename)
                else:
                    with zipped.open(item) as source, temporary.open("xb") as sink:
                        shutil.copyfileobj(source, sink, 1024 * 1024)
                    if temporary.stat().st_size != item.file_size:
                        raise ValueError("Incomplete extraction: " + item.filename)
                    temporary.rename(target)
                row["sha256"] = digest(target)
            items.append(row)
    if "required_member" in spec:
        matches = [r for r in items if r["member"] == spec["required_member"]]
        if len(matches) != 1 or matches[0].get("sha256") != spec["required_member_sha256"]:
            raise ValueError("Original metadata chunk identity mismatch")
    return {"artifact_id": spec["artifact_id"], "archive_sha256": spec["archive_sha256"],
            "status": "EXACT_ARCHIVE_AND_SELECTED_MEMBERS_VERIFIED", "members": items}


def recover_one(spec: dict, out: Path) -> dict:
    folder = out / str(spec["artifact_id"])
    folder.mkdir(parents=True, exist_ok=True)
    archive = folder / "source.zip"
    if not archive.exists():
        temporary = folder / "source.zip.downloading"
        if temporary.exists():
            raise ValueError("Partial archive exists; inspect before retrying")
        print(json.dumps({"artifact": spec["artifact_id"], "stage": "retrieving_archive"}), flush=True)
        download(spec["artifact_id"], temporary, spec["archive_sha256"])
        temporary.rename(archive)
    print(json.dumps({"artifact": spec["artifact_id"], "stage": "verifying_members"}), flush=True)
    report = extract_verified(archive, folder / "members", spec)
    (folder / "recovery_receipt.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"artifact": spec["artifact_id"], "stage": "verified"}), flush=True)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=ROOT / "analysis/v3/upstream_sources.json")
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=2, choices=[1, 2])
    args = parser.parse_args()
    out = args.out_dir.resolve()
    if out == ROOT or (ROOT in out.parents and not any(out.is_relative_to(ROOT / p) for p in ("local_data", "outputs"))):
        raise ValueError("Raw archives must be external or in ignored local_data/outputs")
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    out.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        reports = list(pool.map(lambda spec: recover_one(spec, out), manifest["archives"]))
    summary = {"status": "PINNED_UPSTREAM_ARCHIVES_RECOVERED", "archives": reports,
               "new_source_photo_downloads": 0, "ecological_models_executed": False}
    (out / "upstream_recovery_report.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
