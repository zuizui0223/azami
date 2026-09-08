"""Private, content-addressed snapshots and exact restoration of numerical files.

No cloud transport, deletion, encryption or production authorization is provided.
A caller supplies every file and exact hash; no directory is recursively swept.
The final manifest is written only after all new private blobs are verified.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re

from .enriched_source_cohort import pinned_hash
from .recover_native_source_authority import private_directory, write_new_json
from .workflow import digest

# NPZ checkpoints are copied as opaque, pinned numerical bytes, never executed.
EXTENSIONS = {".csv", ".tsv", ".sqlite", ".json", ".jsonl", ".ndjson", ".txt", ".npz"}


def checked_entries(entries: list[dict]) -> list[dict]:
    if not entries:
        raise ValueError("A private snapshot requires explicit files")
    seen = set()
    for entry in entries:
        name = entry["name"]
        if not isinstance(name, str):
            raise ValueError("Unsafe private numerical member name")
        path = PurePosixPath(name)
        if (not isinstance(name, str) or not name or "\\" in name or ":" in name
                or path.is_absolute() or str(path) != name or ".." in path.parts
                or name.casefold() in {"private_restore_report.json", "incomplete_run.json"}
                or any(p.endswith((".", " ")) for p in path.parts)
                or path.suffix.lower() not in EXTENSIONS
                or any(re.fullmatch(r"(?i)(CON|PRN|AUX|NUL|COM[1-9]|LPT[1-9])(?:\..*)?", p) for p in path.parts)):
            raise ValueError("Unsafe private numerical member name or unsupported file type")
        lower = name.casefold()
        if lower in seen or any(lower.startswith(p + "/") or p.startswith(lower + "/") for p in seen):
            raise ValueError("Duplicate or colliding snapshot member names")
        if not re.fullmatch(r"[0-9a-f]{64}", entry.get("sha256", "")):
            raise ValueError("Every numerical input needs an exact SHA-256")
        seen.add(lower)
    return sorted(entries, key=lambda e: e["name"])


def copy_verified(source: Path, destination: Path, expected: str) -> int:
    if source.is_symlink() or not source.is_file():
        raise ValueError("Snapshot source must be a regular file, not a symbolic link")
    h = hashlib.sha256()
    size = 0
    with source.open("rb") as src, destination.open("xb") as dst:
        for block in iter(lambda: src.read(1024 * 1024), b""):
            dst.write(block)
            h.update(block)
            size += len(block)
        dst.flush()
        os.fsync(dst.fileno())
    if h.hexdigest() != expected or digest(destination) != expected:
        raise ValueError("Private copy SHA-256 mismatch; incomplete output must not be used")
    return size


def snapshot(selection: Path, out: Path) -> dict:
    out = private_directory(out)
    if out.exists():
        raise ValueError("Snapshot destination exists; never replace a prior snapshot")
    specification = json.loads(selection.read_text(encoding="utf-8"))
    if specification.get("schema_version") != 1:
        raise ValueError("Unknown snapshot selection schema")
    entries = checked_entries(specification["files"])
    for entry in entries:
        source = Path(entry["path"])
        if not source.is_absolute() or source.is_symlink() or not source.is_file():
            raise ValueError("Explicit absolute regular input files are required")
    out.mkdir(parents=True)
    (out / "blobs").mkdir()
    try:
        sizes = {}
        files = []
        for index, entry in enumerate(entries, 1):
            sha = entry["sha256"]
            source = Path(entry["path"])
            if sha not in sizes:
                sizes[sha] = copy_verified(source, out / "blobs" / sha, sha)
            else:
                pinned_hash(source, sha, "deduplicated snapshot input")
            files.append({"name": entry["name"], "sha256": sha, "bytes": sizes[sha]})
            if index % 100 == 0:
                print(json.dumps({"stage": "private_snapshot", "files_verified": index}), flush=True)
        manifest = {"schema_version": 1, "kind": "private_numerical_snapshot_v1", "files": files,
                    "raw_images_included": False, "encrypted_by_this_tool": False}
        write_new_json(out / "snapshot_manifest.json", manifest)
        report = {"status": "PRIVATE_LOCAL_SNAPSHOT_WRITTEN_AND_HASH_VERIFIED",
                  "selection_manifest_sha256": digest(selection),
                  "snapshot_manifest_sha256": digest(out / "snapshot_manifest.json"),
                  "files": len(files), "unique_blobs": len(sizes), "unique_blob_bytes": sum(sizes.values()),
                  "off_device_private_restore_verified": False,
                  "production_image_execution_authorized": False, "ecological_fitting_authorized": False}
        write_new_json(out / "snapshot_report.json", report)
        return report
    except Exception as error:
        write_new_json(out / "incomplete_run.json", {"status": "INCOMPLETE_DO_NOT_USE", "error_type": type(error).__name__})
        raise


def restore(snapshot_dir: Path, out: Path, *, expected_manifest_sha256: str) -> dict:
    out = private_directory(out)
    if out.exists():
        raise ValueError("Restore destination exists; never overwrite files")
    if (snapshot_dir / "incomplete_run.json").exists():
        raise ValueError("Cannot restore an incomplete snapshot")
    manifest_path = snapshot_dir / "snapshot_manifest.json"
    pin = pinned_hash(manifest_path, expected_manifest_sha256, "private snapshot manifest")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != 1 or manifest.get("kind") != "private_numerical_snapshot_v1":
        raise ValueError("Unknown private snapshot schema")
    entries = checked_entries(manifest["files"])
    out.mkdir(parents=True)
    try:
        total = 0
        for index, entry in enumerate(entries, 1):
            destination = out / entry["name"]
            destination.parent.mkdir(parents=True, exist_ok=True)
            size = copy_verified(snapshot_dir / "blobs" / entry["sha256"], destination, entry["sha256"])
            if size != entry["bytes"]:
                raise ValueError("Restored byte count differs from manifest")
            total += size
            if index % 100 == 0:
                print(json.dumps({"stage": "private_restore", "files_verified": index}), flush=True)
        report = {"status": "PRIVATE_LOCAL_RESTORE_ALL_MEMBERS_BYTE_VERIFIED",
                  "snapshot_manifest_sha256": pin, "files": len(entries), "restored_bytes": total,
                  "off_device_private_restore_verified": False, "source_files_deleted": 0,
                  "production_image_execution_authorized": False, "ecological_fitting_authorized": False,
                  "limits": ["Restoration proves exact file bytes, not an off-device backup or independent scientific validation.",
                             "These private files are not encrypted by this tool and must not be uploaded as public artifacts."]}
        write_new_json(out / "private_restore_report.json", report)
        return report
    except Exception as error:
        write_new_json(out / "incomplete_run.json", {"status": "INCOMPLETE_DO_NOT_USE", "error_type": type(error).__name__})
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    pack = sub.add_parser("snapshot")
    pack.add_argument("--selection", type=Path, required=True)
    pack.add_argument("--out", type=Path, required=True)
    unpack = sub.add_parser("restore")
    unpack.add_argument("--snapshot-dir", type=Path, required=True)
    unpack.add_argument("--out", type=Path, required=True)
    unpack.add_argument("--expected-manifest-sha256", required=True)
    args = vars(parser.parse_args())
    command = args.pop("command")
    print(json.dumps(snapshot(**args) if command == "snapshot" else restore(**args), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
