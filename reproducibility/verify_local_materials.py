#!/usr/bin/env python3
"""Verify a local/synced copy of the frozen Chapter 1 reproducibility archive.

This verifier is intentionally independent of Google Drive APIs. It checks the
bytes that actually reached a local directory (for example a OneDrive-synced
folder) against the SHA-256 values committed to the repository.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DURABLE = ROOT / "reproducibility" / "durable_archive_manifest.json"
AVAILABILITY = ROOT / "reproducibility" / "material_availability.json"

DIRECT_FILENAMES = {
    8066010557: "artifact-8066010557-figure1-metadata.zip",
    8099953404: "artifact-8099953404-figure1-global.zip",
    8225059018: "artifact-8225059018-figure1-continuous.zip",
    8076736948: "artifact-8076736948-detector.zip",
    8227254443: "artifact-8227254443-historical.zip",
    8983877726: "artifact-8983877726-spatial.zip",
    9612943217: "artifact-9612943217-continuous.zip",
    9632715852: "artifact-9632715852-multilevel.zip",
    8521925441: "artifact-8521925441-qc.zip",
    8521926057: "artifact-8521926057-predictions.zip",
}
NUMERICAL_IDS = {9612943217, 9632715852, 8227254443, 8983877726}


def sha256_file(path: Path, block: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(block), b""):
            digest.update(chunk)
    return digest.hexdigest()


def hash_stream(handle, digest: hashlib._Hash, block: int = 8 * 1024 * 1024) -> None:
    for chunk in iter(lambda: handle.read(block), b""):
        digest.update(chunk)


def require_file(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(path)


def verify_direct(archive: Path, scope: str, durable: dict) -> list[dict]:
    rows: list[dict] = []
    for item in durable["direct_files"]:
        artifact_id = int(item["artifact_id"])
        if scope == "numerical" and artifact_id not in NUMERICAL_IDS:
            continue
        filename = DIRECT_FILENAMES[artifact_id]
        path = archive / filename
        require_file(path)
        observed = sha256_file(path)
        expected = item["sha256"]
        if observed != expected:
            raise ValueError(f"SHA mismatch for {filename}: {observed} != {expected}")
        rows.append({"file": filename, "sha256": observed, "status": "ok"})
    return rows


def verify_annotation(archive: Path, durable: dict) -> list[dict]:
    spec = durable["chunked_files"]["8521924881"]
    total = hashlib.sha256()
    rows: list[dict] = []
    for part in spec["parts"]:
        path = archive / part["name"]
        require_file(path)
        if path.stat().st_size != int(part["size_bytes"]):
            raise ValueError(f"Size mismatch for {path.name}")
        observed = sha256_file(path)
        if observed != part["sha256"]:
            raise ValueError(f"SHA mismatch for {path.name}: {observed} != {part['sha256']}")
        with path.open("rb") as handle:
            hash_stream(handle, total)
        rows.append({"file": path.name, "sha256": observed, "status": "ok"})
    if total.hexdigest() != spec["original_sha256"]:
        raise ValueError("Reassembled annotation SHA does not match the frozen source")
    rows.append({"reassembled": "8521924881", "sha256": total.hexdigest(), "status": "ok"})
    return rows


def verify_exhaustive(archive: Path, durable: dict) -> list[dict]:
    spec = durable["chunked_files"]["8269246732"]
    total = hashlib.sha256()
    rows: list[dict] = []
    for part in sorted(spec["parts"], key=lambda row: int(row["part"])):
        index = int(part["part"])
        path = archive / f"exhaustive45-part-{index:02d}.zip"
        require_file(path)
        with zipfile.ZipFile(path) as zf:
            members = [name for name in zf.namelist() if not name.endswith("/")]
            if len(members) != 1:
                raise ValueError(f"Expected one raw chunk in {path.name}; found {members}")
            with zf.open(members[0]) as handle:
                chunk_digest = hashlib.sha256()
                for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
                    chunk_digest.update(chunk)
                    total.update(chunk)
        observed = chunk_digest.hexdigest()
        if observed != part["raw_sha256"]:
            raise ValueError(
                f"Raw chunk SHA mismatch for {path.name}: {observed} != {part['raw_sha256']}"
            )
        rows.append({"file": path.name, "raw_sha256": observed, "status": "ok"})
    if total.hexdigest() != spec["original_sha256"]:
        raise ValueError("Reassembled exhaustive source SHA does not match the frozen source")
    rows.append({"reassembled": "8269246732", "sha256": total.hexdigest(), "status": "ok"})
    return rows


def check_raster_cache(path: Path) -> list[dict]:
    sources = json.loads(
        (ROOT / "analysis" / "ch1" / "chelsa_process_environment_sources.json").read_text(
            encoding="utf-8"
        )
    )["sources"]
    rows: list[dict] = []
    for source in sources:
        filename = source["url"].rsplit("/", 1)[-1]
        raster = path / filename
        require_file(raster)
        rows.append({"file": filename, "size_bytes": raster.stat().st_size, "status": "present"})
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive-dir", required=True, type=Path)
    parser.add_argument("--scope", choices=("numerical", "all"), default="numerical")
    parser.add_argument("--raster-cache-dir", type=Path)
    parser.add_argument("--report", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    archive = args.archive_dir.resolve()
    durable = json.loads(DURABLE.read_text(encoding="utf-8"))
    availability = json.loads(AVAILABILITY.read_text(encoding="utf-8"))
    rows = verify_direct(archive, args.scope, durable)
    if args.scope == "all":
        rows.extend(verify_annotation(archive, durable))
        rows.extend(verify_exhaustive(archive, durable))
    raster_rows: list[dict] = []
    if args.raster_cache_dir is not None:
        raster_rows = check_raster_cache(args.raster_cache_dir.resolve())
    report = {
        "status": "PASS",
        "scope": args.scope,
        "archive_dir": str(archive),
        "files": rows,
        "raster_cache": raster_rows,
        "raw_process_raster_reconstruction_fully_offline": availability["material_audit_conclusion"][
            "raw_process_raster_reconstruction_fully_offline"
        ]
        if args.raster_cache_dir is None
        else True,
        "note": (
            "A numerical archive PASS verifies the durable binary inputs. Full-27 recomputation still "
            "needs the five frozen CHELSA process rasters via network or --raster-cache-dir; the canonical "
            "runner subsequently validates the reconstructed nine-predictor table by frozen SHA-256."
        ),
    }
    text = json.dumps(report, indent=2, ensure_ascii=False) + "\n"
    if args.report is not None:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(text, encoding="utf-8")
    print(text, end="")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, ValueError, zipfile.BadZipFile) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        raise SystemExit(2)
