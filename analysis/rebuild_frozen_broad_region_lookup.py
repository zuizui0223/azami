#!/usr/bin/env python3
"""Rebuild the exact frozen broad-region lookup used by the sampling audit.

This wrapper pins both source inputs that were recovered during the 2026-09-03
reproducibility audit: the strict-spatial observation cohort and the Natural
Earth 1:50m admin-0 archive.  It delegates polygon assignment to
``build_naturalearth_broad_region_lookup.py`` and then requires the exact frozen
CSV bytes before reporting success.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
EXPECTED_OBSERVATION_SHA256 = "2172e3570f684770d0f919ecd81265c8460574e287bc4fb057db4f719cab7bb0"
EXPECTED_NATURAL_EARTH_ZIP_SHA256 = "5fed433373581fa648920435f937d95f2d3c0200e067409c6478dcdf1b853139"
EXPECTED_OUTPUT_SHA256 = "085c4e8d45ceb34d32c6c961675ce74a4f0a33580f6cdd8ecd2ff1800a6364ff"
EXPECTED_ROWS = 46276
EXPECTED_REGION_COUNTS = {
    "Europe": 23337,
    "North America": 19357,
    "Oceania": 1339,
    "Asia": 1161,
    "South America": 648,
    "Africa": 432,
    "UNMAPPED": 2,
}
EXPECTED_ASSIGNMENT_COUNTS = {
    "polygon_within": 44332,
    "nearest_country_fallback": 1942,
    "unmapped": 2,
}
NATURAL_EARTH_ARCHIVE_NAME = "ne_50m_admin_0_countries.zip"
NATURAL_EARTH_SHAPEFILE = "ne_50m_admin_0_countries.shp"
NATURAL_EARTH_SOURCE_URL = (
    "https://naturalearth.s3.amazonaws.com/50m_cultural/"
    "ne_50m_admin_0_countries.zip"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--observations", required=True, type=Path)
    parser.add_argument("--naturalearth-zip", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require_hash(path: Path, expected: str, label: str) -> str:
    actual = sha256_file(path)
    if actual != expected:
        raise SystemExit(f"{label} SHA-256 mismatch: {actual}; expected {expected}")
    return actual


def safe_extract(zip_path: Path, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    root = destination.resolve()
    with zipfile.ZipFile(zip_path) as archive:
        for member in archive.infolist():
            target = (destination / member.filename).resolve()
            if root not in target.parents and target != root:
                raise SystemExit(f"Unsafe ZIP path: {member.filename}")
        archive.extractall(destination)


def main() -> int:
    args = parse_args()
    observation_sha = require_hash(
        args.observations,
        EXPECTED_OBSERVATION_SHA256,
        "strict-spatial observation cohort",
    )
    naturalearth_sha = require_hash(
        args.naturalearth_zip,
        EXPECTED_NATURAL_EARTH_ZIP_SHA256,
        "Natural Earth archive",
    )

    with tempfile.TemporaryDirectory(prefix="azami-naturalearth-") as temp_name:
        temp = Path(temp_name)
        safe_extract(args.naturalearth_zip, temp)
        shapefiles = sorted(temp.rglob(NATURAL_EARTH_SHAPEFILE))
        if len(shapefiles) != 1:
            raise SystemExit(
                f"Expected exactly one {NATURAL_EARTH_SHAPEFILE}; found {len(shapefiles)}"
            )
        builder_report = temp / "broad_region_builder_report.json"
        args.output.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "analysis" / "build_naturalearth_broad_region_lookup.py"),
                "--observations",
                str(args.observations),
                "--countries",
                str(shapefiles[0]),
                "--output",
                str(args.output),
                "--report",
                str(builder_report),
            ],
            cwd=ROOT,
            check=True,
        )
        builder_payload = json.loads(builder_report.read_text(encoding="utf-8"))

    output_sha = sha256_file(args.output)
    frame = pd.read_csv(args.output, low_memory=False)
    region_counts = {
        str(key): int(value)
        for key, value in frame["broad_region"].value_counts(dropna=False).items()
    }
    assignment_counts = {
        str(key): int(value)
        for key, value in frame["assignment_method"].value_counts(dropna=False).items()
    }
    checks = {
        "observation_sha256_matches": observation_sha == EXPECTED_OBSERVATION_SHA256,
        "naturalearth_zip_sha256_matches": naturalearth_sha
        == EXPECTED_NATURAL_EARTH_ZIP_SHA256,
        "output_sha256_matches": output_sha == EXPECTED_OUTPUT_SHA256,
        "n_rows_matches": int(len(frame)) == EXPECTED_ROWS,
        "region_counts_match": region_counts == EXPECTED_REGION_COUNTS,
        "assignment_counts_match": assignment_counts == EXPECTED_ASSIGNMENT_COUNTS,
        "n_unmapped_matches": int(frame["broad_region"].eq("UNMAPPED").sum()) == 2,
    }
    payload = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "source_observations": str(args.observations),
        "source_observation_sha256": observation_sha,
        "naturalearth_archive": str(args.naturalearth_zip),
        "naturalearth_source_url": NATURAL_EARTH_SOURCE_URL,
        "naturalearth_zip_sha256": naturalearth_sha,
        "output": str(args.output),
        "output_sha256": output_sha,
        "n_rows": int(len(frame)),
        "region_counts": region_counts,
        "assignment_counts": assignment_counts,
        "builder_report": builder_payload,
        "checks": checks,
        "claim_boundary": (
            "This reconstructs the frozen sampling-composition region lookup only; "
            "it is not a biogeographic or ecological region inference."
        ),
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    if payload["status"] != "PASS":
        raise SystemExit("Frozen broad-region reconstruction failed closed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
