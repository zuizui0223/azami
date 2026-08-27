#!/usr/bin/env python3
"""Validate the frozen v2 submission package without requiring protected inputs."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from zipfile import ZipFile


ROOT = Path(__file__).resolve().parents[1]
SUBMISSION = ROOT / "submission"
RESULTS = ROOT / "analysis_outputs" / "v2_full27_environment_atlas_2026-08-27"
FIGURES = ROOT / "manuscript" / "figures" / "v2_submission"
TEXT_SUFFIXES = {".csv", ".json", ".md", ".py", ".txt"}


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def package_bytes(path: Path, name: str) -> bytes:
    if Path(name).suffix.lower() in TEXT_SUFFIXES:
        text = path.read_text(encoding="utf-8").replace("\r\n", "\n").replace("\r", "\n")
        return text.encode("utf-8")
    return path.read_bytes()


def source_path(bundle_name: str) -> Path:
    if bundle_name == "README.md":
        return SUBMISSION / "README.md"
    if bundle_name.startswith("figures/"):
        return FIGURES / bundle_name.removeprefix("figures/")
    if bundle_name.startswith("data/"):
        return RESULTS / bundle_name.removeprefix("data/")
    if bundle_name.startswith("source/"):
        return ROOT / bundle_name.removeprefix("source/")
    if bundle_name.startswith("code/"):
        return ROOT / "analysis" / bundle_name.removeprefix("code/")
    return SUBMISSION / bundle_name


def check(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> None:
    manifest_path = SUBMISSION / "SUBMISSION_MANIFEST.json"
    sums_path = SUBMISSION / "SHA256SUMS.txt"
    zip_path = SUBMISSION / "Azami_Chapter1_v2_submission_candidate.zip"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    rows = manifest["files"]
    check(manifest["status"] == "scientific_hold_pending_independent_validation", "HOLD status changed")
    check(manifest["source_lane"] == "full27_full_environment_only", "source lane changed")
    check(len(rows) == len({row["path"] for row in rows}), "duplicate manifest paths")

    expected_hashes: dict[str, str] = {}
    for row in rows:
        name = str(row["path"])
        path = source_path(name)
        check(path.is_file(), f"missing manifest source: {name}")
        data = package_bytes(path, name)
        check(len(data) == int(row["bytes"]), f"byte count drift: {name}")
        digest = sha256_bytes(data)
        check(digest == row["sha256"], f"SHA-256 drift: {name}")
        expected_hashes[name] = digest

    with ZipFile(zip_path) as archive:
        names = set(archive.namelist())
        expected_names = set(expected_hashes) | {"SUBMISSION_MANIFEST.json", "SHA256SUMS.txt"}
        check(names == expected_names, "ZIP inventory differs from manifest")
        for name, digest in expected_hashes.items():
            check(sha256_bytes(archive.read(name)) == digest, f"ZIP payload drift: {name}")
        check(
            archive.read("SUBMISSION_MANIFEST.json")
            == package_bytes(manifest_path, "SUBMISSION_MANIFEST.json"),
            "ZIP manifest differs from external manifest",
        )

    sums = {
        line.split("  ", 1)[1]: line.split("  ", 1)[0]
        for line in sums_path.read_text(encoding="utf-8").splitlines()
        if "  " in line
    }
    check(all(sums.get(name) == digest for name, digest in expected_hashes.items()), "SHA256SUMS payload mismatch")
    check(sums.get(zip_path.name) == sha256_bytes(zip_path.read_bytes()), "ZIP SHA-256 mismatch")

    validation_reports = [
        RESULTS / "v2_full27_environment_validation.json",
        RESULTS / "sampling" / "v2_full27_sampling_composition_validation.json",
        RESULTS / "v2_full27_sensitivities_validation.json",
    ]
    checks = 0
    for path in validation_reports:
        report = json.loads(path.read_text(encoding="utf-8"))
        check(report["status"] == "PASS", f"frozen validation is not PASS: {path.name}")
        check(report["n_checks"] == report["n_passed"], f"failed frozen check: {path.name}")
        checks += int(report["n_checks"])

    forbidden = (
        "continuous_trait_reanalysis_v1",
        "reviewer_bias_control_v1",
        "pr69_",
        "legacy_lability",
        "trait_lability",
    )
    check(
        not any(any(token in name.lower() for token in forbidden) for name in expected_hashes),
        "legacy path leaked into bundle",
    )
    print(json.dumps({
        "status": "PASS",
        "manifest_files": len(rows),
        "frozen_validation_checks": checks,
        "zip_sha256": sha256_bytes(zip_path.read_bytes()),
        "scientific_status": manifest["status"],
    }, indent=2))


if __name__ == "__main__":
    main()
