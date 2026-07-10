#!/usr/bin/env python3
"""Create a checksummed, machine-readable manifest for a validated submission bundle."""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

PACKAGE_DISTRIBUTIONS = (
    "numpy",
    "pandas",
    "statsmodels",
    "opencv-python-headless",
    "matplotlib",
    "Pillow",
    "requests",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bundle-root", required=True)
    parser.add_argument("--validation-report", required=True)
    parser.add_argument(
        "--config",
        default=str(Path(__file__).with_name("submission_config.json")),
    )
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--git-sha", default="")
    parser.add_argument("--r-session-info", default="")
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def package_versions() -> dict[str, str]:
    versions: dict[str, str] = {}
    for distribution in PACKAGE_DISTRIBUTIONS:
        try:
            versions[distribution] = importlib.metadata.version(distribution)
        except importlib.metadata.PackageNotFoundError:
            versions[distribution] = "not_installed"
    return versions


def resolve_git_sha(requested: str) -> str:
    requested = requested.strip()
    if requested:
        return requested
    environment = os.environ.get("GITHUB_SHA", "").strip()
    if environment:
        return environment
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def csv_shape(path: Path) -> tuple[int, int, list[str]]:
    columns = list(pd.read_csv(path, nrows=0).columns)
    n_rows = 0
    for chunk in pd.read_csv(path, chunksize=100_000, low_memory=False):
        n_rows += len(chunk)
    return n_rows, len(columns), columns


def file_record(path: Path, root: Path) -> dict[str, Any]:
    record: dict[str, Any] = {
        "relative_path": path.relative_to(root).as_posix(),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_file(path),
        "suffix": path.suffix.lower(),
    }
    if path.suffix.lower() == ".csv":
        rows, columns, names = csv_shape(path)
        record.update({
            "table_rows": rows,
            "table_columns": columns,
            "column_names": names,
        })
    elif path.suffix.lower() == ".json":
        value = json.loads(path.read_text(encoding="utf-8"))
        record["json_top_level_type"] = type(value).__name__
        if isinstance(value, dict):
            record["json_top_level_keys"] = sorted(value)
    return record


def build_manifest(
    bundle_root: Path,
    validation_report: Path,
    config_path: Path,
    out_dir: Path,
    git_sha: str,
    r_session_info: Path | None = None,
) -> dict[str, Any]:
    bundle_root = bundle_root.resolve()
    out_dir = out_dir.resolve()
    if not bundle_root.is_dir():
        raise ValueError(f"Bundle root is not a directory: {bundle_root}")

    validation = json.loads(validation_report.read_text(encoding="utf-8"))
    if validation.get("status") != "pass":
        raise ValueError("Validation report is not a passing report")
    config = json.loads(config_path.read_text(encoding="utf-8"))
    if validation.get("analysis_version") != config.get("analysis_version"):
        raise ValueError("Validation report and config use different analysis versions")

    out_dir.mkdir(parents=True, exist_ok=True)
    excluded_roots = {out_dir}
    files: list[Path] = []
    for path in bundle_root.rglob("*"):
        if not path.is_file():
            continue
        resolved = path.resolve()
        if any(root == resolved or root in resolved.parents for root in excluded_roots):
            continue
        if "__pycache__" in path.parts or path.name in {".DS_Store", "Thumbs.db"}:
            continue
        files.append(path)
    files.sort(key=lambda path: path.relative_to(bundle_root).as_posix())
    if not files:
        raise ValueError("Submission bundle contains no files")

    records = [file_record(path, bundle_root) for path in files]
    table = pd.DataFrame([
        {
            "relative_path": record["relative_path"],
            "size_bytes": record["size_bytes"],
            "sha256": record["sha256"],
            "suffix": record["suffix"],
            "table_rows": record.get("table_rows"),
            "table_columns": record.get("table_columns"),
        }
        for record in records
    ])
    table.to_csv(out_dir / "submission_file_manifest.csv", index=False, encoding="utf-8-sig")

    environment: dict[str, Any] = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "git_sha": git_sha,
        "python_version": sys.version,
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "packages": package_versions(),
    }
    if r_session_info is not None:
        if not r_session_info.is_file():
            raise ValueError(f"R session-info file does not exist: {r_session_info}")
        copied = out_dir / "R_sessionInfo.txt"
        copied.write_bytes(r_session_info.read_bytes())
        environment["r_session_info_file"] = copied.name
        environment["r_session_info_sha256"] = sha256_file(copied)
    (out_dir / "submission_environment.json").write_text(
        json.dumps(environment, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    manifest = {
        "analysis_id": config["analysis_id"],
        "analysis_version": config["analysis_version"],
        "created_at_utc": environment["created_at_utc"],
        "git_sha": git_sha,
        "bundle_root_name": bundle_root.name,
        "validation_report_sha256": sha256_file(validation_report),
        "config_sha256": sha256_file(config_path),
        "n_files": len(records),
        "total_size_bytes": int(sum(record["size_bytes"] for record in records)),
        "files": records,
        "interpretation_boundary": {
            "primary": "continuous colour, outline and image-referenced orientation",
            "auxiliary": "exploratory high-resolution involucre/spine image proxies",
            "historical": "alternative-tree sensitivity, not a resolved species history",
            "causal": "macroecological association only; no adaptation or selection claim"
        },
    }
    (out_dir / "submission_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    return manifest


def main() -> None:
    args = parse_args()
    manifest = build_manifest(
        bundle_root=Path(args.bundle_root),
        validation_report=Path(args.validation_report),
        config_path=Path(args.config),
        out_dir=Path(args.out_dir),
        git_sha=resolve_git_sha(args.git_sha),
        r_session_info=Path(args.r_session_info) if args.r_session_info else None,
    )
    print(json.dumps({
        "analysis_version": manifest["analysis_version"],
        "n_files": manifest["n_files"],
        "total_size_bytes": manifest["total_size_bytes"],
        "out_dir": args.out_dir,
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
