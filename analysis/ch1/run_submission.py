#!/usr/bin/env python3
"""Canonical command-line entry point for the frozen Chapter 1 submission analysis.

This wrapper does not reimplement scientific calculations. It validates paths and
runs the exact numbered scripts that generated the frozen results, preserving
provenance while giving authors and reviewers one stable entry point.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PIPELINE = Path(__file__).with_name("pipeline.json")


def load_pipeline() -> dict:
    return json.loads(PIPELINE.read_text(encoding="utf-8"))


def require_file(path: str | Path, label: str) -> Path:
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_file():
        raise SystemExit(f"Missing {label}: {resolved}")
    return resolved


def require_dir(path: str | Path, label: str, create: bool = False) -> Path:
    resolved = Path(path).expanduser().resolve()
    if create:
        resolved.mkdir(parents=True, exist_ok=True)
    if not resolved.is_dir():
        raise SystemExit(f"Missing {label}: {resolved}")
    return resolved


def run(command: list[str]) -> None:
    print("+", " ".join(command), flush=True)
    subprocess.run(command, cwd=ROOT, check=True)


def command_check(_: argparse.Namespace) -> None:
    pipeline = load_pipeline()
    missing = []
    for stage in pipeline["stages"].values():
        script = stage.get("script")
        if script and not (ROOT / script).is_file():
            missing.append(script)
    required_docs = [
        "manuscript/final_claims.json",
        "manuscript/FINAL_MANUSCRIPT_STRATEGY.md",
        "manuscript/RUNBOOK.md",
        "ch1_global/v2/submission_config.json",
    ]
    missing.extend(path for path in required_docs if not (ROOT / path).is_file())
    if missing:
        raise SystemExit("Missing canonical paths:\n- " + "\n- ".join(sorted(missing)))
    print(json.dumps({"status": "ok", "checked_stages": sorted(pipeline["stages"])}, indent=2))


def command_lability(args: argparse.Namespace) -> None:
    observations = require_file(args.observations, "CHELSA observation table")
    coefficients = require_file(args.coefficients, "pooled coefficient table")
    out_dir = require_dir(args.out_dir, "output directory", create=True)
    script = ROOT / load_pipeline()["stages"]["lability"]["script"]
    run([
        sys.executable,
        str(script),
        "--observations", str(observations),
        "--coefficients", str(coefficients),
        "--out-dir", str(out_dir),
        "--min-species-n", str(args.min_species_n),
        "--min-traits", str(args.min_traits),
        "--min-predictors", str(args.min_predictors),
    ])


def command_claims(_: argparse.Namespace) -> None:
    script = ROOT / load_pipeline()["stages"]["claims"]["script"]
    run([sys.executable, str(script)])


def command_manifest(args: argparse.Namespace) -> None:
    bundle_root = require_dir(args.bundle_root, "submission bundle")
    validation_report = require_file(args.validation_report, "validation report")
    out_dir = require_dir(args.out_dir, "manifest output directory", create=True)
    stage = load_pipeline()["stages"]["submission_manifest"]
    command = [
        sys.executable,
        str(ROOT / stage["script"]),
        "--bundle-root", str(bundle_root),
        "--validation-report", str(validation_report),
        "--config", str(ROOT / stage["config"]),
        "--out-dir", str(out_dir),
    ]
    if args.git_sha:
        command.extend(["--git-sha", args.git_sha])
    if args.r_session_info:
        command.extend(["--r-session-info", str(require_file(args.r_session_info, "R sessionInfo"))])
    run(command)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    check = sub.add_parser("check", help="verify canonical scripts and manuscript controls")
    check.set_defaults(func=command_check)

    lability = sub.add_parser("lability", help="run the frozen two-axis lability stage")
    lability.add_argument("--observations", required=True)
    lability.add_argument("--coefficients", required=True)
    lability.add_argument("--out-dir", required=True)
    lability.add_argument("--min-species-n", type=int, default=10)
    lability.add_argument("--min-traits", type=int, default=6)
    lability.add_argument("--min-predictors", type=int, default=3)
    lability.set_defaults(func=command_lability)

    claims = sub.add_parser("claims", help="validate frozen manuscript claims")
    claims.set_defaults(func=command_claims)

    manifest = sub.add_parser("manifest", help="build a checksummed submission manifest")
    manifest.add_argument("--bundle-root", required=True)
    manifest.add_argument("--validation-report", required=True)
    manifest.add_argument("--out-dir", required=True)
    manifest.add_argument("--git-sha", default="")
    manifest.add_argument("--r-session-info", default="")
    manifest.set_defaults(func=command_manifest)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
