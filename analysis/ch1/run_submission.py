#!/usr/bin/env python3
"""Submission-facing control entry point for the frozen Chapter 1 analysis."""
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


def command_check(_: argparse.Namespace) -> None:
    pipeline = load_pipeline()
    missing: list[str] = []
    path_keys = (
        "script",
        "selection_script",
        "annotation_app",
        "finalizer",
        "evaluation_script",
        "measurement_script",
        "join_script",
        "residual_export_script",
        "region_script",
        "input_builder",
        "audit_script",
        "summary_exporter",
        "japan38_observation_exporter",
    )
    for stage in pipeline.get("stages", {}).values():
        for key in path_keys:
            value = stage.get(key)
            if value and not (ROOT / value).is_file():
                missing.append(value)
        workflow = stage.get("workflow")
        if workflow and not (ROOT / workflow).is_file():
            missing.append(workflow)

    required = [
        "README.md",
        "manuscript/SUBMISSION_MANUSCRIPT.md",
        "manuscript/final_claims.json",
        "manuscript/COHORT_FLOW_AND_ANALYSIS_LEDGER.md",
        "manuscript/EXTERNAL_COMPLETION_GATES.md",
        "manuscript/FIGURE_TABLE_MAP.md",
        "manuscript/RUNBOOK.md",
    ]
    missing.extend(path for path in required if not (ROOT / path).is_file())
    if missing:
        raise SystemExit("Missing canonical paths:\n- " + "\n- ".join(sorted(set(missing))))

    print(
        json.dumps(
            {
                "status": "ok",
                "analysis_status": pipeline.get("status"),
                "checked_stages": sorted(pipeline.get("stages", {})),
            },
            indent=2,
        )
    )


def command_claims(_: argparse.Namespace) -> None:
    script = ROOT / load_pipeline()["stages"]["claims"]["script"]
    subprocess.run([sys.executable, str(script)], cwd=ROOT, check=True)


def command_summary(_: argparse.Namespace) -> None:
    pipeline = load_pipeline()
    out = {
        "status": pipeline.get("status"),
        "submission_story": pipeline.get("submission_story"),
        "current_submission_gates": pipeline.get("current_submission_gates"),
        "policy": pipeline.get("policy"),
    }
    print(json.dumps(out, ensure_ascii=False, indent=2))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for name, func, help_text in (
        ("check", command_check, "verify the active submission path"),
        ("claims", command_claims, "validate frozen manuscript claims"),
        ("summary", command_summary, "print submission story and remaining gates"),
    ):
        command = sub.add_parser(name, help=help_text)
        command.set_defaults(func=func)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
