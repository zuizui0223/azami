#!/usr/bin/env python3
"""Submission-facing control entry point for the frozen Chapter 1 analysis."""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
PIPELINE = HERE / "pipeline.json"
CONTRACT = HERE / "submission_contract.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def load_pipeline() -> dict:
    return load_json(PIPELINE)


def load_contract() -> dict:
    return load_json(CONTRACT)


def check(_):
    pipeline = load_pipeline()
    contract = load_contract()
    missing: list[str] = []

    for stage in pipeline.get("stages", {}).values():
        for key in (
            "script", "selection_script", "annotation_app", "finalizer",
            "evaluation_script", "measurement_script", "join_script",
            "residual_export_script", "region_script", "input_builder",
            "audit_script", "contract", "summary_exporter", "japan38_observation_exporter",
            "hypothesis_status", "result_ledger",
        ):
            value = stage.get(key)
            if value and not (ROOT / value).is_file():
                missing.append(value)
        workflow = stage.get("workflow")
        if workflow and not (ROOT / workflow).is_file():
            missing.append(workflow)

    for value in contract.get("required_submission_paths", []):
        if not (ROOT / value).exists():
            missing.append(value)

    if missing:
        raise SystemExit("Missing canonical paths:\n- " + "\n- ".join(sorted(set(missing))))

    forbidden_present = [
        value for value in contract.get("forbidden_active_paths", [])
        if (ROOT / value).exists()
    ]
    forbidden_present.extend(
        value for value in contract.get("forbidden_active_workflows", [])
        if (ROOT / value).exists()
    )
    if forbidden_present:
        raise SystemExit(
            "Submission-only tree contains forbidden legacy paths:\n- "
            + "\n- ".join(sorted(set(forbidden_present)))
        )

    active_stage_names = set(pipeline.get("stages", {}))
    if "eazami_trait_handoff" in active_stage_names:
        raise SystemExit("EAzami handoff must not be an active submission stage")
    if any(name.startswith("legacy_lability") for name in active_stage_names):
        raise SystemExit("Withdrawn lability workflow must not be an active submission stage")
    if int(pipeline.get("schema_version", 0)) < 10:
        raise SystemExit("Active submission pipeline schema must be >= 10")

    print(json.dumps({
        "status": "ok",
        "analysis_status": pipeline.get("status"),
        "checked_stages": sorted(active_stage_names),
        "forbidden_paths_present": [],
    }, indent=2))


def claims(_):
    pipeline = load_pipeline()
    scripts = [
        ROOT / pipeline["stages"]["claims"]["script"],
        ROOT / pipeline["stages"]["hypothesis_recovery"]["script"],
    ]
    for script in scripts:
        subprocess.run([sys.executable, str(script)], cwd=ROOT, check=True)


def summary(_):
    pipeline = load_pipeline()
    contract = load_contract()
    hypothesis_status = load_json(ROOT / pipeline["stages"]["hypothesis_recovery"]["hypothesis_status"])
    out = {
        "status": pipeline.get("status"),
        "submission_story": pipeline.get("submission_story"),
        "current_submission_gates": pipeline.get("current_submission_gates"),
        "completion": hypothesis_status.get("completion"),
        "policy": pipeline.get("policy"),
        "withdrawn_from_headline": contract.get("withdrawn_from_headline"),
    }
    print(json.dumps(out, ensure_ascii=False, indent=2))


def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    for name, func in (("check", check), ("claims", claims), ("summary", summary)):
        command = sub.add_parser(name)
        command.set_defaults(func=func)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
