#!/usr/bin/env python3
"""Compile focused blinded trait-audit responses into canonical evidence rows.

Human responses are validated against the public task packet and ontology. This
script does not load model candidates or the private selection key; that merge
is deliberately deferred to a separate evaluation step.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd


TASK_REQUIRED = {"task_id", "annotation_unit_id", "trait_id"}
RESPONSE_REQUIRED = {"task_id", "human_state", "task_complete"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compile validated human trait-audit response files.")
    parser.add_argument("--tasks", required=True, help="Public blinded_trait_audit_tasks.csv")
    parser.add_argument("--ontology", required=True)
    parser.add_argument(
        "--response",
        action="append",
        required=True,
        help="Repeatable annotator_id=path/to/trait_audit_responses.csv",
    )
    parser.add_argument("--out-dir", required=True)
    return parser.parse_args()


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def parse_response_spec(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise ValueError("Each --response must be annotator_id=path")
    annotator_id, path_text = value.split("=", 1)
    annotator_id, path_text = annotator_id.strip(), path_text.strip()
    if not annotator_id:
        raise ValueError("annotator_id may not be empty")
    path = Path(path_text)
    if not path.is_file():
        raise FileNotFoundError(path)
    return annotator_id, path


def ontology_states(path: Path) -> dict[str, dict[str, Any]]:
    table = pd.read_csv(path, dtype=str, keep_default_na=False)
    required = {"trait_id", "allowed_states", "allow_multiple"}
    if missing := required.difference(table.columns):
        raise ValueError(f"Ontology missing columns: {sorted(missing)}")
    result: dict[str, dict[str, Any]] = {}
    for row in table.to_dict("records"):
        trait_id = text(row["trait_id"])
        states = [state for state in text(row["allowed_states"]).split("|") if state]
        if not trait_id or not states:
            raise ValueError("Ontology has a trait with empty ID or state list")
        result[trait_id] = {
            "states": states,
            "allow_multiple": text(row["allow_multiple"]).lower() in {"true", "1", "yes"},
        }
    return result


def normalize_human_state(raw: Any, spec: dict[str, Any], trait_id: str) -> str:
    value = text(raw)
    if not value:
        raise ValueError(f"{trait_id}: completed task has an empty human_state")
    states = spec["states"]
    allowed = set(states)
    parts = [part.strip() for part in value.split("|") if part.strip()]
    if not parts:
        raise ValueError(f"{trait_id}: completed task has no valid state")
    unknown = set(parts).difference(allowed)
    if unknown:
        raise ValueError(f"{trait_id}: unknown human state(s): {sorted(unknown)}")
    if not spec["allow_multiple"]:
        if len(parts) != 1:
            raise ValueError(f"{trait_id}: exactly one state required")
        return parts[0]
    if "unassessable" in parts and len(parts) != 1:
        raise ValueError(f"{trait_id}: unassessable must be exclusive")
    if len(set(parts)) != len(parts):
        raise ValueError(f"{trait_id}: duplicate states are not allowed")
    return "|".join(state for state in states if state in set(parts))


def main() -> None:
    args = parse_args()
    tasks = pd.read_csv(args.tasks, dtype=str, keep_default_na=False)
    if missing := TASK_REQUIRED.difference(tasks.columns):
        raise ValueError(f"Public task packet missing columns: {sorted(missing)}")
    if tasks.empty or tasks["task_id"].duplicated().any():
        raise ValueError("Public task packet must have nonempty unique task_id values")
    ontology = ontology_states(Path(args.ontology))
    task_lookup = tasks.set_index("task_id")
    unknown_traits = sorted(set(tasks["trait_id"]).difference(ontology))
    if unknown_traits:
        raise ValueError(f"Public task packet references unknown ontology traits: {unknown_traits}")

    canonical_rows: list[dict[str, str]] = []
    completion_rows: list[dict[str, object]] = []
    for raw_spec in args.response:
        annotator_id, path = parse_response_spec(raw_spec)
        response = pd.read_csv(path, dtype=str, keep_default_na=False)
        if missing := RESPONSE_REQUIRED.difference(response.columns):
            raise ValueError(f"Response {path} missing columns: {sorted(missing)}")
        if response["task_id"].duplicated().any():
            raise ValueError(f"Response {path} has duplicate task_id values")
        unknown_task_ids = sorted(set(response["task_id"]).difference(task_lookup.index))
        if unknown_task_ids:
            raise ValueError(f"Response {path} includes task IDs not in public packet: {unknown_task_ids[:10]}")
        completed = response.loc[response["task_complete"].map(text).str.lower().eq("yes")].copy()
        for row in completed.to_dict("records"):
            task_id = text(row["task_id"])
            task = task_lookup.loc[task_id]
            trait_id = text(task["trait_id"])
            canonical_rows.append({
                "task_id": task_id,
                "annotation_unit_id": text(task["annotation_unit_id"]),
                "trait_id": trait_id,
                "human_state": normalize_human_state(row.get("human_state", ""), ontology[trait_id], trait_id),
                "notes": text(row.get("notes", "")),
                "annotator_id": annotator_id,
                "response_source_file": str(path.resolve()),
            })
        completion_rows.append({
            "annotator_id": annotator_id,
            "response_source_file": str(path.resolve()),
            "n_task_rows_in_file": int(len(response)),
            "n_completed": int(len(completed)),
            "n_incomplete": int(len(response) - len(completed)),
            "n_expected_tasks": int(len(tasks)),
        })

    canonical = pd.DataFrame(canonical_rows)
    if canonical.empty:
        raise ValueError("No completed audit tasks were found in supplied response files")
    if canonical.duplicated(["task_id", "annotator_id"]).any():
        raise ValueError("The same annotator has multiple completed responses for one task")
    output = Path(args.out_dir)
    output.mkdir(parents=True, exist_ok=True)
    canonical = canonical.sort_values(["trait_id", "task_id", "annotator_id"]).reset_index(drop=True)
    canonical.to_csv(output / "trait_audit_annotations_canonical.csv", index=False, encoding="utf-8-sig")
    completion = pd.DataFrame(completion_rows).sort_values("annotator_id")
    completion.to_csv(output / "trait_audit_completion_by_annotator.csv", index=False, encoding="utf-8-sig")
    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "n_public_tasks": int(len(tasks)),
        "n_completed_annotation_records": int(len(canonical)),
        "n_unique_completed_tasks": int(canonical["task_id"].nunique()),
        "n_annotators": int(canonical["annotator_id"].nunique()),
        "completed_records_by_trait": canonical["trait_id"].value_counts().sort_index().to_dict(),
        "rule": "Canonical records are human evidence only. No model candidates, uncertainty margins, taxonomy, or locality are merged at this stage.",
    }
    (output / "trait_audit_compile_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
