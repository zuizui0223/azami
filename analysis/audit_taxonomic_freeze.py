#!/usr/bin/env python3
"""Build or validate the Chapter 1 taxonomic freeze table.

The script collects exact names from one or more frozen analysis tables. In
``--initialize`` mode it writes a complete review template without making any
taxonomic judgement. In validation mode it fails closed when names are missing,
duplicated, unresolved, excluded while still active, or insufficiently sourced.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

REQUIRED_DECISION_COLUMNS = {
    "input_name",
    "accepted_name",
    "decision",
    "authority_source",
    "authority_record",
    "checked_date",
    "notes",
}
ALLOWED_DECISIONS = {"accepted", "synonym", "provisional", "excluded"}
NAME_COLUMNS = (
    "taxon_name",
    "scientific_name",
    "species",
    "species_name",
    "accepted_name",
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return [
            {str(k).strip(): (v or "").strip() for k, v in row.items()}
            for row in csv.DictReader(handle)
        ]


def detect_name_column(rows: list[dict[str, str]], explicit: str | None) -> str:
    if not rows:
        raise ValueError("Input taxon table is empty")
    columns = set(rows[0])
    if explicit:
        if explicit not in columns:
            raise ValueError(f"Requested name column {explicit!r} not found")
        return explicit
    for candidate in NAME_COLUMNS:
        if candidate in columns:
            return candidate
    raise ValueError(
        "Could not identify a taxon-name column; pass --name-column explicitly"
    )


def collect_names(
    paths: Iterable[Path], explicit_column: str | None
) -> tuple[set[str], dict[str, int], dict[str, list[str]]]:
    names: set[str] = set()
    counts: Counter[str] = Counter()
    sources: dict[str, list[str]] = defaultdict(list)
    for path in paths:
        rows = read_csv(path)
        column = detect_name_column(rows, explicit_column)
        seen_here: set[str] = set()
        for row in rows:
            name = row.get(column, "").strip()
            if not name:
                continue
            names.add(name)
            counts[name] += 1
            seen_here.add(name)
        for name in sorted(seen_here):
            sources[name].append(str(path))
    return names, dict(counts), dict(sources)


def write_initial_template(
    output_path: Path,
    observed_names: set[str],
    counts: dict[str, int],
    sources: dict[str, list[str]],
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "input_name",
        "accepted_name",
        "decision",
        "authority_source",
        "authority_record",
        "checked_date",
        "analysis_row_count",
        "source_tables",
        "notes",
    ]
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for name in sorted(observed_names):
            writer.writerow(
                {
                    "input_name": name,
                    "accepted_name": name,
                    "decision": "provisional",
                    "authority_source": "",
                    "authority_record": "",
                    "checked_date": "",
                    "analysis_row_count": counts.get(name, 0),
                    "source_tables": ";".join(sources.get(name, [])),
                    "notes": "",
                }
            )


def validate_decisions(
    rows: list[dict[str, str]],
) -> tuple[dict[str, dict[str, str]], list[str]]:
    errors: list[str] = []
    if not rows:
        return {}, ["Decision table is empty"]

    missing_columns = REQUIRED_DECISION_COLUMNS - set(rows[0])
    if missing_columns:
        return {}, [f"Missing decision columns: {sorted(missing_columns)}"]

    by_input: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row_number, row in enumerate(rows, start=2):
        input_name = row["input_name"].strip()
        accepted_name = row["accepted_name"].strip()
        decision = row["decision"].strip().lower()
        if not input_name:
            errors.append(f"row {row_number}: input_name is blank")
            continue
        by_input[input_name].append(row)
        if decision not in ALLOWED_DECISIONS:
            errors.append(
                f"row {row_number}: decision {decision!r} is not one of "
                f"{sorted(ALLOWED_DECISIONS)}"
            )
        if decision != "excluded" and not accepted_name:
            errors.append(f"row {row_number}: accepted_name is required for {decision}")
        if not row["authority_source"].strip():
            errors.append(f"row {row_number}: authority_source is blank")
        if not row["authority_record"].strip():
            errors.append(f"row {row_number}: authority_record is blank")
        if not row["checked_date"].strip():
            errors.append(f"row {row_number}: checked_date is blank")

    for input_name, matches in by_input.items():
        if len(matches) != 1:
            errors.append(
                f"{input_name!r}: {len(matches)} decision rows; exactly one required"
            )

    unique = {name: matches[0] for name, matches in by_input.items() if len(matches) == 1}
    return unique, errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--taxon-table", action="append", required=True, type=Path)
    parser.add_argument("--decisions", type=Path)
    parser.add_argument("--initialize", action="store_true")
    parser.add_argument("--name-column")
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    observed_names, observation_counts, sources = collect_names(
        args.taxon_table, args.name_column
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)

    if args.initialize:
        template_path = args.output_dir / "taxonomic_decisions_to_review.csv"
        write_initial_template(
            template_path, observed_names, observation_counts, sources
        )
        summary = {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "mode": "initialize",
            "input_tables": [str(path) for path in args.taxon_table],
            "n_observed_names": len(observed_names),
            "template": str(template_path),
            "status": "review_required",
        }
        (args.output_dir / "taxonomic_freeze_summary.json").write_text(
            json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        print(f"Wrote review template for {len(observed_names)} names: {template_path}")
        return 0

    if args.decisions is None:
        parser.error("--decisions is required unless --initialize is used")

    decisions, errors = validate_decisions(read_csv(args.decisions))
    missing = sorted(observed_names - set(decisions))
    unused = sorted(set(decisions) - observed_names)
    unresolved = sorted(
        name
        for name in observed_names & set(decisions)
        if decisions[name]["decision"].strip().lower() == "provisional"
    )
    excluded_present = sorted(
        name
        for name in observed_names & set(decisions)
        if decisions[name]["decision"].strip().lower() == "excluded"
    )

    if missing:
        errors.append(f"Observed names without decisions: {missing}")
    if unresolved:
        errors.append(f"Provisional decisions remain: {unresolved}")
    if excluded_present:
        errors.append(f"Excluded names remain in analysis inputs: {excluded_present}")

    accepted_to_inputs: dict[str, list[str]] = defaultdict(list)
    for name in sorted(observed_names & set(decisions)):
        row = decisions[name]
        if row["decision"].strip().lower() != "excluded":
            accepted_to_inputs[row["accepted_name"].strip()].append(name)

    output_rows = []
    for name in sorted(observed_names):
        row = decisions.get(name, {})
        output_rows.append(
            {
                "input_name": name,
                "accepted_name": row.get("accepted_name", ""),
                "decision": row.get("decision", "missing"),
                "authority_source": row.get("authority_source", ""),
                "authority_record": row.get("authority_record", ""),
                "checked_date": row.get("checked_date", ""),
                "analysis_row_count": observation_counts.get(name, 0),
                "source_tables": ";".join(sources.get(name, [])),
                "notes": row.get("notes", ""),
            }
        )

    table_path = args.output_dir / "taxonomic_freeze_decisions.csv"
    with table_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(output_rows[0]) if output_rows else ["input_name"],
        )
        writer.writeheader()
        writer.writerows(output_rows)

    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "mode": "validate",
        "input_tables": [str(path) for path in args.taxon_table],
        "decision_table": str(args.decisions),
        "n_observed_names": len(observed_names),
        "n_accepted_names": len(accepted_to_inputs),
        "n_synonym_collapses": sum(len(v) > 1 for v in accepted_to_inputs.values()),
        "missing_names": missing,
        "unused_decisions": unused,
        "provisional_names": unresolved,
        "excluded_names_present": excluded_present,
        "accepted_name_to_input_names": dict(sorted(accepted_to_inputs.items())),
        "errors": errors,
        "status": "pass" if not errors else "fail",
    }
    (args.output_dir / "taxonomic_freeze_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1
    print(
        f"Taxonomic freeze passed: {len(observed_names)} input names -> "
        f"{len(accepted_to_inputs)} accepted names"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
