#!/usr/bin/env python3
"""Compile reviewed literature evidence into taxon-by-trait summaries.

Conflicting descriptions are retained as conflicts. This script never chooses a
"true" state by majority vote and never exports literature states as image labels.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

REQUIRED = {"evidence_id", "accepted_taxon", "trait_id", "trait_state", "source_id", "review_status"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compile reviewed Chapter 1 literature evidence.")
    parser.add_argument("--evidence", required=True, help="literature_evidence_valid.csv")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--include-review-statuses", default="accepted,adjudicated")
    return parser.parse_args()


def text(value: object) -> str:
    return "" if pd.isna(value) else str(value).strip()


def join_unique(values: pd.Series) -> str:
    return "|".join(sorted({text(value) for value in values if text(value)}))


def main() -> None:
    args = parse_args()
    statuses = {value.strip().lower() for value in args.include_review_statuses.split(",") if value.strip()}
    if not statuses:
        raise ValueError("--include-review-statuses must contain at least one status")
    evidence = pd.read_csv(args.evidence, dtype=str, keep_default_na=False)
    missing = REQUIRED.difference(evidence.columns)
    if missing:
        raise ValueError(f"Evidence table missing required columns: {sorted(missing)}")
    evidence["review_status"] = evidence["review_status"].map(text).str.lower()
    used = evidence.loc[evidence["review_status"].isin(statuses)].copy()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if used.empty:
        summary = pd.DataFrame(columns=[
            "accepted_taxon", "trait_id", "n_evidence", "n_sources", "states_observed",
            "trait_state_consensus", "conflict_flag", "interpretation", "source_ids", "evidence_ids",
        ])
    else:
        grouped = used.groupby(["accepted_taxon", "trait_id"], dropna=False)
        summary = grouped.agg(
            n_evidence=("evidence_id", "size"),
            n_sources=("source_id", lambda values: len({text(value) for value in values if text(value)})),
            states_observed=("trait_state", join_unique),
            source_ids=("source_id", join_unique),
            evidence_ids=("evidence_id", join_unique),
        ).reset_index()
        summary["conflict_flag"] = summary["states_observed"].str.contains(r"\|", regex=True)
        summary["trait_state_consensus"] = summary.apply(
            lambda row: "" if row["conflict_flag"] else row["states_observed"], axis=1
        )
        summary["interpretation"] = summary["conflict_flag"].map(
            {True: "conflicting_reviewed_evidence", False: "single_reviewed_state"}
        )
        summary = summary[[
            "accepted_taxon", "trait_id", "n_evidence", "n_sources", "states_observed",
            "trait_state_consensus", "conflict_flag", "interpretation", "source_ids", "evidence_ids",
        ]].sort_values(["trait_id", "accepted_taxon"])

    conflicts = summary.loc[summary["conflict_flag"]].copy() if not summary.empty else summary.copy()
    used.to_csv(out_dir / "literature_evidence_used.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(out_dir / "literature_trait_summary.csv", index=False, encoding="utf-8-sig")
    conflicts.to_csv(out_dir / "literature_trait_conflicts.csv", index=False, encoding="utf-8-sig")
    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "included_review_statuses": sorted(statuses),
        "n_valid_evidence_rows_input": int(len(evidence)),
        "n_evidence_rows_used": int(len(used)),
        "n_taxon_trait_summaries": int(len(summary)),
        "n_conflicting_taxon_trait_summaries": int(len(conflicts)),
        "note": "A single reviewed literature state is a source summary, not a label for an individual photograph.",
    }
    (out_dir / "literature_trait_summary_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
