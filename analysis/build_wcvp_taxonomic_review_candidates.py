#!/usr/bin/env python3
"""Build authority-backed taxonomic review candidates from WCVP via GBIF.

This script is deliberately not a taxonomic decision engine. It reduces the manual
freeze workload by matching every frozen source name against the WCVP checklist
served by GBIF, while leaving every final decision provisional until reviewed.
"""

from __future__ import annotations

import argparse
import csv
import json
import time
import urllib.parse
import urllib.request
from datetime import date
from pathlib import Path
from typing import Any

import pandas as pd

WCVP_DATASET_KEY = "f382f0ce-323a-4091-bb9f-add557f3a9a2"
GBIF_SEARCH = "https://api.gbif.org/v1/species/search"
NAME_COLUMNS = ("taxon_name", "scientific_name", "species", "species_name", "accepted_name")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--taxon-table", action="append", required=True, type=Path)
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--name-column", default=None)
    p.add_argument("--sleep-seconds", type=float, default=0.08)
    p.add_argument("--timeout-seconds", type=float, default=30.0)
    p.add_argument("--max-retries", type=int, default=4)
    return p.parse_args()


def detect_name_column(frame: pd.DataFrame, explicit: str | None) -> str:
    if explicit:
        if explicit not in frame.columns:
            raise ValueError(f"Requested name column {explicit!r} is absent")
        return explicit
    for column in NAME_COLUMNS:
        if column in frame.columns:
            return column
    raise ValueError(f"No recognized taxon-name column among {list(frame.columns)}")


def normalized(value: Any) -> str:
    return " ".join(str(value).strip().split())


def canonical(result: dict[str, Any]) -> str:
    value = result.get("canonicalName") or result.get("scientificName") or ""
    return normalized(value)


def query_wcvp(name: str, timeout: float, retries: int) -> dict[str, Any]:
    params = urllib.parse.urlencode(
        {
            "dataset_key": WCVP_DATASET_KEY,
            "q": name,
            "limit": 100,
        }
    )
    url = f"{GBIF_SEARCH}?{params}"
    request = urllib.request.Request(
        url,
        headers={
            "Accept": "application/json",
            "User-Agent": "azami-ch1-taxonomic-freeze/1.0",
        },
    )
    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return json.loads(response.read().decode("utf-8"))
        except Exception as error:  # network errors are retried, never converted to absence
            last_error = error
            time.sleep(min(2 ** attempt, 8))
    raise RuntimeError(f"WCVP query failed for {name!r}: {last_error}")


def candidate_accepted_name(result: dict[str, Any]) -> str:
    for key in (
        "accepted",
        "acceptedScientificName",
        "acceptedName",
    ):
        if result.get(key):
            return normalized(result[key])
    status = normalized(result.get("taxonomicStatus", "")).upper()
    synonym = bool(result.get("synonym", False))
    if status == "ACCEPTED" and not synonym:
        return normalized(result.get("scientificName") or result.get("canonicalName") or "")
    return ""


def summarize_match(name: str, payload: dict[str, Any]) -> dict[str, Any]:
    results = list(payload.get("results") or [])
    target = name.casefold()
    exact = [r for r in results if canonical(r).casefold() == target]

    if len(exact) == 1:
        chosen = exact[0]
        quality = "exact_unique"
    elif len(exact) > 1:
        # Prefer one accepted record when that resolves duplicate search hits, but retain
        # the multiplicity flag so a reviewer sees that alternatives existed.
        accepted = [
            r for r in exact
            if normalized(r.get("taxonomicStatus", "")).upper() == "ACCEPTED"
            and not bool(r.get("synonym", False))
        ]
        chosen = accepted[0] if len(accepted) == 1 else exact[0]
        quality = "exact_multiple"
    elif results:
        chosen = results[0]
        quality = "search_only_no_exact"
    else:
        chosen = {}
        quality = "unmatched"

    status = normalized(chosen.get("taxonomicStatus", ""))
    synonym = bool(chosen.get("synonym", False)) if chosen else False
    accepted_name = candidate_accepted_name(chosen)
    scientific_name = normalized(chosen.get("scientificName", ""))
    key = chosen.get("key", "")

    if quality == "exact_unique" and status.upper() == "ACCEPTED" and not synonym:
        recommendation = "accepted"
        review_priority = "routine"
    elif quality.startswith("exact") and (synonym or status.upper() == "SYNONYM"):
        recommendation = "synonym"
        review_priority = "high"
    elif quality == "unmatched":
        recommendation = "unresolved"
        review_priority = "high"
    else:
        recommendation = "manual_review"
        review_priority = "high"

    return {
        "input_name": name,
        "accepted_name": accepted_name,
        "decision": "provisional",
        "authority_source": "WCVP via GBIF checklist",
        "authority_record": f"https://www.gbif.org/species/{key}" if key else "",
        "checked_date": date.today().isoformat(),
        "notes": "",
        "wcvp_match_quality": quality,
        "wcvp_taxonomic_status": status,
        "wcvp_synonym": synonym,
        "wcvp_scientific_name": scientific_name,
        "wcvp_accepted_name_candidate": accepted_name,
        "wcvp_key": key,
        "recommended_decision": recommendation,
        "review_priority": review_priority,
        "n_exact_candidates": len(exact),
        "n_search_results": len(results),
    }


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    names: set[str] = set()
    source_rows: list[dict[str, Any]] = []
    for path in args.taxon_table:
        frame = pd.read_csv(path, low_memory=False)
        column = detect_name_column(frame, args.name_column)
        values = sorted({normalized(v) for v in frame[column].dropna() if normalized(v)})
        names.update(values)
        source_rows.append(
            {
                "path": str(path),
                "name_column": column,
                "n_rows": int(len(frame)),
                "n_unique_names": int(len(values)),
            }
        )

    rows: list[dict[str, Any]] = []
    raw_dir = args.output_dir / "wcvp_raw"
    raw_dir.mkdir(exist_ok=True)
    for index, name in enumerate(sorted(names), start=1):
        payload = query_wcvp(name, args.timeout_seconds, args.max_retries)
        safe = "".join(ch if ch.isalnum() else "_" for ch in name).strip("_")
        (raw_dir / f"{safe}.json").write_text(
            json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
        rows.append(summarize_match(name, payload))
        if index < len(names):
            time.sleep(args.sleep_seconds)

    review = pd.DataFrame(rows)
    review.to_csv(
        args.output_dir / "taxonomic_decisions_wcvp_review.csv",
        index=False,
        encoding="utf-8-sig",
        quoting=csv.QUOTE_MINIMAL,
    )
    high = review.loc[review["review_priority"].eq("high")].copy()
    high.to_csv(
        args.output_dir / "taxonomic_high_priority_review.csv",
        index=False,
        encoding="utf-8-sig",
    )
    pd.DataFrame(source_rows).to_csv(
        args.output_dir / "taxonomic_source_tables.csv",
        index=False,
    )

    summary = {
        "authority": "World Checklist of Vascular Plants (WCVP), accessed through the GBIF checklist API",
        "wcvp_dataset_key": WCVP_DATASET_KEY,
        "checked_date": date.today().isoformat(),
        "n_input_names": int(len(review)),
        "match_quality": review["wcvp_match_quality"].value_counts().to_dict(),
        "recommended_decisions": review["recommended_decision"].value_counts().to_dict(),
        "n_high_priority_review": int(review["review_priority"].eq("high").sum()),
        "decision_boundary": (
            "All decision values remain provisional. A botanist must review the table and replace "
            "provisional rows before analysis/audit_taxonomic_freeze.py can pass."
        ),
    }
    (args.output_dir / "taxonomic_wcvp_candidate_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
