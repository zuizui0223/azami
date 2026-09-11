#!/usr/bin/env python3
"""Validate literature evidence rows against source passages and the trait ontology.

A row is valid only when its exact evidence quote occurs in the supplied source
passage and its state is permitted for the ontology trait. The validator never
uses taxon identity to fill missing evidence.
"""

from __future__ import annotations

import argparse
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

PASSAGE_REQUIRED = {"passage_id", "source_id", "accepted_taxon", "source_taxon_name", "passage_text"}
EVIDENCE_REQUIRED = {
    "evidence_id", "passage_id", "source_id", "accepted_taxon", "source_taxon_name",
    "trait_id", "trait_state", "scope", "life_stage", "geographic_scope",
    "evidence_quote", "evidence_translation", "llm_confidence", "extraction_method",
    "review_status", "reviewer_id", "notes",
}
VALID_REVIEW_STATUS = {"pending", "accepted", "rejected", "adjudicated"}
VALID_CONFIDENCE = {"high", "medium", "low", "na"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate evidence-grounded literature trait extraction.")
    parser.add_argument("--ontology", required=True)
    parser.add_argument("--passages", required=True)
    parser.add_argument("--evidence", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--strict", action="store_true", help="Stop on any invalid row instead of writing rejects.")
    return parser.parse_args()


def text(value: object) -> str:
    return "" if pd.isna(value) else str(value).strip()


def as_bool(value: object) -> bool:
    return text(value).lower() in {"true", "1", "yes"}


def normalize_for_quote(value: object) -> str:
    raw = text(value).casefold()
    return re.sub(r"\s+", " ", raw).strip()


def ontology_lookup(ontology: pd.DataFrame) -> dict[str, dict[str, object]]:
    required = {"trait_id", "allowed_states", "allow_multiple", "literature_extractable"}
    missing = required.difference(ontology.columns)
    if missing:
        raise ValueError(f"Ontology missing columns: {sorted(missing)}")
    lookup: dict[str, dict[str, object]] = {}
    for row in ontology.itertuples(index=False):
        trait_id = text(getattr(row, "trait_id"))
        lookup[trait_id] = {
            "states": set(state for state in text(getattr(row, "allowed_states")).split("|") if state),
            "allow_multiple": as_bool(getattr(row, "allow_multiple")),
            "literature_extractable": as_bool(getattr(row, "literature_extractable")),
        }
    return lookup


def validate_row(row: pd.Series, passages: pd.DataFrame, traits: dict[str, dict[str, object]]) -> str | None:
    passage_id = text(row["passage_id"])
    if passage_id not in passages.index:
        return f"unknown passage_id={passage_id!r}"
    passage = passages.loc[passage_id]
    if text(row["source_id"]) != text(passage["source_id"]):
        return "source_id does not match referenced passage"
    if text(row["accepted_taxon"]) != text(passage["accepted_taxon"]):
        return "accepted_taxon does not match referenced passage"
    trait_id = text(row["trait_id"])
    if trait_id not in traits:
        return f"unknown trait_id={trait_id!r}"
    trait = traits[trait_id]
    if not trait["literature_extractable"]:
        return f"trait_id={trait_id!r} is not literature_extractable"
    state = text(row["trait_state"])
    if state not in trait["states"]:
        return f"trait_state={state!r} is not allowed for trait_id={trait_id!r}"
    quote = normalize_for_quote(row["evidence_quote"])
    if not quote:
        return "evidence_quote is blank"
    if quote not in normalize_for_quote(passage["passage_text"]):
        return "evidence_quote is not an exact contiguous quote from the supplied passage"
    status = text(row["review_status"]).lower()
    if status not in VALID_REVIEW_STATUS:
        return f"unsupported review_status={status!r}"
    confidence = text(row["llm_confidence"]).lower()
    if confidence not in VALID_CONFIDENCE:
        return f"unsupported llm_confidence={confidence!r}"
    if status in {"accepted", "adjudicated"} and not text(row["reviewer_id"]):
        return f"review_status={status!r} requires reviewer_id"
    if not text(row["extraction_method"]):
        return "extraction_method is blank"
    return None


def main() -> None:
    args = parse_args()
    ontology = pd.read_csv(args.ontology, dtype=str, keep_default_na=False)
    passages = pd.read_csv(args.passages, dtype=str, keep_default_na=False)
    evidence = pd.read_csv(args.evidence, dtype=str, keep_default_na=False)
    missing = PASSAGE_REQUIRED.difference(passages.columns)
    if missing:
        raise ValueError(f"Passage table missing required columns: {sorted(missing)}")
    missing = EVIDENCE_REQUIRED.difference(evidence.columns)
    if missing:
        raise ValueError(f"Evidence table missing required columns: {sorted(missing)}")
    if passages["passage_id"].astype(str).duplicated().any():
        raise ValueError("passage_id must be unique")
    if evidence["evidence_id"].astype(str).eq("").any() or evidence["evidence_id"].astype(str).duplicated().any():
        raise ValueError("evidence_id must be nonblank and unique")

    passages = passages.set_index("passage_id", drop=False)
    traits = ontology_lookup(ontology)
    valid_rows: list[dict[str, object]] = []
    reject_rows: list[dict[str, object]] = []
    for _, row in evidence.iterrows():
        reason = validate_row(row, passages, traits)
        record = row.to_dict()
        if reason:
            record["validation_error"] = reason
            reject_rows.append(record)
        else:
            valid_rows.append(record)
    if args.strict and reject_rows:
        examples = [row["validation_error"] for row in reject_rows[:5]]
        raise ValueError(f"Found {len(reject_rows)} invalid evidence rows: {examples}")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    valid = pd.DataFrame(valid_rows, columns=list(evidence.columns))
    rejects = pd.DataFrame(reject_rows, columns=list(evidence.columns) + ["validation_error"])
    valid.to_csv(out_dir / "literature_evidence_valid.csv", index=False, encoding="utf-8-sig")
    rejects.to_csv(out_dir / "literature_evidence_rejected.csv", index=False, encoding="utf-8-sig")
    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "n_input_rows": int(len(evidence)),
        "n_valid_rows": int(len(valid_rows)),
        "n_rejected_rows": int(len(reject_rows)),
        "note": "Valid evidence is still not image ground truth. Accepted/adjudicated review status is required before species-level synthesis.",
    }
    (out_dir / "literature_evidence_validation_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
