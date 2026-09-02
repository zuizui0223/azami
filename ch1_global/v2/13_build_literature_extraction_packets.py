#!/usr/bin/env python3
"""Build evidence-grounded LLM extraction packets from curated literature passages.

This tool does not call an LLM. It prepares deterministic JSONL tasks from
passages that the researcher has lawfully obtained and reviewed. Each task
requires an exact quote from the supplied passage, so an LLM cannot silently
turn species-level descriptions into unsupported image labels.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

PASSAGE_REQUIRED = {
    "passage_id", "source_id", "accepted_taxon", "source_taxon_name",
    "page_or_figure", "language", "passage_text", "rights_checked",
}
EVIDENCE_COLUMNS = [
    "evidence_id", "passage_id", "source_id", "accepted_taxon", "source_taxon_name",
    "trait_id", "trait_state", "scope", "life_stage", "geographic_scope",
    "evidence_quote", "evidence_translation", "llm_confidence", "extraction_method",
    "review_status", "reviewer_id", "notes",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build evidence-grounded LLM tasks from curated literature passages.")
    parser.add_argument("--ontology", required=True)
    parser.add_argument("--passages", required=True, help="CSV of curated passages, not URLs or entire unreviewed books.")
    parser.add_argument("--out-dir", required=True)
    return parser.parse_args()


def text(value: object) -> str:
    return "" if pd.isna(value) else str(value).strip()


def as_bool(value: object) -> bool:
    return text(value).lower() in {"true", "1", "yes"}


def task_prompt(allowed_traits: list[dict[str, object]], row: pd.Series) -> str:
    trait_lines = "\n".join(
        f"- {item['trait_id']}: allowed states = {', '.join(item['states'])}; observation unit = {item['observation_unit']}"
        for item in allowed_traits
    )
    return f"""You are extracting species-level botanical evidence from one supplied passage.

Rules:
1. Extract only statements explicitly supported by the passage. Do not use taxonomic knowledge, images, web search, or inference.
2. For each extracted row, `evidence_quote` must be an exact contiguous quote from the passage in its original language.
3. Use only the trait IDs and states listed below. Do not create a row for absence of information.
4. Do not convert a species description into an image-level label. This task records literature evidence only.
5. Return valid JSON only, with keys `passage_id` and `evidence_rows`. Return an empty `evidence_rows` list when no allowed trait is explicitly supported.
6. Set `review_status` to `pending`; a human reviewer decides acceptance later.

Allowed traits:
{trait_lines}

Passage metadata:
- passage_id: {text(row['passage_id'])}
- source_id: {text(row['source_id'])}
- accepted_taxon: {text(row['accepted_taxon'])}
- source_taxon_name: {text(row['source_taxon_name'])}
- page_or_figure: {text(row['page_or_figure'])}
- language: {text(row['language'])}

Passage:
---
{text(row['passage_text'])}
---

Expected JSON shape:
{{
  "passage_id": "{text(row['passage_id'])}",
  "evidence_rows": [
    {{
      "trait_id": "capitulum_orientation",
      "trait_state": "nodding",
      "scope": "typical flowering heads",
      "life_stage": "anthesis",
      "geographic_scope": "as stated in source",
      "evidence_quote": "exact original-language quote",
      "evidence_translation": "optional translation",
      "llm_confidence": "high|medium|low",
      "extraction_method": "llm_assisted",
      "review_status": "pending",
      "notes": ""
    }}
  ]
}}"""


def main() -> None:
    args = parse_args()
    ontology = pd.read_csv(args.ontology, dtype=str, keep_default_na=False)
    passages = pd.read_csv(args.passages, dtype=str, keep_default_na=False)
    missing = PASSAGE_REQUIRED.difference(passages.columns)
    if missing:
        raise ValueError(f"Passage table missing required columns: {sorted(missing)}")
    if passages["passage_id"].astype(str).duplicated().any():
        raise ValueError("passage_id must be unique")
    if not passages["rights_checked"].map(as_bool).all():
        bad = passages.loc[~passages["rights_checked"].map(as_bool), "passage_id"].tolist()
        raise ValueError(f"Only passages with rights_checked=true may be packetized: {bad}")
    if passages["passage_text"].map(text).eq("").any():
        raise ValueError("passage_text must not be blank")

    literature_traits = ontology.loc[ontology["literature_extractable"].map(as_bool)].copy()
    if literature_traits.empty:
        raise ValueError("Ontology contains no literature_extractable traits")
    allowed_traits = [
        {
            "trait_id": text(row.trait_id),
            "states": [state for state in text(row.allowed_states).split("|") if state],
            "observation_unit": text(row.observation_unit),
        }
        for row in literature_traits.itertuples(index=False)
    ]

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows: list[dict[str, object]] = []
    task_path = out_dir / "literature_extraction_tasks.jsonl"
    with task_path.open("w", encoding="utf-8") as handle:
        for index, row in passages.iterrows():
            passage_id = text(row["passage_id"])
            digest = hashlib.sha256(text(row["passage_text"]).encode("utf-8")).hexdigest()
            task_id = f"lit_{index + 1:06d}"
            payload = {
                "task_id": task_id,
                "passage_id": passage_id,
                "source_id": text(row["source_id"]),
                "accepted_taxon": text(row["accepted_taxon"]),
                "source_taxon_name": text(row["source_taxon_name"]),
                "page_or_figure": text(row["page_or_figure"]),
                "language": text(row["language"]),
                "passage_sha256": digest,
                "prompt": task_prompt(allowed_traits, row),
            }
            handle.write(json.dumps(payload, ensure_ascii=False) + "\n")
            manifest_rows.append({
                "task_id": task_id,
                "passage_id": passage_id,
                "source_id": payload["source_id"],
                "accepted_taxon": payload["accepted_taxon"],
                "source_taxon_name": payload["source_taxon_name"],
                "page_or_figure": payload["page_or_figure"],
                "language": payload["language"],
                "passage_sha256": digest,
            })
    pd.DataFrame(manifest_rows).to_csv(out_dir / "literature_extraction_task_manifest.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(columns=EVIDENCE_COLUMNS).to_csv(out_dir / "literature_evidence_template.csv", index=False, encoding="utf-8-sig")
    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "n_passages": int(len(passages)),
        "n_literature_traits": int(len(literature_traits)),
        "task_file": task_path.name,
        "note": "Tasks require exact original-language evidence quotes. LLM output remains pending until human review.",
    }
    (out_dir / "literature_extraction_packet_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
