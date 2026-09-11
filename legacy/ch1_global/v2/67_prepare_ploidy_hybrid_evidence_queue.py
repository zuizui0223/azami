#!/usr/bin/env python3
"""Create a source-backed ploidy/hybrid evidence queue for tree sensitivities."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", required=True)
    parser.add_argument("--out-csv", required=True)
    parser.add_argument("--schema-json", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    species = pd.read_csv(args.species, dtype=str, keep_default_na=False)
    if "taxon_name" not in species.columns:
        raise ValueError("Species table lacks taxon_name")
    taxa = sorted({name.strip() for name in species["taxon_name"] if name.strip()})
    rows = []
    for taxon in taxa:
        rows.append({
            "accepted_taxon": taxon,
            "source_taxon_name": "",
            "chromosome_number_2n": "",
            "base_chromosome_number_x": "",
            "ploidy_level_numeric": "",
            "ploidy_class": "unknown",
            "polyploid_origin": "unknown",
            "hybrid_status": "unknown",
            "putative_parent_1": "",
            "putative_parent_2": "",
            "nuclear_plastid_conflict_reported": "unknown",
            "source_id": "",
            "source_type": "",
            "source_year": "",
            "page_table_figure": "",
            "evidence_quote": "",
            "geographic_scope": "",
            "review_status": "pending",
            "reviewer": "",
            "notes": "",
        })
    output = pd.DataFrame(rows)
    Path(args.out_csv).parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.out_csv, index=False, encoding="utf-8-sig")
    schema = {
        "n_taxa": len(taxa),
        "allowed_ploidy_class": ["diploid", "polyploid", "mixed_within_taxon", "unknown"],
        "allowed_polyploid_origin": ["auto", "allo", "uncertain_polyploid", "not_applicable", "unknown"],
        "allowed_hybrid_status": ["documented_hybrid", "suspected_hybrid", "no_evidence_reviewed", "unknown"],
        "required_for_accepted_state": [
            "accepted_taxon", "source_id", "page_table_figure", "evidence_quote", "review_status=accepted"
        ],
        "absence_rule": "No report is unknown, not diploid and not non-hybrid.",
        "analysis_rule": (
            "Primary sensitivity excludes only source-reviewed documented/suspected hybrids or polyploids. "
            "Unknown taxa remain unknown and are never silently coded as diploid."
        ),
        "llm_role": "Extract chromosome, cytotype and hybrid evidence with quotes; do not infer missing states.",
    }
    Path(args.schema_json).write_text(json.dumps(schema, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(schema, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
