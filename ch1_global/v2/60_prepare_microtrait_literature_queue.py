#!/usr/bin/env python3
"""Create an evidence-first literature queue for hair and mucilage traits.

The output deliberately separates species-level textual evidence from individual
image measurements. An LLM may extract evidence snippets, but a state is not
accepted without a source, page/figure and verbatim supporting text.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

TRAITS = {
    "involucre_indumentum": {
        "allowed_states": "glabrous|webby|tomentose|glandular|mixed|unknown",
        "required_evidence": "description of involucre or phyllary surface",
    },
    "external_mucilage_visible": {
        "allowed_states": "present|absent|unknown",
        "required_evidence": "explicit glandular/sticky/viscid/mucilage statement or macro observation",
    },
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--species", required=True, help="CSV containing taxon_name")
    p.add_argument("--out-csv", required=True)
    p.add_argument("--schema-json", required=True)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    species = pd.read_csv(args.species, dtype=str, keep_default_na=False)
    if "taxon_name" not in species.columns:
        raise ValueError("Species table lacks taxon_name")
    taxa = sorted(value.strip() for value in species["taxon_name"].unique() if value.strip())
    rows = []
    for taxon in taxa:
        for trait, definition in TRAITS.items():
            rows.append({
                "accepted_taxon": taxon,
                "source_taxon_name": "",
                "trait_id": trait,
                "trait_state": "unknown",
                "allowed_states": definition["allowed_states"],
                "evidence_scope": "species_level_literature_or_macro_observation",
                "life_stage": "flowering_head",
                "geographic_scope": "",
                "source_id": "",
                "source_type": "",
                "source_year": "",
                "page_or_figure": "",
                "evidence_quote": "",
                "evidence_paraphrase": "",
                "required_evidence": definition["required_evidence"],
                "extractor": "",
                "llm_confidence": "",
                "review_status": "pending",
                "reviewer": "",
                "conflict_group": "",
                "notes": "",
            })
    output = pd.DataFrame(rows)
    Path(args.out_csv).parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.out_csv, index=False, encoding="utf-8-sig")
    schema = {
        "semantic_boundary": (
            "Rows are source-backed species-level knowledge. They are never individual-image labels "
            "and must not be merged into image measurements without retaining evidence provenance."
        ),
        "traits": TRAITS,
        "required_for_reviewed_state": [
            "accepted_taxon", "trait_id", "trait_state", "source_id",
            "page_or_figure", "evidence_quote", "review_status=accepted",
        ],
        "absence_rule": (
            "Absence of a statement is unknown, not absent. external_mucilage_visible=absent requires "
            "an explicit negative statement or a reviewed macro observation of the relevant surface."
        ),
        "llm_role": "evidence extraction and synonym normalization only; no unsupported completion",
        "n_taxa": len(taxa),
        "n_queue_rows": len(output),
    }
    Path(args.schema_json).write_text(json.dumps(schema, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(schema, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
