#!/usr/bin/env python3
"""Validate the Chapter 1 multi-layer trait ontology.

The ontology is a contract: species-level literature evidence, individual-image
observations, and detector/model outputs must keep their distinct data layers.
This script validates that contract before source passages or annotation tasks
are created.
"""

from __future__ import annotations

import argparse
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

REQUIRED = {
    "trait_id", "display_name", "layer", "observation_unit", "allowed_states",
    "allow_multiple", "required_visual_context", "assessability_rule",
    "literature_extractable", "image_annotatable", "intended_use",
    "do_not_infer_from", "notes",
}
ALLOWED_LAYERS = {"detector", "image", "literature", "both"}
TRUE_VALUES = {"true", "1", "yes"}
FALSE_VALUES = {"false", "0", "no"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate a Chapter 1 trait ontology CSV.")
    parser.add_argument("--ontology", required=True)
    parser.add_argument("--out-dir", required=True)
    return parser.parse_args()


def text(value: object) -> str:
    return "" if pd.isna(value) else str(value).strip()


def parse_bool(value: object, field: str, trait_id: str) -> bool:
    normalized = text(value).lower()
    if normalized in TRUE_VALUES:
        return True
    if normalized in FALSE_VALUES:
        return False
    raise ValueError(f"trait_id={trait_id}: {field} must be true/false, got {value!r}")


def validate(ontology: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    missing = REQUIRED.difference(ontology.columns)
    if missing:
        raise ValueError(f"Ontology missing required columns: {sorted(missing)}")
    if ontology.empty:
        raise ValueError("Ontology is empty")

    work = ontology.copy()
    work["trait_id"] = work["trait_id"].map(text)
    if work["trait_id"].eq("").any():
        raise ValueError("Ontology contains blank trait_id")
    if work["trait_id"].duplicated().any():
        duplicate = work.loc[work["trait_id"].duplicated(keep=False), "trait_id"].tolist()
        raise ValueError(f"Ontology contains duplicate trait_id values: {duplicate}")
    invalid_ids = [value for value in work["trait_id"] if not re.fullmatch(r"[a-z][a-z0-9_]*", value)]
    if invalid_ids:
        raise ValueError(f"trait_id values must be lowercase snake_case: {invalid_ids}")

    messages: list[str] = []
    state_counts: list[int] = []
    for index, row in work.iterrows():
        trait_id = row["trait_id"]
        layer = text(row["layer"]).lower()
        if layer not in ALLOWED_LAYERS:
            raise ValueError(f"trait_id={trait_id}: unsupported layer {layer!r}")
        states = [state.strip() for state in text(row["allowed_states"]).split("|") if state.strip()]
        if len(states) < 2:
            raise ValueError(f"trait_id={trait_id}: allowed_states needs at least two pipe-separated states")
        if len(states) != len(set(states)):
            raise ValueError(f"trait_id={trait_id}: duplicate allowed state(s): {states}")
        bad_states = [state for state in states if not re.fullmatch(r"[a-z][a-z0-9_]*", state)]
        if bad_states:
            raise ValueError(f"trait_id={trait_id}: states must be lowercase snake_case: {bad_states}")
        allow_multiple = parse_bool(row["allow_multiple"], "allow_multiple", trait_id)
        literature = parse_bool(row["literature_extractable"], "literature_extractable", trait_id)
        image = parse_bool(row["image_annotatable"], "image_annotatable", trait_id)
        if layer == "detector" and literature:
            raise ValueError(f"trait_id={trait_id}: detector-only trait cannot be literature_extractable")
        if layer == "literature" and image:
            raise ValueError(f"trait_id={trait_id}: literature-only trait cannot be image_annotatable")
        if layer == "both" and not (literature and image):
            raise ValueError(f"trait_id={trait_id}: layer=both requires both literature_extractable and image_annotatable")
        if not text(row["assessability_rule"]):
            raise ValueError(f"trait_id={trait_id}: assessability_rule must not be blank")
        if not text(row["do_not_infer_from"]):
            raise ValueError(f"trait_id={trait_id}: do_not_infer_from must not be blank")
        state_counts.append(len(states))
        work.at[index, "allowed_states"] = "|".join(states)
        work.at[index, "allow_multiple"] = allow_multiple
        work.at[index, "literature_extractable"] = literature
        work.at[index, "image_annotatable"] = image
        if allow_multiple:
            messages.append(f"{trait_id}: multiple states are permitted; downstream evidence must use one row per state.")

    work["n_allowed_states"] = state_counts
    return work, messages


def main() -> None:
    args = parse_args()
    ontology_path = Path(args.ontology)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    ontology = pd.read_csv(ontology_path, dtype=str, keep_default_na=False)
    validated, messages = validate(ontology)
    validated.to_csv(out_dir / "trait_ontology_validated.csv", index=False, encoding="utf-8-sig")
    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "ontology": str(ontology_path.resolve()),
        "n_traits": int(len(validated)),
        "n_literature_extractable": int(validated["literature_extractable"].astype(bool).sum()),
        "n_image_annotatable": int(validated["image_annotatable"].astype(bool).sum()),
        "messages": messages,
    }
    (out_dir / "trait_ontology_validation.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
