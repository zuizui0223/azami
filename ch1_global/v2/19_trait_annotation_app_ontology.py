#!/usr/bin/env python3
"""Ontology-driven local Streamlit app for blinded Chapter 1 trait annotation.

It displays tight head, optional head-plus-peduncle context, and source images
without taxonomy, locality, or model output. Trait states come directly from
the versioned ontology rather than a hard-coded v1 label list.
"""

from __future__ import annotations

import argparse
import json
import os
import tempfile
from pathlib import Path
from typing import Any

import pandas as pd

TASK_REQUIRED = {"annotation_unit_id", "source_image", "crop_path"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--packet", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--ontology", default=str(Path(__file__).with_name("ontology") / "ch1_trait_ontology.csv"))
    parser.add_argument("--annotator", default="")
    parser.add_argument("--start-at", type=int, default=1)
    args, _ = parser.parse_known_args()
    return args


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def as_bool(value: Any) -> bool:
    return text(value).lower() in {"true", "1", "yes"}


def load_traits(path: Path) -> list[dict[str, Any]]:
    ontology = pd.read_csv(path, dtype=str, keep_default_na=False)
    required = {"trait_id", "display_name", "layer", "allowed_states", "image_annotatable", "assessability_rule"}
    missing = required.difference(ontology.columns)
    if missing:
        raise ValueError(f"Ontology is missing columns: {sorted(missing)}")
    traits = []
    for row in ontology.itertuples(index=False):
        if not as_bool(getattr(row, "image_annotatable")):
            continue
        if text(getattr(row, "layer")).lower() == "detector":
            continue
        trait_id = text(getattr(row, "trait_id"))
        traits.append({
            "trait_id": trait_id,
            "field": f"trait__{trait_id}",
            "display_name": text(getattr(row, "display_name")),
            "states": [state for state in text(getattr(row, "allowed_states")).split("|") if state],
            "assessability_rule": text(getattr(row, "assessability_rule")),
            "required_visual_context": text(getattr(row, "required_visual_context", "")),
        })
    if not traits:
        raise ValueError("Ontology has no image-annotatable non-detector traits")
    return traits


def safe_image(path_text: Any):
    path = Path(text(path_text))
    if not path.exists():
        return None
    try:
        from PIL import Image
        return Image.open(path).convert("RGB")
    except Exception:
        return None


def atomic_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8-sig", newline="", delete=False, dir=path.parent, suffix=".tmp") as handle:
        temporary = Path(handle.name)
        frame.to_csv(handle, index=False)
    os.replace(temporary, path)


def load_progress(packet: pd.DataFrame, path: Path, response_fields: list[str]) -> pd.DataFrame:
    work = packet.copy()
    for field in response_fields:
        if field not in work.columns:
            work[field] = ""
    if not path.exists():
        return work
    existing = pd.read_csv(path, dtype=str, keep_default_na=False)
    if "annotation_unit_id" not in existing.columns or existing["annotation_unit_id"].duplicated().any():
        raise ValueError("Existing response file lacks unique annotation_unit_id")
    existing = existing.set_index("annotation_unit_id")
    for field in response_fields:
        if field in existing.columns:
            work[field] = work["annotation_unit_id"].map(existing[field]).fillna(work[field]).map(text)
    return work


def ready(row: pd.Series, trait_fields: list[str]) -> bool:
    return all(text(row.get(field, "")) for field in trait_fields) and text(row.get("annotation_complete", "")).lower() == "yes"


def main() -> None:
    args = parse_args()
    try:
        import streamlit as st
    except ImportError as error:
        raise SystemExit("Streamlit is required: python -m pip install -r ch1_global/v2/requirements_annotation_app.txt") from error
    packet = pd.read_csv(args.packet, dtype=str, keep_default_na=False)
    missing = TASK_REQUIRED.difference(packet.columns)
    if missing:
        raise ValueError(f"Packet is missing columns: {sorted(missing)}")
    if packet["annotation_unit_id"].duplicated().any():
        raise ValueError("Packet has duplicate annotation_unit_id")
    traits = load_traits(Path(args.ontology))
    trait_fields = [trait["field"] for trait in traits]
    response_fields = [*trait_fields, "annotation_complete", "notes"]
    out_dir = Path(args.out_dir)
    response_path = out_dir / "ontology_annotation_responses.csv"
    status_path = out_dir / "ontology_annotation_status.json"
    data = load_progress(packet, response_path, response_fields)

    st.set_page_config(page_title="Cirsium Ch.1 ontology annotation", layout="wide")
    st.title("Cirsium Ch.1 — blinded ontology annotation")
    st.caption("Score only visible evidence. Taxonomy, location, source prose, and model predictions are hidden.")
    if "index" not in st.session_state:
        st.session_state.index = max(0, min(len(data) - 1, args.start_at - 1))
    index = max(0, min(len(data) - 1, st.session_state.index))
    completed = data.apply(lambda row: ready(row, trait_fields), axis=1)
    st.sidebar.metric("Completed", f"{int(completed.sum())} / {len(data)}")
    st.sidebar.progress(float(completed.mean()) if len(data) else 0.0)
    if st.sidebar.button("Jump to next incomplete") and (~completed).any():
        st.session_state.index = int((~completed).to_numpy().nonzero()[0][0])
        st.rerun()

    row = data.iloc[index]
    st.subheader(f"Item {index + 1} of {len(data)}")
    columns = st.columns(3)
    panels = [
        ("Head crop", "crop_path"),
        ("Head + peduncle context", "context_crop_path"),
        ("Source image", "source_image"),
    ]
    for column, (caption, field) in zip(columns, panels):
        with column:
            image = safe_image(row.get(field, ""))
            if image is not None:
                st.image(image, caption=caption, use_container_width=True)
            elif field == "context_crop_path":
                st.info("Context crop unavailable; use source image for orientation only when the peduncle is visible.")
            else:
                st.warning(f"{caption} unavailable")

    with st.form(key=f"trait_form_{row['annotation_unit_id']}", clear_on_submit=False):
        values: dict[str, str] = {}
        form_columns = st.columns(3)
        for position, trait in enumerate(traits):
            with form_columns[position % 3]:
                values[trait["field"]] = st.selectbox(
                    trait["display_name"],
                    options=["", *trait["states"]],
                    index=(["", *trait["states"]].index(text(row.get(trait["field"], ""))) if text(row.get(trait["field"], "")) in ["", *trait["states"]] else 0),
                    help=f"Required view: {trait['required_visual_context']}\n\nAssessability: {trait['assessability_rule']}",
                    key=f"{row['annotation_unit_id']}_{trait['field']}",
                )
        notes = st.text_area("Notes", value=text(row.get("notes", "")))
        previous, save, next_button = st.columns(3)
        go_previous = previous.form_submit_button("Save and previous")
        save_here = save.form_submit_button("Save")
        go_next = next_button.form_submit_button("Save and next")
    if go_previous or save_here or go_next:
        for field, value in values.items():
            data.at[data.index[index], field] = value
        data.at[data.index[index], "notes"] = notes
        data.at[data.index[index], "annotation_complete"] = "yes" if all(values.values()) else "no"
        atomic_csv(data, response_path)
        status_path.parent.mkdir(parents=True, exist_ok=True)
        status_path.write_text(json.dumps({
            "packet": str(Path(args.packet).resolve()), "ontology": str(Path(args.ontology).resolve()),
            "annotator": args.annotator, "n_total": len(data),
            "n_complete": int(data.apply(lambda row: ready(row, trait_fields), axis=1).sum()),
            "last_annotation_unit_id": text(row["annotation_unit_id"]),
        }, ensure_ascii=False, indent=2), encoding="utf-8")
        if go_previous:
            st.session_state.index = max(0, index - 1)
        elif go_next:
            st.session_state.index = min(len(data) - 1, index + 1)
        st.rerun()

    nav_left, nav_right = st.columns(2)
    if nav_left.button("Previous without saving", disabled=index == 0):
        st.session_state.index = index - 1
        st.rerun()
    if nav_right.button("Next without saving", disabled=index >= len(data) - 1):
        st.session_state.index = index + 1
        st.rerun()


if __name__ == "__main__":
    main()
