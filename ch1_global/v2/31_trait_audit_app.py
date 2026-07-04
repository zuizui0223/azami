#!/usr/bin/env python3
"""Focused, blinded Streamlit review app for trait-specific audit tasks.

The app reads only the public task packet. It never loads CLIP candidates,
uncertainty margins, taxonomy, locality, or the private selection key.
"""

from __future__ import annotations

import argparse
import json
import os
import tempfile
from pathlib import Path
from typing import Any

import pandas as pd


REQUIRED = {"task_id", "annotation_unit_id", "trait_id", "source_image", "crop_path", "context_crop_path"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--tasks", required=True)
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


def resolve_path(value: Any, packet_root: Path) -> Path:
    path = Path(text(value))
    return path if path.is_absolute() else packet_root / path


def safe_image(value: Any, packet_root: Path):
    path = resolve_path(value, packet_root)
    if not path.is_file():
        return None
    try:
        from PIL import Image
        with Image.open(path) as opened:
            return opened.convert("RGB")
    except Exception:
        return None


def atomic_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8-sig", newline="", delete=False, dir=path.parent, suffix=".tmp") as handle:
        temporary = Path(handle.name)
        frame.to_csv(handle, index=False)
    os.replace(temporary, path)


def load_ontology(path: Path) -> dict[str, dict[str, Any]]:
    table = pd.read_csv(path, dtype=str, keep_default_na=False)
    required = {"trait_id", "display_name", "allowed_states", "allow_multiple", "required_visual_context", "assessability_rule"}
    if missing := required.difference(table.columns):
        raise ValueError(f"Ontology missing columns: {sorted(missing)}")
    result: dict[str, dict[str, Any]] = {}
    for row in table.to_dict("records"):
        trait_id = text(row["trait_id"])
        result[trait_id] = {
            "display_name": text(row["display_name"]),
            "states": [state for state in text(row["allowed_states"]).split("|") if state],
            "allow_multiple": as_bool(row["allow_multiple"]),
            "required_visual_context": text(row["required_visual_context"]),
            "assessability_rule": text(row["assessability_rule"]),
        }
    return result


def load_progress(tasks: pd.DataFrame, response_path: Path) -> pd.DataFrame:
    work = tasks.copy()
    for field in ("human_state", "notes", "task_complete"):
        if field not in work.columns:
            work[field] = ""
    if not response_path.exists():
        return work
    existing = pd.read_csv(response_path, dtype=str, keep_default_na=False)
    if "task_id" not in existing.columns or existing["task_id"].duplicated().any():
        raise ValueError("Existing response file lacks unique task_id")
    existing = existing.set_index("task_id")
    for field in ("human_state", "notes", "task_complete"):
        if field in existing.columns:
            work[field] = work["task_id"].map(existing[field]).fillna(work[field]).map(text)
    return work


def complete(row: pd.Series) -> bool:
    return text(row.get("human_state", "")) != "" and text(row.get("task_complete", "")).lower() == "yes"


def main() -> None:
    args = parse_args()
    try:
        import streamlit as st
    except ImportError as error:
        raise SystemExit("Streamlit is required: python -m pip install -r ch1_global/v2/requirements_annotation_app.txt") from error
    task_path = Path(args.tasks).resolve()
    packet_root = task_path.parent
    tasks = pd.read_csv(task_path, dtype=str, keep_default_na=False)
    if missing := REQUIRED.difference(tasks.columns):
        raise ValueError(f"Task packet missing columns: {sorted(missing)}")
    if tasks.empty or tasks["task_id"].duplicated().any():
        raise ValueError("Task IDs must be nonempty and unique")
    ontology = load_ontology(Path(args.ontology))
    unknown = sorted(set(tasks["trait_id"]).difference(ontology))
    if unknown:
        raise ValueError(f"Tasks refer to unknown ontology traits: {unknown}")
    response_path = Path(args.out_dir) / "trait_audit_responses.csv"
    status_path = Path(args.out_dir) / "trait_audit_status.json"
    data = load_progress(tasks, response_path)

    st.set_page_config(page_title="Cirsium Ch.1 focused trait audit", layout="wide")
    st.title("Cirsium Ch.1 — blinded focused trait audit")
    st.caption("Score visible evidence only. AI suggestions, species names, localities, and model confidence are hidden.")
    if "index" not in st.session_state:
        st.session_state.index = max(0, min(len(data) - 1, args.start_at - 1))
    index = max(0, min(len(data) - 1, st.session_state.index))
    completed = data.apply(complete, axis=1)
    st.sidebar.metric("Completed", f"{int(completed.sum())} / {len(data)}")
    st.sidebar.progress(float(completed.mean()))
    if st.sidebar.button("Jump to next incomplete") and (~completed).any():
        st.session_state.index = int((~completed).to_numpy().nonzero()[0][0])
        st.rerun()

    row = data.iloc[index]
    trait = ontology[text(row["trait_id"])]
    st.subheader(f"Task {index + 1} of {len(data)} — {trait['display_name']}")
    st.caption(f"Required view: {trait['required_visual_context']}  |  Assessability: {trait['assessability_rule']}")
    columns = st.columns(3)
    for column, (caption, field) in zip(columns, (("Head crop", "crop_path"), ("Head + peduncle context", "context_crop_path"), ("Source image", "source_image"))):
        with column:
            image = safe_image(row[field], packet_root)
            if image is not None:
                st.image(image, caption=caption, use_container_width=True)
            else:
                st.warning(f"{caption} unavailable")

    states = trait["states"]
    with st.form(key=f"audit_{row['task_id']}", clear_on_submit=False):
        if trait["allow_multiple"]:
            non_assessable = [state for state in states if state != "unassessable"]
            previous = [state for state in text(row.get("human_state", "")).split("|") if state]
            use_unassessable = st.checkbox("unassessable", value=previous == ["unassessable"])
            selected = st.multiselect(
                trait["display_name"],
                options=non_assessable,
                default=[] if use_unassessable else [state for state in previous if state in non_assessable],
                disabled=use_unassessable,
            )
            human_state = "unassessable" if use_unassessable else "|".join(selected)
        else:
            human_state = st.selectbox(
                trait["display_name"],
                options=["", *states],
                index=(["", *states].index(text(row.get("human_state", ""))) if text(row.get("human_state", "")) in ["", *states] else 0),
            )
        notes = st.text_area("Notes", value=text(row.get("notes", "")))
        previous_button, save_button, next_button = st.columns(3)
        go_previous = previous_button.form_submit_button("Save and previous")
        save = save_button.form_submit_button("Save")
        go_next = next_button.form_submit_button("Save and next")
    if go_previous or save or go_next:
        data.at[data.index[index], "human_state"] = human_state
        data.at[data.index[index], "notes"] = notes
        data.at[data.index[index], "task_complete"] = "yes" if human_state else "no"
        atomic_csv(data, response_path)
        status_path.parent.mkdir(parents=True, exist_ok=True)
        status_path.write_text(json.dumps({
            "tasks": str(task_path), "ontology": str(Path(args.ontology).resolve()), "annotator": args.annotator,
            "n_tasks": len(data), "n_complete": int(data.apply(complete, axis=1).sum()),
            "last_task_id": text(row["task_id"]),
        }, ensure_ascii=False, indent=2), encoding="utf-8")
        if go_previous:
            st.session_state.index = max(0, index - 1)
        elif go_next:
            st.session_state.index = min(len(data) - 1, index + 1)
        st.rerun()

    left, right = st.columns(2)
    if left.button("Previous without saving", disabled=index == 0):
        st.session_state.index = index - 1
        st.rerun()
    if right.button("Next without saving", disabled=index >= len(data) - 1):
        st.session_state.index = index + 1
        st.rerun()


if __name__ == "__main__":
    main()
