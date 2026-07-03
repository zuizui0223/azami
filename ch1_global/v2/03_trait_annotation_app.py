#!/usr/bin/env python3
"""Local Streamlit app for blinded Chapter 1 trait annotation.

This app is deliberately limited to human-scored image traits and image-quality
flags. Detector boxes, keypoints, and segmentation masks are separate tasks.
The app reads a blinded packet and writes only a blinded response CSV.

Install once
------------
python -m pip install -r ch1_global/v2/requirements_annotation_app.txt

Run
---
streamlit run ch1_global/v2/03_trait_annotation_app.py -- \
  --packet "C:\\Users\\zuizui\\cirsium_inat\\ch1_v2_annotation\\packets\\annotator_1_packet.csv" \
  --out-dir "C:\\Users\\zuizui\\cirsium_inat\\ch1_v2_annotation\\responses\\annotator_1"

The app never displays taxonomy, coordinates, model predictions, or source
metadata. It stores an atomic progress CSV after every Save action, so closing
the browser does not lose completed rows.
"""

from __future__ import annotations

import argparse
import json
import os
import tempfile
from pathlib import Path
from typing import Any

import pandas as pd


TASK_COLUMNS = {
    "annotation_unit_id",
    "source_image",
    "crop_path",
    "needs_trait_labels",
}

RESPONSE_FIELDS: dict[str, list[str]] = {
    "flowering_qc": [
        "",
        "usable",
        "not_flowering",
        "too_small_or_blurred",
        "occluded",
        "wrong_taxon_or_non-Cirsium",
        "uncertain",
    ],
    "camera_tilt_flag": [
        "",
        "vertical_reference_clear",
        "vertical_reference_uncertain",
        "tilted_or_rotated",
        "not_assessable",
    ],
    "view_angle_class": ["", "lateral", "oblique", "frontal", "not_assessable"],
    "head_orientation_class": ["", "upright", "inclined", "nodding", "not_assessable"],
    "flower_colour_class": ["", "white", "pale", "dark", "not_assessable"],
    "colour_qc": ["", "usable", "illumination_or_white_balance_unusable", "not_assessable"],
    "involucral_cover_score_0_to_4": ["", "0", "1", "2", "3", "4", "not_assessable"],
    "head_shape_class": ["", "globose", "depressed", "cylindrical", "not_assessable"],
    "corolla_display_score_0_to_4": ["", "0", "1", "2", "3", "4", "not_assessable"],
}

RESPONSE_COLUMNS = [*RESPONSE_FIELDS, "annotation_complete", "notes"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--packet", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--annotator", default="")
    parser.add_argument("--start-at", type=int, default=1)
    args, _ = parser.parse_known_args()
    return args


def normalize(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def validate_packet(packet: pd.DataFrame) -> None:
    missing = TASK_COLUMNS.difference(packet.columns)
    if missing:
        raise ValueError(f"Packet is missing required columns: {sorted(missing)}")
    if packet["annotation_unit_id"].duplicated().any():
        raise ValueError("Packet has duplicated annotation_unit_id values")


def merge_existing(packet: pd.DataFrame, existing_path: Path) -> pd.DataFrame:
    """Load existing blinded responses by ID without changing packet row order."""
    work = packet.copy()
    for column in RESPONSE_COLUMNS:
        if column not in work.columns:
            work[column] = ""
    if not existing_path.exists():
        return work

    existing = pd.read_csv(existing_path, dtype=str, keep_default_na=False)
    if "annotation_unit_id" not in existing.columns:
        raise ValueError(f"Existing progress file lacks annotation_unit_id: {existing_path}")
    if existing["annotation_unit_id"].duplicated().any():
        raise ValueError(f"Existing progress file has duplicate IDs: {existing_path}")

    existing = existing.set_index("annotation_unit_id")
    for column in RESPONSE_COLUMNS:
        if column in existing.columns:
            work[column] = work["annotation_unit_id"].map(existing[column]).fillna(work[column]).map(normalize)
    return work


def is_complete(row: pd.Series) -> bool:
    return normalize(row.get("annotation_complete", "")).lower() == "yes"


def response_ready(values: dict[str, str]) -> bool:
    """Require all human-scored fields to be considered before completion."""
    return all(normalize(values.get(field, "")) for field in RESPONSE_FIELDS)


def atomic_write_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8-sig", newline="", delete=False, dir=path.parent, suffix=".tmp") as handle:
        temp_path = Path(handle.name)
        df.to_csv(handle, index=False)
    os.replace(temp_path, path)


def atomic_write_json(payload: dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", delete=False, dir=path.parent, suffix=".tmp") as handle:
        temp_path = Path(handle.name)
        json.dump(payload, handle, ensure_ascii=False, indent=2)
    os.replace(temp_path, path)


def safe_image(path_text: Any):
    path = Path(normalize(path_text))
    if not path.exists():
        return None
    try:
        from PIL import Image

        return Image.open(path).convert("RGB")
    except Exception:
        return None


def option_index(options: list[str], current: Any) -> int:
    text = normalize(current)
    return options.index(text) if text in options else 0


def main() -> None:
    args = parse_args()
    try:
        import streamlit as st
    except ImportError as error:
        raise SystemExit(
            "Streamlit is required. Install with: python -m pip install -r "
            "ch1_global/v2/requirements_annotation_app.txt"
        ) from error

    packet_path = Path(args.packet)
    out_dir = Path(args.out_dir)
    response_path = out_dir / "annotation_responses.csv"
    status_path = out_dir / "annotation_status.json"

    st.set_page_config(page_title="Cirsium Ch.1 blinded annotation", layout="wide")
    st.title("Cirsium Ch.1 — blinded trait annotation")
    st.caption("Judge visible traits only. Species, locality, and model predictions are deliberately hidden.")

    try:
        packet = pd.read_csv(packet_path, dtype=str, keep_default_na=False)
        validate_packet(packet)
        data = merge_existing(packet, response_path)
    except Exception as error:
        st.error(str(error))
        st.stop()

    if "annotation_index" not in st.session_state:
        st.session_state.annotation_index = max(0, min(len(data) - 1, args.start_at - 1))
    if st.session_state.annotation_index >= len(data):
        st.session_state.annotation_index = max(0, len(data) - 1)

    completed = data.apply(is_complete, axis=1)
    n_complete = int(completed.sum())
    st.sidebar.metric("Completed", f"{n_complete} / {len(data)}")
    st.sidebar.progress(n_complete / len(data) if len(data) else 0.0)
    st.sidebar.caption("Responses save atomically after every Save action.")

    jump_candidates = [i + 1 for i, done in enumerate(completed.tolist()) if not done]
    if jump_candidates and st.sidebar.button("Jump to next incomplete"):
        st.session_state.annotation_index = jump_candidates[0] - 1
        st.rerun()

    index = st.session_state.annotation_index
    row = data.iloc[index]
    st.subheader(f"Item {index + 1} of {len(data)}")

    left, right = st.columns(2)
    with left:
        crop = safe_image(row.get("crop_path", ""))
        if crop is not None:
            st.image(crop, caption="Head crop", use_container_width=True)
        else:
            st.warning("Head crop could not be opened. Mark image quality as appropriate.")
    with right:
        source = safe_image(row.get("source_image", ""))
        if source is not None:
            st.image(source, caption="Source image", use_container_width=True)
        else:
            st.info("Source image unavailable; score only what is visible in the crop.")

    st.markdown("**Definitions:** Orientation uses a visible biological vertical reference, not the image frame. Cover/display are primary only for lateral views.")
    with st.form(key=f"annotation_form_{row['annotation_unit_id']}", clear_on_submit=False):
        response: dict[str, str] = {}
        col1, col2, col3 = st.columns(3)
        field_columns = [col1, col2, col3]
        for position, (field, options) in enumerate(RESPONSE_FIELDS.items()):
            with field_columns[position % 3]:
                response[field] = st.selectbox(
                    field.replace("_", " "),
                    options=options,
                    index=option_index(options, row.get(field, "")),
                    key=f"{row['annotation_unit_id']}_{field}",
                )
        notes = st.text_area("Notes (optional)", value=normalize(row.get("notes", "")))
        previous, save, save_next = st.columns(3)
        go_previous = previous.form_submit_button("Save and previous")
        save_here = save.form_submit_button("Save")
        go_next = save_next.form_submit_button("Save and next")

    if go_previous or save_here or go_next:
        response["notes"] = notes
        response["annotation_complete"] = "yes" if response_ready(response) else "no"
        for column, value in response.items():
            data.at[data.index[index], column] = value
        atomic_write_csv(data, response_path)
        atomic_write_json(
            {
                "packet": str(packet_path.resolve()),
                "response_file": str(response_path.resolve()),
                "annotator": args.annotator,
                "n_total": len(data),
                "n_complete": int(data.apply(is_complete, axis=1).sum()),
                "last_annotation_unit_id": normalize(row["annotation_unit_id"]),
                "last_index_1_based": index + 1,
            },
            status_path,
        )
        if go_previous:
            st.session_state.annotation_index = max(0, index - 1)
        elif go_next:
            st.session_state.annotation_index = min(len(data) - 1, index + 1)
        st.rerun()

    navigation_left, navigation_right = st.columns(2)
    if navigation_left.button("Previous without saving", disabled=index == 0):
        st.session_state.annotation_index = index - 1
        st.rerun()
    if navigation_right.button("Next without saving", disabled=index >= len(data) - 1):
        st.session_state.annotation_index = index + 1
        st.rerun()


if __name__ == "__main__":
    main()
