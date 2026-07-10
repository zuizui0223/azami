#!/usr/bin/env python3
"""Focused, blinded Streamlit review app for trait-specific audit tasks.

Local mode reads an on-disk public task packet and writes progress atomically::

    streamlit run ch1_global/v2/31_trait_audit_app.py -- \
      --tasks /path/to/blinded_trait_audit_tasks.csv \
      --out-dir /path/to/responses/annotator

When launched without ``--tasks`` and ``--out-dir`` (for example on Streamlit
Community Cloud), the app switches to upload mode. The annotator uploads only
the public blinded packet ZIP and may optionally upload a previous response CSV
to resume. Cloud storage is temporary, so the current response CSV is always
available as a browser download.

The app never loads CLIP candidates, uncertainty margins, taxonomy, locality,
or the private selection key.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import re
import shutil
import stat
import tempfile
import uuid
import zipfile
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Any

import pandas as pd


REQUIRED = {"task_id", "annotation_unit_id", "trait_id", "source_image", "crop_path", "context_crop_path"}
FORBIDDEN_PRIVATE_FIELDS = {
    "ai_candidate_state",
    "runner_up_state",
    "similarity_margin",
    "selection_stratum",
    "double_label",
}
FORBIDDEN_ARCHIVE_TOKENS = {
    "private_audit_key",
    "trait_audit_selection_key",
    "representative_trait_holdout_key",
}
PREFERRED_TASK_FILENAMES = (
    "blinded_trait_audit_tasks.csv",
    "representative_trait_holdout_tasks.csv",
)
DEFAULT_ONTOLOGY = Path(__file__).with_name("ontology") / "ch1_trait_ontology.csv"
MAX_UNCOMPRESSED_PACKET_BYTES = 2 * 1024 * 1024 * 1024


@dataclass(frozen=True)
class AppConfig:
    task_path: Path
    out_dir: Path
    ontology_path: Path
    annotator: str
    start_at: int
    cloud_mode: bool


@dataclass(frozen=True)
class PublicPacket:
    task_path: Path
    ontology_path: Path | None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--tasks")
    parser.add_argument("--out-dir")
    parser.add_argument("--ontology", default=str(DEFAULT_ONTOLOGY))
    parser.add_argument("--annotator", default="")
    parser.add_argument("--start-at", type=int, default=1)
    args, _ = parser.parse_known_args()
    return args


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def as_bool(value: Any) -> bool:
    return text(value).lower() in {"true", "1", "yes"}


def safe_component(value: str, fallback: str = "annotator") -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip()).strip("._-")
    return cleaned[:80] or fallback


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
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8-sig", newline="", delete=False, dir=path.parent, suffix=".tmp"
    ) as handle:
        temporary = Path(handle.name)
        frame.to_csv(handle, index=False)
    os.replace(temporary, path)


def load_ontology(path: Path) -> dict[str, dict[str, Any]]:
    table = pd.read_csv(path, dtype=str, keep_default_na=False)
    required = {
        "trait_id",
        "display_name",
        "allowed_states",
        "allow_multiple",
        "required_visual_context",
        "assessability_rule",
    }
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


def _safe_member_path(name: str) -> PurePosixPath:
    normalized = name.replace("\\", "/")
    path = PurePosixPath(normalized)
    if path.is_absolute() or not path.parts or any(part in {"", ".", ".."} for part in path.parts):
        raise ValueError(f"Unsafe path in uploaded ZIP: {name!r}")
    lowered = "/".join(path.parts).lower()
    if any(token in lowered for token in FORBIDDEN_ARCHIVE_TOKENS):
        raise ValueError("Uploaded ZIP contains a private audit/holdout key. Upload only the public blinded packet.")
    return path


def extract_public_packet_zip(payload: bytes, destination: Path) -> None:
    """Safely extract a public packet ZIP without path traversal or private keys."""
    destination.mkdir(parents=True, exist_ok=True)
    root = destination.resolve()
    try:
        archive = zipfile.ZipFile(io.BytesIO(payload))
    except zipfile.BadZipFile as error:
        raise ValueError("The uploaded file is not a readable ZIP archive") from error

    with archive:
        infos = archive.infolist()
        if not infos:
            raise ValueError("The uploaded ZIP is empty")
        total_size = sum(max(0, info.file_size) for info in infos)
        if total_size > MAX_UNCOMPRESSED_PACKET_BYTES:
            raise ValueError("The uploaded packet expands beyond the 2 GiB safety limit")

        checked: list[tuple[zipfile.ZipInfo, PurePosixPath]] = []
        for info in infos:
            member = _safe_member_path(info.filename)
            if info.flag_bits & 0x1:
                raise ValueError("Encrypted ZIP members are not supported")
            unix_mode = info.external_attr >> 16
            if unix_mode and stat.S_ISLNK(unix_mode):
                raise ValueError("Symbolic links are not allowed in uploaded ZIP packets")
            target = (destination / Path(*member.parts)).resolve()
            try:
                target.relative_to(root)
            except ValueError as error:
                raise ValueError(f"Unsafe path in uploaded ZIP: {info.filename!r}") from error
            checked.append((info, member))

        for info, member in checked:
            target = destination / Path(*member.parts)
            if info.is_dir():
                target.mkdir(parents=True, exist_ok=True)
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            with archive.open(info) as source, target.open("wb") as sink:
                shutil.copyfileobj(source, sink)


def _csv_columns(path: Path) -> set[str] | None:
    try:
        return set(pd.read_csv(path, dtype=str, keep_default_na=False, nrows=0).columns)
    except Exception:
        return None


def discover_public_packet(root: Path) -> PublicPacket:
    """Locate one blinded task table and an optional bundled ontology."""
    csv_files = sorted(path for path in root.rglob("*.csv") if path.is_file())
    if not csv_files:
        raise ValueError("Uploaded ZIP contains no CSV files")

    task_candidates: list[Path] = []
    for path in csv_files:
        columns = _csv_columns(path)
        if columns is None:
            continue
        leaked = FORBIDDEN_PRIVATE_FIELDS.intersection(columns)
        if leaked:
            relative = path.relative_to(root)
            raise ValueError(
                f"Uploaded ZIP is not blinded: {relative} contains private fields {sorted(leaked)}"
            )
        if REQUIRED.issubset(columns):
            task_candidates.append(path)

    if not task_candidates:
        raise ValueError(
            "No app-ready task CSV was found. Expected columns: " + ", ".join(sorted(REQUIRED))
        )

    preferred = [
        path
        for filename in PREFERRED_TASK_FILENAMES
        for path in task_candidates
        if path.name == filename
    ]
    choices = preferred or task_candidates
    if len(choices) != 1:
        names = ", ".join(str(path.relative_to(root)) for path in choices)
        raise ValueError(f"Uploaded ZIP contains multiple possible task CSVs: {names}")
    task_path = choices[0]

    ontology_candidates = [
        path for path in csv_files if path.name == "ch1_trait_ontology.csv" and path != task_path
    ]
    ontology_path = ontology_candidates[0] if len(ontology_candidates) == 1 else None
    return PublicPacket(task_path=task_path, ontology_path=ontology_path)


def validate_resume_csv(payload: bytes, task_ids: set[str]) -> pd.DataFrame:
    try:
        frame = pd.read_csv(io.BytesIO(payload), dtype=str, keep_default_na=False)
    except Exception as error:
        raise ValueError("The resume file is not a readable CSV") from error
    if "task_id" not in frame.columns or frame["task_id"].duplicated().any():
        raise ValueError("Resume CSV must contain one unique row per task_id")
    unknown = sorted(set(frame["task_id"]).difference(task_ids))
    if unknown:
        preview = ", ".join(unknown[:5])
        raise ValueError(f"Resume CSV contains task IDs absent from this packet: {preview}")
    return frame


def _clear_cloud_session(st, session_root: Path) -> None:
    shutil.rmtree(session_root, ignore_errors=True)
    for key in list(st.session_state):
        if key.startswith("trait_audit_") or key == "index":
            del st.session_state[key]
    st.rerun()


def cloud_config(st, args: argparse.Namespace) -> AppConfig:
    st.info(
        "Cloud upload mode: enter an annotator ID and upload the public blinded audit packet ZIP. "
        "Never upload the private selection/holdout key."
    )
    annotator = st.text_input("Annotator ID", value=text(st.session_state.get("trait_audit_annotator", "")))
    st.session_state.trait_audit_annotator = annotator
    packet_upload = st.file_uploader(
        "Public blinded audit packet ZIP",
        type=["zip"],
        help="Use the public packet or packaged audit-app artifact, never the private key directory.",
        key="trait_audit_packet_upload",
    )
    resume_upload = st.file_uploader(
        "Previous trait_audit_responses.csv (optional, to resume)",
        type=["csv"],
        key="trait_audit_resume_upload",
    )

    if not annotator.strip():
        st.warning("Enter an annotator ID to keep response files attributable and separate.")
        st.stop()
    if packet_upload is None:
        st.stop()

    if "trait_audit_session_id" not in st.session_state:
        st.session_state.trait_audit_session_id = uuid.uuid4().hex
    session_root = Path(tempfile.gettempdir()) / "cirsium_trait_audit_cloud" / st.session_state.trait_audit_session_id
    if st.sidebar.button("Reset uploaded audit session"):
        _clear_cloud_session(st, session_root)

    payload = packet_upload.getvalue()
    digest = hashlib.sha256(payload).hexdigest()
    packet_dir = session_root / "packets" / digest
    packet_state_key = f"{digest}:{len(payload)}"
    if st.session_state.get("trait_audit_packet_state_key") != packet_state_key:
        shutil.rmtree(packet_dir, ignore_errors=True)
        try:
            extract_public_packet_zip(payload, packet_dir)
            packet = discover_public_packet(packet_dir)
        except Exception:
            shutil.rmtree(packet_dir, ignore_errors=True)
            raise
        st.session_state.trait_audit_packet_state_key = packet_state_key
        st.session_state.trait_audit_task_path = str(packet.task_path)
        st.session_state.trait_audit_bundled_ontology = str(packet.ontology_path) if packet.ontology_path else ""
        st.session_state.index = max(0, args.start_at - 1)

    task_path = Path(st.session_state.trait_audit_task_path)
    bundled_ontology = text(st.session_state.get("trait_audit_bundled_ontology", ""))
    ontology_path = Path(bundled_ontology) if bundled_ontology else Path(args.ontology).expanduser().resolve()
    out_dir = session_root / "responses" / safe_component(annotator)
    response_path = out_dir / "trait_audit_responses.csv"

    runtime_key = f"{digest}:{annotator.strip()}"
    if st.session_state.get("trait_audit_runtime_key") != runtime_key:
        st.session_state.trait_audit_runtime_key = runtime_key
        st.session_state.index = max(0, args.start_at - 1)

    if resume_upload is not None:
        resume_payload = resume_upload.getvalue()
        resume_key = f"{runtime_key}:{hashlib.sha256(resume_payload).hexdigest()}"
        if st.session_state.get("trait_audit_resume_key") != resume_key:
            tasks = pd.read_csv(task_path, dtype=str, keep_default_na=False)
            resumed = validate_resume_csv(resume_payload, set(tasks["task_id"].map(text)))
            atomic_csv(resumed, response_path)
            st.session_state.trait_audit_resume_key = resume_key
            st.session_state.index = max(0, args.start_at - 1)

    return AppConfig(
        task_path=task_path,
        out_dir=out_dir,
        ontology_path=ontology_path,
        annotator=annotator.strip(),
        start_at=args.start_at,
        cloud_mode=True,
    )


def local_config(args: argparse.Namespace) -> AppConfig:
    if bool(args.tasks) != bool(args.out_dir):
        raise ValueError("Local mode requires both --tasks and --out-dir. Omit both to use Cloud upload mode.")
    return AppConfig(
        task_path=Path(args.tasks).expanduser().resolve(),
        out_dir=Path(args.out_dir).expanduser().resolve(),
        ontology_path=Path(args.ontology).expanduser().resolve(),
        annotator=text(args.annotator),
        start_at=args.start_at,
        cloud_mode=False,
    )


def csv_download_bytes(frame: pd.DataFrame) -> bytes:
    return frame.to_csv(index=False).encode("utf-8-sig")


def main() -> None:
    args = parse_args()
    try:
        import streamlit as st
    except ImportError as error:
        raise SystemExit(
            "Streamlit is required: python -m pip install -r ch1_global/v2/requirements_trait_audit_app.txt"
        ) from error

    st.set_page_config(page_title="Cirsium Ch.1 focused trait audit", layout="wide")
    st.title("Cirsium Ch.1 — blinded focused trait audit")
    st.caption("Score visible evidence only. AI suggestions, species names, localities, and model confidence are hidden.")

    try:
        config = local_config(args) if args.tasks or args.out_dir else cloud_config(st, args)
        task_path = config.task_path
        packet_root = task_path.parent
        tasks = pd.read_csv(task_path, dtype=str, keep_default_na=False)
        if missing := REQUIRED.difference(tasks.columns):
            raise ValueError(f"Task packet missing columns: {sorted(missing)}")
        if tasks.empty or tasks["task_id"].duplicated().any():
            raise ValueError("Task IDs must be nonempty and unique")
        ontology = load_ontology(config.ontology_path)
        unknown = sorted(set(tasks["trait_id"]).difference(ontology))
        if unknown:
            raise ValueError(f"Tasks refer to unknown ontology traits: {unknown}")
        response_path = config.out_dir / "trait_audit_responses.csv"
        status_path = config.out_dir / "trait_audit_status.json"
        data = load_progress(tasks, response_path)
    except Exception as error:
        st.error(str(error))
        st.stop()

    if "index" not in st.session_state:
        st.session_state.index = max(0, min(len(data) - 1, config.start_at - 1))
    index = max(0, min(len(data) - 1, st.session_state.index))
    completed = data.apply(complete, axis=1)
    st.sidebar.metric("Completed", f"{int(completed.sum())} / {len(data)}")
    st.sidebar.progress(float(completed.mean()))
    if config.cloud_mode:
        st.sidebar.warning("Cloud storage is temporary. Download the response CSV regularly and before closing the session.")
    else:
        st.sidebar.caption(f"Responses: {response_path}")
    st.sidebar.download_button(
        "Download current response CSV",
        data=csv_download_bytes(data),
        file_name=f"{safe_component(config.annotator)}_trait_audit_responses.csv",
        mime="text/csv",
    )
    if st.sidebar.button("Jump to next incomplete") and (~completed).any():
        st.session_state.index = int((~completed).to_numpy().nonzero()[0][0])
        st.rerun()

    row = data.iloc[index]
    trait = ontology[text(row["trait_id"])]
    st.subheader(f"Task {index + 1} of {len(data)} — {trait['display_name']}")
    st.caption(
        f"Required view: {trait['required_visual_context']}  |  Assessability: {trait['assessability_rule']}"
    )
    columns = st.columns(3)
    views = (
        ("Head crop", "crop_path"),
        ("Head + peduncle context", "context_crop_path"),
        ("Source image", "source_image"),
    )
    for column, (caption, field) in zip(columns, views):
        with column:
            image = safe_image(row[field], packet_root)
            if image is not None:
                st.image(image, caption=caption, width="stretch")
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
            options = ["", *states]
            current = text(row.get("human_state", ""))
            human_state = st.selectbox(
                trait["display_name"],
                options=options,
                index=options.index(current) if current in options else 0,
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
        status_path.write_text(
            json.dumps(
                {
                    "tasks": str(task_path),
                    "ontology": str(config.ontology_path.resolve()),
                    "annotator": config.annotator,
                    "n_tasks": len(data),
                    "n_complete": int(data.apply(complete, axis=1).sum()),
                    "last_task_id": text(row["task_id"]),
                },
                ensure_ascii=False,
                indent=2,
            ),
            encoding="utf-8",
        )
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
