"""Shared validation helpers for submission-facing tabular audits."""

from __future__ import annotations

from typing import Iterable

import pandas as pd


def require_columns(frame: pd.DataFrame, columns: Iterable[str], label: str) -> None:
    """Require named columns and report all missing fields together."""
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")


def assert_unique(frame: pd.DataFrame, keys: list[str], label: str) -> None:
    """Fail with representative duplicate keys instead of silently de-duplicating."""
    duplicated = frame.duplicated(keys, keep=False)
    if duplicated.any():
        examples = frame.loc[duplicated, keys].head(10).to_dict("records")
        raise ValueError(f"{label} has duplicate keys {keys}: {examples}")


def require_complete_text(frame: pd.DataFrame, column: str, label: str | None = None) -> None:
    """Require a non-empty text value in every row of a column."""
    require_columns(frame, [column], label or "table")
    missing = frame[column].isna() | frame[column].astype(str).str.strip().eq("")
    if missing.any():
        name = label or column
        raise ValueError(f"{int(missing.sum())} rows lack {name}")
