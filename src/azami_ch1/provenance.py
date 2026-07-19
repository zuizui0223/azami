"""Shared provenance helpers for submission-facing Chapter 1 utilities."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping


def utc_now_iso() -> str:
    """Return an explicit UTC ISO-8601 timestamp."""
    return datetime.now(timezone.utc).isoformat()


def write_json(
    path: str | Path,
    payload: Mapping[str, Any],
    *,
    include_generated_utc: bool = False,
) -> Path:
    """Write deterministic, UTF-8 JSON with a final newline.

    ``generated_utc`` is inserted only when requested and not already supplied.
    Parent directories are created automatically.
    """
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    output = dict(payload)
    if include_generated_utc:
        output.setdefault("generated_utc", utc_now_iso())
    destination.write_text(
        json.dumps(output, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return destination
