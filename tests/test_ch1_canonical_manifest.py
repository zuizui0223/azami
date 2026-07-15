from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_canonical_manifest_paths_exist() -> None:
    pipeline = ROOT / "analysis" / "ch1" / "pipeline.json"
    payload = json.loads(pipeline.read_text(encoding="utf-8"))
    assert payload["schema_version"] == 1
    assert payload["policy"]["move_numbered_scripts"] is False
    for stage in payload["stages"].values():
        script = stage.get("script")
        if script:
            assert (ROOT / script).is_file(), script


def test_manuscript_writing_controls_exist() -> None:
    required = (
        "manuscript/MANUSCRIPT_OUTLINE.md",
        "manuscript/FIGURE_TABLE_MAP.md",
        "manuscript/FINAL_MANUSCRIPT_STRATEGY.md",
        "manuscript/final_claims.json",
        "manuscript/RUNBOOK.md",
    )
    for relative in required:
        assert (ROOT / relative).is_file(), relative
