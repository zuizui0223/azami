from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "ch1_global" / "v2" / "ANALYSIS_MANIFEST.tsv"


def load_rows() -> list[dict[str, str]]:
    with MANIFEST.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_manifest_has_only_live_repository_paths() -> None:
    rows = load_rows()
    assert rows
    for row in rows:
        relative = row["file"]
        assert relative, row
        assert (ROOT / relative).exists(), relative


def test_manifest_contains_current_reproducibility_entrypoints() -> None:
    rows = load_rows()
    by_stage = {row["stage"]: row for row in rows}
    required = {
        "full27_atlas": "analysis/run_geb_v2_full27_environment_atlas.py",
        "full27_validation": "analysis/validate_geb_v2_full27_environment_atlas.py",
        "figure_rebuild": ".github/workflows/rebuild-frozen-figures.yml",
        "integrity_ci": ".github/workflows/reproducibility-integrity.yml",
        "source_artifact_catalog": "reproducibility/actions_artifact_catalog.json",
        "frozen_outputs": "analysis_outputs/README.md",
    }
    for stage, path in required.items():
        assert by_stage[stage]["file"] == path


def test_removed_submission_workflows_are_not_marked_active() -> None:
    text = MANIFEST.read_text(encoding="utf-8")
    stale = (
        ".github/workflows/ch1-lightweight-spde-inla.yml",
        ".github/workflows/ch1-submission-ci.yml",
        ".github/workflows/ch1-build-submission-bundle.yml",
    )
    for path in stale:
        assert path not in text
