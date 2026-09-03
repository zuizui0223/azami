from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "ch1_global" / "v2" / "ANALYSIS_MANIFEST.tsv"
CATALOG = ROOT / "reproducibility" / "actions_artifact_catalog.json"


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
        "canonical_rebuild_runner": "analysis/rebuild_frozen_analysis.py",
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


def test_artifact_catalog_retains_reconstruction_checkpoints() -> None:
    payload = json.loads(CATALOG.read_text(encoding="utf-8"))
    by_id = {int(row["artifact_id"]): row for row in payload["artifacts"]}
    required_ids = {
        8076736948,
        8269246732,
        8227254443,
        9612943217,
        9632715852,
        8521924881,
        8521926057,
    }
    assert required_ids.issubset(by_id)
    assert by_id[8076736948]["model_weight_sha256"] == (
        "4078e0510532852681b65ee529cd82237b649ec99b17c4ca5f1da460a62d2bed"
    )
    for artifact_id in (9612943217, 9632715852):
        row = by_id[artifact_id]
        assert row["github_digest"].removeprefix("sha256:") == row["local_archive_sha256"]
    historical = by_id[8227254443]
    assert historical["github_digest"] == (
        "sha256:499061e7a49f9455cf8c367fe26e313b7e0e33b2280d2354717e61a90ea8c6bc"
    )
    assert historical["tree_files"]["historical_trees/gbotb_lcvp_scenario2_randomized.trees"] == (
        "sha256:82655f79297e44a6630a599d8b0a1dc6f85e792812b8b01f2234e567a478e3af"
    )
    assert payload["historical_recovery"]["immutable_tag"] == "azami-ch1-v2-2026-08-27"
