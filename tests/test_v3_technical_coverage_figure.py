"""Public aggregate figure checks; no raw images or identifiers needed."""
import csv
import json
from pathlib import Path

import numpy as np
import pytest

from analysis.v3 import plot_technical_coverage as figure
from analysis.v3.workflow import digest

ROOT = Path(__file__).resolve().parents[1]
SUMMARY = ROOT / "analysis_outputs/v3/technical_sensitivity_summary_20260907.csv"
VERIFY = ROOT / "reproducibility/v3_technical_summary_verification_20260907.json"


def test_public_aggregate_counts_and_full_grid():
    endpoints, n, counts, loss = figure.load_data(SUMMARY, VERIFY)
    assert n == 2853 and len(endpoints) == 27 and loss.shape == (27, 13)
    assert counts.tolist() == [995]+[1250]*9+[1136]*4+[160]*9+[53]*4
    assert np.all(loss[-4:, 5] == 100)
    labels = dict(zip(endpoints, figure.LABELS))
    assert labels["bract_projection_roughness"] == "Projection roughness"
    assert labels["bract_spread_fraction"] == "Bract spread fraction"
    assert labels["bract_projection_maximum"] == "Projection maximum"


def modified_summary(tmp_path, change):
    with SUMMARY.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames
        rows = list(reader)
    change(rows)
    target = tmp_path / "summary.csv"
    with target.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    receipt = json.loads(VERIFY.read_text(encoding="utf-8"))
    receipt["summary_sha256"] = digest(target)
    verification = tmp_path / "receipt.json"
    verification.write_text(json.dumps(receipt), encoding="utf-8")
    return target, verification


def test_missing_loss_is_not_zero(tmp_path):
    def change(rows):
        next(r for r in rows if r["condition"] == "bbox_left_5pct" and r["exposure_stratum"] == "all_cached")["component_weighted_loss_fraction"] = ""
    args = modified_summary(tmp_path, change)
    assert np.isnan(figure.load_data(*args)[3][0, 0])


def test_changed_hash_rejected(tmp_path):
    receipt = json.loads(VERIFY.read_text(encoding="utf-8"))
    receipt["summary_sha256"] = "0"*64
    target = tmp_path / "receipt.json"
    target.write_text(json.dumps(receipt), encoding="utf-8")
    with pytest.raises(ValueError, match="hash differs"):
        figure.load_data(SUMMARY, target)


def test_duplicate_grid_rejected(tmp_path):
    def change(rows):
        rows[3] = rows[0].copy()
    with pytest.raises(ValueError, match="Complete, unique"):
        figure.load_data(*modified_summary(tmp_path, change))


def test_smoke_export_and_no_overwrite(tmp_path):
    out = tmp_path / "figure"
    report = figure.render(SUMMARY, VERIFY, out)
    assert report["ecological_models_executed"] is False
    assert (out / "technical_coverage.pdf").stat().st_size > 1000
    with pytest.raises(FileExistsError):
        figure.render(SUMMARY, VERIFY, out)
