"""T3 source audit tests: anatomy is not the same as evolutionary origins."""
import csv
from pathlib import Path

import pytest

from analysis.v4 import t3_source_evidence_gate as gate

REAL = Path(__file__).resolve().parents[1] / "analysis" / "v4" / "results" / "20261023_T3_SENECIONEAE_POSTURE_SOURCE_CANDIDATES.csv"

def write_rows(path, rows):
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=gate.SCHEMA)
        w.writeheader()
        w.writerows(rows)

def test_real_nine_independent_botanical_taxon_descriptions_not_mk_transitions():
    rows = gate.read_evidence(REAL)
    report = gate.summarize(rows)
    assert len(rows) == 9
    assert report["nodding_species_descriptions_not_independent_gains"] == 6
    assert report["erect_species_descriptions_not_ancestral_states"] == 1
    assert report["qualified_uncertain_species"] == 1
    assert report["ontogenetically_variable_species"] == 1
    assert set(report["source_genus_counts"]) == {
        "Cremanthodium", "Ligularia", "Senecio", "Xenophyllum"
    }
    assert report["strict_reviewed_same_tip_phenotypes"] == 0
    assert not report["historical_repeated_nodding_origins_estimated"]
    assert not report["rain_adaptation_supported_by_this_screen"]

def test_developmental_posture_change_and_rather_nodding_are_not_D():
    rows = {r["taxon"]: r for r in gate.read_evidence(REAL)}
    assert rows["Senecio humbertii"]["source_posture"] == "V"
    assert rows["Senecio humbertii"]["life_stage_scope"] == "developmentally_changes"
    assert rows["Xenophyllum staffordiae"]["source_posture"] == "?"
    assert rows["Xenophyllum esquilachense"]["source_posture"] == "U"

def test_missing_review_and_tree_prohibits_promotion_to_T3(tmp_path):
    rows = gate.read_evidence(REAL)
    rows[0]["strict_T3_admissible"] = "yes"
    fp = tmp_path / "incorrect_promotion.csv"
    write_rows(fp, rows)
    with pytest.raises(ValueError, match="Prematurely promoted"):
        gate.read_evidence(fp)

def test_ontogenetic_change_requires_explicit_stage(tmp_path):
    rows = gate.read_evidence(REAL)
    rows[-1]["life_stage_scope"] = "stage_unspecified"
    fp = tmp_path / "bad_development.csv"
    write_rows(fp, rows)
    with pytest.raises(ValueError, match="Developmental posture"):
        gate.read_evidence(fp)

def test_duplicate_taxon_cannot_be_counted_twice(tmp_path):
    rows = gate.read_evidence(REAL)
    rows.append(dict(rows[0], record_id="new_run"))
    fp = tmp_path / "duplicate.csv"
    write_rows(fp, rows)
    with pytest.raises(ValueError, match="Duplicated source species"):
        gate.read_evidence(fp)
