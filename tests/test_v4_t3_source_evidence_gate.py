"""T3 source audit tests: anatomy is not the same as evolutionary origins."""
import csv
import json
from pathlib import Path

import pytest

from analysis.v4 import t3_source_evidence_gate as gate

REAL = Path(__file__).resolve().parents[1] / "analysis" / "v4" / "results" / "20261023_T3_SENECIONEAE_POSTURE_SOURCE_CANDIDATES.csv"

def write_rows(path, rows):
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=gate.SCHEMA)
        w.writeheader()
        w.writerows(rows)

def test_real_ten_botanical_records_not_mk_transitions():
    rows = gate.read_evidence(REAL)
    report = gate.summarize(rows)
    assert len(rows) == 10
    assert report["nodding_species_descriptions_not_independent_gains"] == 6
    assert report["erect_species_descriptions_not_ancestral_states"] == 1
    assert report["qualified_uncertain_species"] == 1
    assert report["ontogenetically_variable_species"] == 2
    assert report["qualified_age_dependent_nonbinary_species"] == 1
    assert report["source_taxon_concept_contests"] == 2
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
    assert rows["Xenophyllum staffordiae"]["life_stage_scope"] == "age_dependent_tendency"
    assert rows["Senecio nigrocinctus"]["source_posture"] == "V"
    assert rows["Senecio nigrocinctus"]["life_stage_scope"] == "developmentally_changes"
    assert rows["Senecio eliseae"]["taxonomic_status"] == "contested_2018_2022"
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

def test_conflicting_taxonomic_concepts_not_promoted_even_if_other_evidence_flags_yes(tmp_path):
    rows = gate.read_evidence(REAL)
    indexed = {r["taxon"]: r for r in rows}
    eliseae = indexed["Senecio eliseae"]
    eliseae["life_stage_scope"] = "anthesis_independently_verified"
    eliseae["taxon_concept_independently_reviewed"] = "yes"
    eliseae["blind_double_annotation_completed"] = "yes"
    eliseae["same_tip_rooted_nuclear_tree_verified"] = "yes"
    eliseae["strict_T3_admissible"] = "yes"
    fp = tmp_path / "taxon_synonym_collision.csv"
    write_rows(fp, rows)
    with pytest.raises(ValueError, match="Contested author taxon concept"):
        gate.read_evidence(fp)

def test_age_tendency_is_not_misclassified_as_fixed_nodding(tmp_path):
    rows = gate.read_evidence(REAL)
    staffordiae = next(r for r in rows if r["taxon"] == "Xenophyllum staffordiae")
    staffordiae["source_posture"] = "D"
    fp = tmp_path / "wrong_age_posture.csv"
    write_rows(fp, rows)
    with pytest.raises(ValueError, match="Age-related tendency"):
        gate.read_evidence(fp)

def test_cremanthodium_ligularia_competing_2016_2026_concepts_stay_separate(tmp_path):
    rows = gate.read_evidence(REAL)
    target = next(r for r in rows if r["taxon"] == "Cremanthodium spathulifolium")
    assert target["source_posture"] == "D"
    assert target["taxonomic_status"] == "contested_2016_2026"
    assert "Ligularia lagotis" in target["coding_limit"]
    target["life_stage_scope"] = "anthesis_independently_verified"
    target["taxon_concept_independently_reviewed"] = "yes"
    target["blind_double_annotation_completed"] = "yes"
    target["same_tip_rooted_nuclear_tree_verified"] = "yes"
    target["strict_T3_admissible"] = "yes"
    outfile = tmp_path / "taxonomy_conflict.csv"
    write_rows(outfile, rows)
    with pytest.raises(ValueError, match="Contested author taxon concept"):
        gate.read_evidence(outfile)

def test_external_phylogenomics_public_source_selection_does_not_impute_tips():
    manifest = json.loads((
        REAL.parent / "20261023_T3_PUBLIC_SENECIONEAE_NUCLEAR_SOURCE_SELECTION.json"
    ).read_text(encoding="utf-8"))
    a = manifest["source_evidence"]["Ren2021"]
    b = manifest["source_evidence"]["Ren2024"]
    assert a["reported_Ligularia_Cremanthodium_Parasenecio_samples"] == 27
    assert b["reported_Ligularia_Cremanthodium_Parasenecio_representatives"] == 5
    assert b["retained_nuclear_loci"] == 727
    assert b["published_dataset_zip_sha256"] == (
        "962495cf8131fadd3a1cdf1ea88204432f2aaa81c9e6c0a97435756f0ad88139"
    )
    assert b["published_dataset_zip_size_bytes"] == 1583779312
    assert b["published_tree_member"].endswith("2_astral_analysis/astral.tre")
    assert not a["sample_names_and_exact_T3_tip_overlap_verified"]
    assert not b["sample_names_and_exact_T3_tip_overlap_verified"]
    assert not b["direct_tree_ingested"]
    assert not manifest["selection"]["joint_Cardueae_Senecioneae_independently_rooted_species_tree_available"]
    assert not manifest["selection"]["independent_nodding_gains_inferred"]
    assert not manifest["selection"]["rain_by_transition_model_fitted"]
