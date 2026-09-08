"""Offline source/calendar/dependence linkage tests; no traits or network."""
import csv
import json
from pathlib import Path
import sqlite3

import pytest

from analysis.v3 import enriched_source_cohort as enrich
from analysis.v3.prepare_observation_annotations import annotate, api_fields


def make_fixture(tmp_path, prior_components=True):
    annotation_path = tmp_path / "observation_annotations.sqlite"
    reconciliation_path = tmp_path / "source_reconciliation.sqlite"
    cohort_path = tmp_path / "ecological_source_cohort.csv"
    inputs = {}
    for obs, lat in (("1", 40), ("2", 40), ("3", -40), ("4", 0)):
        fields = api_fields({"id": int(obs), "taxon": {"id": 100, "name": "Cirsium fixture", "rank": "species"},
                             "observed_on": "2024-02-29", "obscured": False, "captive": False,
                             "geojson": {"coordinates": [20, lat]}, "positional_accuracy": 12,
                             "user": {"id": 731}, "quality_grade": "research", "annotations": []})
        derived = annotate(fields, [], "api")
        inputs[obs] = {"obs_id": obs, "preferred_source_kind": "api", "selected_fields_json": json.dumps(fields),
                       "conflicted_fields_json": "[]", **derived}
        if prior_components:
            inputs[obs]["component_id"] = "prior_a" if obs in {"1", "4"} else "prior_" + obs
    with sqlite3.connect(annotation_path) as db:
        numeric = {"analysis_latitude", "analysis_longitude", "position_accuracy_m", "observed_year", "doy", "days_in_year",
                   "sin_doy", "cos_doy", "south_indicator", "south_sin", "south_cos"}
        columns = list(inputs["1"])
        db.execute("CREATE TABLE annotations(" + ",".join(name + (" REAL" if name in numeric else " TEXT") for name in columns) + ")")
        db.executemany("INSERT INTO annotations VALUES(" + ",".join("?" for _ in columns) + ")", [[row[name] for name in columns] for row in inputs.values()])
    with sqlite3.connect(reconciliation_path) as db:
        db.execute("CREATE TABLE source_observations(obs_id TEXT PRIMARY KEY)")
        db.executemany("INSERT INTO source_observations VALUES(?)", [(obs,) for obs in inputs])
        db.execute("CREATE TABLE source_links(obs_id TEXT,photo_id TEXT,PRIMARY KEY(obs_id,photo_id))")
        db.executemany("INSERT INTO source_links VALUES(?,?)", [("1", "11"), ("2", "11"), ("2", "22"), ("3", "22"), ("4", "44")])
    rows = []
    for obs in ("1", "3", "4"):
        a = inputs[obs]
        rows.append({"obs_id": obs, "accepted_key": "100", "accepted_name": "Cirsium fixture",
                     "source_taxon_name": a["source_taxon_name"], "source_taxon_rank": "species",
                     "analysis_latitude": a["analysis_latitude"], "analysis_longitude": a["analysis_longitude"],
                     "observation_month": 2, "position_accuracy_status": a["position_accuracy_status"],
                     "position_accuracy_m": a["position_accuracy_m"], "quarter_degree_cell": "fixture",
                     "native_range_status": "native"})
    with cohort_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=enrich.COHORT_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    return cohort_path, annotation_path, reconciliation_path


def pins(paths):
    return {"expected_" + role + "_sha256": enrich.digest(path) for role, path in zip(("cohort", "annotations", "reconciliation"), paths)}


def run(paths, out):
    return enrich.build(*paths, out, **pins(paths))


def read_rows(out):
    with (out / "enriched_ecological_source_cohort_private.csv").open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def test_real_source_probe_cannot_claim_a_native_cohort(tmp_path):
    paths = make_fixture(tmp_path)
    report = enrich.probe(paths[1], paths[2], expected_annotations_sha256=enrich.digest(paths[1]),
                          expected_reconciliation_sha256=enrich.digest(paths[2]))
    assert report["dependence"]["full_source_observations"] == 4
    assert report["dependence"]["known_dependence_components"] == 1
    assert report["native_cohort_verified"] is False
    assert report["ecological_fitting_authorized"] is False
    assert "obs_id" not in json.dumps(report)


def test_enrichment_preserves_membership_calendar_source_metadata_and_prior_links(tmp_path):
    paths = make_fixture(tmp_path)
    original_hashes = pins(paths)
    out = tmp_path / "enriched"
    report = run(paths, out)
    rows = read_rows(out)
    assert {row["obs_id"] for row in rows} == {"1", "3", "4"}
    assert len({row["dependence_component_id"] for row in rows}) == 1
    assert {row["dependence_component_source_observations"] for row in rows} == {"4"}
    assert {row["doy"] for row in rows} == {"60"}
    assert {row["days_in_year"] for row in rows} == {"366"}
    south = next(row for row in rows if row["obs_id"] == "3")
    north = next(row for row in rows if row["obs_id"] == "1")
    assert south["south_sin"] == south["sin_doy"]
    assert float(north["south_sin"]) == 0
    assert south["source_user_id"] == "731" and south["source_quality_grade"] == "research"
    assert all(row["image_quality_covariates_status"] == "not_joined_requires_endpoint_matched_measurements" for row in rows)
    assert report["cohort_rows"] == 3 and report["source_rows_deleted"] == 0
    assert report["dependence"]["full_source_observations"] == 4
    assert report["dependence"]["reconciled_source_photo_links"] == 5
    assert report["calendar"]["hemisphere_counts"] == {"north": 1, "south": 1, "equatorial": 1}
    assert report["output_csv_sha256"] == enrich.digest(out / "enriched_ecological_source_cohort_private.csv")
    assert report["ecological_fitting_authorized"] is False and report["trait_files_read"] == 0
    assert pins(paths) == original_hashes
    public_text = json.dumps(report)
    assert all(token not in public_text for token in ('"obs_id"', 'Cirsium fixture', '"source_user_id"', '"analysis_latitude"', str(tmp_path)))


def test_absent_prior_components_use_full_source_bridges_before_cohort_subset(tmp_path):
    paths = make_fixture(tmp_path, prior_components=False)
    out = tmp_path / "enriched"
    report = run(paths, out)
    by_id = {row["obs_id"]: row for row in read_rows(out)}
    # Observation 2 is not in the cohort, but links observations 1 and 3.
    assert by_id["1"]["dependence_component_id"] == by_id["3"]["dependence_component_id"]
    assert by_id["4"]["dependence_component_id"] != by_id["1"]["dependence_component_id"]
    assert report["dependence"]["prior_annotation_components_preserved"] is False


@pytest.mark.parametrize("role", ["cohort", "annotations", "reconciliation"])
def test_wrong_input_pin_fails_before_output(tmp_path, role):
    paths = make_fixture(tmp_path)
    expected = pins(paths)
    expected["expected_" + role + "_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="identity differs"):
        enrich.build(*paths, tmp_path / "out", **expected)
    assert not (tmp_path / "out").exists()


def test_sibling_receipts_bind_source_chain(tmp_path):
    paths = make_fixture(tmp_path)
    hashes = {role: enrich.digest(path) for role, path in zip(("cohort", "annotations", "reconciliation"), paths)}
    (tmp_path / "observation_annotations_report.json").write_text(json.dumps({"output_database_sha256": hashes["annotations"],
        "execution_contract": {"source_reconciliation_sha256": hashes["reconciliation"]}}))
    (tmp_path / "source_reconciliation_report.json").write_text(json.dumps({"ledger_sha256": hashes["reconciliation"]}))
    (tmp_path / "ecological_source_cohort_report.json").write_text(json.dumps({"cohort_sha256": hashes["cohort"]}))
    report = run(paths, tmp_path / "good")
    assert all(row["present"] for row in report["checked_upstream_receipts"].values())
    (tmp_path / "observation_annotations_report.json").write_text(json.dumps({"output_database_sha256": hashes["annotations"],
        "execution_contract": {"source_reconciliation_sha256": "0" * 64}}))
    with pytest.raises(ValueError, match="different reconciliation"):
        run(paths, tmp_path / "wrong")


@pytest.mark.parametrize("field,value,match", [("doy", 61, "calendar"), ("analysis_latitude", 41, "coordinates"),
                                             ("source_taxon_name", "changed", "state"), ("captive_state", "unknown", "state")])
def test_saved_annotation_mismatch_is_not_silently_repaired(tmp_path, field, value, match):
    paths = make_fixture(tmp_path)
    with sqlite3.connect(paths[1]) as db:
        db.execute("UPDATE annotations SET " + field + "=? WHERE obs_id='1'", (value,))
    with pytest.raises(ValueError, match=match):
        run(paths, tmp_path / "out")
    assert not (tmp_path / "out" / "enriched_source_cohort_report.json").exists()


def test_conflicted_date_blocks_existing_eligible_row(tmp_path):
    paths = make_fixture(tmp_path)
    with sqlite3.connect(paths[1]) as db:
        db.execute("UPDATE annotations SET conflicted_fields_json='[\"observed_on\"]' WHERE obs_id='1'")
    with pytest.raises(ValueError, match="conflicts with eligible"):
        run(paths, tmp_path / "out")
    assert json.loads((tmp_path / "out" / "incomplete_run.json").read_text())["status"] == "INCOMPLETE_DO_NOT_USE"


def test_missing_annotation_even_outside_cohort_blocks_full_membership(tmp_path):
    paths = make_fixture(tmp_path)
    with sqlite3.connect(paths[1]) as db:
        db.execute("DELETE FROM annotations WHERE obs_id='2'")
    with pytest.raises(ValueError, match="membership differ"):
        run(paths, tmp_path / "out")


def test_trait_column_not_admitted_and_existing_output_not_overwritten(tmp_path):
    paths = make_fixture(tmp_path)
    run(paths, tmp_path / "out")
    with pytest.raises(ValueError, match="preserve every previous"):
        run(paths, tmp_path / "out")
    with paths[0].open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    with paths[0].open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[*enrich.COHORT_FIELDS, "corolla_lab_chroma"])
        writer.writeheader()
        writer.writerows({**row, "corolla_lab_chroma": 10} for row in rows)
    with pytest.raises(ValueError, match="do not admit trait"):
        run(paths, tmp_path / "blocked")


def test_component_ids_stable_under_cohort_reordering(tmp_path):
    paths = make_fixture(tmp_path)
    run(paths, tmp_path / "first")
    with paths[0].open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    with paths[0].open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=enrich.COHORT_FIELDS)
        writer.writeheader()
        writer.writerows(reversed(rows))
    run(paths, tmp_path / "second")
    assert read_rows(tmp_path / "first") == read_rows(tmp_path / "second")
