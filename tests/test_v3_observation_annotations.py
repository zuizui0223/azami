"""Synthetic calendar/privacy/provenance tests; no ecological model fitting."""
import importlib
import json
import math
from pathlib import Path
import sqlite3
import sys

import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
prepare=importlib.import_module("analysis.v3.prepare_observation_annotations")


def fields(**kwargs):
    result=prepare.api_fields({"id":1,"observed_on":"2024-02-29","obscured":False,"geojson":{"coordinates":[20,40]}})
    result.update(kwargs)
    return result


def test_leap_year_calendar_and_hemisphere_interactions():
    north=prepare.annotate(fields(),[],"api")
    south=prepare.annotate(fields(latitude="-40"),[],"api")
    assert north["doy"]==60 and north["days_in_year"]==366
    assert north["sin_doy"]==pytest.approx(math.sin(2*math.pi*59/366))
    assert north["sin_doy"]==south["sin_doy"]
    assert north["south_sin"]==0 and south["south_sin"]==south["sin_doy"]
    assert north["hemisphere"]=="north" and south["hemisphere"]=="south"


@pytest.mark.parametrize("raw",["2023-02-29","2024-02","20240229","2024-02-29T12:00:00Z","not a date"])
def test_coarse_invalid_or_timestamp_dates_not_silently_repaired(raw):
    result=prepare.annotate(fields(observed_on=raw),[],"api")
    assert result["date_status"]=="invalid_or_not_exact_day"
    assert result["doy"] is None and result["sin_doy"] is None


@pytest.mark.parametrize("change",[{"geoprivacy":"private"},{"geoprivacy":"obscured"},{"obscured":"true"}])
def test_restricted_coordinates_do_not_produce_analysis_coordinates(change):
    result=prepare.annotate(fields(**change),[],"api")
    assert result["coordinate_status"]=="restricted"
    assert result["analysis_latitude"] is None and result["south_indicator"] is None
    assert result["date_status"]=="exact_day"  # Calendar can remain useful separately.


def test_unknown_privacy_is_not_false_and_private_fields_never_read():
    data=prepare.api_fields({"id":1,"private_location":"45,120","private_geojson":{"coordinates":[120,45]}})
    assert data["latitude"]=="" and data["obscured"]=="unknown"
    assert prepare.annotate(data,[],"api")["coordinate_status"]=="privacy_unknown"


@pytest.mark.parametrize("lat,lon",[("91","0"),("0","181"),("NaN","0"),("Inf","0")])
def test_invalid_coordinates_remain_unavailable(lat,lon):
    result=prepare.annotate(fields(latitude=lat,longitude=lon),[],"api")
    assert result["analysis_latitude"] is None


def test_api_missing_boolean_not_filled_from_flattened_false():
    api=fields(captive="unknown")
    metadata=fields(captive="false")
    kind,n,chosen,conflicts=prepare.consensus([("metadata",metadata),("api",api)])
    assert kind=="api" and chosen["captive"]=="unknown" and n==1
    assert prepare.annotate(chosen,conflicts,kind)["captive_state"]=="unknown"


def test_conflicts_block_only_affected_fields_without_last_row_choice():
    a,b=fields(),fields(observed_on="2024-03-01")
    first=prepare.consensus([("api",a),("api",b)])
    second=prepare.consensus([("api",b),("api",a)])
    assert first==second
    kind,n,data,conflicts=first
    result=prepare.annotate(data,conflicts,kind)
    assert result["date_status"]=="conflicting_source_versions"
    assert result["analysis_latitude"]==40


def test_coordinate_conflict_blocks_location_and_hemisphere():
    kind,n,data,conflicts=prepare.consensus([("api",fields()),("api",fields(latitude="-40"))])
    result=prepare.annotate(data,conflicts,kind)
    assert result["coordinate_status"]=="conflicting_source_versions"
    assert result["analysis_latitude"] is None and result["hemisphere"]=="unknown"


@pytest.mark.parametrize("value,state",[("","missing"),("-1","invalid"),("nan","invalid"),("0","reported_zero"),("15","reported_positive")])
def test_precision_states_are_not_spatial_resolution_acceptance(value,state):
    result=prepare.annotate(fields(positional_accuracy=value),[],"api")
    assert result["position_accuracy_status"]==state
    assert result["native_status"]=="not_assessed_full_source"


def test_complete_archive_to_annotation_pipeline_preserves_sources(tmp_path):
    import test_v3_source_reconciliation as recovery
    archives,merged,manifest,contract=recovery.fixture(tmp_path,compressed=True)
    rec=tmp_path/"reconciliation"
    recovery.audit.reconcile(archives,merged,rec,manifest,contract)
    groups=tmp_path/"groups"
    groups.mkdir()
    path=groups/"dependence_groups.sqlite"
    with sqlite3.connect(path) as db:
        db.execute("CREATE TABLE observations(obs_id TEXT PRIMARY KEY,component_id TEXT,fold INTEGER)")
        db.executemany("INSERT INTO observations VALUES(?,?,?)",[(str(i),"group"+str(i),i%5) for i in range(1,7)])
    (groups/"dependence_groups_report.json").write_text(json.dumps({"output_database_sha256":prepare.digest(path)}))
    out=tmp_path/"annotations"
    result=prepare.build(archives,rec,groups,out,manifest)
    assert result["counts"]["observations_retained"]==result["counts"]["annotation_rows"]==6
    assert result["counts"]["source_metadata_rows_retained"]==3
    assert result["counts"]["source_api_rows_retained"]==5
    assert result["states"]["preferred_source_kind"]=={"api":5,"metadata":1}
    assert not any(result["integrity_checks"].values())
    with sqlite3.connect(out/"observation_annotations.sqlite") as db:
        assert db.execute("SELECT COUNT(*) FROM annotations WHERE native_status='not_assessed_full_source'").fetchone()[0]==6
    assert prepare.digest(path)==json.loads((groups/"dependence_groups_report.json").read_text())["output_database_sha256"]
    with pytest.raises(ValueError,match="Preserve"):
        prepare.build(archives,rec,groups,out,manifest)
