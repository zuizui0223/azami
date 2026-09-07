import copy
import importlib
import json
from pathlib import Path
import sqlite3
import sys

import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
aggregate=importlib.import_module("analysis.v3.build_observation_measurements")


def head(value=10):
    endpoints={r["endpoint_id"]:{"value":None,"status":"missing"} for r in aggregate.registry()}
    for key,v in [("orientation_image_vertical_angle",value),("corolla_lab_chroma",20),("corolla_lab_lightness",50)]:
        endpoints[key]={"value":v,"status":"usable"}
    flower={"n_pixels":300,"lab_chroma":18,"lab_lightness":45}
    background={"support_status":"available","lab_chroma":8,"lab_lightness":40}
    diagnostics={"head_min_dimension_px":200,"head_laplacian_variance":100,"paired_colour":{"floral_union":flower,"non_head_context":background,"green_non_head_context":background}}
    return {"endpoints":endpoints,"diagnostics":diagnostics,"status":"measured"}


def test_raw_qc_failed_values_retained_but_not_in_eligible_mean():
    a,b=head(10),head(100)
    b["endpoints"]["orientation_image_vertical_angle"]["status"]="bad_qc"
    table,pairs=aggregate.aggregate_heads([a,b])
    row=next(r for r in table if r["endpoint_id"]=="orientation_image_vertical_angle")
    assert row["raw_mean"]==55 and row["eligible_mean"]==10
    assert row["n_heads"]==2 and row["n_eligible_heads"]==1


def test_hue_and_composition_keep_joint_head_sets():
    sample=head()
    for key in aggregate.COMPOSITION:
        sample["endpoints"][key]={"value":.25,"status":"usable"}
    sample["endpoints"][aggregate.COMPOSITION[0]]["value"]=.5
    sample["endpoints"][aggregate.HUE[0]]={"value":0,"status":"usable"}
    ok=aggregate.eligible_endpoints(sample["endpoints"])
    assert not any(ok[k] for k in aggregate.COMPOSITION+aggregate.HUE)


def test_context_pair_uses_same_heads_and_retains_lost_support():
    a,b=head(),head()
    b["diagnostics"]["paired_colour"]["green_non_head_context"]={"support_status":"insufficient_pixels","lab_chroma":1000,"lab_lightness":1000}
    endpoints,pairs=aggregate.aggregate_heads([a,b])
    row=next(r for r in pairs if r["context_kind"]=="green_non_head_context" and r["statistic"]=="lab_chroma")
    assert row["n_eligible_colour_heads"]==2 and row["n_pairs"]==1
    assert row["floral_mean"]==18 and row["context_mean"]==8 and row["contrast_mean"]==10
    assert row["legacy_mean"]==20


def test_missing_quality_does_not_change_covariate_denominator():
    a,b=head(),head()
    del b["diagnostics"]["head_min_dimension_px"]
    table,pairs=aggregate.aggregate_heads([a,b])
    row=next(r for r in table if r["endpoint_id"]=="orientation_image_vertical_angle")
    assert row["eligible_mean"]==10 and row["head_min_dimension_mean_px"] is None


def test_empty_image_is_missing_not_zero():
    table,pairs=aggregate.aggregate_heads([])
    assert len(table)==27 and len(pairs)==4
    assert all(r["raw_mean"] is None and r["eligible_mean"] is None and r["n_heads"]==0 for r in table)


def test_identical_pixels_cannot_choose_favourable_processing_result():
    aggregate.equivalent_pixel_summaries([{"mean":1},{"mean":1}])
    with pytest.raises(ValueError,match="inconsistent"):
        aggregate.equivalent_pixel_summaries([{"mean":1},{"mean":2}])


def fixture(tmp_path):
    roots=[tmp_path/name for name in ("workspace","detection","measurement","dependence")]
    for path in roots:
        path.mkdir()
    paths=[p/name for p,name in zip(roots,["image_workspace.sqlite","detection.sqlite","measurements.sqlite","dependence_groups.sqlite"])]
    with sqlite3.connect(paths[0]) as db:
        db.executescript("CREATE TABLE objects(sha256,decoded_rgb_sha256); CREATE TABLE photo_versions(photo_id,sha256); CREATE TABLE source_links(obs_id,photo_id);")
        db.executemany("INSERT INTO objects VALUES(?,?)",[(s,"pixels"+s) for s in "ABCDE"])
        db.executemany("INSERT INTO photo_versions VALUES(?,?)",[("p1","A"),("alias","A"),("p2","B"),("p3","C"),("p5","D"),("p5","E")])
        db.executemany("INSERT INTO source_links VALUES(?,?)",[("1","p1"),("1","alias"),("1","p2"),("2","p1"),("3","p3"),("4","p4"),("5","p5")])
    jobs=[(f"A{i}","A",10) for i in range(100)]+[("B1","B",0)]
    with sqlite3.connect(paths[1]) as db:
        db.executescript("CREATE TABLE jobs(sha256,status); CREATE TABLE detections(head_id,sha256);")
        db.executemany("INSERT INTO jobs VALUES(?,?)",[("A","detected"),("B","detected"),("C","no_detection"),("D","no_detection"),("E","error")])
        db.executemany("INSERT INTO detections VALUES(?,?)",[(h,s) for h,s,v in jobs])
    with sqlite3.connect(paths[2]) as db:
        db.executescript("CREATE TABLE jobs(head_id TEXT PRIMARY KEY,sha256,status); CREATE TABLE endpoints(head_id,endpoint_id,value,status); CREATE INDEX head_endpoint ON endpoints(head_id); CREATE TABLE details(head_id TEXT PRIMARY KEY,raw_and_diagnostics_json);")
        for h,s,v in jobs:
            data=head(v)
            db.execute("INSERT INTO jobs VALUES(?,?,?)",(h,s,"measured"))
            db.execute("INSERT INTO details VALUES(?,?)",(h,json.dumps({"diagnostics":data["diagnostics"]})))
            db.executemany("INSERT INTO endpoints VALUES(?,?,?,?)",[(h,key,row["value"],row["status"]) for key,row in data["endpoints"].items()])
    with sqlite3.connect(paths[3]) as db:
        db.execute("CREATE TABLE observations(obs_id,component_id,fold)")
        db.executemany("INSERT INTO observations VALUES(?,?,?)",[(str(i),"group"+str(max(2,i)),i%5) for i in range(1,6)])
    hashes=[aggregate.digest(p) for p in paths]
    reports=[{"workspace_database_sha256":hashes[0]},
             {"detection_database_sha256":hashes[1],"job_states":{"detected":2,"no_detection":2,"error":1},"execution_contract":{"source_workspace_sha256":hashes[0]}},
             {"measurement_database_sha256":hashes[2],"job_states":{"measured":101},"execution_contract":{"detector_database_sha256":hashes[1],"feature_specification":{"source_sha256_text_lf":{aggregate.REGISTRY.relative_to(aggregate.ROOT).as_posix():aggregate.text_digest(aggregate.REGISTRY)}}}},
             {"output_database_sha256":hashes[3],"source_workspace_sha256":hashes[0],"counts":{"observations_retained":5,"source_links_retained":7}}]
    for p,name,report in zip(roots,["image_workspace_report.json","detection_report.json","measurement_report.json","dependence_groups_report.json"],reports):
        (p/name).write_text(json.dumps(report))
    return roots,paths,reports


def test_full_pipeline_preserves_unmeasured_and_ambiguous_sources(tmp_path):
    roots,paths,reports=fixture(tmp_path)
    result=aggregate.build(*roots,tmp_path/"output")
    assert result["counts"]["source_observations_retained"]==5
    assert result["counts"]["source_photo_links_referenced"]==7
    assert result["counts"]["photo_ids_with_unresolved_pixel_version_choice"]==1
    with sqlite3.connect(tmp_path/"output/observation_measurements.sqlite") as db:
        assert db.execute("SELECT COUNT(*) FROM observation_endpoint_inventory").fetchone()[0]==135
        rows=dict(db.execute("SELECT obs_id,value FROM observation_endpoint_inventory WHERE endpoint_id='orientation_image_vertical_angle'"))
        assert rows=={"1":5,"2":10,"3":None,"4":None,"5":None}
        states=dict(db.execute("SELECT obs_id,measurement_support FROM observation_endpoint_inventory WHERE endpoint_id='orientation_image_vertical_angle'"))
        assert states["3"]=="detector_negative_images_only"
        assert states["4"]==states["5"]=="no_selected_cached_image"
        assert db.execute("SELECT n_shared_pixel_images FROM observation_support WHERE obs_id='1'").fetchone()[0]==1
    assert result["endpoint_coverage"]["orientation_image_vertical_angle"]=={"observations":2,"observation_image_links":3,"unique_pixel_images":2}
    for path,sha in zip(paths,result["execution_contract"]["input_database_sha256"].values()):
        assert aggregate.digest(path)==sha


def test_incomplete_measurement_schedule_stops(tmp_path):
    roots,paths,reports=fixture(tmp_path)
    reports[2]["job_states"]["pending"]=1
    (roots[2]/"measurement_report.json").write_text(json.dumps(reports[2]))
    with pytest.raises(ValueError,match="Complete"):
        aggregate.build(*roots,tmp_path/"output")
    assert not (tmp_path/"output").exists()
