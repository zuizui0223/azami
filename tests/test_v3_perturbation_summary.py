import importlib
import json
import math
from pathlib import Path
import sqlite3
import sys

import pandas as pd
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
summary=importlib.import_module("analysis.v3.summarize_perturbations")


def frame(rows):
    result=pd.DataFrame(rows,columns=["head_id","sha256","component_id","baseline_value","perturbed_value","baseline_status","perturbed_status"])
    result["development_exposed"]=0
    return summary.scalar_frame(result)


def test_equal_image_then_equal_component_weight():
    rows=[(f"a{i}","imageA","groupA",0,10,"usable","usable") for i in range(100)]
    rows += [("b","imageB","groupA",0,0,"usable","usable"),("c","imageC","groupB",0,1,"usable","usable")]
    result=summary.summarize_frame(frame(rows))
    assert result["mean_absolute_change"]==3  # ((10+0)/2 + 1)/2, not 1001/102.
    assert result["paired_images"]==3 and result["paired_components"]==2
    assert result["spearman_paired_component_means"] is None


def test_attrition_includes_groups_without_survivors():
    data=frame([("a","A","a",0,10,"usable","usable"),("b","B","b",1,None,"usable","low_quality")])
    result=summary.summarize_frame(data)
    assert result["mean_absolute_change"]==10
    assert result["component_weighted_loss_fraction"]==.5
    assert result["lost_usable_heads"]==1 and result["paired_components"]==1


def test_all_missing_is_unavailable_not_zero():
    result=summary.summarize_frame(frame([("a","A","a",None,None,"missing","missing")]))
    assert result["mean_absolute_change"] is None
    assert result["component_weighted_loss_fraction"] is None
    assert result["neither_usable_heads"]==1
    json.dumps(result,allow_nan=False)


def test_qc_transitions_retain_denominator():
    rows=[(str(i),str(i),str(i),1,2,a,b) for i,(a,b) in enumerate([("usable","usable"),("usable","bad"),("bad","usable"),("bad","bad")])]
    result=summary.summarize_frame(frame(rows))
    assert result["scheduled_heads"]==4
    assert result["paired_usable_heads"]==result["lost_usable_heads"]==result["gained_usable_heads"]==result["neither_usable_heads"]==1


def test_spearman_can_be_perfect_with_large_systematic_shift():
    rows=[(str(i),str(i),str(i),i,i+100,"usable","usable") for i in range(4)]
    result=summary.summarize_frame(frame(rows))
    assert result["spearman_paired_component_means"]==pytest.approx(1)
    assert result["mean_absolute_change"]==100
    assert result["mean_signed_change"]==100


def joint_table(values):
    rows=[]
    for key,a,b in values:
        rows.append(dict(head_id="h",sha256="image",component_id="c",development_exposed=0,endpoint_id=key,
                         baseline_value=a,perturbed_value=b,baseline_status="usable",perturbed_status="usable"))
    return pd.DataFrame(rows)


def test_hue_wraps_at_zero_instead_of_358_degree_error():
    ids=["corolla_hue_sin","corolla_hue_cos"]
    table=joint_table([(ids[0],math.sin(math.radians(359)),math.sin(math.radians(1))),
                       (ids[1],math.cos(math.radians(359)),math.cos(math.radians(1)))])
    result=summary.summarize_frame(summary.joint_frame(table,ids,"hue"),joint=True)
    assert result["mean_absolute_change"]==pytest.approx(2)
    assert result["mean_signed_change"] is None
    assert result["spearman_paired_component_means"] is None


def test_composition_total_variation_and_closure():
    ids=["a","b","c","d"]
    table=joint_table(list(zip(ids,[1.,0.,0.,0.],[0.,1.,0.,0.])))
    result=summary.summarize_frame(summary.joint_frame(table,ids,"composition"),joint=True)
    assert result["mean_absolute_change"]==1
    table.loc[0,"perturbed_value"]=.5
    result=summary.summarize_frame(summary.joint_frame(table,ids,"composition"),joint=True)
    assert result["lost_usable_heads"]==1 and result["mean_absolute_change"] is None


def test_joint_missing_part_is_not_a_zero():
    ids=["a","b","c","d"]
    table=joint_table(list(zip(ids,[1,0,0,0],[None,None,None,None])))
    result=summary.summarize_frame(summary.joint_frame(table,ids,"composition"),joint=True)
    assert result["lost_usable_heads"]==1
    assert result["mean_absolute_change"] is None


def test_partial_run_cannot_be_summarized_as_complete(tmp_path):
    (tmp_path/"perturbation_report.json").write_text(json.dumps({"status":"TECHNICAL_PERTURBATION_PARTIAL","job_states":{"pending":1}}))
    with pytest.raises(ValueError,match="Incomplete"):
        summary.run(tmp_path,tmp_path,tmp_path,tmp_path/"out")


def complete_inputs(tmp_path):
    roots=[tmp_path/name for name in ("probes","baseline","groups")]
    for path in roots:
        path.mkdir()
    probes,baseline,groups=roots
    paths=[probes/"perturbations.sqlite",baseline/"measurements.sqlite",groups/"dependence_groups.sqlite"]
    conditions=[{"id":"condition_"+str(i)} for i in range(14)]
    with sqlite3.connect(paths[0]) as db:
        db.executescript("CREATE TABLE jobs(head_id,sha256,component_id); CREATE TABLE endpoints(head_id,condition,endpoint_id,value,status);")
        db.executemany("INSERT INTO jobs VALUES(?,?,?)",[("h1","image1","group1"),("h2","image2","group2")])
        db.executemany("INSERT INTO endpoints VALUES(?,?,?,?,?)",[
            (head,condition["id"],endpoint["endpoint_id"],None,"missing")
            for head in ("h1","h2") for condition in conditions for endpoint in summary.registry()])
    with sqlite3.connect(paths[1]) as db:
        db.execute("CREATE TABLE endpoints(head_id,endpoint_id,value,status)")
        db.executemany("INSERT INTO endpoints VALUES(?,?,?,?)",[(head,e["endpoint_id"],None,"missing") for head in ("h1","h2") for e in summary.registry()])
    with sqlite3.connect(paths[2]) as db:
        db.execute("CREATE TABLE components(component_id,development_pool_exposed,n_cached_objects)")
        db.executemany("INSERT INTO components VALUES(?,?,?)",[("group1",1,1),("group2",0,1)])
    receipt={"status":"TECHNICAL_PERTURBATION_COMPLETED","job_states":{"completed":2},
             "output_database_sha256":summary.digest(paths[0]),
             "execution_contract":{"input_database_sha256":{"measurement":summary.digest(paths[1]),"dependence_groups":summary.digest(paths[2])},"contract":{"conditions":conditions}}}
    (probes/"perturbation_report.json").write_text(json.dumps(receipt))
    return roots,paths,receipt


def test_complete_summary_preserves_missing_denominators(tmp_path):
    roots,paths,receipt=complete_inputs(tmp_path)
    result=summary.run(*roots,tmp_path/"out")
    table=pd.read_csv(tmp_path/"out/technical_sensitivity_summary.csv")
    assert result["summary_rows"]==len(table)==1218
    assert table.mean_absolute_change.isna().all()
    assert table[table.exposure_stratum.eq("all_cached")].scheduled_heads.eq(2).all()
    assert table[table.exposure_stratum.ne("all_cached")].scheduled_heads.eq(1).all()
    assert not result["independent_accuracy_estimated"]


def test_changed_input_stops_before_output(tmp_path):
    roots,paths,receipt=complete_inputs(tmp_path)
    with sqlite3.connect(paths[0]) as db:
        db.execute("DELETE FROM endpoints WHERE rowid=1")
    with pytest.raises(ValueError,match="identity changed"):
        summary.run(*roots,tmp_path/"out")
    assert not (tmp_path/"out").exists()


def test_missing_condition_row_stops_even_with_matching_receipt_hash(tmp_path):
    roots,paths,receipt=complete_inputs(tmp_path)
    with sqlite3.connect(paths[0]) as db:
        db.execute("DELETE FROM endpoints WHERE rowid=1")
    receipt["output_database_sha256"]=summary.digest(paths[0])
    (roots[0]/"perturbation_report.json").write_text(json.dumps(receipt))
    with pytest.raises(ValueError,match="denominator"):
        summary.run(*roots,tmp_path/"out")
    assert not (tmp_path/"out").exists()
