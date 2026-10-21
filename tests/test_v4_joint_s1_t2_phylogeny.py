"""Frozen 244/194 S1 evidence and independent T2 morphology join tests."""
import json
from pathlib import Path
import pandas as pd
import pytest
from analysis.v4 import joint_s1_t2_phylogeny as j

ROOT=Path(__file__).resolve().parents[1]

def test_frozen_exposure_cohorts_and_existing_phi_fit():
    p=pd.read_csv(ROOT/j.S1_PATH)
    s=pd.read_csv(ROOT/j.PGLS_PATH)
    models=json.loads((ROOT/j.MODEL_PATH).read_text())
    cover=json.loads((ROOT/j.COVERAGE_PATH).read_text())[j.ENDPOINT]
    rows=j.s1_table(p,s,models,cover)
    assert len(rows)==10
    assert {r["exposure"] for r in rows}==set(j.EXPOSURES)
    assert {r["S1_full_accepted_keys"] for r in rows}=={244}
    assert {r["S1_tree_covered_keys"] for r in rows}=={194}
    assert {r["placement_scenario_count"] for r in rows}=={52}
    b12=next(r for r in rows if r["term"]=="chelsa_bio12__between")
    b15=next(r for r in rows if r["term"]=="chelsa_bio15__between")
    assert b12["S1_effect_image_degrees"]==pytest.approx(.4381991922885207)
    assert b12["supported_across_52_pointwise"]=="all_positive"
    assert b15["placement_95pct_model_interval_envelope"][0]<0
    assert b15["supported_across_52_pointwise"]=="no_uniform_pointwise_support"
    models["selection"]["selected"][0]="chelsa_vpd_mean"
    with pytest.raises(ValueError,match="five-exposure"):
        j.s1_table(p,s,models,cover)

def fake_states():
    ea=[{"accepted_taxon":f"Cirsium East{i}","orientation_state":"D" if i<5 else "U",
         "orientation_resolved":"true" if i<17 else "false"} for i in range(20)]
    na=[{"Species":f"Cirsium_NA{i}","HeadPosition":"1" if i<5 else "0"} for i in range(63)]
    return ea,na

def test_state_source_independence_and_WCVP_collision_blocks():
    ea,na=fake_states()
    states=j.orientation_state_rows(ea,na)
    assert len(states)==80 and sum(r["source"].startswith("EAzami") for r in states)==17
    assert all("image" not in row for row in states)
    source={
        "1":{"source_names":["Cirsium East0"],"n":50,"median_image_angle_deg":98.1},
        "2":{"source_names":["Cirsium East1","Cirsium East2"],"n":90,"median_image_angle_deg":60.0},
        "3":{"source_names":["Cirsium NA0"],"n":3,"median_image_angle_deg":80.0}
    }
    tree_idx={"1":{"tree_matches":{"Cirsium_East0"}}}
    joined={r["source_taxon"]:r for r in j.match_sources(states,source,tree_idx)}
    assert joined["Cirsium East0"]["join_status"]=="exact_S1_and_tree"
    assert joined["Cirsium East1"]["join_status"]=="blocked_WCVP_source_name_collision"
    assert joined["Cirsium NA0"]["join_status"]=="exact_S1_only_tree_uncovered"
    assert joined["Cirsium East5"]["join_status"]=="not_in_S1_orientation_cohort"
    ea[1]["accepted_taxon"]=ea[0]["accepted_taxon"]
    with pytest.raises(ValueError,match="Duplicate East Asian"):
        j.orientation_state_rows(ea,na)

def test_frozen_244_requires_native_eligible_population():
    models=json.loads((ROOT/j.MODEL_PATH).read_text())
    cover=json.loads((ROOT/j.COVERAGE_PATH).read_text())[j.ENDPOINT]
    frame=pd.DataFrame({"obs_id":[1,2],"accepted_key":["a","a"],
       "taxon_name":["Cirsium x","Cirsium y"],
       j.ENDPOINT:[90.,100.],j.ENDPOINT+"__eligible":[True,True],
       **{v:[0.,0.] for v in j.EXPOSURES}})
    with pytest.raises(ValueError,match="Native 244/170066"):
        j.actual_native_roster(frame,models,cover)
    frame.loc[1,"obs_id"]=1
    with pytest.raises(ValueError,match="IDs not unique"):
        j.actual_native_roster(frame,models,cover)
