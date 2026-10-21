"""Synthetic and source-boundary tests for ONE A→E workflow.

No synthetic output is presented as an empirical ecological conclusion.
"""
from __future__ import annotations
import json
import numpy as np
import pandas as pd
import pytest

from analysis.v4 import native_gates_ae as ae
from analysis.v4.native_gates_bcd import (
    fixed_sample, gate_c, gate_d_same_sample, wasserstein_parts,
)


def synthetic_native(n_taxa=5, n_each=100):
    rows=[]
    for k in range(n_taxa):
        for j in range(n_each):
            row={
                "obs_id":f"obs-{k:02d}-{j:03d}",
                "accepted_key":str(k+101),
                "taxon_name":f"Cirsium_synthetic{k}",
                "latitude":.01 + (j%12)*.001,
                "longitude":1.+(j%12)*5.1,
                "observed_on":"2026-04-10" if j%2 else "2026-09-10",
                "quality_grade":"research" if j%2 else "needs_id",
                "chelsa_bio12":700+j*5+k*10,
                "chelsa_rsds_mean":15+np.sin(j*.2),
            }
            for i,name in enumerate(ae.CORE9):
                if name=="corolla_hue_sin": v=np.sin(j*.087 +k*.1)
                elif name=="corolla_hue_cos": v=np.cos(j*.087 +k*.1)
                else:v=20+i+k*.2+np.sin(j*.09+i)+j*.005
                row[name]=v
                row[name+"__eligible"]=True
            rows.append(row)
    frame=pd.DataFrame(rows)
    obs=pd.DataFrame({
        "obs_id":frame.obs_id,"user_id":[f"u-{i//n_each}-{i%30}" for i in range(len(frame))],
        "quality_grade":frame.quality_grade,
    })
    mapping=[
        {"accepted_key":str(i+101),
         "source_names":["Cirsium_brevicaule","Cirsium_irumtiense"] if i==0 else
           ["Cirsium_kawakamii"] if i==1 else [f"Cirsium_synthetic{i}"],
         "tree_matches":[] if i==0 else [f"Cirsium_synthetic{i}"],
         "status":"ambiguous" if i==0 else "matched"}
        for i in range(n_taxa)
    ]
    return frame,obs,mapping


def test_gate_a_excludes_shared_accepted_key_without_discarding_original_sources(monkeypatch):
    frame,obs,mapping=synthetic_native()
    monkeypatch.setattr(ae,"EXPECTED_N",len(frame))
    monkeypatch.setattr(ae,"EXPECTED_TAXA",5)
    monkeypatch.setattr(ae,"EXPECTED_CORE9",(len(frame),5))
    result,selected=ae.gate_a(frame,mapping,obs,min_cells=5,min_observers=10,
                              min_balanced=50,max_taxa=4)
    assert result["status"]=="observer_cell_balanced_selection_constructed"
    assert result["matched_nine_output_taxa"]==5
    assert result["repeated_name_accepted_key_collisions"]==1
    assert result["matched_rows_in_collided_accepted_keys"]==100
    assert len(selected)>0 and "101" not in set(selected.accepted_key)
    assert selected.duplicated(["accepted_key","_cell","user_id"]).sum()==0
    assert result["source_original_observers"]>0
    assert result["observation_month"]["status"].startswith("source_vs")
    e=ae.gate_e(mapping)
    assert e["ryukyu_pair_share_wcvp_key"]
    assert e["eazami_source_concepts"][2]["status"].startswith("COLLISION")
    assert e["eazami_source_concepts"][3]["status"].startswith("COLLISION")
    assert e["eazami_source_concepts"][4]["status"]=="single_source_name_exact"
    assert e["eazami_source_concepts"][0]["status"]=="not_in_native_source_as_exact_concept"
    # DIL and GitHub artifact summaries must not expose private user or obs IDs.
    public=json.dumps(result)
    assert "obs-00" not in public and "u-0-" not in public
    # Actual shuffled-source invariance: positional order never selects an ID.
    a2,s2=ae.gate_a(frame.sample(frac=1,random_state=13),mapping,
         obs.sample(frac=1,random_state=71),min_cells=5,min_observers=10,
         min_balanced=50,max_taxa=4)
    assert a2==result
    assert set(s2.obs_id)==set(selected.obs_id)


def test_gate_a_original_observers_fail_closed_instead_of_photo_independence(monkeypatch):
    frame,obs,mapping=synthetic_native(5,20)
    monkeypatch.setattr(ae,"EXPECTED_N",len(frame))
    monkeypatch.setattr(ae,"EXPECTED_TAXA",5)
    monkeypatch.setattr(ae,"EXPECTED_CORE9",(len(frame),5))
    result,selected=ae.gate_a(frame,mapping,None)
    assert result["status"]=="blocked_missing_original_observers"
    assert result["predeclared_selection"]["selected_taxa"]==0
    assert selected.empty
    incomplete=obs.iloc[:-1]
    with pytest.raises(ValueError,match="observer roster"):
        ae.gate_a(frame,mapping,incomplete)
    wrong=obs.copy()
    wrong.loc[0,"obs_id"]=wrong.loc[1,"obs_id"]
    with pytest.raises(ValueError,match="observer roster"):
        ae.gate_a(frame,mapping,wrong)


def test_true_wasserstein_and_common_module_kde_not_parametric_ellipsoid():
    a=np.array([[0.],[1.],[2.]])
    w,centre,shape=wasserstein_parts(a,a+2)
    assert w==pytest.approx(2)
    assert centre==pytest.approx(2)
    assert shape==pytest.approx(0,abs=1e-6)
    frame,obs,_=synthetic_native(6,110)
    frame["user_id"]=obs.user_id
    # Already one-per-cell×observer, row grain remains source observation.
    frame["_cell"]=pd.Series(np.arange(len(frame))%12,dtype="int64")
    result=gate_c(frame,replicates=2,n=10,kde_draws=32,max_taxa=6)
    assert result["status"]=="completed_native_density_screen_not_functional_hypervolume"
    assert result["shared_taxon_sample_every_module"]
    for module in ("presentation","colour","outline"):
        v=result["modules"][module]
        assert v["taxa"]==6 and v["taxon_pairs"]==15
        assert v["median_empirical_wasserstein_w2"]>=0
        assert 0<=v["median_density_overlap"]<=1
        assert v["median_90pct_KDE_probability_volume"]>0
    assert result["modules"]["colour"]["dimensions"]==3
    assert result["modules"]["presentation"]["dimensions"]==1


def test_environment_comparison_preserves_full_native_median_and_same_sample():
    frame,obs,_=synthetic_native(6,100)
    frame["user_id"]=obs.user_id
    frame["_cell"]=pd.Series(np.arange(len(frame))%12,dtype="int64")
    d=gate_d_same_sample(frame,frame,{},n=10,minimum_side=3)
    assert d["same_candidate_observer_cell_frame_as_gates_B_C_D"]
    assert set(d["predictors"])=={"chelsa_bio12","chelsa_rsds_mean"}
    for mod in ("presentation","colour","outline"):
        assert "median_low_high_empirical_w2" in d["predictors"]["chelsa_bio12"][mod]


def test_unified_main_input_requires_same_verified_native_source_and_models(tmp_path):
    from analysis.v4.run_native_abcde import validate_inputs
    fake=tmp_path/"native.csv.gz"; fake.write_bytes(b"synthetic-not-scientific-input")
    paths=[tmp_path/k for k in ("receipt.json","models.json","integration.json","spatial.json")]
    for path in paths:path.write_text("{}")
    with pytest.raises(ValueError,match="corrected native source bytes"):
        validate_inputs(fake,*paths)
