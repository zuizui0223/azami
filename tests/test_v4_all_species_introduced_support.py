"""All-source Cirsium introduced-status taxon support and leakage checks.

Source names are evaluated BEFORE reading any photographed head-angle magnitude,
except legitimate image-output availability, which is documented selection.
No specimen coordinates, observer IDs or head images are exported.
"""
from __future__ import annotations

import pandas as pd
import pytest

from analysis.v4 import native_introduced_opportunity as old
from analysis.v4 import all_species_introduced_support as study


def miniature_source(n=55,collision=False):
    rows=[]
    for name,key in (("Cirsium arvense", "1001"),("Cirsium vulgare","1002")):
        for state in ("native","introduced"):
            for i in range(n):
                if name=="Cirsium vulgare" and state=="introduced" and i>=35:
                    # Original missing annotations are never biological zero
                    flag=False
                else:
                    flag=True
                rows.append({
                    "obs_id":f"{name}_{state}_{i}",
                    "taxon_name":name,
                    "taxon_rank":"species",
                    "accepted_key":int(key) if not collision or name=="Cirsium arvense" else 1001,
                    "native_primary":state=="native",
                    "eligibility_reason":"include" if state=="native" else "range_introduced",
                    "user_id":f"observer_{name}_{state}_{i%33}",
                    "latitude":float(-38+(i%17)*3.1),
                    "longitude":float(1+(i%16)*11),
                    "tdwg1":1 if state=="native" else 7,
                    "orientation_eligible":flag,
                    "colour_eligible":True,
                    "has_archived_measurement_row":flag,
                    "quality_grade":"research" if flag else "needs_id",
                    "source_photo_count":2
                })
    data=pd.DataFrame(rows)
    return data


def _run(monkeypatch,data):
    monkeypatch.setattr(study,"SOURCE_N",len(data))
    monkeypatch.setattr(study,"NATIVE_N",len(data)//2)
    return study.records_by_source_species(data)


def test_screen_all_species_not_just_two_preselected_names(monkeypatch):
    data=miniature_source(n=120)
    data=pd.concat([data,pd.DataFrame([
        {**data.iloc[0].to_dict(),"obs_id":"new_taxon_native",
         "taxon_name":"Cirsium specimenx", "accepted_key":1003,
         "native_primary":True,"eligibility_reason":"include"}
    ])],ignore_index=True)
    monkeypatch.setattr(study,"SOURCE_N",len(data))
    monkeypatch.setattr(study,"NATIVE_N",int(data.native_primary.sum()))
    result=study.records_by_source_species(data)
    assert result["source_species_rank_names"]==3
    assert result["source_species_with_both_wild_regions"]==2
    assert result["species_passing_basic_observer_cell_measurement_support"]==2
    assert result["species_passing_stronger_observer_cell_measurement_support"]==1
    byname={x["source_taxon_name"]:x for x in result["contiguous_source_species_counts"]}
    assert byname["Cirsium specimenx"]["both_ranges_pass_basic_measurement_support"] is False
    assert byname["Cirsium vulgare"]["range_support"]["introduced"]["orientation_eligible_observations"]==35
    assert byname["Cirsium vulgare"]["range_support"]["introduced"]["original_source_observations"]==120
    assert "obs_id" not in str(result)
    assert "observer_Cirsium" not in str(result)
    assert "observer_Cirsium_arvense" not in str(result)
    assert "latitude" not in str(result)


def test_source_WCVP_key_name_collision_must_disqualify_both_names(monkeypatch):
    d=miniature_source(n=55,collision=True)
    a=_run(monkeypatch,d)
    assert a["mixed_accepted_key_count_in_all_wild_source_species"]==1
    assert a["species_passing_basic_observer_cell_measurement_support"]==0
    assert all(row["taxon_concept"]=="ACCEPTED_KEY_COLLISION" for row in a["contiguous_source_species_counts"])


def test_inactive_or_no_photo_response_cannot_appear_as_two_range_candidate(monkeypatch):
    d=miniature_source(n=55)
    d.loc[d.taxon_name.eq("Cirsium arvense")&d.eligibility_reason.eq("range_introduced"),
          "orientation_eligible"]=False
    a=_run(monkeypatch,d)
    matched={x["source_taxon_name"]:x for x in a["contiguous_source_species_counts"]}
    assert matched["Cirsium arvense"]["both_ranges_pass_basic_measurement_support"] is False
    assert a["species_passing_basic_observer_cell_measurement_support"]==1


def test_eligibility_does_not_read_measured_angle_or_rain(monkeypatch):
    d=miniature_source(n=55)
    a=_run(monkeypatch,d)
    d["chelsa_bio12"]=range(len(d))
    d["orientation_angle_degrees_median"]=range(len(d))
    d["orientation_angle_degrees_median"]=180-d["orientation_angle_degrees_median"]
    b=_run(monkeypatch,d)
    assert a==b


def test_inconsistent_source_census_or_wild_range_fails_closed(monkeypatch):
    d=miniature_source(n=55)
    monkeypatch.setattr(study,"SOURCE_N",220)
    monkeypatch.setattr(study,"NATIVE_N",110)
    with pytest.raises(ValueError,match="observation census"):
        study.records_by_source_species(d.iloc[:-1])
    d.loc[0,"native_primary"]=False
    with pytest.raises(ValueError,match="partition"):
        study.records_by_source_species(d)
