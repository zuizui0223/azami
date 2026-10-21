"""Exact introduced-range source status and phenotype availability guards.

Only synthetic rows; no raw geographic points or human IDs in tests.
"""
from __future__ import annotations

import io
import json
import zipfile

import numpy as np
import pandas as pd
import pytest

from analysis.v4 import native_introduced_opportunity as m


def _test_photo_rows():
    return pd.DataFrame([
        {"obs_id":"o1","photo_id":"p1","user_id":"u1",
         "taxon_name":"Cirsium vulgare","taxon_rank":"species",
         "latitude":45.,"longitude":-100.,
         "quality_grade":"research","observed_on":"2025-05-01",
         "source_photo_count":1,
         "native_range_status":"introduced","native_primary":False,
         "eligibility_reason":"range_introduced","tdwg_level3_code":"ABT","tdwg1":7},
        {"obs_id":"o2","photo_id":"p2","user_id":"u2",
         "taxon_name":"Cirsium vulgare","taxon_rank":"species",
         "latitude":45.,"longitude":-100.,
         "quality_grade":"needs_id","observed_on":"2025-05-01",
         "source_photo_count":1,
         "native_range_status":"native","native_primary":True,
         "eligibility_reason":"include","tdwg_level3_code":"ABT","tdwg1":7},
        {"obs_id":"o3","photo_id":"p3","user_id":"u3",
         "taxon_name":"Cirsium arvense","taxon_rank":"species",
         "latitude":45.,"longitude":-100.,
         "quality_grade":"casual","observed_on":"2025-05-01",
         "source_photo_count":1,
         "native_range_status":"introduced","native_primary":False,
         "eligibility_reason":"captive","tdwg_level3_code":"ABT","tdwg1":7},
        {"obs_id":"o4","photo_id":"p4","user_id":"u4",
         "taxon_name":"Cirsium arvense","taxon_rank":"species",
         "latitude":45.,"longitude":-100.,
         "quality_grade":"research","observed_on":"2025-05-01",
         "source_photo_count":1,
         "native_range_status":"unlisted","native_primary":False,
         "eligibility_reason":"range_unlisted","tdwg_level3_code":"ABT","tdwg1":7},
    ])


def test_2000s_TDWG_continent_source_level_and_native_intro_not_complement():
    frame=_test_photo_rows()
    statuses=m.range_labels(frame)
    assert statuses.tolist()==[
        "WCVP_introduced_TDWG","WCVP_native_TDWG",
        "captive","WCVP_uncertain_or_unlisted"
    ]
    assert int(frame.eligibility_reason.eq("range_introduced").sum())==1
    assert int(frame.native_primary.sum())==1
    assert int((~frame.native_primary).sum())==3
    assert sum(m.CONTINENTS[k] in ("North America","South America")
               for k in (7,8))==2
    cells=m.angular_cells(frame.latitude,frame.longitude)
    assert cells.nunique()==1
    with pytest.raises(ValueError,match="valid coordinates"):
        m.angular_cells(pd.Series([91.]),pd.Series([0.]))


def test_historical_observation_measurements_only_one_value_per_original_obs(monkeypatch):
    mrows=_test_photo_rows().iloc[:2].copy()
    monkeypatch.setattr(m,"SOURCE_N",2)
    payload=pd.DataFrame([
       {"obs_id":"o1","orientation_angle_degrees_median":85.,
        "n_usable_heads_orientation":1,
        "corolla_lab_chroma_median":5.,"n_usable_heads_colour":2},
       {"obs_id":"o2","orientation_angle_degrees_median":85.,
        "n_usable_heads_orientation":0,
        "corolla_lab_chroma_median":5.,"n_usable_heads_colour":0}
    ]).to_csv(index=False).encode()
    out=m.append_observation_measures(mrows,payload)
    assert out.orientation_eligible.tolist()==[True,False]
    assert out.colour_eligible.tolist()==[True,False]
    assert out.has_archived_measurement_row.tolist()==[True,True]
    dupe=pd.read_csv(io.BytesIO(payload))
    dupe=pd.concat([dupe,dupe.iloc[[0]]],ignore_index=True)
    with pytest.raises(ValueError,match="duplicate observation"):
        m.append_observation_measures(mrows,dupe.to_csv(index=False).encode())


def test_unmeasured_introduced_observation_is_not_zero_orientation(monkeypatch):
    src=_test_photo_rows().iloc[:2].copy()
    monkeypatch.setattr(m,"SOURCE_N",2)
    out=m.append_observation_measures(src,None)
    assert not out.orientation_eligible.any()
    assert not out.has_archived_measurement_row.any()


def test_authority_blob_requires_original_source_artifact_hash():
    with pytest.raises(ValueError,match="authority artifact ZIP"):
        m.authority_source(b"fake-archive")


def test_native_plus_introduced_plus_other_partition_excludes_captive_and_unknown(monkeypatch):
    frame=_test_photo_rows().copy()
    frame["has_archived_measurement_row"]=True
    frame["orientation_eligible"]=frame.obs_id.eq("o1")
    frame["colour_eligible"]=frame.obs_id.eq("o1")
    monkeypatch.setattr(m,"SOURCE_N",4)
    monkeypatch.setattr(m,"SOURCE_PHOTOS",4)
    monkeypatch.setattr(m,"NATIVE_N",1)
    monkeypatch.setattr(m,"NATIVE_ACCEPTED",1)
    monkeypatch.setattr(m,"NATIVE_MEASURED_N",1)
    monkeypatch.setattr(m,"TARGET",("Cirsium vulgare",))
    frame["accepted_key"]=[101,101,102,102]
    report=m.aggregate(frame)
    o=report["overall"]
    assert (o["native_primary_wild_observations"],o["introduced_primary_wild_observations"],
           o["excluded_or_undetermined_observations"])==(1,1,2)
    assert o["status_partition_conserved"]
    assert o["source_measurement_coverage"]["introduced_orientation_eligible_obs"]==1
    assert report["focal_native_intro_orientation_support"][0]["both_range_orientation_supported"] is False
    assert report["focal_source_species_by_TDWG_lvl1"][1]["TDWG_lvl1_group"]=="North America"


def test_decrypt_verification_requires_source_sha_even_before_key():
    with pytest.raises(ValueError,match="Encrypted original source identity changed"):
        m.decrypt_checked(b"untrusted","",m.ORIGINAL_SEALED_SHA,m.ORIGINAL_SHA,m.ORIGINAL_AAD)
