"""Full S2 environmental × phenotype landscapes, source integrity and privacy."""
import json
import numpy as np
import pandas as pd
import pytest

from analysis.v4 import s2_joint_rain_angle_landscape as s2
from analysis.v4.native_spatial import cells


def source_pair(monkeypatch):
    rng=np.random.default_rng(20261021)
    rows=[]
    for name in s2.SPECIES:
        for status in ("native","introduced"):
            for i in range(1200):
                lat=-58.+(i%18)*5.3
                lon=-168.+((i//18)%20)*12.5
                rain=400+(i%91)*15+110*(status=="introduced")
                rain+=rng.normal(0,0.2)
                angle=75+0.008*(rain-600)+rng.normal(0,18)
                region=1 if status=="native" else (7 if i%3==0 else 8 if i%3==1 else 5)
                rows.append({
                    "obs_id":f"ob-{name}-{status}-{i}",
                    "user_id":f"photographer-{name}-{status}-{i}",
                    "taxon_name":name,"range_state":status,"tdwg1":region,
                    "latitude":float(lat),"longitude":float(lon),
                    "_cell":int(cells([lat],[lon])[0]),
                    "observed_on":"2024-07-12",
                    "chelsa_bio01":9+rng.normal(),
                    "chelsa_bio04":810+rng.normal(0,50),
                    "chelsa_bio12":rain,
                    "chelsa_bio15":35+rng.normal(0,5),
                    "chelsa_rsds_mean":13+rng.normal(0,0.5),
                    "orientation_angle_degrees_median":float(np.clip(angle,0,180)),
                })
    photo=pd.DataFrame(rows)
    source=photo.drop(columns=["orientation_angle_degrees_median"]).copy()
    monkeypatch.setattr(s2,"EXPECTED_SOURCE",len(source))
    monkeypatch.setattr(s2,"EXPECTED_PHOTO",len(photo))
    return source,photo


def test_photo_level_joint_and_all_introduced_mixture(monkeypatch):
    wild,photo=source_pair(monkeypatch)
    report=s2.evaluate(wild,photo)
    assert report["status"]=="ORIGINAL_S2_JOINT_RAIN_AND_IMAGE_ANGLE_ALL_INTRODUCED"
    assert report["visualization_only_no_trait_width_or_tail_indices_selected"]
    for name,case in report["species"].items():
        nr=len(case["display_climate_BIO12_edges_mm_year"])-1
        na=len(case["display_angle_edges_image_degrees"])-1
        assert nr>10 and na==36
        for status in ("native","all_introduced"):
            part=case[status]
            prob=np.asarray(part["original_photo_joint_mass"])
            assert prob.shape==(nr,na)
            assert np.isclose(prob.sum(),1)
            assert np.all(prob>=0)
            assert len(part["original_wild_source_rain_bin_mass"])==nr
            conditional=part["photo_angle_conditional_on_rain_or_null"]
            assert len(conditional)==nr
            observed=np.array(part["photo_n_per_display_rain_bin"])
            for i,row in enumerate(conditional):
                if row[0] is None:
                    assert (observed[i]<s2.MIN_DRAW_ROWS or
                            part["photo_cells_per_display_rain_bin"][i]<s2.MIN_DRAW_CELLS)
                else:
                    assert np.isclose(sum(row),1)
        mix=case["introduced_continents_supplementary"]
        assert mix["all_introduced_joint_mass_recovered_exactly"]
        assert np.isclose(sum(g["fraction_global_cell_weight"]
                              for g in mix["groups"].values()),1)
        total=np.zeros((nr,na),float)
        for v in mix["groups"].values():
            total+=v["fraction_global_cell_weight"]*np.asarray(v["joint_mass_within_region"])
        assert np.allclose(total,case["all_introduced"]["original_photo_joint_mass"])


def test_climate_plotting_grid_cannot_change_with_orientation(monkeypatch):
    wild,photo=source_pair(monkeypatch)
    result=s2.evaluate(wild,photo)
    inverted=photo.copy()
    inverted["orientation_angle_degrees_median"]=180-inverted["orientation_angle_degrees_median"]
    alt=s2.evaluate(wild,inverted)
    for name in s2.SPECIES:
        a=result["species"][name]
        b=alt["species"][name]
        assert a["display_climate_BIO12_edges_mm_year"]==b["display_climate_BIO12_edges_mm_year"]
        assert a["native"]["original_wild_source_rain_bin_mass"]==(
            b["native"]["original_wild_source_rain_bin_mass"])
        assert a["native"]["original_photo_joint_mass"]!=(
            b["native"]["original_photo_joint_mass"])


def test_source_missing_cannot_pass_and_output_aggregate_only(monkeypatch,tmp_path):
    wild,photo=source_pair(monkeypatch)
    with pytest.raises(ValueError,match="metadata"):
        s2.evaluate(wild.drop(columns=["tdwg1"]),photo)
    saved=s2.save(wild,photo,tmp_path)
    raw=(tmp_path/"S2_full_joint_rain_angle_landscape.json").read_text()
    assert json.loads(raw)["status"]==saved["status"]
    for token in ('"obs_id"','"user_id"','"latitude"','"longitude"'):
        assert token not in raw
    assert not any(saved["privacy"].values())
