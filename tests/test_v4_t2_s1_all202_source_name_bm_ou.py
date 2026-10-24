"""Unit tests that all202 tips are exactly source named, not 194 alias keys."""
import numpy as np
import pandas as pd
import pytest
from analysis.v4 import t2_s1_all202_source_name_bm_ou as a


def fake_source(tmp_path):
    rng=np.random.default_rng(20261022)
    rows=[]
    for i in range(244):
        collided=i<7
        for j in range(2 if collided else 1):
            included=(j<2 if i<4 else j==0) if collided else i<198
            # First seven collided; 191 clean matching rows are indices 7:198.
            rows.append({
                "accepted_key":str(206000000+i),
                "original_source_name":f"Cirsium_plant_{i:03d}_{j}",
                "original_photo_count":677,
                "mean_image_angle_deg":float(rng.uniform(40,140)),
                "median_image_angle_deg":float(rng.uniform(40,140)),
                "WCVP_key_collides_source_names":collided,
                "exact_source_name_has_frozen_placement_tip":included,
                **{e:float(rng.normal(100,20)) for e in a.EXPOSURES},
            })
    assert len(rows)==251
    rows[0]["original_photo_count"]+=170066-677*251
    df=pd.DataFrame(rows)
    output=tmp_path/"demo.csv"
    df.to_csv(output,index=False)
    return output


def test_full_251_source_inventories_still_include_202_exact_named_original_tips(tmp_path):
    path=fake_source(tmp_path)
    groups,sha=a.read_source(path)
    assert len(groups["matched202"])==202
    assert len(groups["clean191"])==191
    assert len(groups["n_ge_5"])==202
    assert len(groups["n_ge_10"])==202
    assert len({x["tip"] for x in groups["matched202"]})==202
    assert len(a.designs(groups))==6
    assert sum(x["collision"] for x in groups["matched202"])==11
    assert len(sha)==64


def test_frozen_five_covariates_and_photo_angle_not_nodding():
    assert a.MATCHED_TIPS==202
    assert a.ORIGINAL_ALIAS_CLEAN==191
    assert a.SOURCE_UNITS==251
    assert a.ACCEPTED_KEYS==244
    assert a.TOTAL_NATIVE==170066
    assert a.EXPOSURES[2]=="chelsa_bio12"


def test_missing_one_exact_placement_fails_closed(tmp_path):
    f=fake_source(tmp_path)
    d=pd.read_csv(f)
    d.loc[d.exact_source_name_has_frozen_placement_tip].iloc[0]
    index=d.index[d.exact_source_name_has_frozen_placement_tip][0]
    d.loc[index,"exact_source_name_has_frozen_placement_tip"]=False
    d.to_csv(f,index=False)
    with pytest.raises(ValueError,match="202 exact original terminal"):
        a.read_source(f)


def test_source_operational_alias_collision_is_not_promotion_to_true_species(tmp_path):
    f=fake_source(tmp_path)
    df=pd.read_csv(f,dtype={"accepted_key":str})
    assert df.groupby("accepted_key").size().gt(1).sum()==7
    assert len(df)-df.accepted_key.nunique()==7
    assert set(a.CASE_NAMES)=={t[0] for t in a.designs(None)}
