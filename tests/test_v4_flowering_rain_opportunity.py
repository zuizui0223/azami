"""Source/season specification tests; no ecological effect is invented here."""
import copy
import numpy as np
import pandas as pd
import pytest

from analysis.v4.flowering_rain_opportunity import (
    flowering_windows,date_columns, SEASON_WINDOW, MIN_EVENTS,
)


def source_example():
    rows=[]
    # All seven frozen multi-name accepted-key conflicts are still represented
    # in the source before eligible phenotypes, not accepted as source species.
    collision=[str(1000+i) for i in range(7)]
    mapping=[
        {"accepted_key":k,"source_names":["Cirsium_a","Cirsium_b"],
         "tree_matches":[],"status":"ambiguous"}
        for k in collision
    ]
    mapping.append({"accepted_key":"42","source_names":["Cirsium_arvense"],
                    "tree_matches":["Cirsium_arvense"],"status":"matched"})
    for key in collision:
        rows.append({
            "accepted_key":key,"obs_id":"conflict-"+key,
            "user_id":"conflict-user-"+key, "taxon_name":"Cirsium_a",
            "observed_on":"2023-06-02","latitude":-24.,
            "longitude":10.,"corolla_lab_chroma__eligible":True,
            "corolla_lab_chroma":20.,"orientation_image_vertical_angle__eligible":True,
            "orientation_image_vertical_angle":68.
        })
    # Exactly one species×20deg latitude×90deg longitude region, but three 5deg cells.
    for i in range(96):
        month=[12,1,2][i%3]
        year=[2021,2022,2023][i%3]
        cell=i%4
        rows.append({
            "accepted_key":"42","obs_id":f"native-{i}",
            "user_id":f"u-{i%32}","taxon_name":"Cirsium_arvense",
            "observed_on":f"{year}-{month:02d}-11",
            "latitude":-24.,"longitude":8+cell*5.5,
            "corolla_lab_chroma__eligible":True,"corolla_lab_chroma":25+i*.1,
            "orientation_image_vertical_angle__eligible":True,
            "orientation_image_vertical_angle":70.+i*.01
        })
    return pd.DataFrame(rows),mapping


def test_phenology_three_month_circular_window_with_taxon_concept_quarantine():
    f,mapping=source_example()
    a,private=flowering_windows(f,mapping)
    assert a["reference_month_count"]==3
    assert a["taxon_concept_collision_keys_excluded"]==7
    assert a["taxon_region_combinations_examined"]==1
    assert a["phenology_qualified_taxon_region_combinations"]==1
    assert a["phenology_qualified_unique_accepted_keys"]==1
    assert a["source_native_observations"]==103
    assert a["valid_native_dates"]==103
    assert len(private)==1
    assert list(private.iloc[0][["season_month_1","season_month_2","season_month_3"]]) == [12,1,2]
    assert private.iloc[0]["source_distinct_years"]>=2
    assert a["photo_floret_indicator"].startswith("corolla_lab_chroma__eligible")
    assert not a["seasonal_precipitation_from_monthly_CHELSA_extracted"]
    assert not a["per_photo_head_identity_verified"]


def test_season_choice_ignores_orientation_magnitude_colour_value_and_precipitation():
    f,mapping=source_example()
    a,b=flowering_windows(f,mapping)
    g=f.sample(frac=1,random_state=19).reset_index(drop=True)
    g["orientation_image_vertical_angle"]=np.linspace(0,180,len(g))
    g["corolla_lab_chroma"]=100.
    g["chelsa_bio12"]=np.linspace(0,1000,len(g))
    new,c=flowering_windows(g,mapping)
    assert a["phenology_qualified_taxon_region_combinations"]==new["phenology_qualified_taxon_region_combinations"]
    assert b.sort_values("accepted_key").reset_index(drop=True).equals(
        c.sort_values("accepted_key").reset_index(drop=True)
    )


def test_date_assessability_and_observer_gate_fail_closed():
    f,mapping=source_example()
    bad=f.copy()
    bad.loc[7,"observed_on"]="not-a-date"
    out=date_columns(bad)
    assert out._month.isna().sum()==1
    bad=f.drop(columns=["user_id"])
    with pytest.raises(ValueError,match="Original observer identities"):
        flowering_windows(bad,mapping)
    bad=f.copy()
    bad.loc[7,"obs_id"]=bad.loc[8,"obs_id"]
    with pytest.raises(ValueError,match="Source observations not unique"):
        flowering_windows(bad,mapping)
    bad=f.copy()
    bad.loc[7,"latitude"]=94.
    with pytest.raises(ValueError,match="Earth"):
        flowering_windows(bad,mapping)


def test_disjunct_source_region_does_not_conflate_hemispheres():
    frame,mapping=source_example()
    subset=frame.loc[frame.accepted_key.eq("42")].copy()
    southern=subset.copy()
    southern["obs_id"]="south-"+southern.obs_id
    southern["user_id"]="s-"+southern.user_id
    southern["latitude"]=-24
    northern=subset.copy()
    northern["obs_id"]="north-"+northern.obs_id
    northern["user_id"]="n-"+northern.user_id
    northern["latitude"]=25
    # In the northern hemisphere, local blooming is June–August.
    old=pd.to_datetime(northern.observed_on)
    northern["observed_on"]=[
        f"{year}-{mon:02d}-11" for year,mon in
        zip(old.dt.year,old.dt.month.map({12:6,1:7,2:8}))
    ]
    rest=frame.loc[~frame.accepted_key.eq("42")]
    both=pd.concat([rest,southern,northern],ignore_index=True)
    report,wins=flowering_windows(both,mapping)
    assert len(wins)==2
    seasons={tuple(int(y) for y in row[["season_month_1","season_month_2","season_month_3"]])
             for _,row in wins.iterrows()}
    assert seasons=={(12,1,2),(6,7,8)}
    assert report["phenology_qualified_taxon_region_combinations"]==2


def test_absent_inferred_phenology_not_misreported_as_zero_flowering():
    f,mapping=source_example()
    idx=f.accepted_key.eq("42")
    f.loc[idx,"corolla_lab_chroma__eligible"]=False
    report,win=flowering_windows(f,mapping)
    assert win.empty
    assert report["phenology_qualified_taxon_region_combinations"]==0
    assert report["season_coverage_fraction_distribution"]["median"] is None
    assert report["source_floret_proxy_eligible_dated"]==7
