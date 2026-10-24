"""Test full-S1 source-name disaggregation without fictional 251 accepted species."""
import numpy as np
import pandas as pd
import pytest
from analysis.v4 import t2_s1_source_name_resolved_all251 as mod


def synthetic():
    records=[]
    five=dict(zip(mod.EXPOSURES,[10.,200.,1200.,50.,140.]))
    for key,name,values in (
        ("a","Cirsium_brevicaule",[60.,62.]),
        ("a","Cirsium_irumtiense",[75.,80.,81.]),
        ("b","Cirsium_vulgare",[67.,69.])):
        for v in values:
            records.append({"accepted_key":key,"taxon_name":name,
                mod.ENDPOINT:float(v),mod.ENDPOINT+"__eligible":True,
                **{k:float(value)+float(v)/100
                   for k,value in five.items()}})
    df=pd.DataFrame(records)
    mapper={
        "a":{"source_names":["Cirsium_brevicaule","Cirsium_irumtiense"],
             "tree_matches":["Cirsium_brevicaule","Cirsium_irumtiense"]},
        "b":{"source_names":["Cirsium_vulgare"],
             "tree_matches":["Cirsium_vulgare"]},
    }
    return df,mapper


def test_real_original_names_maintained_within_same_accepted_key():
    frame,mapx=synthetic()
    valid,unit=mod.source_rows(frame,mapx)
    assert len(valid)==7
    assert len(unit)==3
    assert len(unit.loc[unit.accepted_key=="a"])==2
    a=unit.set_index("original_source_name")
    assert np.isclose(a.loc["Cirsium_brevicaule","mean_image_angle_deg"],61)
    assert np.isclose(a.loc["Cirsium_irumtiense","mean_image_angle_deg"],(75+80+81)/3)
    assert a.loc["Cirsium_brevicaule","exact_source_name_has_frozen_placement_tip"]
    assert a.loc["Cirsium_irumtiense","exact_source_name_has_frozen_placement_tip"]
    assert (unit.biological_species_taxonomic_revision_verified==False).all()
    assert (unit.botanical_gravity_posture_inferred_from_photo_angle==False).all()


def test_distinct_source_names_cannot_ever_be_swapped_inside_accepted_key():
    frame,mapper=synthetic()
    mapper["a"]["source_names"]=["Cirsium_brevicaule"]
    with pytest.raises(ValueError,match="absent"):
        mod.source_rows(frame,mapper)


def test_invalid_photo_axis_excluded_or_fail():
    frame,mapper=synthetic()
    frame.loc[0,mod.ENDPOINT]=195.
    with pytest.raises(ValueError,match="outside EXIF"):
        mod.source_rows(frame,mapper)


def test_source_row_finite_original_five_climates_required():
    frame,mapper=synthetic()
    frame.loc[0,mod.EXPOSURES[2]]=np.nan
    valid,unit=mod.source_rows(frame,mapper)
    assert len(valid)==6
    assert unit.loc[unit.original_source_name=="Cirsium_brevicaule",
                    "original_photo_count"].iloc[0]==1


def test_original_full_scope_is_244_accepted_not_251_valid_species():
    assert mod.EXPECTED==(170066,244,251,7)
    assert mod.FOCAL_KEY=="206745753"
    assert mod.FOCAL=={"Cirsium_brevicaule","Cirsium_irumtiense"}
