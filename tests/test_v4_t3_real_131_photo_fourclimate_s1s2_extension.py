"""T3 source-only 131 real ITS-tip photo candidates may NOT imply phenotypes."""
import numpy as np
import pytest
import pandas as pd
from analysis.v4 import t3_real_131_photo_fourclimate_s1s2_extension as m

def test_same_S1_S2_original_four_climate_axes():
    assert m.FOUR==("chelsa_bio18","chelsa_bio05","chelsa_rsds_mean","elevation_m")
    assert m.MAX_PHOTO_RECORDS==30
    assert m.MIN_INDEPENDENT_PHOTO_SUPPORT==10
    assert m.MIN_INDEPENDENT_CELLS==3
    assert len(m.ALLOWED)==3

def test_only_exact_source_species_with_original_georeferenced_stillimage_is_candidate():
    r={"species":"Cremanthodium campanulatum","decimalLatitude":28.1,
       "decimalLongitude":102.3,"key":111,"license":"http://creativecommons.org/licenses/by/4.0/legalcode",
       "media":[{"type":"StillImage","identifier":"https://example.org/photo.jpg"}]}
    ok=m.photo_entry(r,"Cremanthodium campanulatum",123)
    assert ok["_latitude"]==28.1 and ok["_longitude"]==102.3
    assert ok["source_taxon"]=="Cremanthodium campanulatum"
    assert "orientation" not in ok
    assert m.photo_entry(r,"Cremanthodium decaisnei",123) is None
    no_photo={**r,"media":[]}
    assert m.photo_entry(no_photo,"Cremanthodium campanulatum",123) is None
    no_coords={**r,"decimalLatitude":None}
    assert m.photo_entry(no_coords,"Cremanthodium campanulatum",123) is None

def test_all131_census_is_taxon_source_not_binary_posture_labels():
    census=[{"original_GenBank_ITS_taxon":"Cremanthodium campanulatum",
             "GBIF_strict_species_key":123,
             "n_photo_observation_valid":3},
            {"original_GenBank_ITS_taxon":"Cremanthodium decaisnei",
             "GBIF_strict_species_key":None,"n_photo_observation_valid":0}]
    frame=[
      {"source_taxon":"Cremanthodium campanulatum","_latitude":28.0,"_longitude":102.0,
       "gbif_occurrence_key":"1",**{z:10.0 for z in m.FOUR}},
      {"source_taxon":"Cremanthodium campanulatum","_latitude":29.0,"_longitude":102.1,
       "gbif_occurrence_key":"2",**{z:11.0 for z in m.FOUR}},
      {"source_taxon":"Cremanthodium campanulatum","_latitude":30.0,"_longitude":102.2,
       "gbif_occurrence_key":"3",**{z:12.0 for z in m.FOUR}}]
    # aggregate conserves source census when no images found for a species.
    out=m.aggregate(frame,census,"original_roster_sha",{})
    assert out["n_georeferenced_photo_records_before_all4_climate"]==3
    assert out["n_actual_qualified_T3_S1_equivalent_continuous_head_angle_measurements"]==0
    assert out["n_species_qualify_for_MANUAL_photo_angle_audit"]==0
    assert out["no_cross_tribe_angle_environment_effect_fitted"]
    assert out["by_genus"]["Cremanthodium"]["n_original_ITS_source_taxa"]==2
