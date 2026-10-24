"""Four new climates use original EAzami CTMC and exact count-preserving maps."""
import numpy as np
import pandas as pd
from analysis.v4 import t2_eazami_original_ctmc_fourclimate_transition as m

def test_frozen_EAzami_original_binary_CTMC_nuclear_contract():
    assert m.AXES==("chelsa_bio18","chelsa_bio05","chelsa_rsds_mean","elevation_m")
    assert m.STATE_PANELS=={"10":(9,5,4,126),"5":(12,7,5,792),"3":(13,7,6,1716)}
    assert m.TREE_BLOB=="d7f8c0b96ee9c257c815db06a5a0aa9fdae4841e"
    assert m.SCRIPT_BLOB=="09731b563907b3ec2778ea04086d729771376f7e"
    assert m.CROSSWALK_BLOB=="f655544e77f1659323c1e915567dacbb51b01cc3"

def test_four_climates_require_real_source_occurrences_at_each_original_panel():
    taxa=["Cirsium alpicola","Cirsium brevicaule"]
    data=[]
    for t,n in zip(taxa,(4,5)):
        for k in range(n):
            data.append({"scientific_name_query":t,
               "chelsa_bio18":500+k*10,
               "chelsa_bio05":20+k*.5,
               "chelsa_rsds_mean":12+k*.2,
               "elevation_m":100+k*22})
    frame=pd.DataFrame(data)
    ref=type("FrozenRef",(),{"zscore":staticmethod(
        lambda x: (np.asarray(x)-np.asarray(x).mean())/np.asarray(x).std(ddof=1))})
    env,gate=m.panel_data(ref,frame,taxa,5)
    assert env is None
    assert gate["missing"]["Cirsium alpicola"]==4
    frame.loc[frame.scientific_name_query=="Cirsium alpicola","chelsa_bio18"]=np.nan
    env,gate=m.panel_data(ref,frame,taxa,1)
    assert env is None and gate["status"].startswith("STOP_")

def test_EAzami_source_metadata_and_not_original_flora_angle_substitution():
    assert m.EXPECTED_RSDS_SHA=="cc1ad14a70a1a940be3d01d918d4c3f7e35db95b86b12af9be5487e1bbc1fc7d"
    assert m.N_ORIG_OCC=={"Japan":693,"Taiwan":55}
    assert set(m.EXPECTED_RASTER_SHA)=={"chelsa_bio18","chelsa_bio05","elevation_m"}
    assert m.LEGACY_BIO15==[.267907,.269417]
    assert m.LEGACY_BIO1==[-.199306,-.199151]
