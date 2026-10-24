"""S1/S2 four-axis extension retains public source, original image-angle grain."""
import numpy as np
import pandas as pd
from analysis.v4 import s1_s2_real_bio18_bio05_rsds_elevation as m

def synthetic_frame(n=205,photos=20):
    rng=np.random.default_rng(20261024)
    v=[]
    for sp in range(n):
        for im in range(photos):
            a=rng.normal(0,1,4)
            v.append({"accepted_key":str(sp),"source_name":f"Cirsium_source{sp}",
                 "obs_id":f"{sp}_{im}",m.S1_ENDPOINT:70+a[0]*.5-a[1]+a[2],
                 "chelsa_bio18":1000+150*a[0]+sp*.1,
                 "chelsa_bio05":30+3*a[1]+sp*.01,
                 "chelsa_rsds_mean":14+1.5*a[2],
                 "elevation_m":400+180*a[3]+sp*.1})
    return pd.DataFrame(v)

def test_four_variables_are_really_requested_and_distinct():
    assert m.EXPOSURES==("chelsa_bio18","chelsa_bio05","chelsa_rsds_mean","elevation_m")
    assert m.BIO5_URL.endswith("CHELSA_bio05_1981-2010_V.2.1.tif")
    assert m.BIO18_URL.endswith("CHELSA_bio18_1981-2010_V.2.1.tif")
    assert m.ELEV_URL.endswith("wc2.1_30s_elev.zip")
    assert m.S1_N==170066 and m.S2_PHOTO_N==173737
    assert m.S2_WILD_N==297622 and m.OUTPUT_SPECIES==251
    assert m.ORIGINAL_251_SHA=="47b32972b20564ade163907ef6fdf3eb2a69f8b42c8d87b8cf4ea75ca1682dc1"

def test_original_source_head_photo_continuous_angle_4axis_HC3():
    f=synthetic_frame()
    X,sd=m.four_design(f)
    y=f[m.S1_ENDPOINT].to_numpy(float)
    fit=m.robust_hc3(y,X,["intercept",*m.EXPOSURES])
    assert len(fit["coefficients"])==5
    assert all(len(v)==2 for v in fit["HC3_95CI"].values())
    assert fit["coefficients"]["chelsa_bio18"]>0
    assert fit["coefficients"]["chelsa_bio05"]<0
    assert sd["elevation_m"]["SD"]>0

def test_taxon_means_and_within_groups_never_dichotomized(monkeypatch):
    f=synthetic_frame()
    monkeypatch.setattr(m,"S1_N",len(f))
    among,desc=m.s1_between(f)
    assert len(among)==205
    assert desc["min_photo_10"]["n_source_names"]==205
    assert desc["min_photo_1"]["n_original_photo_observations"]==len(f)
    within=m.s1_within(f)
    assert within["n_source_taxon_clusters"]==205
    assert "chelsa_bio05" in within["coefficients_degrees_per_full_S1_photo_SD"]
    assert "chelsa_bio18" in within["coefficients_degrees_per_full_S1_photo_SD"]

def test_no_extra_gbif_site_classification_imputed():
    assert m.NEW_DOWNLOADS["chelsa_bio18"][1]==549240986
    assert m.NEW_DOWNLOADS["chelsa_bio05"][1]==165253968
    assert m.NEW_DOWNLOADS["worldclim_elevation"][1]==338188707
