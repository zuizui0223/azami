"""Original 244 accepted-key versus 251 source-name BIO12 sensitivity gates."""
import numpy as np
import pandas as pd
import pytest
from analysis.v4 import t2_s1_251_vs244_bio12_taxonomy_sensitivity as a


def fake(monkeypatch):
    monkeypatch.setattr(a,"N_KEYS",55)
    monkeypatch.setattr(a,"N_UNITS",62)
    rng=np.random.default_rng(20261022)
    records=[]
    for key in range(55):
        source_count=2 if key<7 else 1
        for m in range(source_count):
            b=1000+600*rng.random()
            p=25+int(rng.integers(0,20))
            original=f"Cirsium_source_{key}_{m}"
            five={
                a.EXPOSURES[0]:5+20*rng.random(),
                a.EXPOSURES[1]:250+700*rng.random(),
                a.EXPOSURES[2]:b,
                a.EXPOSURES[3]:20+70*rng.random(),
                a.EXPOSURES[4]:10+20*rng.random(),
            }
            y=65+.004*b+1.3*m+rng.normal()
            records.append({"accepted_key":str(key),"original_source_name":original,
                "original_photo_count":p,"mean_image_angle_deg":y,**five})
    frame=pd.DataFrame(records)
    monkeypatch.setattr(a,"N_PHOTOS",int(frame.original_photo_count.sum()))
    return frame


def test_reversible_key_coarsening_preserves_all_original_counts(monkeypatch):
    x=fake(monkeypatch)
    w=a.make_accepted(x)
    assert len(x)==62 and len(w)==55
    assert w.original_photo_count.sum()==x.original_photo_count.sum()
    assert x.groupby("accepted_key").size().gt(1).sum()==7


def test_cluster_robust_per100mm_and_paired_bootstrap(monkeypatch):
    x=fake(monkeypatch)
    w=a.make_accepted(x)
    pair=a.bootstrap_difference(x,w,draws=101,seed=42)
    assert pair["paired_resamples"]==101
    assert len(pair["source_name_minus_accepted_key_BIO12_slope_95_percentile"])==2
    assert pair["number_of_simulated_new_species"]==0
    for d in (x,w):
        fit=a.slope_fit(d)
        assert np.isclose(a.point_slope(d),fit["BIO12_image_angle_deg_per_100mm"])
        assert fit["spatial_and_true_nuclear_phylogeny_not_adjusted"]
        assert len(fit["WCVP_key_cluster_robust_95_interval"])==2


def test_wrong_source_data_does_not_pass_full_census_contract(tmp_path,monkeypatch):
    x=fake(monkeypatch)
    x.loc[0,"original_photo_count"]+=1
    p=tmp_path/"source.csv"
    x.to_csv(p,index=False)
    with pytest.raises(ValueError,match="identity contract"):
        a.read_units(p)


def test_exact_biological_history_boundaries():
    assert a.SEED==20261022
    assert a.EXPOSURES[2]=="chelsa_bio12"
    assert a.N_KEYS==244
    assert a.N_UNITS==251
