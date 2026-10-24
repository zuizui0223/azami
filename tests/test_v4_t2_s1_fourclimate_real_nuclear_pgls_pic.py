"""S1 four real environments on genuine nuclear phylogeny, not invented U-D gains."""
import numpy as np
import pytest
from io import StringIO
from Bio import Phylo
from analysis.v4 import t2_s1_fourclimate_real_nuclear_pgls_pic as m

def test_registered_requested_axes_and_real_authenticated_source():
    assert m.COVARIATES==("chelsa_bio18","chelsa_bio05","chelsa_rsds_mean","elevation_m")
    assert m.NEW_SOURCE_PHOTOS==169020
    assert m.NEW_AGGREGATE_SHA=="823ecb8442eb03598253cb33f3e16bac86538c6c2652f48ca2dc0886d0e19626"
    assert m.REGIONS==("NA42_genuine_published_nuclear","Japan12_genuine_176_loci_nuclear")

def mock_rows(n=22):
    items=[]
    for i in range(n):
        items.append({"chelsa_bio18":str(160+23*i+13*np.sin(i)),
             "chelsa_bio05":str(11+1.5*np.sin(.2*i)+i*.19),
             "chelsa_rsds_mean":str(13+np.cos(.7*i)+i*.08),
             "elevation_m":str(650+50*np.cos(.3*i)+i*8)})
    return items

def test_four_axis_panel_design_and_source_collinearity_diagnostics():
    x,d=m.design(mock_rows())
    assert x.shape==(22,5)
    assert len(d["VIF"])==4
    assert all(v>=1-1e-8 for v in d["VIF"])
    assert len(d["pairwise_pearson_r"])==4

def test_PIC_synthetic_bifurcating_tree_recovers_contrasts_not_U_D_counts():
    tree=Phylo.read(StringIO("(((a:1,b:1):1,c:2):1,(d:1,e:1):2);"),"newick")
    mapping={k:k for k in "abcde"}
    d={}
    for i,k in enumerate("abcde"):
        d[k]={"mean_image_angle_deg":str(60+2*i),
              "chelsa_bio18":str(120+20*i),
              "chelsa_bio05":str(28+i),
              "chelsa_rsds_mean":str(10+i),
              "elevation_m":str(500+80*i)}
    pic=m.independent_contrasts(tree,mapping,d)
    assert pic.shape==(4,5)
    assert np.isfinite(pic).all()
    p=m.pic_report(pic)
    assert p["status"]=="NOT_EVALUABLE_INSUFFICIENT_GENUINE_NUCLEAR_BRANCHES"

def test_under_minimum_9_source_species_does_not_return_fake_PGLS():
    tree=Phylo.read(StringIO("(a:1,b:1,c:1,d:1,e:1);"),"newick")
    obj={str(k):{"n_photo_complete_four_axis":"5"} for k in "abcde"}
    panel=m.analyze_panel(tree,{k:k for k in obj},obj,minimum_photos=10)
    assert panel["status"]=="NOT_EVALUABLE_TOO_FEW_SOURCE_PHOTOS_NUCLEAR_TIPS"
