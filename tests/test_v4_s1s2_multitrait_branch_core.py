"""No fake biodiversity evidence: synthetic unit tests of scientific method guards."""
from __future__ import annotations
from io import StringIO

import numpy as np
import pytest
from Bio import Phylo

from analysis.v4 import s1s2_multitrait_branch_core as m

TRAITS=(
    ("orientation_image_vertical_angle","presentation",True),
    ("corolla_lab_lightness","flower_colour",True),
    ("corolla_lab_chroma","flower_colour",True),
    ("corolla_hue_sin","flower_colour",True),
    ("corolla_hue_cos","flower_colour",True),
    ("visible_floret_fraction","floral_display",False),
    ("corolla_white_pixel_fraction","colour_composition",False),
    ("corolla_redmagenta_pixel_fraction","colour_composition",False),
    ("corolla_purple_pixel_fraction","colour_composition",False),
    ("corolla_yellow_pixel_fraction","colour_composition",False),
    ("capitulum_outline_aspect_ratio","whole_head_shape",True),
    ("capitulum_outline_circularity","whole_head_shape",True),
    ("capitulum_outline_solidity","whole_head_shape",True),
    ("capitulum_width_profile_cv","whole_head_shape",True),
    ("involucre_length_width_ratio","involucre_body",True),
    ("involucre_apical_taper_ratio","involucre_body",True),
    ("involucre_basal_taper_ratio","involucre_body",True),
    ("bract_projection_roughness","phyllary_contour_proxies",True),
    ("bract_projection_p95","phyllary_contour_proxies",True),
    ("bract_projection_maximum","phyllary_contour_proxies",True),
    ("bract_spread_fraction","phyllary_contour_proxies",True),
    ("bract_projection_peak_density","phyllary_contour_proxies",True),
    ("bract_projection_asymmetry","phyllary_contour_proxies",True),
    ("involucre_surface_edge_density","involucre_surface_texture",True),
    ("involucre_surface_lbp_entropy","involucre_surface_texture",True),
    ("involucre_surface_high_frequency_energy","involucre_surface_texture",True),
    ("involucre_surface_specular_fraction","involucre_specularity",True),
)


def fixture_registry():
    roles={
     "involucre_surface_texture":"technical_validation_only",
     "involucre_specularity":"technical_validation_only_not_stickiness",
    }
    groups=sorted({group for _,group,_ in TRAITS})
    return {
      "endpoints":[{"id":name,"mapping_group":group,
              "source_v2_status":"measured" if measured else "unexecuted_no_measurement"}
              for name,group,measured in TRAITS],
      "dependent_measurement_groups":[{"group":group,"role":roles.get(group,"primary_image_geometry"),
              "n_endpoints":sum(g==group for _,g,_ in TRAITS)} for group in groups],
    }


def toy_tree():
    t="(((T1:1,T2:1):1,(T3:1,T4:1):1):1,((T5:1,T6:1):1,((T7:1,T8:1):1,(T9:1,T10:1):1):1):1);"
    return Phylo.read(StringIO(t),"newick")


def sample():
    n=10
    T=np.linspace(-2,2,n)
    E=np.column_stack([T+0.03*np.sin(T),-T+.03*np.cos(T),np.sin(T*2),np.cos(T*2)])
    return T,E


def test_full_measured_inventory_preserves_unmeasured_and_technical():
    r=m.selection_registry(fixture_registry())
    assert (len(r["measured"]),len(r["unmeasured"]),len(r["technical_only"]),len(r["inferential"]))==(22,5,4,18)
    x=m.tests_from_registry(r)
    assert "block:flower_colour_joint" in x
    assert "block:corolla_hue_circular" in x
    assert "endpoint:corolla_hue_sin" not in x
    assert "endpoint:involucre_surface_specular_fraction" not in x
    assert "endpoint:bract_projection_p95" in x
    assert set(k for v in x.values() for k in v)==set(r["inferential"])


def test_fail_closed_on_changed_floral_ontology():
    r=fixture_registry()
    r["endpoints"][3]["source_v2_status"]="unexecuted_no_measurement"
    with pytest.raises(ValueError,match="22 / unavailable 5"):
        m.selection_registry(r)


def test_brownian_operator_constant_and_reconstructed_branch_filter():
    t=toy_tree();names=[f"T{i}" for i in range(1,11)]
    D,mask=m.branch_operator(t,names)
    assert D.shape[1]==10 and len(mask)==len(D)
    assert mask.sum()>=4 and (~mask).sum()==10
    assert np.allclose(D @ np.ones(10),0,atol=1e-8)
    with pytest.raises(ValueError,match="eight"):
        m.branch_operator(t,names[:6])
    with pytest.raises(ValueError,match="do not match"):
        m.branch_operator(t,names[:-1]+["T99"])


def test_continuous_climate_signed_and_internal_only_explicitly_separate():
    T,E=sample();D,internal=m.branch_operator(toy_tree(),[f"T{i}" for i in range(1,11)])
    all_r=m.alignment_test(D,E,T,nperm=999,seed=100)
    int_r=m.alignment_test(D[internal],E,T,nperm=999,seed=100)
    assert all_r["chelsa_bio18"]["score"]>0
    assert int_r["chelsa_bio18"]["score"]>0
    for r in (all_r,int_r):
        assert set(r)==set(m.FOUR)
        assert all(x["score_type"]=="signed_scalar" for x in r.values())
        assert all(0.001<=x["conditional_tip_shuffle_fraction"]<=1 for x in r.values())


def test_flower_hue_vector_joint_and_nominal_score_is_not_directional():
    T,E=sample();D,_=m.branch_operator(toy_tree(),[f"T{i}" for i in range(1,11)])
    hue=np.column_stack((np.sin(T),np.cos(T)))
    result=m.alignment_test(D,E,hue,nperm=999,seed=200)
    assert all(x["score_type"]=="multivariate_nondirectional_norm" for x in result.values())
    assert all(x["score"]>=0 for x in result.values())
    assert result["chelsa_bio18"]["score"]<1.5


def test_module_coshift_and_null_conditional_on_one_tree():
    T,_=sample();D,internal=m.branch_operator(toy_tree(),[f"T{i}" for i in range(1,11)])
    a=np.column_stack([T,np.sin(T)])
    out=m.co_shift_test(D[internal],a,a.copy(),nperm=999,seed=300)
    assert np.isclose(out["co_shift_cosine"],1,atol=1e-8)
    assert 0.001<=out["conditional_tip_shuffle_fraction"]<=1
    with pytest.raises(ValueError,match="At least 999"):
        m.co_shift_test(D,a,a,nperm=9)


def test_holm_fwer_and_seed_determinism():
    p=m.holm_family([("A",0.01),("B",0.04),("C",0.002)])
    assert np.isclose(p["C"],0.006)
    assert np.isclose(p["A"],0.02)
    assert np.isclose(p["B"],0.04)
    assert m.stable_seed("colour+shape") == m.stable_seed("colour+shape")
    assert m.stable_seed("colour+shape") != m.stable_seed("shape+colour")


def test_no_8_tips_no_constant_env_no_optimistic_permutation():
    D,_=m.branch_operator(toy_tree(),[f"T{i}" for i in range(1,11)])
    t,e=sample()
    with pytest.raises(ValueError,match="At least 999"):
        m.alignment_test(D,e,t,nperm=99)
    with pytest.raises(ValueError,match="constant"):
        m.alignment_test(D,np.ones((10,4)),t,nperm=999)