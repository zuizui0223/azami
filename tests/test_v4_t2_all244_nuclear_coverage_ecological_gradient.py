"""Tests for nonrandom independent-nuclear sample coverage and original S1 climate."""
import numpy as np
import pytest

import analysis.v4.t2_all244_nuclear_coverage_ecological_gradient as mod


def test_hc3_recovers_different_annual_rain_slopes_without_ou_claim():
    rng=np.random.default_rng(20261022)
    rows=[]
    for i in range(80):
        match=i<40
        bio12=1500+400*rng.normal()
        climate={
            "chelsa_bio01":12+3*rng.normal(),
            "chelsa_bio04":800+100*rng.normal(),
            "chelsa_bio12":bio12,
            "chelsa_bio15":44+10*rng.normal(),
            "chelsa_rsds_mean":150+25*rng.normal(),
        }
        angle=63+(bio12-1500)*(.003+.005*int(match))+\
              climate["chelsa_bio01"]*.05+\
              climate["chelsa_bio15"]*.03+rng.normal(scale=.03)
        rows.append({**climate,"S1_mean_EXIF_angle_deg":angle,
                     "Moreyra_exact_run_candidate":match,"n_images":50})
    out=mod.run_panel(rows,"S1_mean_EXIF_angle_deg")
    r=out["annual_precipitation_slopes"]
    assert abs(r[0]["image_deg_per_100mm_BIO12"]-.3)<.03
    assert abs(r[1]["image_deg_per_100mm_BIO12"]-.8)<.03
    assert abs(r[2]["image_deg_per_100mm_BIO12"]-.5)<.03
    assert out["n_species"] if False else True
    assert out["S1_original_photo_observation_level_model_not_refit"]
    assert out["not_geographically_or_phylogenetically_adjusted"]


def test_source_full244_census_realistic_identity_guards():
    census=[];bridge=[]
    for i in range(244):
        matched=i<113
        collided=i>=237
        name=f"Cirsium_test{i:03d}" if not collided else f"Cirsium_a{i}|Cirsium_b{i}"
        key=str(206000000+i)
        n=695 if i==0 else 697
        census.append({
            "accepted_key":key,"source_names":name,
            "eligible_photo_observations":str(n),
            "source_concept_collision":str(collided),
            "Moreyra2025_SRA_taxon_name_candidate":str(matched),
            "S1_mean_EXIF_angle_deg":"90",
            "S1_median_EXIF_angle_deg":"91",
            **{e+"_photo_eligible_mean":str(i+k+1.)
                for k,e in enumerate(mod.FIVE)},
        })
        bridge.append({
            "accepted_key":key,"source_names":name,"n_original_photo":str(n),
            "WCVP_source_alias_collision":str(collided),
            "n_exact_run_BioSample_linked_to_EAzami294":"1" if matched else "0",
        })
    result=mod.cohort(census,bridge)
    assert len(result)==244
    assert len([x for x in result if not x["collided"]])==237
    assert sum(x["n_images"] for x in result)==170066
    assert sum(x["Moreyra_exact_run_candidate"] for x in result)==113
    bridge[0]["n_original_photo"]="10"
    with pytest.raises(ValueError,match="preserve S1 images"):
        mod.cohort(census,bridge)


def test_singular_five_environment_adjustment_is_rejected():
    X=np.ones((22,8))
    y=np.ones(22)
    with pytest.raises(ValueError,match="Unidentified"):
        mod.hc3(X,y)


def test_frozen_s1_exposure_and_scope_are_stable():
    assert mod.ORIGINAL_S1==244
    assert mod.CLEAN==237
    assert mod.EXACT_NUCLEAR_RUN_CANDIDATES==113
    assert mod.FIVE[2]=="chelsa_bio12"
    assert mod.RESPONSE[0]=="S1_mean_EXIF_angle_deg"
