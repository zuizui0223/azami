"""Likelihood checks for S1 194/191/171/143 continuous photo-angle signal."""
from io import StringIO
import json

from Bio import Phylo
import numpy as np
import pytest

from analysis.v4 import s1_continuous_phylogenetic_signal as sig
from analysis.v4.joint_s1_t2_phylogeny import EXPOSURES


def make_tree(n=28):
    left=[f"t{i}:0.3" for i in range(n//2)]
    right=[f"t{i}:0.3" for i in range(n//2,n)]
    return Phylo.read(StringIO("(("+",".join(left)+"):5.0,("+
                                   ",".join(right)+"):5.0);"),"newick")


def test_pagel_profile_invariant_to_angle_origin_and_units():
    tree=make_tree()
    tips=sorted(t.name for t in tree.get_terminals())
    rows=[{"tip":name} for name in tips]
    C=sig.covariance_for(tree,rows)
    x=np.ones((len(rows),1))
    y=np.array([55.+(0 if int(name[1:])<14 else 25)+int(name[1:])*.04
                for name in tips])
    fitted=sig.LambdaProfile(C,x,nugget=.05)
    result=fitted.fit(y)
    translated=fitted.fit(y+27)
    multiplied=fitted.fit(y*1.7)
    assert result["n_taxa"]==28
    assert 0<=result["lambda_ML"]<=1
    assert result["likelihood_ratio_against_lambda0"]>=0
    assert result["lambda_ML"]==pytest.approx(translated["lambda_ML"],abs=2e-4)
    assert result["lambda_ML"]==pytest.approx(multiplied["lambda_ML"],abs=2e-4)
    assert result["lambda_ML"]>0.1


def test_lrt_null_bootstrap_is_an_actual_refitted_null_not_chisq():
    tree=make_tree()
    tips=sorted(t.name for t in tree.get_terminals())
    C=sig.covariance_for(tree,[{"tip":name} for name in tips])
    rng=np.random.default_rng(144)
    climate=rng.normal(size=(len(tips),5))
    X=np.column_stack([np.ones(len(tips)),climate])
    model=sig.LambdaProfile(C,X,nugget=.25)
    y=45+climate[:,0]*6+rng.normal(size=len(tips))
    result=model.bootstrap_lambda0(y,reps=25,seed=13)
    assert result["requested_draws"]==25
    assert 0<result["boundary_calibrated_bootstrap_p"]<=1
    assert result["null"]=="parametric_bootstrap_phylogenetically_independent_tips_given_original_five_climates"
    assert result["null_lrt_q50_q95_q99"][0]<=result["null_lrt_q50_q95_q99"][1]


def test_five_climate_model_cannot_be_replaced_by_binary_U_D_or_new_predictor():
    rng=np.random.default_rng(16)
    rows=[]
    for i in range(12):
        climate=(rng.normal(size=5)*[3.,12.,100.,30.,2.]).tolist()
        rows.append({"tip":f"Cirsium_{i}","S1_five_climate_means":
                    dict(zip(EXPOSURES,climate))})
    ordered,X=sig.model_design(rows,True)
    assert X.shape==(12,6)
    assert set(ordered[0]["S1_five_climate_means"])==set(EXPOSURES)
    _,X0=sig.model_design(rows,False)
    assert X0.shape==(12,1)
    rows[0]["S1_five_climate_means"]["chelsa_bio15"]=1e6
    with pytest.raises(ValueError,match="Constant"):
        for r in rows:
            r["S1_five_climate_means"]["chelsa_bio01"]=1.
        sig.model_design(rows,True)


def test_s1_194_original_file_needs_climate_but_does_not_download_any_images(tmp_path):
    rows=[]
    for i in range(194):
        tip=(["Cirsium_parryi","Cirsium_setidens","Cirsium_yezoense"][i]
             if i<3 else f"Cirsium_species_{i}")
        rows.append({
            "tip":tip,"accepted_key":f"k{i}",
            "n_eligible_images":1 if i<193 else 169912-193,
            "source_concept_collision":i<3,
            "image_mean_deg":50+(.1*i),"image_median_deg":48+(.1*i),
            "S1_five_climate_means":{
              "chelsa_bio01":10.,"chelsa_bio04":200.,
              "chelsa_bio12":900.,"chelsa_bio15":40.,
              "chelsa_rsds_mean":9.}
        })
    raw={"source_audit":{"source_sha256":sig.SHA,
        "native_observation_source_authenticated":True,
        "nonindependent_52_tree_ensemble":True},
         "summary":{"n_tree_taxon_keys_with_original_continuous_orientation":194,
                    "n_original_images_in_S1":170066,
                    "S1_original_five_exposures":list(EXPOSURES)},
         "source_cohort_194_species":rows}
    path=tmp_path/"S1_CONTINUOUS_194.json"
    path.write_text(json.dumps(raw))
    assert len(sig.verify_source(path))==194
    assert {k:len(v) for k,v in sig.taxon_sets(rows).items()}["collision_free191"]==191
    rows[3].pop("S1_five_climate_means")
    path.write_text(json.dumps(raw))
    with pytest.raises(ValueError,match="five climate"):
        sig.verify_source(path)
