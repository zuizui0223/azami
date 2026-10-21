"""Tests for same-source/phylogeny T2 five-variable sharp-null diagnostic."""
import json
from copy import deepcopy
from io import StringIO
from pathlib import Path

from Bio import Phylo
import numpy as np
import pytest

from analysis.v4 import t2_source_blocked_five_climate as m
from analysis.v4.joint_s1_t2_phylogeny import EXPOSURES

def fake_rows():
    np_rng=np.random.default_rng(51)
    rows=[]
    for i in range(46):
        na=i<42
        source=m.NA if na else m.EA
        nodding=(i in (0,9,25,42,43,44))
        vals=np_rng.normal(loc=0,scale=1,size=5)
        vals += np.array([2, 4, 900, 45, 8])
        rows.append({
            "accepted_key":f"k{i}",
            "source_taxon":f"Cirsium t{i:02}",
            "join_status":"exact_S1_and_tree",
            "s1_phylogeny_194_key_match":True,
            "s1_observations":20,
            "state":"nodding" if nodding else "erect",
            "source":source,
            "s1_five_climate_centroids":dict(zip(EXPOSURES,vals.tolist()))
        })
    return rows


def test_46_42_explicit_source_counts_and_complete_fixed_five():
    rows=fake_rows()
    whole=m.cohort(rows)
    na=m.cohort([r for r in rows if r["source"]==m.NA])
    assert whole["n"]==46 and whole["n_nodding"]==6
    assert whole["NA"]=={"n":42,"D":3}
    assert whole["EA"]=={"n":4,"D":3}
    assert whole["X"].shape==(46,7)
    assert na["X"].shape==(42,6), "Never retain all-zero extra source coefficient"
    assert na["n_nodding"]==3


def test_exact_label_universe_conditioned_on_source_frequency_and_p():
    full=m.cohort(fake_rows())
    labels,count=m.exact_source_maps(full)
    assert count==45920
    assert np.all(labels[:,full["lane"]==0].sum(axis=1)==3)
    assert np.all(labels[:,full["lane"]==1].sum(axis=1)==3)
    outcome=m.score(full,labels)
    assert outcome["null_maps"]==45920
    assert outcome["all_five_covariates_in_model"]==list(EXPOSURES)
    for raw,adjusted in zip(outcome["source_stratified_exact_sharp_null_two_sided_p"],
                            outcome["two_axis_maxT_source_stratified_exact_p"]):
        assert 0<raw<=adjusted<=1


def test_original_north_american_42_tree_same_tip_and_regional_joins():
    rows=fake_rows()
    na=m.cohort([r for r in rows if r["source"]==m.NA])
    labels,n=m.exact_source_maps(na)
    assert n==11480
    newick="("+",".join(f"{z}:1" for z in na["names"])+");"
    tree=Phylo.read(StringIO(newick),"newick")
    covariance=m.tree_covariance_for(tree,na["names"])
    outcome=m.score(na,labels,covariance=covariance)
    assert outcome["n_taxa"]==42
    assert outcome["model"].startswith("Brownian")
    assert outcome["null_maps"]==11480
    with pytest.raises(ValueError,match="absent from tree"):
        m.tree_covariance_for(tree,na["names"]+["Cirsium_missing"])


def test_source_concepts_not_included_without_five_frozen_covariates():
    receipt=json.loads(m.RECEIPT.read_text(encoding="utf-8"))
    assert receipt["results"]["model_scope"]["orientation_eligible_accepted_keys"]==244
    assert receipt["results"]["T2_PGLS_52_exploratory"]["n_nodding"]==6
    data={"status":"empirical_full_244_native_join_source_sha_verified",
          "fixed_covariates":list(EXPOSURES),
          "model_scope":{"orientation_eligible_accepted_keys":244,
             "orientation_eligible_observations":170066,
             "placement_tree_covered_keys":194,
             "placement_tree_covered_observations":169912,
             "grafted_placement_scenarios":52},
          "provenance":{"S1_model_source_sha256":receipt["source_provenance"]["native_sha256"],
             "EAzami_source_git_blob":receipt["source_provenance"]["EAzami_fdt4_git_blob"],
             "NA_source_git_blob":receipt["source_provenance"]["NA_head_position_git_blob"]},
          "T2_independent_state_inventory":{"taxon_source_audit":fake_rows()}}
    rows=m.authenticate_input(data,receipt)
    assert len(rows)==46
    data["fixed_covariates"]=["chelsa_bio12"]
    with pytest.raises(ValueError,match="contract changed"):
        m.authenticate_input(data,receipt)
    data["fixed_covariates"]=list(EXPOSURES)
    data["T2_independent_state_inventory"]["taxon_source_audit"][0]["s1_five_climate_centroids"].pop("chelsa_bio12")
    with pytest.raises(ValueError,match="five-climate"):
        m.authenticate_input(data,receipt)
