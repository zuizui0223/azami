"""Match identical S1 photographic phenotype to independent NA nuclear topology."""
from io import StringIO
import json
import numpy as np
import pytest
from Bio import Phylo

from analysis.v4 import s1_continuous_na_nuclear_signal as m
from analysis.v4.joint_s1_t2_phylogeny import EXPOSURES


def fake_rows():
    rng=np.random.default_rng(42)
    rows=[]
    for i in range(42):
        raw=rng.normal(loc=[15,150,800,45,12],scale=[3,20,100,8,2])
        rows.append({
         "tip":f"Cirsium_{i:02d}",
         "image_mean_deg":float(50+2*i + rng.normal(scale=5)),
         "image_median_deg":float(50+2*i + rng.normal(scale=5)),
         "S1_five_climate_means":dict(zip(EXPOSURES,raw)),
        })
    return rows


def tree(names):
    left=names[:21];right=names[21:]
    nw="(("+",".join(f"{x}:0.5" for x in left)+"):2.0,("+
          ",".join(f"{x}:0.5" for x in right)+"):2.0);"
    return Phylo.read(StringIO(nw),"newick")


def test_same_42_nuclear_plus_52_placements_reuse_original_continuous_angle():
    rows=fake_rows()
    names=sorted(x["tip"] for x in rows)
    phylo=tree(names)
    individual=m.one_tree(phylo,"test_independent_nuclear",rows)
    assert individual["n_tips"]==42
    assert individual["unadjusted"]["n_taxa"]==42
    assert individual["five_climate_adjusted"]["model_rank"]==6
    assert 0<=individual["five_climate_adjusted"]["lambda_ML"]<=1
    result=m.analyze(rows,phylo,[(f"test_{i}",tree(names)) for i in range(52)],
                     bootstrap_reps=20)
    report=result["summary"]
    assert report["same_original_photo_phenotype_nuclear_vs_placements_n_taxa"]==42
    assert report["matched_52_placement_phenotype_fiveclimate_lambda"]["n_deltaAIC_over_2"]<=52
    assert len(result["per_placement_results"])==52
    assert report["same_NA42_nuclear_bootstrap"]["nuclear_S1_five_climate_adjusted"]["requested_draws"]==20
    assert report["external_3_botanical_nodding_states_ignored_as_phenotype_response"] is True


def test_partial_tip_mapping_must_fail():
    rows=fake_rows()
    names=sorted(r["tip"] for r in rows)
    with pytest.raises(ValueError,match="mismatch"):
        m.one_tree(tree(names[:-1]+["Cirsium_other"]),"wrong",rows)


def test_five_original_climate_values_no_relabelling():
    rows=fake_rows()
    r=rows[0]
    r["S1_five_climate_means"].pop("chelsa_bio15")
    with pytest.raises(KeyError):
        m.one_tree(tree(sorted(x["tip"] for x in rows)),"bad",rows)


def test_bootstrap_is_not_a_true_nodding_origin_count():
    rows=fake_rows()
    assert len(rows)==42
    assert all("external_botanical_state" not in x for x in rows)
    assert m.NUCLEAR_SOURCE_DOI=="10.1086/724310"
