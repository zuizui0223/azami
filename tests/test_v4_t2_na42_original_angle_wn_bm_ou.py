"""Synthetic-only correctness gates for NA42 independent published nuclear benchmark."""
from io import StringIO
from pathlib import Path

from Bio import Phylo
import numpy as np
import pytest

import analysis.v4.t2_na42_original_angle_wn_bm_ou as mod
from analysis.v4.joint_s1_t2_phylogeny import EXPOSURES


def test_true_nuclear_model_receipt_preserves_full_s1_scope(monkeypatch):
    rng=np.random.default_rng(22091)
    n=42
    cl=rng.normal(size=(n,5))
    rows=[]
    for j in range(n):
        rows.append({
            "tip":f"Cirsium_t{j:02d}",
            "image_mean_deg":float(60+2*cl[j,2]+rng.normal()),
            "image_median_deg":float(61+2*cl[j,2]+rng.normal()),
            "S1_five_climate_means":{
                k:float(cl[j,i]) for i,k in enumerate(EXPOSURES)}
        })
    newick="("+",".join(f"Cirsium_t{j:02d}:1.0" for j in range(n))+");"
    tree=Phylo.read(StringIO(newick),"newick")
    monkeypatch.setattr(mod,"source_rows",lambda path:rows)
    monkeypatch.setattr(mod,"author_nuclear_tree",lambda:tree)
    out=mod.run(Path("unused-source-stub"))
    assert out["matched_independent_nuclear_taxa"]==42
    assert out["source_original_S1_taxa"]==244
    assert len(out["fits"])==4
    assert all(r["best_model_AICc"] in {"WN","BM","OU1"} for r in out["fits"])
    assert out["claim_stop"]["not_global_Cirsium_nuclear_tree"]
    assert out["all_published_biological_nodding_scores_used_as_outcomes"] is False


def test_source_42_tip_contract_blocks_missing_taxon(monkeypatch):
    monkeypatch.setattr(mod,"source_rows",lambda path:[{"tip":"Cirsium_test"}])
    monkeypatch.setattr(mod,"author_nuclear_tree",lambda:None)
    with pytest.raises(ValueError,match="42-tip"):
        mod.prepare(Path("no-source-stub"))
