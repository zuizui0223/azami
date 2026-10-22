"""Synthetic paired-run and fixed author241 overlap safeguards."""
from collections import defaultdict
import json
import pytest
from analysis.v4 import t2_113taxa_author241_shared_locus_screen as m


def make_mock(monkeypatch):
    monkeypatch.setattr(m,"N_LOCI",3)
    monkeypatch.setattr(m,"N_RUNS",9)
    monkeypatch.setattr(m,"N_TAXA",7)
    names=sorted(m.TARGET_NAMES)
    measured={}
    paired=[]
    for idx,name in enumerate(names):
        for rep in range(2 if idx<2 else 1):
            run=f"SRR{idx+1}0{rep}"
            key=str(206700000+idx)
            locus=[]
            for j in range(3):
                # all specimens share L0 and L1 at 25% depth3;
                # species idx6 fails the L2 locus
                dp3=.70 if j<2 else (.1 if idx==6 else .7)
                locus.append({
                   "locus":f"L{j}","dp3_breadth":dp3,
                   "dp1_breadth":.85,
                   "reference_choice_maximizes_coverage_not_orthology":True,
                   "best_original_reference":f"saff-L{j}",
                })
            measured[run]={
                "run":run,"biosample":f"SAMN{idx}0{rep}","accepted_key":key,
                "S1_source_taxon":name,
                "n_loci_dp3_ge50pct":sum(z["dp3_breadth"]>=.5 for z in locus),
                "loci":locus,
            }
            paired.append({
               "run":run,"biosample":f"SAMN{idx}0{rep}",
               "accepted_key":key,"source_30k_to_120k_valid":"True",
               "first30k_n_loci_dp3_ge50pct":str(measured[run]["n_loci_dp3_ge50pct"]),
               "first120k_n_loci_dp3_ge50pct":str(100+idx),
            })
    return measured,paired


def test_all113_style_cohort_is_taxon_de_duplicated_and_one_run_per_tip(monkeypatch):
    actual,paired=make_mock(monkeypatch)
    loci,chosen,summary=m.compute(actual,paired)
    assert len(loci)==3 and len(chosen)==7
    assert summary["original_30k_loci_all7_ge25pct"]==2
    assert summary["original_30k_loci_at_least6of7_ge25pct"]==3
    assert summary["original_30k_loci_all7_ge25pct_same_ref"]==2
    assert summary["n_seven_deep_source_taxa"]==7
    assert all(c["source_name_candidate_voucher_unverified"] for c in chosen)
    assert summary["all7_real120k_locus_overlap_not_calculated"]
    assert not summary["nuclear_tree_fit_allowed_from_breadth_only"]


def test_original_run_to_biosample_cannot_change_between_30k_and_120k(monkeypatch):
    actual,paired=make_mock(monkeypatch)
    paired[2]["biosample"]="SAMN_WRONG"
    with pytest.raises(ValueError,match="BioSample"):
        m.compute(actual,paired)


def test_locus_presence_cannot_be_extrapolated_from_120k_aggregate_counts(monkeypatch):
    actual,paired=make_mock(monkeypatch)
    # Same-run 30k named-locus matrix no longer agrees with validated
    # independent 30k summary after deliberate tampering.
    paired[0]["first30k_n_loci_dp3_ge50pct"]="0"
    with pytest.raises(ValueError,match="contradict"):
        m.compute(actual,paired)


def test_unknown_ref_copy_is_not_certified_as_ortholog(monkeypatch):
    actual,paired=make_mock(monkeypatch)
    next(iter(actual.values()))["loci"][0]["reference_choice_maximizes_coverage_not_orthology"]=False
    with pytest.raises(ValueError,match="safeguard"):
        m.compute(actual,paired)


def test_full_run_set_and_source_species_are_necessary(monkeypatch):
    actual,paired=make_mock(monkeypatch)
    del actual[next(iter(actual))]
    with pytest.raises(ValueError,match="Expected full"):
        m.compute(actual,paired)


def test_missing_any_one_of_ten_author241_empirical_shards_fails(tmp_path):
    (tmp_path/"shard_0.json").write_text(json.dumps({"shard_id":0}))
    with pytest.raises(ValueError,match="Ten original"):
        m.read_perrun(tmp_path)
