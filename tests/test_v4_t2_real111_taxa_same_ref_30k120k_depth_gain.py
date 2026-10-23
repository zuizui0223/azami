"""Synthetic tests of actual same-specimen 30k/120k source-locus occupancy logic."""
import pytest

from analysis.v4 import t2_real111_taxa_same_ref_30k120k_depth_gain as g


def setup(monkeypatch):
    monkeypatch.setattr(g,"N_SOURCE_RUNS",3)
    monkeypatch.setattr(g,"N_SOURCE_KEYS",3)
    monkeypatch.setattr(g,"N_MEASURED",1)
    monkeypatch.setattr(g,"N_LOCI",2)
    monkeypatch.setattr(g,"SHORT",{"SRR2","SRR3"})
    source=[{
        "run":f"SRR{i}","biosample":f"SAMN{i}",
        "accepted_key":f"K{i}","S1_source_taxon":(
           "Cirsium conspicuum" if i==1 else f"Cirsium test{i}")
    } for i in range(1,4)]
    original={}
    newer={}
    for s in source:
        run=s["run"]
        a={
            "run":run,"biosample":s["biosample"],"accepted_key":s["accepted_key"],
        }
        if run=="SRR1":
            base=[
                {"locus":"L1","best_original_reference":"lett-L1",
                 "dp1_breadth":.8,"dp3_breadth":.25,"dp5_breadth":.1,
                 "dp3_callable_sites":75,
                 "reference_choice_maximizes_coverage_not_orthology":True},
                {"locus":"L2","best_original_reference":"saff-L2",
                 "dp1_breadth":.3,"dp3_breadth":.05,"dp5_breadth":0,
                 "dp3_callable_sites":10,
                 "reference_choice_maximizes_coverage_not_orthology":True},
            ]
            new=[{**base[0],"dp1_breadth":1.0,
                  "dp3_breadth":.75,"dp5_breadth":.5,"dp3_callable_sites":225},
                 {**base[1],"best_original_reference":"lett-L2",
                  "dp1_breadth":.75,"dp3_breadth":.55,"dp5_breadth":.25,
                  "dp3_callable_sites":110}]
            original[run]={**a,"loci":base}
            newer[run]={**a,"per_locus_241":new,
                        "status":"real_ENA_120k_same_run_author241_locus_Q20",
                        "n_loci_dp3_ge50pct":2,"n_loci_dp3_ge25pct":2}
        else:
            original[run]={**a,"loci":[]}
            newer[run]={**a,"status":"original_ENA_prefix_below_120k_pairs"}
    return source,original,newer


def test_actual_copy_and_taxon_preserving_same_run_gain(monkeypatch):
    src,old,new=setup(monkeypatch)
    chosen,missing,matrix,loci,summary=g.summarize(src,old,new)
    assert len(chosen)==1 and chosen[0]["run"]=="SRR1"
    assert {x for z in missing for x in z["source_run_ids"].split("|")}=={"SRR2","SRR3"}
    assert len(matrix)==2 and len(loci)==2
    assert matrix[0]["same_reference_copy_at_both_depths"] is True
    assert matrix[1]["same_reference_copy_at_both_depths"] is False
    assert summary["real_111x241_taxon_locus_comparisons"]==2
    assert summary["across_taxa_occupancy"]["30k_dp3_25pct_atleast80pct_taxa"]["n_loci_any_reference"]==1
    assert summary["across_taxa_occupancy"]["120k_dp3_25pct_atleast80pct_taxa"]["n_loci_any_reference"]==2
    assert summary["increase_on_same_111_taxa_and_241_targets"]["120k_dp3_25pct_atleast80pct_taxa"]["n_loci_any_ref_gain_120k_minus30k"]==1
    assert summary["phylogenetic_tree_fitted_from_this_matrix"] is False


def test_missing_taxon_source_is_not_read_as_zero(monkeypatch):
    src,old,new=setup(monkeypatch)
    del new["SRR3"]
    with pytest.raises(ValueError,match="Complete source"):
        g.summarize(src,old,new)


def test_exact_public_run_biosample_mismatch_fails(monkeypatch):
    src,old,new=setup(monkeypatch)
    new["SRR1"]["biosample"]="SAMNnotmatch"
    with pytest.raises(ValueError,match="identity disagreement"):
        g.summarize(src,old,new)


def test_reference_copy_flags_cannot_be_claimed_as_orthology(monkeypatch):
    src,old,new=setup(monkeypatch)
    new["SRR1"]["per_locus_241"][0]["reference_choice_maximizes_coverage_not_orthology"]=False
    with pytest.raises(ValueError,match="Unverified reference-copy"):
        g.summarize(src,old,new)


def test_original_111_by_241_contract():
    assert g.N_SOURCE_RUNS==125 and g.N_SOURCE_KEYS==113
    assert g.N_MEASURED==111 and g.N_LOCI==241
    assert g.TAXON_PCT==(.5,.65,.8,.9)
