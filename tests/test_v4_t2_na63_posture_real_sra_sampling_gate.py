"""Offline guards for published botanical HeadPosition × real 125 SRA map overlap."""
import pytest
from analysis.v4 import t2_na63_posture_real_sra_sampling_gate as g


def synthetic():
    records=[]
    original={}
    score={}
    sources=[(f"Cirsium_erect{i:02d}",0) for i in range(24)]
    sources += [("Cirsium_conspicuum",1),("Cirsium_scopulorum",1)]
    assert len(sources)==26
    idx=0
    for name,state in sources:
        score[name]=state
        n=2 if name=="Cirsium_scopulorum" else 1
        if name=="Cirsium_erect00":n=6
        for j in range(n):
            run=f"SRR9{idx:07d}";idx+=1
            loci=(589 if j==0 else 965) if name=="Cirsium_scopulorum" else 900
            records.append({"run":run,"pairs":30000,
                "reference_loci_at_least_one_alignment":loci})
            original[run]=[{"tree_code":name.replace("_"," ")}]
    assert idx==32
    for i in range(93):
        run=f"SRR8{i:07d}"
        records.append({"run":run,"pairs":30000,
            "reference_loci_at_least_one_alignment":850})
        original[run]=[{"tree_code":f"Cirsium unscored{i:03d}"}]
    receipt={
       "read_results":{"mapped_run_count":125,"failed_run_count":0},
       "workflow_conclusion":"success",
       "per_original_public_run":records,
       "source_moreyra_name_matched_taxon_keys":113,
    }
    return receipt,original,score


def test_original_author_literature_names_and_replicates_kept_distinct():
    assert g.normalize_author_name("Cirsium scopulorum 2")=="Cirsium_scopulorum"
    assert g.normalize_author_name("Cirsium japonicum var. fukienense 2")==\
            "Cirsium_japonicum_var._fukienense"
    assert g.normalize_author_name("Cirsium japonicum")!=\
           g.normalize_author_name("Cirsium japonicum var. fukienense")
    assert g.TRAIT_BLOB_SHA=="fa127f0587b676e6e004ea121424d26069586bcb"


def test_strict_nodding_two_species_vs_24_erect_is_a_read_mapping_gate_only():
    src,orig,states=synthetic()
    d=g.evaluate(src,orig,states)
    assert d["n_original_na_taxa_with_exact_named_S1_Moreyra_mapped_runs"]==26
    assert d["n_original_na_nodding_taxa_with_exact_named_runs"]==2
    assert d["n_exact_mapped_runs_with_author_botanical_state"]==32
    assert d["n_runs_nodding"]==3
    assert not d.get("independent_ancestral_event_count_estimated",False)
    assert d["missing_author_scored_nodding_taxa_are_not_reclassified_as_erect"]
    scop=next(x for x in d["species_details"] if x["taxon"]=="Cirsium_scopulorum")
    assert scop["reference_locus_min"]==589
    assert scop["reference_locus_max"]==965


def test_ambiguous_same_run_taxon_name_blocks_matching():
    src,orig,states=synthetic()
    orig[src["per_original_public_run"][0]["run"]].append(
       {"tree_code":"Cirsium totally_different"})
    with pytest.raises(ValueError,match="ambiguous"):
        g.evaluate(src,orig,states)


def test_unverified_or_incomplete_original_read_data_cannot_pass():
    src,orig,states=synthetic()
    src["read_results"]["mapped_run_count"]=124
    with pytest.raises(ValueError,match="125-run empirical"):
        g.evaluate(src,orig,states)


def test_unscored_ambiguous_SRA_is_retained_as_conflict_not_falsely_assigned():
    src,orig,states=synthetic()
    last=src["per_original_public_run"][-1]["run"]
    orig[last].append({"tree_code":"Cirsium unscored_different"})
    result=g.evaluate(src,orig,states)
    assert result["n_unassignable_source_conflict_runs"]==1
    assert result["multiple_source_tree_code_conflict_runs"][0]["run"]==last
    assert result["multiple_source_tree_code_conflict_runs"][0]["read_mapping_success_preserved"]
    assert result["n_original_na_nodding_taxa_with_exact_named_runs"]==2
