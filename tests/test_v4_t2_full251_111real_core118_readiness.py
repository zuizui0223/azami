"""Synthetic full-S1-name retention tests; empirical run uses 251 originals."""
import pytest

from analysis.v4 import t2_full251_111real_core118_readiness as mod


def sample_data(monkeypatch):
    values={"N_SOURCE":6,"N_ACCEPTED":5,"N_PHOTOS":6,
            "N_COLLISIONS":1,"N_REAL_NAME_CANDIDATES":4,
            "N_MEASURED":2,"N_LOCI":1,
            "N_MATCHED_CORE":1,"N_MISSING_SHORT":2}
    for key,val in values.items(): monkeypatch.setattr(mod,key,val)
    names=("Cirsium arvense","Cirsium vulgare",
           "Cirsium conspicuum","Cirsium scopulorum",
           "Cirsium brevicaule","Cirsium irumtiense")
    source=[]
    for i,name in enumerate(names):
        source.append({
            "accepted_key":str(min(i+1,5)),
            "S1_original_source_name":name,
            "original_S1_photo_count":"1",
            "source_same_biosample_run_120k_Q20_successes":"1" if i<2 else "0",
            "Q20_120k_read_state":(
                "measured_public_prefix120k_Q20_unverified_orthologs" if i<2 else
                "no120k_due_to_short_public_fastq" if i<4 else
                "no_Moreyra125_SRA_in_source_name_candidate_panel"),
        })
    chosen=[{"accepted_key":str(i+1),"name":names[i],
             "run":f"SRR{i+1}","biosample":f"SAMN{i+1}"} for i in range(2)]
    matrix=[
        {"accepted_key":str(i+1),"author241_locus":"L001",
         "public_run":f"SRR{i+1}","biosample":f"SAMN{i+1}",
         "120k_dp3_breadth":".75" if i==0 else ".25",
         "120k_best_reference":"lett-L001"}
        for i in range(2)]
    loci=[{"locus":"L001","taxa_120k_dp3_ge25pct_same_ref":"2",
           "modal_120k_dp3_ge25pct_ref":"lett-L001"}]
    return source,chosen,matrix,loci


def test_original251_and_two_name_one_accepted_key_preserved(monkeypatch):
    data=sample_data(monkeypatch)
    rows,summary=mod.compile_coverage(*data)
    assert len(rows)==6
    assert summary["original_WCVP_accepted_keys"]==5
    assert summary["original_source_name_units"]==6
    assert summary["S1_units_with_real_Q20_120k_NA_source"]==2
    assert summary["S1_units_short_original_SRA_archives"]==2
    assert summary["S1_units_without_original_Moreyra125_source"]==2
    assert summary["n_original_source_photo_units_dropped"]==0
    assert summary["source_focus"]["Cirsium arvense"][
        "n_original118_same_ref_loci_dp3_ge25pct"]==1
    assert summary["source_focus"]["Cirsium vulgare"][
        "n_original118_same_ref_loci_dp3_ge25pct"]==1
    assert summary["source_focus"]["Cirsium brevicaule"][
        "n_original118_same_ref_loci_dp3_ge25pct"] is None
    assert summary["reference_compatibility_not_ortholog_certificate"]


def test_duplicate_or_wrong_biosample_blocks_full_source_join(monkeypatch):
    source,chosen,matrix,loci=sample_data(monkeypatch)
    chosen[0]["biosample"]="SAMN_STALE"
    with pytest.raises(ValueError,match="same verified SRA"):
        mod.compile_coverage(source,chosen,matrix,loci)


def test_unverified_focal_taxon_lumped_into_other_species_blocks(monkeypatch):
    source,chosen,matrix,loci=sample_data(monkeypatch)
    chosen[0]["accepted_key"]="5"
    with pytest.raises(ValueError):
        mod.compile_coverage(source,chosen,matrix,loci)


def test_80percent_copy_gate_enforces_no_reference_imputation(monkeypatch):
    source,chosen,matrix,loci=sample_data(monkeypatch)
    loci[0]["taxa_120k_dp3_ge25pct_same_ref"]="1"
    with pytest.raises(ValueError,match="118 source loci changed"):
        mod.compile_coverage(source,chosen,matrix,loci)


def test_original_empirical_global_contract():
    assert (mod.N_SOURCE,mod.N_ACCEPTED,mod.N_MEASURED,mod.N_LOCI)==(251,244,111,241)
    assert mod.N_MATCHED_CORE==118
