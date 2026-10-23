"""Technical independent 111-tip IQ-TREE input stops before falsely inferred adaptation."""
from pathlib import Path

import pytest
from analysis.v4 import t2_real111_q20_technical_ml_input as m


def synthetic(monkeypatch):
    monkeypatch.setattr(m,"SOURCE_TAXA",4)
    monkeypatch.setattr(m,"SOURCE_TARGETS",2)
    monkeypatch.setattr(m,"N_QUORUM",3)
    monkeypatch.setattr(m,"MIN_GENE_LOCI",1)
    monkeypatch.setattr(m,"MIN_DNA_SITES",4)
    monkeypatch.setattr(m,"MIN_INFORMATIVE_SITES",1)
    monkeypatch.setattr(m,"MIN_SOURCE_TAXON_CALLS",3)
    sample=[{"phylogenetic_tip":f"WCVP{i}","taxon_name":name}
            for i,name in enumerate(("Cirsium flodmanii","Cirsium scopulorum",
                                      "Cirsium tracyi","Cirsium vulgare"))]
    genes=[{
        "locus":f"L{i}","original_reference_copy":f"saff-L{i}",
        "n_sites_same_reference_called_at_least80pct_taxa":"4",
        "retain_min100_Q20_positions":"True",
    } for i in range(2)]
    alignment={"WCVP0":"AAAACCCC","WCVP1":"AAAAGCCC",
               "WCVP2":"GGGGCCCC","WCVP3":"GGGGNCCC"}
    summary={
        "original_S1_WCVP_keys_retained":244,
        "original_S1_source_name_operational_units_retained":251,
        "all_original_111_S1_biological_samples":4,
        "n_original_source_candidate_loci":2,
        "predeclared_site_occupancy_taxa":3,
        "technical_ML_alignment_source_gate_passed":True,
        "strict_single_copy_orthology_gate_passed":False,
        "orientation_evolution_rainfall_test_permitted":False,
        "n_loci_with_atleast100_site_positions_in_80pct_taxa":2,
        "n_full111_alignable_positions_called_by_80pct_taxa":8,
    }
    return sample,genes,alignment,summary


def test_real_missing_N_is_kept_and_no_gene_orthology_claim(monkeypatch):
    samples,genes,seqs,src=synthetic(monkeypatch)
    partitions,out=m.validate(src,genes,samples,seqs)
    assert len(partitions)==2
    assert out["all111_source_Q20_80pct_reference_positions"]==8
    assert out["parsimony_informative_Q20_positions"]>=1
    assert out["single_copy_and_homeolog_orthology_verified"] is False
    assert out["adaptive_radiation_or_rain_transition_supported"] is False
    assert len(out["source_independent_NA63_quartet_exact_name_overlap_if_available"])==4


def test_missing_source_bases_cannot_be_filled_by_reference(monkeypatch):
    samples,genes,seqs,src=synthetic(monkeypatch)
    seqs["WCVP3"]=seqs["WCVP3"].replace("N","-")
    with pytest.raises(ValueError,match="source.*gate|Q20|Invalid"):
        m.validate(src,genes,samples,seqs)


def test_missing_original_sample_rejected(monkeypatch):
    samples,genes,seqs,src=synthetic(monkeypatch)
    del seqs["WCVP0"]
    with pytest.raises(ValueError,match="all 111 specimen"):
        m.validate(src,genes,samples,seqs)


def test_fasta_duplicate_taxon_and_reference_imputed_input_fail(tmp_path):
    f=tmp_path/"bad.fa"
    f.write_text(">WCVP1\nACGT\n>WCVP1\nACGT\n")
    with pytest.raises(ValueError,match="one original"):
        m.fasta_read(f)
    f.write_text(">WCVP1\nACGR\n")
    with pytest.raises(ValueError,match="Invalid original"):
        m.fasta_read(f)


def test_nuclear_118_is_not_original_moreyra350():
    assert m.SOURCE_TAXA==111
    assert m.SOURCE_TARGETS==118
    assert m.N_QUORUM==89
    assert m.MIN_SOURCE_TAXON_CALLS==1000
