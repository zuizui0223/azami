"""Predeclared original 111 public BioSample / modal reference readiness tests."""
import math
import pytest
from analysis.v4 import t2_real111_original241_phylo_preflight as p


def fixture(monkeypatch):
    monkeypatch.setattr(p,"MEASURED",6)
    monkeypatch.setattr(p,"MOREYRA_CANDIDATES",8)
    monkeypatch.setattr(p,"MISSING_SHORT",2)
    monkeypatch.setattr(p,"AUTHOR_LOCI",4)
    monkeypatch.setattr(p,"EXACT_REF_QUORUM",.80)
    taxa=[]
    for i,name in enumerate([
        "Cirsium conspicuum","Cirsium scopulorum","Cirsium vulgare",
        "Cirsium flodmanii","Cirsium italicum","Cirsium tracyi",
        "Cirsium shortone","Cirsium shorttwo"]):
        taxa.append({
           "accepted_key":str(206700000+i),
           "original_source_taxon":name,
           "selected_single_real_run":("SRR"+str(i+1000)) if i<6 else "",
           "selected_single_biosample":("SAMN"+str(i+1000)) if i<6 else "",
           "real120k_depth_measured":"true" if i<6 else "false",
        })
    loci=[]
    for i,count in enumerate([6,5,4,5]):
        name=f"L{i}"
        loci.append({
           "locus":name,
           "n_taxa_Q20_dp3_ge25pct_same_original_reference_copy":str(count),
           "n_taxa_Q20_dp3_ge25pct":str(count+1 if count<6 else 6),
           "n_taxa_Q20_dp3_ge50pct":"4",
           "modal_Q20_dp3_ge25pct_original_ref":f"saff-{name}",
        })
    return taxa,loci


def test_main_quorum_predeclared_and_nodding_preserved(monkeypatch):
    taxa,loci=fixture(monkeypatch)
    observed,missing,candidates,summary=p.preflight(taxa,loci)
    assert len(observed)==6 and len(missing)==2
    assert len(candidates)==3
    assert summary["source_taxon_number_for_gate"]==5
    assert not summary["technical_ML_source_locus_gate_passed"]
    assert set(summary["source_botanical_nodding_author_names_present"])=={
        "Cirsium conspicuum","Cirsium scopulorum"}
    assert set(x["phylogenetic_tip"] for x in observed)=={
        "WCVP"+str(206700000+i) for i in range(6)}
    assert all(x["single_original_reference_copy_is_not_validated_ortholog"]
               for x in candidates)
    assert all(x["not_measured_120k_not_missing_gene"] for x in missing)


def test_no_mixing_original_reference_copy_or_fake_taxon_identity(monkeypatch):
    taxa,loci=fixture(monkeypatch)
    loci[0]["modal_Q20_dp3_ge25pct_original_ref"]="OTHER-L0"
    with pytest.raises(ValueError,match="unambiguous original"):
        p.preflight(taxa,loci)
    taxa,loci=fixture(monkeypatch)
    taxa[0]["selected_single_biosample"]="SAMN1001"
    with pytest.raises(ValueError,match="111 real BioSamples"):
        p.preflight(taxa,loci)
    taxa,loci=fixture(monkeypatch)
    taxa[0]["accepted_key"]="NotWCVP"
    with pytest.raises(ValueError,match="not numeric"):
        p.preflight(taxa,loci)


def test_missing_short_sra_not_retrospectively_selected(monkeypatch):
    taxa,loci=fixture(monkeypatch)
    taxa[-1]["selected_single_real_run"]="SRRFAKE"
    with pytest.raises(ValueError,match="Unavailable 120k"):
        p.preflight(taxa,loci)


def test_source_density_cannot_be_stretched_as_orthology(monkeypatch):
    taxa,loci=fixture(monkeypatch)
    loci[0]["n_taxa_Q20_dp3_ge25pct_same_original_reference_copy"]="7"
    with pytest.raises(ValueError,match="Invalid nuclear Q20"):
        p.preflight(taxa,loci)
    assert p.MASTER_ACCEPTED_KEYS==244
    assert p.ORIGINAL_SOURCE_NAMES==251
    assert p.MIN_EXPECTED_PHYLO_CANDIDATE_LOCI==30
