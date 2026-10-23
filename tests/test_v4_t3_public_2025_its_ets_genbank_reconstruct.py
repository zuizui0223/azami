"""Published 163 accession+154 name cohort; no paralog-free nuclear claims."""
from copy import deepcopy
from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from analysis.v4 import t3_public_2025_its_ets_genbank_reconstruct as m
import pytest

REAL=Path(__file__).resolve().parents[1]/"analysis"/"v4"/"results"/"20261023_T3_TAIWANIA2025_ORIGINAL163_ITS_ETS_GENBANK_ACCESSIONS.csv"

def gb(name="Cremanthodium campanulatum",marker="ITS",acc="KY979098"):
    result=SeqRecord(Seq("ACTG"*130),id=acc+".1",
          description=name+" "+("internal" if marker=="ITS" else "external")+" transcribed spacer, rDNA")
    result.annotations["molecule_type"]="DNA"
    result.annotations["organism"]=name
    return result

def test_original_2025_published_163_rows_retain_nine_duplicate_taxon_entries():
    data=m.read_csv(REAL)
    accession=m.validate_original(data)
    assert len(data)==163
    assert len({r["original_published_scientific_name"] for r in data})==154
    assert len(accession)>=190
    assert sum(bool(r["ETS_GenBank_accession"]) for r in data)==90
    assert {r["original_published_scientific_name"] for r in data if r["ITS_GenBank_accession"]=="KY979098"}=={"Cremanthodium campanulatum"}
    assert "Parasenecio roborowskii" in {r["original_published_scientific_name"] for r in data}
    assert sum(r["original_published_scientific_name"]=="Parasenecio roborowskii" for r in data)==2
    assert sum(r["original_published_scientific_name"]=="Parasenecio cyclotus" for r in data)==2

def test_public_genbank_ITS_not_unlinked_orthologous_ETS_gene():
    r=m.check_record(gb(),"Cremanthodium campanulatum","ITS","KY979098")
    assert r[1] is None and r[0]["n_original_bases"]==520
    assert r[0]["source_sample_same_in_ITS_and_ETS_not_verified"]
    failed=m.check_record(gb(name="Cremanthodium humile"),"Cremanthodium campanulatum","ITS","KY979098")
    assert failed[0] is None and failed[1]=="organism_name_differs_from_published_source_name"
    plastid=gb()
    plastid.description="Cremanthodium campanulatum complete plastid ITS"
    assert m.check_record(plastid,"Cremanthodium campanulatum","ITS","KY979098")[1]=="nonnuclear_title"

def test_original_suppl_invalid_name_or_tampered_provenance_is_rejected():
    original=m.read_csv(REAL)
    x=deepcopy(original)
    x[0]["ITS_GenBank_accession"]="SPAM"
    with pytest.raises(ValueError,match="Invalid GenBank"):
        m.validate_original(x)
    x=deepcopy(original)
    x[2]["ITS_and_ETS_independent_unlinked_gene_loci"]="yes"
    with pytest.raises(ValueError,match="unlinked"):
        m.validate_original(x)

class StubNCBI:
    n_calls=1
    n_remote_bytes=525
    def fetch_batch(self,accessions):
        return [gb(name="Cremanthodium campanulatum",marker="ITS",acc="KY979098")
                for _ in accessions if "KY979098" in accessions][:1]

def test_exact_species_NCBI_allocation_refuses_missing_as_zero_locus():
    data=m.read_csv(REAL)
    # Only the first 3 table records are considered in this deterministic mock.
    scored,taxa,stats=m.collect_original(data,StubNCBI(),taxa_limit=3)
    assert stats["source_rows_queried"]==3
    assert stats["n_unique_taxa_with_exact_published_organism_ITS"]==1
    assert "Cremanthodium campanulatum" in taxa
    assert all(z["ITS_source_status"] in {
        "original_GenBank_exact_species_confirmed",
        "NCBI_unavailable_not_biological_locus_absence"} for z in scored)
    assert stats["original_ITS_ETS_rDNA_spacers_physically_linked_not_independent_nuclear_orthologs"]
    assert not stats["source_failure_batches_not_ITS_absence"]
