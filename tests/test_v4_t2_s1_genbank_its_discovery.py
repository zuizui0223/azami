"""Source-level GenBank ITS accessions: strict exact original S1 species provenance."""
from collections import Counter
from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature

from analysis.v4 import t2_s1_genbank_its_discovery as m


SOURCE=Path(__file__).resolve().parents[1]/"analysis"/"v4"/"results"/m.S1_SOURCE_BASENAME


def real():
    return m.read_s1(SOURCE)


def candidate(name="Cirsium arvense",title="Cirsium arvense internal transcribed spacer 1, 5.8S and ITS2"):
    r=SeqRecord(Seq("ACGT"*151),id="AY123456.2",description=title)
    r.annotations["organism"]=name
    r.annotations["molecule_type"]="DNA"
    r.features=[SeqFeature(type="source",qualifiers={
        "organism":[name],"specimen_voucher":["NMNS-X100"]
    })]
    return r


def test_original_source_all251_names_and_244_accepted_keys_not_discarded():
    rows=real()
    assert len(rows)==251
    assert len({r["accepted_key"] for r in rows})==244
    assert sum(int(r["original_photo_count"]) for r in rows)==170066
    assert sum(r["botanical_U_D_state"]=="D" for r in rows)==12
    assert sum(r["botanical_U_D_state"]=="U" for r in rows)==54
    assert sum(r["nuclear_real120k_BioSample"]=="True" for r in rows)==111
    assert all(len(r["original_source_name"].split())==2 for r in rows)


def test_exact_nuclear_ITS_and_voucher_accession_are_admissible_not_orthology():
    r,reason=m.qualifying_its(candidate(),"Cirsium arvense")
    assert not reason
    assert r["GenBank_accession_version"]=="AY123456.2"
    assert r["source_specimen_voucher_or_isolate"]=="NMNS-X100"
    assert r["S1_species_to_GB_record_exact_name"]
    assert not r["ITS_phylogeny_single_copy_orthology_confirmed"]


def test_ncbi_taxonomic_mismatch_and_plastid_never_enter_its_tree():
    x=candidate(name="Cirsium vulgare")
    out,why=m.qualifying_its(x,"Cirsium arvense")
    assert out is None and why=="different_NCBI_source_organism"
    x=candidate()
    x.features[0].qualifiers["organelle"]=["plastid:chloroplast"]
    out,why=m.qualifying_its(x,"Cirsium arvense")
    assert out is None and why=="organelle_feature_excluded"
    x=candidate(title="Cirsium arvense plastid internal transcribed spacer")
    out,why=m.qualifying_its(x,"Cirsium arvense")
    assert out is None and why=="organelle_or_genome_title_excluded"


def test_no_unrecorded_marker_or_unsupported_sequence_accepted():
    x=candidate(title="Cirsium arvense chloroplast matK gene")
    result,why=m.qualifying_its(x,"Cirsium arvense")
    assert result is None
    x=candidate()
    x.seq=Seq("ACGT"*70)
    result,why=m.qualifying_its(x,"Cirsium arvense")
    assert result is None and why=="ITS_length_outside_common_region"
    x=candidate()
    x.seq=Seq("N"*510)
    result,why=m.qualifying_its(x,"Cirsium arvense")
    assert result is None and why=="ITS_too_many_ambiguous_bases"


class EmptyClient:
    calls=4
    def search(self,name):
        return [],0,"exact species query"
    def fetch(self,ids):
        raise AssertionError("Do not fetch if ESearch says zero")


class RealFakeClient:
    calls=4
    def search(self,name):
        return ["123456"],1,"exact source term"
    def fetch(self,ids):
        return [candidate("Cirsium arvense")]


def test_negative_and_taxonomic_not_covered_are_distinct(tmp_path):
    rows=real()[:3]
    ledger,chosen,reasons=m.discover(rows,EmptyClient())
    assert len(ledger)==3 and not chosen
    assert {r["NCBI_ITS_discovery_status"] for r in ledger}=={
        "NCBI_search_no_ITS_candidates"}
    r=m.export(SOURCE,ledger,chosen,reasons,tmp_path,max_species=3,
               request_count=4)
    assert r["scope"]=="S1_PREFIX_PILOT"
    assert not r["all251_queries_completed_without_network_failures"]
    assert (tmp_path/"S1_GENBANK_ITS_EXACT_SPECIES_ACCESSIONS.csv").exists()


def test_retrieval_failure_must_not_be_reported_as_sequence_absence():
    class Broken:
        def search(self,name):
            raise RuntimeError("upstream remote API failed")
    rows=real()[:2]
    ledger,chosen,_=m.discover(rows,Broken())
    assert not chosen and len(ledger)==2
    assert all(z["NCBI_ITS_discovery_status"]==
               "NCBI_access_unresolved_not_biological_absence" for z in ledger)
    assert all("upstream remote API failed" in z["source_access_error"] for z in ledger)


def test_title_based_marker_filter_cannot_claim_full_ITS_only(monkeypatch):
    x=candidate(title="Cirsium arvense ITS2")
    actual,why=m.qualifying_its(x,"Cirsium arvense")
    assert why is None
    assert "ITS" in actual["original_GenBank_definition"]
    assert not actual["ITS_phylogeny_single_copy_orthology_confirmed"]
