"""Source identity and nested public-read Q20 depth scaling contracts (synthetic tests)."""
import gzip
import pytest
from analysis.v4 import t2_four_species_30k_120k_nuclear_depth as m


def short_reads(n=8):
    return b"".join((f"@read{i}\nACGT\n+\nIIII\n".encode()) for i in range(n))


def test_partial_gzip_fastq_prefix_is_complete_records_only():
    raw=gzip.compress(short_reads(8))
    reads,truncated=m.fastq_prefix(raw,limit=7)
    assert not truncated
    assert len(reads)==28
    assert all(len(reads[i+1].rstrip())==len(reads[i+3].rstrip())
               for i in range(0,len(reads),4))
    bad=gzip.compress(b"@bad\nACGT\n+\nIII\n")
    with pytest.raises(ValueError,match="Invalid ENA source FASTQ"):
        m.fastq_prefix(bad,limit=2)


def test_bounded_only_120k_original_gzip_prefix_source():
    assert m.PAIRS==(30000,120000)
    assert m.MIN_LARGE_PAIRS==120000
    assert max(m.MAX_BYTES)==64*1024*1024
    assert len(m.RUNS)==4
    assert {z[2] for z in m.RUNS.values()}=={
        "S2_source_native_introduced","NA63_author_botanical_nodding"}
    assert m.LOCI==241


def fake_source():
    original=[]
    for run,(bio,name,kind) in m.RUNS.items():
        original.append({"run":run,"biosample":bio,
                         "S1_source_taxon":name,
                         "original_tree_code":name.replace("_"," ")+" 1",
                         "accepted_key":f"K_{run}"})
    return original


def test_original_four_run_BioSamples_and_taxon_types(monkeypatch):
    monkeypatch.setattr(m,"SOURCE_N",4)
    wanted=m.validate_sources(fake_source())
    assert len(wanted)==4
    original=fake_source()
    original[0]["biosample"]="SAMN99999999"
    with pytest.raises(ValueError,match="BioSample"):
        m.validate_sources(original)
    original=fake_source()
    original[2]["original_tree_code"]="Cirsium invented"
    with pytest.raises(ValueError,match="voucher label"):
        m.validate_sources(original)


def test_four_taxon_nested_read_depth_only_no_ortholog_proof(monkeypatch,tmp_path):
    monkeypatch.setattr(m,"SOURCE_N",4)
    monkeypatch.setattr(m,"MIN_LARGE_PAIRS",12)
    monkeypatch.setattr(m,"PAIRS",(3,12))
    monkeypatch.setattr(m,"LOCI",2)
    reads=[b"@x\n",b"ACGT\n",b"+\n",b"IIII\n"]*12
    def fetch(row):
        return [reads,reads],[{"mate":1},{"mate":2}]
    def mapper(row,reference,loci,source,path):
        n=len(source[0])//4
        return {"input_read_pairs":n,
                "n_loci_dp3_ge25pct":1 if n==3 else 2,
                "n_loci_dp3_ge50pct":0 if n==3 else 1,
                "n_loci_dp3_100bp":0 if n==3 else 1,
                "total_best_reference_sites_dp3":100 if n==3 else 200}
    ans=m.run(fake_source(),None,["A","B"],tmp_path,acquire=fetch,mapper=mapper)
    assert len(ans["pilot_taxa"])==4
    assert ans["full_EAzami_nuclear_species_tree_inferred"] is False
    assert ans["full_S1_source_names_still_251"] is True
    assert all(z["large_first_120k"]["n_loci_dp3_ge50pct"]==1
               for z in ans["pilot_taxa"])


def test_read_scaling_that_loses_callable_sites_fails(monkeypatch,tmp_path):
    monkeypatch.setattr(m,"SOURCE_N",4)
    monkeypatch.setattr(m,"MIN_LARGE_PAIRS",12)
    monkeypatch.setattr(m,"PAIRS",(3,12))
    monkeypatch.setattr(m,"LOCI",1)
    r=[b"@x\n",b"A\n",b"+\n",b"I\n"]*12
    def fetch(_):return [r,r],[{},{}]
    def mapper(_row,_ref,_loci,reads,_path):
        n=len(reads[0])//4
        return {"input_read_pairs":n,"n_loci_dp3_ge25pct":1 if n==3 else 0,
                "n_loci_dp3_ge50pct":0,"n_loci_dp3_100bp":0,
                "total_best_reference_sites_dp3":1}
    with pytest.raises(ValueError,match="decreased"):
        m.run(fake_source(),None,["A"],tmp_path,acquire=fetch,mapper=mapper)
