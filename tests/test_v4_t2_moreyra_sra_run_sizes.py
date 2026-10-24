"""NCBI SRA RunInfo-size audit tests: no raw FASTQ; no silent missing focal run."""
import pytest
from analysis.v4 import t2_moreyra_sra_run_sizes as mod


def test_runinfo_parses_compressed_archive_units_only():
    csv=("Run,BioSample,size_MB,LibraryStrategy,BioProject,LibraryLayout,spots,bases\n"
         "SRR100,SAMN100,103.5,TARGETED_CAPTURE,PRJNA957074,PAIRED,12,50000\n"
         "SRR200,SAMN200,82.25,OTHER,PRJNA957074,SINGLE,10,70000\n"
         "SRR999,SAMN999,,WGS,PRJNA957074,SINGLE,9,10\n")
    read=mod.parse_runinfo(csv)
    assert sorted(read)==["SRR100","SRR200"]
    assert read["SRR100"]["size_MB"]==103.5
    source=[
        {"run":"SRR100","biosample":"SAMN100","accepted_key":"x","S1_source_taxon":"Cirsium_a",
         "name_reconciliation_priority":"low","tree_code_vs_sra_name":"exact"},
        {"run":"SRR200","biosample":"SAMN200","accepted_key":"y","S1_source_taxon":"Cirsium_b",
         "name_reconciliation_priority":"high",
         "tree_code_vs_sra_name":"different_submitted_or_published_name"}]
    rows,summary=mod.audit(source,read)
    assert len(rows)==2
    assert summary["compressed_SRA_size_total_MB"]==185.75
    assert summary["source_matched_accepted_key_count"]==2
    assert summary["n_runs_with_original_name_disagreement"]==1
    assert summary["all_125_fastq_downloaded"] is False


def test_missing_requested_run_blocks_complete_archive_claim():
    read=mod.parse_runinfo(
        "Run,BioSample,size_MB\nSRR100,SAMN100,100\n")
    source=[{"run":"SRR200","biosample":"SAMN200","accepted_key":"x",
             "S1_source_taxon":"Cirsium_x"}]
    with pytest.raises(ValueError,match="lacks 1 pinned public runs"):
        mod.audit(source,read)


def test_wrong_biosample_stops_even_when_run_size_present():
    read=mod.parse_runinfo("Run,BioSample,size_MB\nSRR100,SAMN100,100\n")
    source=[{"run":"SRR100","biosample":"SAMN101","accepted_key":"x",
             "S1_source_taxon":"Cirsium_x"}]
    with pytest.raises(ValueError,match="BioSample identifier differs"):
        mod.audit(source,read)


def test_non_csv_response_fails_closed():
    with pytest.raises(ValueError,match="no delimited"):
        mod.parse_runinfo("<html>unavailable</html>")
