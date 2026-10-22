"""ENA authoritative file-report/HTTP Range gzip-prefix parser tests."""
import gzip
import hashlib
from pathlib import Path

import pytest

from analysis.v4 import t2_s2_ena_partial_fastq as ena


def simulated_fastq(n=1500):
    lines=[]
    for i in range(n):
        lines.extend((f"@x_{i}/1\n", "ACGTACGT\n", "+\n", "FFFFFFFF\n"))
    return gzip.compress("".join(lines).encode("utf-8"))


def test_authoritative_ena_fastq_report_and_bounded_partial_read_access(tmp_path):
    raw=simulated_fastq(1600)
    digest=hashlib.md5(raw).hexdigest()
    path_base="ftp.sra.ebi.ac.uk/vol1/fastq/SRR252/SRR25265555/SRR25265555"
    files=f"{path_base}_1.fastq.gz;{path_base}_2.fastq.gz"
    tsv=(
        "run_accession\tfastq_ftp\tfastq_md5\tfastq_bytes\n"
        f"SRR25265555\t{files}\t{digest};{digest}\t{len(raw)};{len(raw)}\n"
    ).encode()
    def request(url,range_end=None):
        if "filereport?" in url:
            assert "SRR25265555" in url
            return tsv
        assert url.startswith("https://ftp.sra.ebi.ac.uk/vol1/fastq/")
        return raw[:range_end] if range_end is not None else raw
    evidence=ena.source_fastq_files("SRR25265555",request=request)
    assert len(evidence)==2
    r=ena.retrieve_run("Cirsium arvense","SRR25265555","SAMN34240259",
                       tmp_path,report=lambda run,req:evidence,request=request)
    assert r["paired_fastq_reads"]==1600
    assert r["max_spots"]==30000
    assert r["full_file_MD5_verified"] is False
    assert r["retrieval_route"].startswith("ENA_official")
    assert (tmp_path/"SRR25265555"/"fastq"/"SRR25265555_1.fastq").exists()
    assert (tmp_path/"SRR25265555"/"fastq"/"SRR25265555_2.fastq").exists()


def test_exact_run_path_and_mate_validation_prevents_other_taxon_reads():
    tsv=(
        "run_accession\tfastq_ftp\tfastq_md5\tfastq_bytes\n"
        "SRR25265555\tftp.sra.ebi.ac.uk/vol1/fastq/other/other_1.fastq.gz;"
        "ftp.sra.ebi.ac.uk/vol1/fastq/other/other_2.fastq.gz\t"
        +("0"*32+";"+"1"*32)+"\t10000;10000\n"
    ).encode()
    with pytest.raises(ValueError,match="file path"):
        ena.source_fastq_files("SRR25265555",request=lambda _:tsv)


def test_partial_gzip_is_not_reported_as_complete_archive():
    raw=simulated_fastq(5000)
    parsed,cut=ena.fastq_prefix(raw[:len(raw)//3],5000)
    assert cut
    assert 0<len(parsed)//4<5000


def test_corrupt_gzip_or_fastq_record_fails_recovery():
    with pytest.raises((ValueError, OSError)):
        ena.fastq_prefix(gzip.compress(b"bad_header\nAAAA\n+\nFFFF\n"),10)


def test_exact_original_s2_source_accessions_and_bounded_limits():
    assert ena.RUNS==(
        ("Cirsium arvense","SRR25265555","SAMN34240259"),
        ("Cirsium vulgare","SRR25265662","SAMN34240353"),
    )
    assert ena.MAX_READS==30000
    assert max(ena.HTTP_RANGE_TRIES)==32*1024*1024
