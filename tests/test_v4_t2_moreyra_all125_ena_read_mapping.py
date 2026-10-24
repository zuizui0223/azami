"""Fail-closed checks for 125-run source-pinned ENA read-to-locus mapping."""
import gzip
import hashlib
from pathlib import Path
import pytest

from analysis.v4 import t2_moreyra_all125_ena_read_mapping as m


def fake_source():
    runs=[
       {"run":"SRR25265555","biosample":"SAMN34240259",
        "S1_source_taxon":"Cirsium_arvense","original_tree_code":"Cirsium arvense",
        "accepted_key":"206745077","experiment":"SRX21011661",
        "tree_code_vs_sra_name":"exact"},
       {"run":"SRR25265662","biosample":"SAMN34240353",
        "S1_source_taxon":"Cirsium_vulgare","original_tree_code":"Cirsium vulgare 1",
        "accepted_key":"206739917","experiment":"SRX21011555",
        "tree_code_vs_sra_name":"exact"}]
    for i in range(123):
        idx=i if i<111 else i-111
        runs.append({"run":f"SRR9000000{i:03d}",
          "biosample":f"SAMN9000000{i:03d}",
          "S1_source_taxon":f"Cirsium_species{idx:03d}",
          "original_tree_code":f"Cirsium species{idx:03d}",
          "accepted_key":f"OTHER{idx:03d}",
          "experiment":"SRXTEST",
          "tree_code_vs_sra_name":"different_submitted_or_published_name" if i%5==0 else "exact"})
    return runs


def test_real_source_exact_125_and_113_taxon_groups(tmp_path):
    import csv
    path=tmp_path/"source.csv"
    rows=fake_source()
    with path.open("w",newline="",encoding="utf-8") as out:
        writer=csv.DictWriter(out,fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    source=m.load_source(path)
    assert len(source)==125
    assert len({x["accepted_key"] for x in source})==113
    assert {x["run"] for x in source}=={x["run"] for x in rows}
    with path.open("w",newline="",encoding="utf-8") as out:
        writer=csv.DictWriter(out,fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows[:-1])
    with pytest.raises(ValueError,match="125-run"):
        m.load_source(path)


def test_ena_official_file_report_is_mate_and_original_run_strict():
    run="SRR25265555"
    prefix=f"ftp.sra.ebi.ac.uk/vol1/fastq/SRR252/{run}/{run}"
    digest="0"*32
    report=("run_accession\tfastq_ftp\tfastq_md5\tfastq_bytes\n"
            f"{run}\t{prefix}_1.fastq.gz;{prefix}_2.fastq.gz"
            f"\t{digest};{digest}\t123456;234567\n").encode()
    def request(url,range_end=None):
        assert run in url and "filereport" in url
        return report
    files=m.authoritative_fastq_files(run,request=request)
    assert len(files)==2 and files[0]["mate"]==1 and files[1]["mate"]==2
    assert files[0]["full_gzip_size_bytes"]==123456
    bad=report.replace(f"{run}_2.fastq.gz".encode(),b"wrong_2.fastq.gz")
    with pytest.raises(ValueError,match="FASTQ URL"):
        m.authoritative_fastq_files(run,request=lambda _:bad)


def test_source_prefix_real_gzip_limited_not_random_and_cleanup_local(tmp_path):
    run="SRR25265555"
    fastq=b"".join((f"@x{i}/1\nACGTACGT\n+\nFFFFFFFF\n").encode() for i in range(1400))
    raw=gzip.compress(fastq)
    urls=[{"url":f"https://ftp.sra.ebi.ac.uk/{run}_{i}.fastq.gz",
           "full_gzip_size_bytes":len(raw),
           "archive_full_MD5_unchecked":hashlib.md5(raw).hexdigest(),
           "mate":i} for i in (1,2)]
    def fetch(url,range_end=None):
        return raw[:range_end] if range_end is not None else raw
    rec={"run":run}
    files,summary=m.source_prefix(rec,tmp_path,report=lambda run,request:urls,fetch=fetch)
    assert summary["read_pairs"]==1400
    assert not summary["ENA_file_MD5_full_archive_verified"]
    assert summary["first_HTTP_prefix_not_random_sample"]
    assert sum(1 for _ in files[0].open("rb"))==5600


def test_all_125_assigned_once_across_10_shards():
    source=sorted(fake_source(),key=lambda x:x["run"])
    all_ids=[]
    for shard in range(10):
        def simulate(row,ref,path):
            return {
               "run":row["run"],"biosample":row["biosample"],
               "status":"mapped_real_public_ENA_FASTQ_prefix",
               "mapped_reference_loci_at_least_one_alignment":500,
               "read_pairs":1000,
            }
        result=m.collect_shard(source,shard,10,Path("/unused"),Path("/unused"),
                               analysis=simulate)
        all_ids.extend(z["run"] for z in result["outcomes"])
        assert result["failed_run_count"]==0
        assert result["mapped_run_count"]==len(result["outcomes"])
        assert not result["rooted_global_nuclear_species_tree_inferred"]
    assert len(all_ids)==len(set(all_ids))==125


def test_failed_public_source_never_promoted_as_real_nuclear_tip():
    source=sorted(fake_source(),key=lambda x:x["run"])
    def simulate(row,ref,path):
        raise TimeoutError("ENA fail explicitly retained")
    data=m.collect_shard(source,0,10,"/scratch","/ref",analysis=simulate)
    assert data["failed_run_count"]==13
    assert all(x["status"]=="failed" and not x["tip_can_enter_nuclear_species_tree"]
               for x in data["outcomes"])


def test_maximum_original_source_boundaries():
    assert m.SOURCE_ROWS==125 and m.SOURCE_TAXA==113
    assert m.SHARDS==10
    assert m.MAX_READS==30000
    assert m.MIN_READS==1000
    assert m.REFERENCE_LOCI==1061
    assert m.REFERENCE_GIT_BLOB_SHA=="4f89e234007f367ffa8aa5e2be536bc44f31f445"
