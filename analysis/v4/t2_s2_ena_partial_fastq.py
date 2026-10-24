"""Fetch bounded REAL ENA mirrored FASTQ prefixes for two exact Moreyra S2 runs.

The original NCBI SRA Toolkit returned "Failed to call external services"
on GitHub Actions even though the authenticated Compositae1061 reference and
source BioSample crosswalk both passed. This alternative discovers FASTQ file
URLs/MD5 from ENA's documented public filereport API rather than inventing
FTP filenames. Downloads at most 32 MiB per mate through HTTP Range, then
emits up to 30000 complete FASTQ records per mate. No full FASTQ integrity
claim: partial content cannot verify the archive-published complete-file MD5.

For public target-capture *access validation and locus-overlap only*.
Outputs are never a 241-locus species alignment or phylogenetic tree.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import time
import urllib.error
import urllib.parse
import urllib.request

RUNS=(
    ("Cirsium arvense","SRR25265555","SAMN34240259"),
    ("Cirsium vulgare","SRR25265662","SAMN34240353"),
)
FIELDS="run_accession,fastq_ftp,fastq_md5,fastq_bytes"
API="https://www.ebi.ac.uk/ena/portal/api/filereport"
MAX_READS=30000
HTTP_RANGE_TRIES=(8*1024*1024,16*1024*1024,32*1024*1024)


def request_bytes(url,range_end=None):
    headers={"User-Agent":"azami-EAzami-source-verified-ENA-nuclear-pilot/1.0"}
    if range_end is not None:
        headers["Range"]=f"bytes=0-{range_end-1}"
    for attempt in range(3):
        try:
            req=urllib.request.Request(url,headers=headers)
            with urllib.request.urlopen(req,timeout=85) as response:
                data=response.read(range_end+1 if range_end else 300000)
                return data
        except (urllib.error.HTTPError,urllib.error.URLError,TimeoutError) as err:
            if attempt==2:
                raise RuntimeError("ENA public source access failed after bounded retries") from err
            time.sleep(attempt+1)
    raise AssertionError("unreachable")


def source_fastq_files(run,request=request_bytes):
    if run not in {r[1] for r in RUNS}:
        raise ValueError("Unrecognized public S2 original source accession")
    q=urllib.parse.urlencode({"accession":run,"result":"read_run",
                               "fields":FIELDS,"format":"tsv"})
    raw=request(API+"?"+q).decode("utf-8-sig")
    parsed=list(csv.DictReader(io.StringIO(raw),delimiter="\t"))
    if len(parsed)!=1 or parsed[0].get("run_accession")!=run:
        raise ValueError("Official ENA filereport did not return exact original accession")
    rec=parsed[0]
    paths=rec.get("fastq_ftp","").split(";")
    digests=rec.get("fastq_md5","").split(";")
    sizes=rec.get("fastq_bytes","").split(";")
    if len(paths)!=2 or len(digests)!=2 or len(sizes)!=2:
        raise ValueError("Expected two official paired-end ENA FASTQ files")
    files=[]
    for i,(path,digest,size) in enumerate(zip(paths,digests,sizes),start=1):
        if (not path.startswith("ftp.sra.ebi.ac.uk/vol1/fastq/") or
                not path.endswith(f"_{i}.fastq.gz") or
                not path.rsplit("/",1)[-1].startswith(run) or
                not re.fullmatch(r"[0-9a-f]{32}",digest.lower())):
            raise ValueError("ENA file path or full-file MD5 inconsistent with accession")
        size=int(size)
        if not 100<=size<=100000000000:
            raise ValueError("Invalid ENA archived gzipped FASTQ byte size")
        files.append({"url":"https://"+path,"archive_full_MD5_unchecked":digest,
                      "full_gzip_size_bytes":size,"mate":i})
    return files


def fastq_prefix(raw,max_reads):
    if max_reads<=0 or max_reads>MAX_READS:
        raise ValueError("Unbounded S2 FASTQ extraction prohibited")
    result=[]
    truncated=False
    try:
        fh=gzip.GzipFile(fileobj=io.BytesIO(raw),mode="rb")
        for _ in range(max_reads):
            header=fh.readline()
            if not header:
                break
            seq=fh.readline()
            plus=fh.readline()
            qual=fh.readline()
            if not seq or not plus or not qual:
                raise EOFError("Truncated FASTQ record")
            if (not header.startswith(b"@") or not plus.startswith(b"+") or
                    len(seq.rstrip(b"\r\n"))!=len(qual.rstrip(b"\r\n"))):
                raise ValueError("Malformed public gzipped ENA FASTQ record")
            result.extend((header,seq,plus,qual))
    except (EOFError, OSError) as err:
        # A valid HTTP prefix may end inside one gzip member; retry a bigger
        # Range rather than presenting this as a complete FASTQ/MD5 check.
        truncated=True
    return result,truncated


def fetch_one_file(item,request=request_bytes,max_reads=MAX_READS):
    observed=None
    for max_bytes in HTTP_RANGE_TRIES:
        limit=min(item["full_gzip_size_bytes"],max_bytes)
        raw=request(item["url"],range_end=limit)
        if len(raw)>limit:
            raise ValueError("ENA remote range response exceeded bounded byte request")
        seq,truncated=fastq_prefix(raw,max_reads)
        n=len(seq)//4
        if n==max_reads or (not truncated and len(raw)==item["full_gzip_size_bytes"] and n>0):
            return seq,{"downloaded_prefix_bytes":len(raw),
                         "fastq_records_extracted":n,
                         "full_gzip_size_bytes":item["full_gzip_size_bytes"],
                         "only_partial_gzip_md5_available":True,
                         "full_gzip_MD5_from_ENA_not_cryptographically_checked":
                               item["archive_full_MD5_unchecked"],
                         "HTTP_range_max_bytes":limit,
                         "source_fastq_mate":item["mate"]}
        observed=(n,len(raw),truncated)
    raise ValueError(f"Only insufficient valid FASTQ prefix for mate {item['mate']}: {observed}")


def retrieve_run(taxon,run,bio,outdir,report=source_fastq_files,request=request_bytes):
    files=report(run,request)
    mates=[]
    metrics=[]
    for f in files:
        seq,m=fetch_one_file(f,request)
        mates.append(seq);metrics.append(m)
    n1=len(mates[0])//4;n2=len(mates[1])//4
    if min(n1,n2)<1000:
        raise ValueError("ENA public source prefix too small for meaningful read-access audit")
    usable=min(n1,n2)
    folder=Path(outdir)/run
    (folder/"fastq").mkdir(parents=True,exist_ok=True)
    for i,reads in enumerate(mates,start=1):
        (folder/"fastq"/f"{run}_{i}.fastq").write_bytes(b"".join(reads[:4*usable]))
    receipt={
        "run":run,"taxon":taxon,"biosample":bio,
        "source":"Moreyra2025","max_spots":MAX_READS,
        "forward_fastq_reads":usable,
        "paired_fastq_reads":usable,
        "retrieval_route":"ENA_official_filreports_gzip_HTTP_Range_prefix",
        "ENA_metadata_IS_source_run_file_not_full_fastq_verification":True,
        "full_file_MD5_verified":False,
        "mate_partial_source_details":metrics,
    }
    (folder/"extraction_metrics.json").write_text(
        json.dumps(receipt,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print("T2_REAL_ENA_S2_PUBLIC_FASTQ_READS="+json.dumps({
          "run":run,"paired_read_records":usable,
          "retrieval_route":receipt["retrieval_route"]},sort_keys=True),flush=True)
    return receipt


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out-dir",type=Path,required=True)
    a=p.parse_args()
    for taxon,run,bio in RUNS:
        retrieve_run(taxon,run,bio,a.out_dir)


if __name__=="__main__":
    main()
