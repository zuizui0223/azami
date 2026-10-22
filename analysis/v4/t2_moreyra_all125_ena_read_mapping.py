"""Actual bounded ENA nuclear-target-read evidence across ALL 125 Moreyra S1 runs.

Source: original authenticated full-244 S1 x 113 Moreyra taxon / 125-SRR
public accession crosswalk, reused from EAzami294 BioSample-verified input.
ENA REST filereport supplies exact gzip paired-end FASTQ filenames; HTTP Range
retrieves at most 32 MiB per mate and extracts <=30,000 complete read pairs.
BWA aligns those *actual public reads* to the independently authenticated
original Compositae1061 1,061-locus reference, and samtools provides locus
OVERLAP counts. This does not assemble genes, infer a species tree, or solve
taxon concept equivalence. No read/BAM files are uploaded in artifacts.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import re
import shutil
import subprocess
import time
import urllib.parse

from analysis.v4 import t2_s2_ena_partial_fastq as ena
from analysis.v4.t2_s2_bounded_public_reads import read_idxstats

SOURCE_ROWS=125
SOURCE_TAXA=113
SHARDS=10
MAX_READS=30000
MIN_READS=1000
API=ena.API
FIELDS=ena.FIELDS
REFERENCE_LOCI=1061
REFERENCE_GIT_BLOB_SHA="4f89e234007f367ffa8aa5e2be536bc44f31f445"
S2_REFERENCE={
 "SRR25265555":("SAMN34240259","Cirsium_arvense"),
 "SRR25265662":("SAMN34240353","Cirsium_vulgare"),
}


def load_source(path):
    with Path(path).open(newline="",encoding="utf-8-sig") as handle:
        rows=list(csv.DictReader(handle))
    required={"accepted_key","S1_source_taxon","original_tree_code",
              "biosample","run","experiment","tree_code_vs_sra_name"}
    if len(rows)!=SOURCE_ROWS or not required<=set(rows[0]):
        raise ValueError("Source-authenticated complete 125-run Moreyra registry changed")
    if len({r["run"] for r in rows})!=SOURCE_ROWS:
        raise ValueError("Original Moreyra source contains duplicate SRR accession")
    if len({r["accepted_key"] for r in rows})!=SOURCE_TAXA:
        raise ValueError("Original S1-244 to Moreyra 113 taxon candidate coverage changed")
    for row in rows:
        if (not re.fullmatch(r"SRR[0-9]+",row["run"]) or
            not re.fullmatch(r"SAMN[0-9]+",row["biosample"]) or
            not row["original_tree_code"].startswith("Cirsium ") or
            not row["S1_source_taxon"].startswith("Cirsium_")):
            raise ValueError("Unverified original Moreyra public SRA run identity")
    byrun={r["run"]:r for r in rows}
    for run,(bio,taxon) in S2_REFERENCE.items():
        if byrun[run]["biosample"]!=bio or byrun[run]["S1_source_taxon"]!=taxon:
            raise ValueError("S2 source reference voucher drift")
    return sorted(rows,key=lambda r:r["run"])


def authoritative_fastq_files(run,request=ena.request_bytes):
    """Accept only official ENA paired FASTQ files for one whitelisted source run."""
    if not re.fullmatch(r"SRR[0-9]+",run):
        raise ValueError("Invalid SRA run for ENA mirrored read mapping")
    url=API+"?"+urllib.parse.urlencode({
        "accession":run,"result":"read_run","fields":FIELDS,"format":"tsv"})
    text=request(url).decode("utf-8-sig")
    rows=list(csv.DictReader(io.StringIO(text),delimiter="\t"))
    if len(rows)!=1 or rows[0].get("run_accession")!=run:
        raise ValueError("ENA filereport does not match authenticated run")
    item=rows[0]
    paths=item.get("fastq_ftp","").split(";")
    md5s=item.get("fastq_md5","").split(";")
    sizes=item.get("fastq_bytes","").split(";")
    if not(len(paths)==len(md5s)==len(sizes)==2):
        raise ValueError("ENA source does not provide original two paired FASTQ mates")
    out=[]
    for mate,(path,md5,size) in enumerate(zip(paths,md5s,sizes),1):
        if (not path.startswith("ftp.sra.ebi.ac.uk/vol1/fastq/") or
            not path.endswith(f"_{mate}.fastq.gz") or
            not path.rsplit("/",1)[-1].startswith(run) or
            not re.fullmatch(r"[0-9a-f]{32}",md5.lower())):
            raise ValueError("ENA paired public FASTQ URL/MD5 does not match source")
        length=int(size)
        if not 100<=length<=100000000000:
            raise ValueError("Invalid source FASTQ archive bytes")
        out.append({
            "url":"https://"+path,"archive_full_MD5_unchecked":md5,
            "full_gzip_size_bytes":length,"mate":mate,
        })
    return out


def source_prefix(row,workdir,report=authoritative_fastq_files,fetch=ena.request_bytes):
    """Materialize bounded first-spot FASTQ; never claim random unbiased subsample."""
    run=row["run"]
    files=report(run,fetch)
    mates=[]
    metrics=[]
    for item in files:
        seq,detail=ena.fetch_one_file(item,request=fetch,max_reads=MAX_READS)
        mates.append(seq)
        metrics.append(detail)
    count=min(len(x)//4 for x in mates)
    if count<MIN_READS or count>MAX_READS:
        raise ValueError("Too few verified paired reads in authentic bounded ENA source")
    folder=Path(workdir)
    folder.mkdir(parents=True,exist_ok=True)
    out=[]
    for i,records in enumerate(mates,1):
        dest=folder/f"{run}_{i}.fastq"
        dest.write_bytes(b"".join(records[:4*count]))
        out.append(dest)
    return out,{
       "max_spots":MAX_READS,"read_pairs":count,
       "downloaded_fastq_prefix_bytes":sum(x["downloaded_prefix_bytes"] for x in metrics),
       "ENA_file_MD5_full_archive_verified":False,
       "first_HTTP_prefix_not_random_sample":True,
    }


def run_commands(commands,timeout=240):
    return subprocess.run(commands,check=True,stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,timeout=timeout)


def align_and_count(row,ref,workdir,*,fetch=ena.request_bytes,report=authoritative_fastq_files):
    run=row["run"]
    folder=Path(workdir)/run
    folder.mkdir(parents=True,exist_ok=True)
    try:
        fastq,source=source_prefix(row,folder,report=report,fetch=fetch)
        bam=folder/"mapped.sorted.bam"
        bwa=subprocess.Popen(
            ["bwa","mem","-t","2",str(ref),str(fastq[0]),str(fastq[1])],
            stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        try:
            sort=run_commands(["samtools","sort","-@","2","-m","128M",
                               "-o",str(bam),"-"],timeout=240) if False else None
            # Use one bounded stream: BWA stdout must feed samtools stdin.
            proc=subprocess.Popen(
                ["samtools","sort","-@","2","-m","128M","-o",str(bam),"-"],
                stdin=bwa.stdout,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            bwa.stdout.close()
            _,sort_err=proc.communicate(timeout=240)
            _,bwa_err=bwa.communicate(timeout=240)
            if proc.returncode!=0 or bwa.returncode!=0:
                raise RuntimeError("BWA/samtools alignment failed: "+
                                   (sort_err+bwa_err).decode("utf-8","replace")[-400:])
        finally:
            if bwa.poll() is None: bwa.kill()
        run_commands(["samtools","index",str(bam)])
        indexed=run_commands(["samtools","idxstats",str(bam)]).stdout.decode("utf-8")
        count=run_commands(["samtools","view","-c","-q","20",
                            "-F","0x904",str(bam)]).stdout.decode().strip()
        idx=folder/"idxstats.tsv"
        idx.write_text(indexed,encoding="utf-8")
        loci=read_idxstats(idx)
        if loci["source_unique_loci"]!=REFERENCE_LOCI:
            raise ValueError("Original reference unique locus count changed")
        hq=int(count)
        if hq<0 or hq>2*source["read_pairs"] or hq>loci["mapped_alignments_from_idxstats"]:
            raise ValueError("Impossible mapped-read count from public pilot")
        return {
            "status":"mapped_real_public_ENA_FASTQ_prefix",
            "run":run,"biosample":row["biosample"],
            "accepted_key":row["accepted_key"],
            "S1_source_taxon":row["S1_source_taxon"],
            "original_tree_code":row["original_tree_code"],
            "original_taxon_name_conflict":row["tree_code_vs_sra_name"]!="exact",
            **source,
            "mapq20_primary_alignments":hq,
            **loci,
            "full_241_locus_ortholog_recovery":False,
            "tip_can_enter_nuclear_species_tree":False,
        }
    finally:
        shutil.rmtree(folder,ignore_errors=True)


def collect_shard(source,shard_id,shard_count,workdir,reference,*,analysis=align_and_count):
    if not(0<=shard_id<shard_count<=SHARDS):
        raise ValueError("Invalid bounded run batch/shard ID")
    selected=[x for i,x in enumerate(source) if i%shard_count==shard_id]
    if not selected:
        raise ValueError("No source runs in shard")
    output=[]
    for idx,row in enumerate(selected,1):
        start=time.monotonic()
        try:
            value=analysis(row,reference,workdir)
            if value["run"]!=row["run"] or value["biosample"]!=row["biosample"]:
                raise ValueError("Analysed accession does not match source roster")
            value["elapsed_seconds"]=round(time.monotonic()-start,3)
            output.append(value)
            print("T2_ALL125_SOURCE_NUCLEAR_RUN_OK="+json.dumps({
                "run":row["run"],"mapped_loci":value.get("mapped_reference_loci_at_least_one_alignment"),
                "read_pairs":value.get("read_pairs")}),flush=True)
        except Exception as err:
            failed={"status":"failed","run":row["run"],
                "biosample":row["biosample"],"accepted_key":row["accepted_key"],
                "S1_source_taxon":row["S1_source_taxon"],
                "reason_type":type(err).__name__,
                "reason":str(err)[:350],
                "elapsed_seconds":round(time.monotonic()-start,3),
                "tip_can_enter_nuclear_species_tree":False}
            output.append(failed)
            print("T2_ALL125_SOURCE_NUCLEAR_RUN_FAILED="+json.dumps(failed,
                  sort_keys=True),flush=True)
    if {x["run"] for x in output}!={r["run"] for r in selected}:
        raise ValueError("Lost source accession during bounded public read run")
    return {
       "schema":"azami_v4_all125_Moreyra_public_read_30k_prefix_locus_audit_v1",
       "study_sample_status":"REAL_PUBLIC_SHORT_READ_BOUNDED_MAPPING_ONLY",
       "shard_id":shard_id,"shard_total":shard_count,
       "requested_run_count":len(selected),
       "mapped_run_count":sum(x["status"].startswith("mapped_real") for x in output),
       "failed_run_count":sum(x["status"]=="failed" for x in output),
       "original_reference_git_blob_sha1":REFERENCE_GIT_BLOB_SHA,
       "rooted_global_nuclear_species_tree_inferred":False,
       "sequence_data_uploaded":False,
       "original_350_locus_tree_reproduced":False,
       "outcomes":output,
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--reference",type=Path,required=True)
    p.add_argument("--shard-id",type=int,required=True)
    p.add_argument("--shard-count",type=int,default=SHARDS)
    p.add_argument("--scratch",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    if not a.reference.is_file():
        raise ValueError("Original audited nuclear Compositae1061 reference missing")
    source=load_source(a.source)
    result=collect_shard(source,a.shard_id,a.shard_count,a.scratch,a.reference)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,sort_keys=True,ensure_ascii=False,
                          indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print("T2_ALL125_PUBLIC_SRA_REAL_NUCLEAR_COVERAGE_SHARD="+json.dumps({
       "shard":a.shard_id,"n":result["requested_run_count"],
       "ok":result["mapped_run_count"],
       "failed":result["failed_run_count"]}),flush=True)


if __name__=="__main__":
    main()
