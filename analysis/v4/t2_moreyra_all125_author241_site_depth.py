"""Site-depth audit on all 125 real S1-linked Moreyra 2025 ENA short-read prefixes.

Each run's first <=30k authentic ENA paired FASTQ reads is aligned against the
original authenticated Compositae1061 1061-locus reference with BWA. Using
samtools depth -q20 (base quality) -Q20 (mapping quality), compute observed
>=1/3/5x site breadth for the PUBLIC AUTHOR FILTERED no-warning 241 loci.

Reference copy: for EACH locus, choose the BEST single original Asteraceae
reference sequence by depth3 covered fraction. This is an upper-bound mapping
coverage diagnostic, NEVER an orthology-certified Cirsium gene consensus:
2-3 references per locus are different and cannot be summed positionwise.
Coordinates, FASTQs, BAMs are erased before artifact publication.
"""
from __future__ import annotations
import argparse
from collections import defaultdict, Counter
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import time

from analysis.v4 import t2_moreyra_all125_ena_read_mapping as upstream

RUN_COUNT=125
TAXA_COUNT=113
LOCUS_COUNT=241
EXPECTED_LOCUS_LIST_SHA="d561c6e393b1964fdd4b3acf14fda8b10f2f43923b1074cd35f86bfed07ebf73"
REFERENCE_SHA=upstream.REFERENCE_GIT_BLOB_SHA
BREADTH_THRESHOLDS=(.1,.25,.5,.75)
MIN_CALLABLE_LENGTH=100
DEPTHS=(1,3,5)


def loci_from_pinned_author_list(path):
    raw=Path(path).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=EXPECTED_LOCUS_LIST_SHA:
        raise ValueError("Source author conservatively filtered 241 locus list digest drift")
    loci=[line.strip() for line in raw.decode("utf-8").splitlines() if line.strip()]
    if len(loci)!=LOCUS_COUNT or len(set(loci))!=LOCUS_COUNT:
        raise ValueError("Author 241 unique filtered loci expected")
    return loci


def reference_lengths(idxstats):
    lengths={}
    loci=defaultdict(list)
    for row in idxstats.strip().splitlines():
        parts=row.split("\t")
        if len(parts)!=4:
            raise ValueError("Malformed samtools idxstats record")
        name,size,mapped,_unmapped=parts
        if name=="*":
            continue
        if not name.startswith(("lett-","sunf-","saff-")):
            raise ValueError("Unexpected Asteraceae reference contig")
        if name in lengths or int(size)<50 or int(mapped)<0:
            raise ValueError("Duplicated/bad original reference sequence")
        lengths[name]=int(size)
        loci[name.split("-",1)[1]].append(name)
    if len(loci)!=1061:
        raise ValueError("Original 1061-locus reference IDs no longer agree")
    return lengths,loci


def count_depth_stream(lines,lengths):
    """Read samtools depth output; use only occupied positions (no fabricated 0x calls)."""
    summary=defaultdict(lambda:[0,0,0])
    previous={}
    for line in lines:
        if isinstance(line,bytes):
            line=line.decode("ascii")
        fields=line.rstrip("\n").split("\t")
        if len(fields)!=3:
            raise ValueError("Unrecognized samtools depth output")
        contig,pos,rawdp=fields
        if contig not in lengths:
            raise ValueError("Unrecognized nuclear target contig in depth")
        position=int(pos); depth=int(rawdp)
        if not (1<=position<=lengths[contig] and depth>=1):
            raise ValueError("Invalid aligned nuclear site depth")
        if position<=previous.get(contig,0):
            raise ValueError("Duplicate/out-of-order site depth records")
        previous[contig]=position
        counter=summary[contig]
        for i,min_depth in enumerate(DEPTHS):
            if depth>=min_depth: counter[i]+=1
    return dict(summary)


def summarize_loci(locus_ids,lengths,refs,depths):
    result=[]
    for locus in locus_ids:
        candidates=refs.get(locus,[])
        if not candidates:
            raise ValueError("No homologous reference contig for author-filtered locus")
        def rank(ref):
            d=depths.get(ref,[0,0,0])
            return (d[1]/lengths[ref],d[2]/lengths[ref],
                    d[0]/lengths[ref],ref)
        best=max(candidates,key=rank)
        c=depths.get(best,[0,0,0])
        length=lengths[best]
        if any(z<0 or z>length for z in c) or not c[2]<=c[1]<=c[0]:
            raise ValueError("Reference coverage exceeds sequence length")
        row={
           "locus":locus,"best_original_reference":best,
           "n_original_reference_copies":len(candidates),
           "reference_length_nt":length,
           "dp1_callable_sites":c[0],
           "dp3_callable_sites":c[1],
           "dp5_callable_sites":c[2],
           "dp1_breadth":c[0]/length,
           "dp3_breadth":c[1]/length,
           "dp5_breadth":c[2]/length,
           "dp3_at_least100bp":c[1]>=MIN_CALLABLE_LENGTH,
           "reference_choice_maximizes_coverage_not_orthology":True,
        }
        result.append(row)
    if len(result)!=LOCUS_COUNT:
        raise ValueError("Locus audit did not retain all 241 original candidates")
    return result


def per_run_metrics(rows):
    if len(rows)!=LOCUS_COUNT:
        raise ValueError("All 241 per-run loci required")
    summary={
       "n_author_low_paralog_241_candidates":LOCUS_COUNT,
       "n_loci_dp1_any":sum(z["dp1_callable_sites"]>0 for z in rows),
       "n_loci_dp3_any":sum(z["dp3_callable_sites"]>0 for z in rows),
       "n_loci_dp5_any":sum(z["dp5_callable_sites"]>0 for z in rows),
       "n_loci_dp3_100bp":sum(z["dp3_at_least100bp"] for z in rows),
       "n_loci_dp3_ge10pct":sum(z["dp3_breadth"]>=.1 for z in rows),
       "n_loci_dp3_ge25pct":sum(z["dp3_breadth"]>=.25 for z in rows),
       "n_loci_dp3_ge50pct":sum(z["dp3_breadth"]>=.5 for z in rows),
       "n_loci_dp3_ge75pct":sum(z["dp3_breadth"]>=.75 for z in rows),
       "total_best_reference_bp":sum(z["reference_length_nt"] for z in rows),
       "total_best_reference_sites_dp1":sum(z["dp1_callable_sites"] for z in rows),
       "total_best_reference_sites_dp3":sum(z["dp3_callable_sites"] for z in rows),
       "total_best_reference_sites_dp5":sum(z["dp5_callable_sites"] for z in rows),
    }
    if any(summary[key]>LOCUS_COUNT for key in summary if key.startswith("n_loci_")):
        raise ValueError("Impossible source filtered locus coverage count")
    return summary


def real_depth_run(source,reference,locus_ids,scratch):
    run=source["run"]
    folder=Path(scratch)/run
    folder.mkdir(parents=True,exist_ok=True)
    try:
        fastq,proof=upstream.source_prefix(source,folder)
        bam=folder/"original_ref_reads.sorted.bam"
        bwa=subprocess.Popen(
           ["bwa","mem","-t","2",str(reference),str(fastq[0]),str(fastq[1])],
           stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        try:
            sorter=subprocess.Popen(
                ["samtools","sort","-@","2","-m","128M","-o",str(bam),"-"],
                stdin=bwa.stdout,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            bwa.stdout.close()
            _stdout,sort_err=sorter.communicate(timeout=240)
            _stdout,bwa_err=bwa.communicate(timeout=240)
            if sorter.returncode!=0 or bwa.returncode!=0:
                raise RuntimeError("BWA/samtools real source alignment failed "+
                    (sort_err+bwa_err).decode("utf-8","replace")[-200:])
        finally:
            if bwa.poll() is None:bwa.kill()
        upstream.run_commands(["samtools","index",str(bam)])
        ix=upstream.run_commands(["samtools","idxstats",str(bam)]).stdout.decode()
        lengths,refs=reference_lengths(ix)
        depth_cmd=["samtools","depth","-q","20","-Q","20",str(bam)]
        proc=subprocess.Popen(depth_cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                              text=True,bufsize=4096)
        try:
            raw=count_depth_stream(proc.stdout,lengths)
            proc.stdout.close()
            err=proc.stderr.read()
            proc.wait(timeout=180)
            if proc.returncode!=0:
                raise RuntimeError("samtools verified mapping/base quality depth failed "+err[-500:])
        finally:
            if proc.poll() is None:
                proc.kill()
                proc.wait()
        locus_rows=summarize_loci(locus_ids,lengths,refs,raw)
        counts=per_run_metrics(locus_rows)
        return {
           "status":"real_ENA_reference_site_depth_q20_Q20",
           "run":run,"biosample":source["biosample"],"accepted_key":source["accepted_key"],
           "S1_source_taxon":source["S1_source_taxon"],
           "original_tree_code":source["original_tree_code"],
           "real_ena_read_pairs":proof["read_pairs"],
           "ENA_file_MD5_full_archive_verified":False,
           "SRA_prefix_sample_not_random":True,
           **counts,
           "loci":locus_rows,
           "reference_species_tree_inferred":False,
           "independently_verified_orthologous_full_gene_sequences":False,
        }
    finally:
        shutil.rmtree(folder,ignore_errors=True)


def collect(source,locus_ids,reference,scratch,shard_id,shards=10,runner=real_depth_run):
    if not (0<=shard_id<shards<=10):
        raise ValueError("Shard must be in 0..9")
    if len(source)!=RUN_COUNT or len({z["accepted_key"] for z in source})!=TAXA_COUNT:
        raise ValueError("Original 125 source run/113 taxon data lost")
    chosen=[r for i,r in enumerate(source) if i%shards==shard_id]
    result=[]
    for row in chosen:
        start=time.monotonic()
        try:
            v=runner(row,reference,locus_ids,scratch)
            if v["run"]!=row["run"] or v["biosample"]!=row["biosample"]:
                raise ValueError("Source authenticity changed during real sequencing data acquisition")
            v["elapsed_seconds"]=round(time.monotonic()-start,2)
            result.append(v)
            print("T2_REAL_241_LOCUS_DEPTH_OK="+json.dumps({k:v[k] for k in (
                "run","real_ena_read_pairs","n_loci_dp3_ge25pct","n_loci_dp3_ge50pct")}),
                flush=True)
        except Exception as err:
            result.append({
               "run":row["run"],"accepted_key":row["accepted_key"],
               "biosample":row["biosample"],"status":"failed",
               "error_class":type(err).__name__,"error":str(err)[:350],
               "elapsed_seconds":round(time.monotonic()-start,2),
            })
            print("T2_REAL_241_LOCUS_DEPTH_FAILED="+json.dumps(result[-1]),flush=True)
    if {x["run"] for x in result}!={x["run"] for x in chosen}:
        raise ValueError("Missing source run after read depth collection")
    return {
      "schema":"azami_v4_author241_site_depth_prefix_reads_v1",
      "scope":"original S1 113 taxa x 125 runs, no full FASTQ assembly or phylogeny",
      "shard_id":shard_id,"n_shards":shards,
      "requested":len(chosen),
      "successful":sum(x["status"].startswith("real_ENA") for x in result),
      "failures":sum(x["status"]=="failed" for x in result),
      "S1_full244_taxon_registry_retained_elsewhere":True,
      "read_data_in_artifact":False,"BAM_in_artifact":False,
      "public_original_reference_sha":REFERENCE_SHA,
      "author_conservative_241_SHA256":EXPECTED_LOCUS_LIST_SHA,
      "outcomes":result,
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--loci",type=Path,required=True)
    p.add_argument("--reference",type=Path,required=True)
    p.add_argument("--scratch",type=Path,required=True)
    p.add_argument("--shard-id",type=int,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    source=upstream.load_source(a.source)
    loci=loci_from_pinned_author_list(a.loci)
    if not a.reference.is_file():
        raise ValueError("Original verified Compositae1061 DNA source absent")
    result=collect(source,loci,a.reference,a.scratch,a.shard_id)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,
                               sort_keys=True,allow_nan=False)+"\n",encoding="utf-8")
    print("T2_REAL_241_DEPTH_SHARD_COMPLETE="+json.dumps({
       "id":result["shard_id"],"run_count":result["requested"],
       "succeeded":result["successful"],"failed":result["failures"]}),
       flush=True)


if __name__=="__main__":
    main()
