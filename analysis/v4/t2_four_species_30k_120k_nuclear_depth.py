"""Matched real 30k-vs-120k public ENA read-depth test on four Cirsium.

Two genuine source S2 focal taxa (C. arvense/C. vulgare), two independently
author-scored nodding taxa (C. conspicuum/C. scopulorum). These four are a
depth-information diagnostic; NOT the full S1 population or an independent
fitness/selection comparison. The full 125-run 30k census remains primary.

Read the SAME run's first 120k pairs and run the original Q20 site-depth
protocol independently on first 30k and 120k subsets. No reads or BAMs
are archived. Four runs x original author 241 loci are source-verified.
"""
from __future__ import annotations

import argparse
import gzip
import io
import json
from pathlib import Path
import shutil
import subprocess
import time

from analysis.v4 import t2_moreyra_all125_ena_read_mapping as upstream
from analysis.v4 import t2_moreyra_all125_author241_site_depth as dp
from analysis.v4 import t2_s2_ena_partial_fastq as ena

PAIRS=(30000,120000)
MAX_BYTES=(16*1024*1024,32*1024*1024,64*1024*1024)
RUNS={
 "SRR25265555":("SAMN34240259","Cirsium_arvense","S2_source_native_introduced"),
 "SRR25265662":("SAMN34240353","Cirsium_vulgare","S2_source_native_introduced"),
 "SRR30887263":("SAMN44017833","Cirsium_conspicuum","NA63_author_botanical_nodding"),
 "SRR30887312":("SAMN44017912","Cirsium_scopulorum","NA63_author_botanical_nodding"),
}
MIN_LARGE_PAIRS=120000
SOURCE_N=125
LOCI=241


def fastq_prefix(raw,limit=MIN_LARGE_PAIRS):
    """Read complete valid FASTQ records only; tolerate incomplete gzip tail."""
    if not 1<=limit<=MIN_LARGE_PAIRS:
        raise ValueError("Unbounded public FASTQ prefix request")
    out=[]
    truncated=False
    with gzip.GzipFile(fileobj=io.BytesIO(raw),mode="rb") as file:
        try:
            for _ in range(limit):
                a=file.readline()
                if not a:break
                b=file.readline()
                c=file.readline()
                d=file.readline()
                if not all((b,c,d)):
                    raise EOFError("Incomplete FASTQ read")
                if not a.startswith(b"@") or not c.startswith(b"+") or \
                    len(b.rstrip(b"\r\n"))!=len(d.rstrip(b"\r\n")):
                    raise ValueError("Invalid ENA source FASTQ format")
                out.extend((a,b,c,d))
        except (EOFError,OSError):
            truncated=True
    return out,truncated


def acquire_larger_prefix(author_row,request=ena.request_bytes,report=upstream.authoritative_fastq_files):
    """Only official paired ENA endpoints, bounded to <=64MiB per mate."""
    meta=report(author_row["run"],request)
    mates=[]
    details=[]
    for sample in meta:
        got=None
        for limit in MAX_BYTES:
            end=min(limit,sample["full_gzip_size_bytes"])
            raw=request(sample["url"],range_end=end)
            if len(raw)>end:
                raise ValueError("Invalid ENA HTTP Range > requested maximum")
            reads,truncated=fastq_prefix(raw)
            found=len(reads)//4
            if found==MIN_LARGE_PAIRS:
                got=reads
                details.append({"mate":sample["mate"],
                    "read_pairs":found,"downloaded_prefix_bytes":len(raw),
                    "first_prefix_not_random":True,"full_file_MD5_verified":False})
                break
            if not truncated and len(raw)==sample["full_gzip_size_bytes"]:
                # Truly short public FASTQ: do not extrapolate to 120k.
                raise ValueError("ENA source full FASTQ has fewer than 120k paired reads")
        if got is None:
            raise ValueError("Official ENA bounded 64MiB prefix cannot yield 120k complete reads")
        mates.append(got)
    if len(mates)!=2 or any(len(x)!=4*MIN_LARGE_PAIRS for x in mates):
        raise ValueError("Official 120k mate pair counts not identical")
    return mates,details


def validate_sources(original):
    if len(original)!=SOURCE_N or len({r["run"] for r in original})!=SOURCE_N:
        raise ValueError("Source original 125-run census missing")
    selected={row["run"]:row for row in original if row["run"] in RUNS}
    if set(selected)!=set(RUNS):
        raise ValueError("Four source-verified Cirsium runs not all present")
    for run,(biosample,name,category) in RUNS.items():
        row=selected[run]
        original_name="_".join(row["S1_source_taxon"].replace(" ","_").split())
        if row["biosample"]!=biosample or original_name!=name:
            raise ValueError("Original source SRA BioSample or taxon identity changed")
        tree_label=row["original_tree_code"].replace(" ","_").split("_")
        if len(tree_label)<2 or "_".join(tree_label[:2])!=name:
            raise ValueError("Original SRA submitted voucher label no longer matches")
    return selected


def map_read_pairs(row,reference,loci,reads,path):
    """Reuse the same original 241-locus site depth functionality."""
    path.mkdir(parents=True,exist_ok=True)
    fastqs=[]
    for i,contents in enumerate(reads,start=1):
        filepath=path/(row["run"]+f"_{i}.fastq")
        filepath.write_bytes(b"".join(contents))
        fastqs.append(filepath)
    bam=path/"aligned.sorted.bam"
    bwa=subprocess.Popen(["bwa","mem","-t","2",str(reference),*(str(x) for x in fastqs)],
                         stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    try:
        sorter=subprocess.Popen(["samtools","sort","-@","2","-m","128M","-o",str(bam),"-"],
                                stdin=bwa.stdout,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        bwa.stdout.close()
        _o,sort_error=sorter.communicate(timeout=240)
        _o,bwa_error=bwa.communicate(timeout=240)
        if sorter.returncode!=0 or bwa.returncode!=0:
            raise RuntimeError("Source real read mapping failed: "+
                (sort_error+bwa_error).decode("utf-8","replace")[-320:])
    finally:
        if bwa.poll() is None:bwa.kill()
    dp.upstream.run_commands(["samtools","index",str(bam)])
    idx=dp.upstream.run_commands(["samtools","idxstats",str(bam)]).stdout.decode()
    lengths,refs=dp.reference_lengths(idx)
    process=subprocess.Popen(["samtools","depth","-q","20","-Q","20",str(bam)],
                             stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    try:
        depth=dp.count_depth_stream(process.stdout,lengths)
        process.stdout.close()
        errors=process.stderr.read()
        process.wait(timeout=180)
        if process.returncode:
            raise RuntimeError("Samtools real read depth failed: "+errors[-320:])
    finally:
        if process.poll() is None:
            process.kill()
            process.wait()
    detail=dp.summarize_loci(loci,lengths,refs,depth)
    return {"input_read_pairs":len(reads[0])//4,
            **dp.per_run_metrics(detail),
            "read_depth_not_orthology":True,
            "individual_locus_reference_choice_is_optimistic":True}


def run(source,ref,loci,scratch,acquire=acquire_larger_prefix,mapper=map_read_pairs):
    if len(loci)!=LOCI or len(set(loci))!=LOCI:
        raise ValueError("Expected original author-verified 241-locus set")
    selected=validate_sources(source)
    results=[]
    for run in RUNS:
        row=selected[run]
        root=Path(scratch)/run
        root.mkdir(parents=True,exist_ok=True)
        try:
            mates,meta=acquire(row)
            if len(mates)!=2 or any(len(z)!=4*MIN_LARGE_PAIRS for z in mates):
                raise ValueError("Bounded ENA paired reads incomplete")
            small=[records[:4*PAIRS[0]] for records in mates]
            small_stats=mapper(row,ref,loci,small,root/"30k")
            large_stats=mapper(row,ref,loci,mates,root/"120k")
            if small_stats["input_read_pairs"]!=PAIRS[0] or \
               large_stats["input_read_pairs"]!=PAIRS[1]:
                raise ValueError("Empirical paired depth comparison not exact 30k/120k")
            for key in ("n_loci_dp3_ge25pct","n_loci_dp3_ge50pct",
                        "n_loci_dp3_100bp","total_best_reference_sites_dp3"):
                if large_stats[key]<small_stats[key]:
                    raise ValueError("Quality-filtered site-depth metric decreased on nested source reads: "+key)
            results.append({
              "run":run,"biosample":row["biosample"],
              "S1_accepted_key":row["accepted_key"],
              "source_taxon":row["S1_source_taxon"],
              "original_head_state_source_context":RUNS[run][2],
              "ENA_partial_gzip_source_report":meta,
              "small_first_30k":small_stats,
              "large_first_120k":large_stats,
              "depth_not_complete_genes":True,
            })
            print("T2_120K_SAME_RUN_DEPTH_PASS="+json.dumps({
              "run":run,"Q20_dp3_50pct_30k":small_stats["n_loci_dp3_ge50pct"],
              "Q20_dp3_50pct_120k":large_stats["n_loci_dp3_ge50pct"],
              "Q20_dp3_25pct_30k":small_stats["n_loci_dp3_ge25pct"],
              "Q20_dp3_25pct_120k":large_stats["n_loci_dp3_ge25pct"]}),flush=True)
        finally:
            shutil.rmtree(root,ignore_errors=True)
    return {
      "status":"REAL_MORPHOLOGICALLY_INFORMATIVE_FOUR_RUNS_WITHIN_RUN_30K_120K_DEPTH_EXPERIMENT",
      "full_S1_accepted_key_census_still_244":True,
      "full_S1_source_names_still_251":True,
      "full_EAzami_nuclear_species_tree_inferred":False,
      "no_direct_head_adaptation_inference":True,
      "quality_filter":"baseQ20,mapQ20,author241,max-one-original-reference-per-locus",
      "same_run_nesting":True,
      "read_counts_compared":list(PAIRS),
      "maximum_gzip_prefix_bytes_per_mate":max(MAX_BYTES),
      "full_fastq_MD5_verified":False,
      "all_source_reads_deleted_before_artifact":True,
      "pilot_taxa":results,
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source",type=Path,required=True)
    parser.add_argument("--loci",type=Path,required=True)
    parser.add_argument("--reference",type=Path,required=True)
    parser.add_argument("--scratch",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    a=parser.parse_args()
    source=upstream.load_source(a.source)
    loci=dp.loci_from_pinned_author_list(a.loci)
    out=run(source,a.reference,loci,a.scratch)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(out,ensure_ascii=False,indent=2,
                                  sort_keys=True,allow_nan=False)+"\n",encoding="utf-8")
    print("T2_SOURCE_ACTUAL_4X_DEPTH_SCALING="+json.dumps({
        "runs":len(out["pilot_taxa"]),
        "quality50_30k":[x["small_first_30k"]["n_loci_dp3_ge50pct"] for x in out["pilot_taxa"]],
        "quality50_120k":[x["large_first_120k"]["n_loci_dp3_ge50pct"] for x in out["pilot_taxa"]],
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
