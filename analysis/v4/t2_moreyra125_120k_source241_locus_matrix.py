"""Actual all-125 source-run, original-241 locus-specific 120k Q20 audit.

Scale the source-verifiable successful SAME 123 public SRA runs of the prior
125-run 30k-vs-120k Q20 experiment. The two documented source archives with
fewer than 120k real paired reads are recorded as NOT_MEASURED, never 0 genes.

For successful runs, retrieve EXACT same original ENA first120k paired FASTQs,
use the original Compositae1061 reference (Git blob 4f89e...) and original
author no-warning 241 list (SHA d561...), and BWA/samtools baseQ20 MAPQ20.
Emit 241 *per-locus* depth-breadth and chosen single reference-copy labels,
not putative independently certified orthologs and NOT a nuclear species tree.
No sequence read/BAM or coordinate files are archived.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import shutil
import time

from analysis.v4 import t2_moreyra_all125_ena_read_mapping as original
from analysis.v4 import t2_moreyra_all125_author241_site_depth as author
from analysis.v4 import t2_four_species_30k_120k_nuclear_depth as larger

RUNS=125
TAXA=113
LOCI=241
SHARDS=10
PAIRS=120000
PRIOR_UNAVAILABLE={"SRR30887185","SRR30887248"}
METRICS=("n_loci_dp3_ge25pct","n_loci_dp3_ge50pct",
         "n_loci_dp3_100bp","total_best_reference_sites_dp3")


def read_csv(path):
    with Path(path).open(newline="",encoding="utf-8-sig") as f:
        rows=list(csv.DictReader(f))
    if not rows:
        raise ValueError("Missing original source 120k empirical CSV")
    return rows


def certify_source(source,original_120k,loci):
    if len(source)!=RUNS or len({r["run"] for r in source})!=RUNS or \
       len({r["accepted_key"] for r in source})!=TAXA:
        raise ValueError("Original 125-run / 113 accepted-key source lost")
    if len(original_120k)!=RUNS or len({r["run"] for r in original_120k})!=RUNS:
        raise ValueError("Original completed 125-run 120k source audit missing")
    if len(loci)!=LOCI or len(set(loci))!=LOCI:
        raise ValueError("Source-pinned author 241 locus IDs changed")
    byrun={r["run"]:r for r in original_120k}
    if set(byrun)!={r["run"] for r in source}:
        raise ValueError("Original 30k/120k run inventory not identical")
    failed=set()
    for r in source:
        old=byrun[r["run"]]
        if any(r[k]!=old[k] for k in ("biosample","accepted_key")):
            raise ValueError("Public SRA BioSample or taxon identity changed")
        if old["full_120k_comparison_status"]!="real_ENA_120k_quality20_depth":
            if old["full_120k_comparison_status"]!="failed":
                raise ValueError("Unknown prior source 120k sample status")
            failed.add(r["run"])
        else:
            if old["source_30k_to_120k_valid"].lower()!="true":
                raise ValueError("Published 120k mapping claimed valid against original false gate")
            for k in METRICS:
                x=int(old["first120k_"+k])
                if x<0 or (k.startswith("n_loci_") and x>LOCI):
                    raise ValueError("Original measured 120k nuclear breadth out of range")
    if failed!=PRIOR_UNAVAILABLE:
        raise ValueError("Two original short ENA source runs changed: "+str(sorted(failed)))
    return byrun


def measure_one(row,old,loci,reference,workdir,
                acquisition=larger.acquire_larger_prefix,
                mapper=larger.map_read_pairs):
    folder=Path(workdir)/row["run"]
    folder.mkdir(parents=True,exist_ok=True)
    try:
        mates,_metadata=acquisition(row)
        if len(mates)!=2 or any(len(x)!=PAIRS*4 for x in mates):
            raise ValueError("Missing original matched 120k pair reads")
        result=mapper(row,reference,loci,mates,folder/"align",include_per_locus=True)
        if int(result["input_read_pairs"])!=PAIRS:
            raise ValueError("Real 120k original accession depth uses wrong read-pair count")
        for field in METRICS:
            if int(result[field])!=int(old["first120k_"+field]):
                raise ValueError("Same source 120k Q20 depth does not replicate earlier verified audit: "+field)
        detailed=result.get("per_locus_241")
        if (not isinstance(detailed,list) or len(detailed)!=LOCI or
            {x["locus"] for x in detailed}!=set(loci)):
            raise ValueError("Not all original author241 per-locus true counts recovered")
        for x in detailed:
            if (x["reference_choice_maximizes_coverage_not_orthology"] is not True
                or not 0<=x["dp5_breadth"]<=x["dp3_breadth"]<=x["dp1_breadth"]<=1):
                raise ValueError("Invalid source genomic depth quality or reference copy")
        return {
            "run":row["run"],"biosample":row["biosample"],"accepted_key":row["accepted_key"],
            "source_name":row["S1_source_taxon"],
            "status":"real_ENA_120k_same_run_author241_locus_Q20",
            "original_120k_pairs":PAIRS,
            "depth_quality":"samtools baseQ>=20, mappingQ>=20",
            **{k:int(result[k]) for k in METRICS},
            "per_locus_241":detailed,
            "individual_gene_consensus_reconstructed":False,
            "true_phylogenetic_orthology_certified":False,
        }
    finally:
        shutil.rmtree(folder,ignore_errors=True)


def execute(source,prior,loci,reference,workdir,shard,nshards=SHARDS,runner=measure_one):
    if not 0<=shard<nshards==SHARDS:
        raise ValueError("Original all125 full ten-shard schedule required")
    old=certify_source(source,prior,loci)
    selected=[r for i,r in enumerate(source) if i%nshards==shard]
    output=[]
    for r in selected:
        before=old[r["run"]]
        if r["run"] in PRIOR_UNAVAILABLE:
            output.append({
                "run":r["run"],"biosample":r["biosample"],
                "accepted_key":r["accepted_key"],"source_name":r["S1_source_taxon"],
                "status":"original_ENA_prefix_below_120k_pairs",
                "original_120k_pairs":None,
                "not_zero_loci_or_gene_loss":True,
                "prior_explanation":before["source_run_unavailability_reason"],
            })
            continue
        time0=time.monotonic()
        try:
            record=runner(r,before,loci,reference,workdir)
            if any(record.get(k)!=r[k] for k in ("run","biosample","accepted_key")):
                raise ValueError("Actual source specimen failed accession identity check")
            record["elapsed_seconds"]=round(time.monotonic()-time0,2)
            output.append(record)
            print("T2_ALL125_120K_REAL_AUTHOR241_LOCI_OK="+json.dumps({
                "run":r["run"],"n_Q20_dp3_half":record["n_loci_dp3_ge50pct"]}),flush=True)
        except Exception as err:
            output.append({
                "run":r["run"],"biosample":r["biosample"],"accepted_key":r["accepted_key"],
                "status":"failed","error":type(err).__name__+": "+str(err)[:350],
            })
            print("T2_ALL125_120K_REAL_AUTHOR241_LOCI_FAILED="+json.dumps(output[-1]),flush=True)
    if len(output)!=len(selected) or {r["run"] for r in output}!={r["run"] for r in selected}:
        raise ValueError("Missing original run in 120k shared-locus matrix")
    return {
        "schema":"azami_v4_Moreyra125_120k_real241_per_locus_Q20_v1",
        "shard":shard,"total_shards":nshards,"n_original_125_run_source":RUNS,
        "n_original_113_photo_taxon_candidates":TAXA,
        "n_author_no_warning_loci":LOCI,
        "requested":len(selected),
        "mapped":sum(x["status"]=="real_ENA_120k_same_run_author241_locus_Q20" for x in output),
        "original_short_runs":sum(x["status"]=="original_ENA_prefix_below_120k_pairs" for x in output),
        "unexpected_failures":sum(x["status"]=="failed" for x in output),
        "missing_is_not_genetic_absence":True,
        "published_full_FASTQ_MD5_not_verified":True,
        "only_bounded_first120k_not_random":True,
        "raw_FASTQ_or_BAM_in_artifact":False,
        "species_tree_inferred":False,
        "outcomes":output,
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source",required=True,type=Path)
    p.add_argument("--original120k",required=True,type=Path)
    p.add_argument("--loci",required=True,type=Path)
    p.add_argument("--reference",required=True,type=Path)
    p.add_argument("--shard",required=True,type=int)
    p.add_argument("--scratch",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    rows=original.load_source(a.source)
    old=read_csv(a.original120k)
    loci=author.loci_from_pinned_author_list(a.loci)
    if not a.reference.is_file():raise ValueError("Original Compositae1061 reference not recovered")
    report=execute(rows,old,loci,a.reference,a.scratch,a.shard)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(report,sort_keys=True,indent=2,ensure_ascii=False,
                                 allow_nan=False)+"\n",encoding="utf-8")
    print("T2_ALL125_120K_SITE_LOCI_SHARD="+json.dumps({
        "shard":a.shard,"requested":report["requested"],"mapped":report["mapped"],
        "short":report["original_short_runs"],"failed":report["unexpected_failures"]},
        sort_keys=True),flush=True)


if __name__=="__main__":
    main()
