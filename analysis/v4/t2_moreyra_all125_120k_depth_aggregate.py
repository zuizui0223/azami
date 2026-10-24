"""Source-verified all125 original 30k-vs-120k 241-locus Q20 depth census.

All 125 run accessions are retained even if a 120k public gzip prefix is
temporarily unavailable. Never discard failed/missing reads as zero loci or
assert all 113 taxa have a species tree. The same originally documented 30k
run is used for each paired comparison (no pooling independent BioSamples).
"""
from __future__ import annotations

import argparse
from collections import defaultdict,Counter
import csv
import json
from pathlib import Path
import statistics

RUNS=125
TAXA=113
SHARDS=10
MINIMUM_QUALITY_LOCUS_COUNT=100
DEPTH_KEYS=("n_loci_dp3_ge25pct","n_loci_dp3_ge50pct",
            "n_loci_dp3_100bp","n_loci_dp3_ge75pct",
            "total_best_reference_sites_dp3")


def read_csv(path):
    with Path(path).open(newline="",encoding="utf-8-sig") as handle:
        rows=list(csv.DictReader(handle))
    if not rows:
        raise ValueError("Required source-authentic per-run 241-locus data missing")
    return rows


def read_shards(path):
    paths=sorted(Path(path).rglob("shard_*.json"))
    if len(paths)!=SHARDS:
        raise ValueError("Full ten completed 120k source shards are required")
    reports=[json.loads(p.read_text(encoding="utf-8")) for p in paths]
    if {x["shard"] for x in reports}!=set(range(SHARDS)):
        raise ValueError("Missing/duplicate source run shard")
    for x in reports:
        if x["shard_count"]!=SHARDS or x["raw_FASTQ_or_BAM_saved_in_artifact"]:
            raise ValueError("Invalid 125-run SRA prefix or raw sequence artifact proof")
        if x["requested"]!=len(x["outcomes"]):
            raise ValueError("Original source SRA run lost inside 120k shard")
    return reports


def reconcile(original30k,receipts):
    if (len(original30k)!=RUNS or
        len({z["run"] for z in original30k})!=RUNS or
        len({z["accepted_key"] for z in original30k})!=TAXA):
        raise ValueError("Original 125/113 30k actual site-depth cohort changed")
    frozen={x["run"]:x for x in original30k}
    changed={}
    for shard in receipts:
        for row in shard["outcomes"]:
            run=row["run"]
            if run not in frozen or run in changed:
                raise ValueError("Unknown/duplicated 125 source run in expanded 120k sample")
            previous=frozen[run]
            if (row["biosample"]!=previous["biosample"] or
                row["accepted_key"]!=previous["accepted_key"]):
                raise ValueError("Different voucher/BioSample/accepted taxon in 120k test")
            if row["status"] not in {"real_ENA_120k_quality20_depth","failed"}:
                raise ValueError("Invalid source-run quality status")
            if row["status"].startswith("real_"):
                if row["read_pairs"]!=120000:
                    raise ValueError("Non-120k run misreported as comparable")
                for key in DEPTH_KEYS:
                    now=int(row[key]);before=int(previous[key])
                    if now<0 or now>241 and key.startswith("n_loci_"):
                        raise ValueError("Invalid quality depth count for original241")
                    if now<before:
                        raise ValueError("Nested source run lost Q20 covered sites: "+run+" "+key)
            changed[run]=row
    if set(changed)!=set(frozen):
        raise ValueError("125 original source SRA run subset incomplete")
    outcomes=[]
    for accession,old in sorted(frozen.items()):
        now=changed[accession]
        record={"run":accession,"biosample":old["biosample"],
                "accepted_key":old["accepted_key"],
                "source_taxon":old["S1_source_taxon"],
                "original_first30k_read_pairs":int(old["real_ena_read_pairs"]),
                "full_120k_comparison_status":now["status"],
                "source_30k_to_120k_valid":now["status"].startswith("real_"),
                "source_run_unavailability_reason":now.get("error_message",""),
                "no_independent_sample_pooling":True}
        for key in DEPTH_KEYS:
            record["first30k_"+key]=int(old[key])
            record["first120k_"+key]=(int(now[key])
                if now["status"].startswith("real_") else None)
        outcomes.append(record)
    return outcomes


def summarize(runs):
    if len(runs)!=RUNS:
        raise ValueError("Original complete S1 Moreyra 125 source run roster missing")
    completed=[x for x in runs if x["source_30k_to_120k_valid"]]
    failed=[x for x in runs if not x["source_30k_to_120k_valid"]]
    bytax=defaultdict(list)
    for x in runs:
        bytax[x["accepted_key"]].append(x)
    if len(bytax)!=TAXA:
        raise ValueError("All 113 S1 original public genotype name candidates required")
    taxon=[]
    for key,items in sorted(bytax.items()):
        usable=[x for x in items if x["source_30k_to_120k_valid"]]
        source_name={x["source_taxon"] for x in items}
        if len(source_name)!=1:
            raise ValueError("One WCVP source key collapsed two different original names")
        taxon.append({
            "accepted_key":key,
            "source_name":next(iter(source_name)),
            "public_original_run_count":len(items),
            "successfully_depth120k_mapped_runs":len(usable),
            "best_30k_single_run_dp3_ge50pct":max(x["first30k_n_loci_dp3_ge50pct"] for x in items),
            "best_120k_single_run_dp3_ge50pct":(
                max(x["first120k_n_loci_dp3_ge50pct"] for x in usable) if usable else None),
            "best_30k_single_run_dp3_ge25pct":max(x["first30k_n_loci_dp3_ge25pct"] for x in items),
            "best_120k_single_run_dp3_ge25pct":(
                max(x["first120k_n_loci_dp3_ge25pct"] for x in usable) if usable else None),
            "minimum_shared_gene_recovery_not_automatically_certified":True,
        })
    def q(a):
        if not a:return None
        vals=sorted(a)
        return {"min":vals[0],"q25":vals[int((len(vals)-1)*.25)],
                "median":statistics.median(vals),
                "q75":vals[int((len(vals)-1)*.75)],"max":vals[-1]}
    before=[z["first30k_n_loci_dp3_ge50pct"] for z in completed]
    after=[z["first120k_n_loci_dp3_ge50pct"] for z in completed]
    summary={
        "schema":"azami_v4_full125_nested_original_Q20_author241_30k120k_v1",
        "original_S1_WCVP_keys":244,
        "original_S1_source_name_units":251,
        "original_S1_nuclear_candidate_accepted_keys":TAXA,
        "original_same_source_runs_all_retained":RUNS,
        "n_runs_with_recovered_120k_read_depth":len(completed),
        "n_runs_without_complete_120k_prefix":len(failed),
        "technical_failure_accessions":[{"run":x["run"],"reason":x["source_run_unavailability_reason"]} for x in failed],
        "n_source_taxa_with_any_120k_test":sum(x["successfully_depth120k_mapped_runs"]>0 for x in taxon),
        "n_source_taxa_with_best_single_120k_dp3_ge50pct_ge100":sum(
             x["best_120k_single_run_dp3_ge50pct"] is not None and
             x["best_120k_single_run_dp3_ge50pct"]>=MINIMUM_QUALITY_LOCUS_COUNT for x in taxon),
        "paired_same_run_30k_Q20_dp3_ge50pct_loci":q(before),
        "paired_same_run_120k_Q20_dp3_ge50pct_loci":q(after),
        "n_paired_successes_with_higher_depth3_half_locus_coverage":sum(a>b for a,b in zip(after,before)),
        "quality_coverage_not_orthology":True,
        "read_prefix_first_N_not_random_sampling":True,
        "biological_species_tree_completed":False,
        "strict_multiple_phylogenetic_ortholog_alignments_completed":False,
    }
    return taxon,summary


def write_csv(path,rows):
    with path.open("w",newline="",encoding="utf-8") as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]),extrasaction="raise")
        writer.writeheader();writer.writerows(rows)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--original30k",type=Path,required=True)
    p.add_argument("--real120k-shards",type=Path,required=True)
    p.add_argument("--out-dir",type=Path,required=True)
    a=p.parse_args()
    old=read_csv(a.original30k)
    receipts=read_shards(a.real120k_shards)
    merged=reconcile(old,receipts)
    taxa,summary=summarize(merged)
    a.out_dir.mkdir(parents=True,exist_ok=True)
    write_csv(a.out_dir/"T2_ALL125_NESTED_PUBLIC_30K_VS_120K_Q20.csv",merged)
    write_csv(a.out_dir/"T2_ALL113_TAXON_SINGLE_RUN_30K_VS_120K.csv",taxa)
    (a.out_dir/"T2_125_SOURCE_30K120K_Q20_SUMMARY.json").write_text(
        json.dumps(summary,indent=2,ensure_ascii=False,sort_keys=True)+"\n",encoding="utf-8")
    print("T2_ALL125_REAL_NESTED_30K120K_AUTHOR241_COMPARISON="+json.dumps(
        summary,ensure_ascii=False,sort_keys=True),flush=True)


if __name__=="__main__":
    main()
