"""Build capacity-aware full-target SRA accession shards from authentic RunInfo.

A shard limit is a *compressed SRA scheduling ceiling*, NOT a guarantee of
FASTQ expansion, peak RAM or locus inference. All 113 candidate S1 species
remain in the registry and all 125 runs are accounted for. Oversize groups
are explicitly routed to HPC, never silently excluded.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import json
import math
from pathlib import Path

RUN_N=125
TAXON_N=113
SUM_ARCHIVE_MB=57630.0
DEFAULT_GROUP_CAP_MB=4000.0


def read_source(path):
    with Path(path).open(newline="",encoding="utf-8-sig") as f:
        records=list(csv.DictReader(f))
    required={"S1_accepted_key","source_taxon","SRA_run","BioSample",
              "size_MB_compressed_SRA_metadata","source_taxon_name_disagreement",
              "voucher_accepted_concept_verified"}
    if len(records)!=RUN_N or not required<=set(records[0]):
        raise ValueError("125-source-run sizing input missing/changed")
    keys=set()
    records2=[]
    for r in records:
        run=r["SRA_run"]
        if run in keys or not run.startswith("SRR"):
            raise ValueError("Duplicate or invalid SRA run accession")
        keys.add(run)
        v=float(r["size_MB_compressed_SRA_metadata"])
        if not math.isfinite(v) or v<=0:
            raise ValueError("Source NCBI archived compressed size invalid")
        if r["voucher_accepted_concept_verified"] not in {"False","false","0"}:
            raise ValueError("Voucher-to-accepted-concept not yet verified: stop")
        records2.append({**r,"archive_MB":v})
    if (len({r["S1_accepted_key"] for r in records2})!=TAXON_N
            or abs(sum(r["archive_MB"] for r in records2)-SUM_ARCHIVE_MB)>1e-5):
        raise ValueError("Frozen authenticated 113 candidate taxa / 57.63GB changed")
    return records2


def schedule(rows,cap=DEFAULT_GROUP_CAP_MB):
    if cap<=0:
        raise ValueError("Invalid compressed SRA group limit")
    bytax=defaultdict(list)
    for r in rows:
        bytax[r["S1_accepted_key"]].append(r)
    tasks=[]
    for accepted_key,runs in bytax.items():
        total=sum(float(r["archive_MB"]) for r in runs)
        if len({r["source_taxon"] for r in runs})!=1:
            raise ValueError("Accepted-key SRA run group combines different source names")
        tasks.append({"accepted_key":accepted_key,
                      "source_taxon":runs[0]["source_taxon"],
                      "n_runs":len(runs),
                      "archive_MB":total,
                      "needs_HPC_oversize":total>cap,
                      "source_name_conflict_any":
                         any(str(r["source_taxon_name_disagreement"])=="True"
                             for r in runs)})
    tasks.sort(key=lambda x:(-x["archive_MB"],x["accepted_key"]))
    bins=[]
    over=0
    for task in tasks:
        if task["needs_HPC_oversize"]:
            over+=1
            task["assigned_shard"]=f"HPC_oversize_{over:03d}"
            task["note"]="large public SRA compressed taxon group; require verified disk/streaming and voucher check"
            continue
        ix=next((i for i,b in enumerate(bins)
                 if b["sum_MB"]+task["archive_MB"]<=cap),None)
        if ix is None:
            bins.append({"sum_MB":0.0,"taxa":[]})
            ix=len(bins)-1
        bins[ix]["sum_MB"]+=task["archive_MB"]
        bins[ix]["taxa"].append(task["accepted_key"])
        task["assigned_shard"]=f"standard_compressed_cap_{ix+1:03d}"
        task["note"]="capacity-planning-only: compressed size does not guarantee FASTQ peak disk"
    if len(tasks)!=TAXON_N:
        raise ValueError("113 S1 candidate taxa not retained in shard plan")
    map_key={x["accepted_key"]:x for x in tasks}
    runplan=[]
    for r in rows:
        task=map_key[r["S1_accepted_key"]]
        runplan.append({
            "accepted_key":r["S1_accepted_key"],
            "source_taxon":r["source_taxon"],
            "SRA_run":r["SRA_run"],
            "BioSample":r["BioSample"],
            "size_MB_compressed_SRA_metadata":r["archive_MB"],
            "shard":task["assigned_shard"],
            "needs_large_scratch":task["needs_HPC_oversize"],
            "verified_accepted_taxonomic_concept":False,
            "FASTQ_processed":False,
            "phylogeny_constructed":False,
        })
    if len(runplan)!=RUN_N or len({x["SRA_run"] for x in runplan})!=RUN_N:
        raise ValueError("Lost SRA accession during capacity planning")
    summary={
        "status":"REAL_125_PUBLIC_RUN_CAPACITY_PLAN_ONLY_NO_FASTQ_OR_TREE",
        "number_source_S1_taxa":TAXON_N,
        "number_original_public_SRA_runs":RUN_N,
        "compressed_SRA_size_total_MB":sum(t["archive_MB"] for t in tasks),
        "target_max_compressed_MB_per_standard_shard":cap,
        "n_normal_capacity_shards":len(bins),
        "n_oversized_species_groups":over,
        "oversized_groups":[{"source_taxon":t["source_taxon"],
                             "accepted_key":t["accepted_key"],
                             "compressed_MB":t["archive_MB"],
                             "n_runs":t["n_runs"]}
                            for t in tasks if t["needs_HPC_oversize"]],
        "standard_shard_archive_MB":[round(b["sum_MB"],3) for b in bins],
        "minimum_scratch_or_FASTQ_peak_storage_estimable_from_SRA_compressed_MB":False,
        "source_taxonomy_voucher_verified":False,
        "original_public_350_locus_tree_reproduced":False,
        "full_244_independent_tree_completed":False,
        "source_interpretation":"113 source name candidates, 125 original public runs; not 113 verified nuclear tips",
    }
    return tasks,runplan,summary


def write_csv(dest,rows):
    with dest.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),extrasaction="raise")
        w.writeheader()
        w.writerows(rows)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ncbi-sizes",type=Path,required=True)
    p.add_argument("--out-dir",type=Path,required=True)
    p.add_argument("--group-compressed-cap-mb",type=float,default=DEFAULT_GROUP_CAP_MB)
    a=p.parse_args()
    src=read_source(a.ncbi_sizes)
    groups,shards,receipt=schedule(src,a.group_compressed_cap_mb)
    a.out_dir.mkdir(parents=True,exist_ok=True)
    write_csv(a.out_dir/"T2_MOREYRA_113_SPECIES_COMPRESSION_SHARDS.csv",groups)
    write_csv(a.out_dir/"T2_MOREYRA_125_SRA_RUN_EXECUTION_QUEUE.csv",shards)
    (a.out_dir/"T2_MOREYRA_SHARD_PLAN_RECEIPT.json").write_text(
        json.dumps(receipt,ensure_ascii=False,sort_keys=True,indent=2)+"\n",
        encoding="utf-8")
    print("T2_MOREYRA_CAPACITY_SHARDS="+json.dumps(receipt,sort_keys=True,
                                                   ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
