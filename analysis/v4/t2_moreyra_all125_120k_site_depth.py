"""Scale source-verified 120k real ENA Q20-author241 depth to ALL 125 S1 runs.

Proof question: which originally linked 113 taxon candidates acquire enough
callable original low-paralog loci once read depth quadruples relative to
the independently successful original 30k all125 source experiment?

This is 125 full *source accessions*, not whole 294 EAzami samples and not
orthologous consensus sequences. Save only per-run numeric metrics; erase
FASTQ/BAM for every run. Preserve technical failures and biological gaps.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import shutil
import time

from analysis.v4 import t2_moreyra_all125_ena_read_mapping as source
from analysis.v4 import t2_moreyra_all125_author241_site_depth as depth
from analysis.v4 import t2_four_species_30k_120k_nuclear_depth as more

RUNS=125
TAXA=113
SHARDS=10
TARGET_PAIRS=120000
METRICS=("n_loci_dp1_any","n_loci_dp3_any","n_loci_dp5_any",
         "n_loci_dp3_100bp","n_loci_dp3_ge10pct","n_loci_dp3_ge25pct",
         "n_loci_dp3_ge50pct","n_loci_dp3_ge75pct",
         "total_best_reference_bp","total_best_reference_sites_dp1",
         "total_best_reference_sites_dp3","total_best_reference_sites_dp5")


def execute(source_rows,reference,loci,scratch,shard,*,nshards=10,
            acquire=more.acquire_larger_prefix,mapper=more.map_read_pairs):
    if not 0<=shard<nshards==SHARDS:
        raise ValueError("Full exactly 10-shard schedule mandatory")
    if len(source_rows)!=RUNS or len({x["run"] for x in source_rows})!=RUNS or \
       len({x["accepted_key"] for x in source_rows})!=TAXA:
        raise ValueError("Source 125 run / 113 accepted source keys no longer intact")
    if len(loci)!=depth.LOCUS_COUNT or len(set(loci))!=len(loci):
        raise ValueError("Source-author original 241 unique candidate loci required")
    chosen=[r for i,r in enumerate(source_rows) if i%nshards==shard]
    results=[]
    for row in chosen:
        t0=time.monotonic()
        folder=Path(scratch)/row["run"]
        folder.mkdir(parents=True,exist_ok=True)
        try:
            mates,_meta=acquire(row)
            if len(mates)!=2 or any(len(z)!=TARGET_PAIRS*4 for z in mates):
                raise ValueError("120k complete paired ENA original reads are not available")
            info=mapper(row,reference,loci,mates,folder/"align")
            if info["input_read_pairs"]!=TARGET_PAIRS or \
               any(int(info[k])<0 for k in METRICS):
                raise ValueError("Invalid bounded 120k Q20 source mapping result")
            rowout={
                "run":row["run"],"biosample":row["biosample"],
                "accepted_key":row["accepted_key"],
                "S1_source_taxon":row["S1_source_taxon"],
                "original_tree_code":row["original_tree_code"],
                "status":"real_ENA_120k_quality20_depth",
                "read_pairs":info["input_read_pairs"],
                **{k:int(info[k]) for k in METRICS},
                "quality_is_baseQ20_and_mappingQ20":True,
                "source_first_prefix_not_random":True,
                "original_241_source_loci_are_not_completed_orthologs":True,
                "single_specimen_copy_ambiguity_not_resolved":True,
            }
            results.append(rowout)
            print("T2_120K_ALL125_REAL_Q20_RUN_OK="+json.dumps({
                "run":rowout["run"],
                "Q20_dp3_ge25pct":rowout["n_loci_dp3_ge25pct"],
                "Q20_dp3_ge50pct":rowout["n_loci_dp3_ge50pct"],
            }),flush=True)
        except Exception as err:
            result={"run":row["run"],"biosample":row["biosample"],
                    "accepted_key":row["accepted_key"],
                    "status":"failed",
                    "error_type":type(err).__name__,
                    "error_message":str(err)[:500]}
            results.append(result)
            print("T2_120K_ALL125_REAL_Q20_RUN_UNAVAILABLE="+json.dumps(result),flush=True)
        finally:
            shutil.rmtree(folder,ignore_errors=True)
    if {r["run"] for r in results}!={r["run"] for r in chosen}:
        raise ValueError("Source run IDs lost during full S1 nuclear depth scale-up")
    return {
        "schema":"azami_v4_all125_original_S1_run_120k_Q20_author241_site_depth_v1",
        "shard":shard,"shard_count":nshards,
        "source_original_S1_accepted_keys":244,
        "source_S1_source_name_units":251,
        "source_linked_real_sra_runs":RUNS,
        "source_name_matched_taxa":TAXA,
        "requested":len(chosen),
        "successes":sum(x["status"].startswith("real_ENA_") for x in results),
        "failures":sum(x["status"]=="failed" for x in results),
        "read_pairs_per_succeeded_run":TARGET_PAIRS,
        "original_author_241_loci":depth.LOCUS_COUNT,
        "quality":"samtools depth baseQ20 mapQ20; best original reference copy per locus",
        "nuclear_tree_generated":False,
        "raw_FASTQ_or_BAM_saved_in_artifact":False,
        "outcomes":results,
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--loci",type=Path,required=True)
    p.add_argument("--reference",type=Path,required=True)
    p.add_argument("--scratch",type=Path,required=True)
    p.add_argument("--shard",type=int,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    src=source.load_source(a.source)
    loci=depth.loci_from_pinned_author_list(a.loci)
    if not a.reference.is_file():
        raise ValueError("Source pinned original Compositae1061 reference absent")
    report=execute(src,a.reference,loci,a.scratch,a.shard)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(report,indent=2,sort_keys=True,
                                ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    print("T2_FULL125_ACTUAL_120K_DEPTH_SHARD="+json.dumps({
        k:report[k] for k in ("shard","requested","successes","failures")
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
