"""Audit original 241 common nuclear-target loci across all 113 S1 taxa and 7 deep taxa.

Inputs: complete 125 exact run x 241 per-locus Q20/30k empirical shards and
same-run actual 120k depth results, NOT reconstructed sequences or nuclear tree.
At any taxon, multiple public runs contribute only taxon-level max-depth
PRESENCE, not sample pooling or extra phylogenetic tips.

Shared named Compositae1061 targets are a candidate alignment *universe* only:
the best original reference copy may differ between taxa and taxonomic identity
remains unverified. No inference from 30k sample loci to 120k joint occupancy.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
from collections import defaultdict,Counter

N_RUNS=125
N_TAXA=113
N_LOCI=241
N_SCENARIOS=10
TARGET_N=7
TARGET_NAMES={
 "Cirsium_vulgare","Cirsium_tracyi","Cirsium_pulcherrimum",
 "Cirsium_italicum","Cirsium_flodmanii","Cirsium_latifolium",
 "Cirsium_scopulorum"
}
DEPTHS=(.25,.5)
def norm(x):
    return "_".join(str(x).strip().split()).replace(" ","_")

def read_csv(path):
    with Path(path).open(newline="",encoding="utf-8-sig") as fp:
        rows=list(csv.DictReader(fp))
    if not rows:
        raise ValueError("Required source-verified Q20 original archive missing")
    return rows

def read_perrun(root):
    paths=sorted(Path(root).rglob("shard_*.json"))
    if len(paths)!=N_SCENARIOS:
        raise ValueError("Ten original empirically measured 30k read run shards required")
    reports=[json.loads(p.read_text(encoding="utf-8")) for p in paths]
    if {x["shard_id"] for x in reports}!=set(range(N_SCENARIOS)):
        raise ValueError("Duplicate/missing empirical 30k shard ID")
    by_run={}
    for report in reports:
        if report["read_data_in_artifact"] or report["BAM_in_artifact"] or report["n_shards"]!=10:
            raise ValueError("Artifact improperly includes raw genomic reads or wrong source")
        for r in report["outcomes"]:
            run=r["run"]
            if run in by_run or r["status"]!="real_ENA_reference_site_depth_q20_Q20":
                raise ValueError("Missing/duplicate mapped original 30k run")
            if len(r["loci"])!=N_LOCI:
                raise ValueError("Author 241 data absent from source")
            if len({x["locus"] for x in r["loci"]})!=N_LOCI:
                raise ValueError("Author241 locus duplicates within a run")
            by_run[run]=r
    if len(by_run)!=N_RUNS or len({x["accepted_key"] for x in by_run.values()})!=N_TAXA:
        raise ValueError("Original 125-run / 113-taxon numerical source incomplete")
    return by_run

def compute(runs,joined120):
    if len(runs)!=N_RUNS or len(joined120)!=N_RUNS:
        raise ValueError("Expected full 125 empirical paired-source runs")
    bysource={r["run"]:r for r in joined120}
    if set(bysource)!=set(runs):
        raise ValueError("30k per-locus and 120k same-source run IDs are not identical")
    taxa=defaultdict(list)
    for accession,row in runs.items():
        actual=bysource[accession]
        if (row["accepted_key"]!=actual["accepted_key"] or
            row["biosample"]!=actual["biosample"]):
            raise ValueError("Identical original BioSample accession pairing not verified")
        if int(actual["first30k_n_loci_dp3_ge50pct"])!=int(row["n_loci_dp3_ge50pct"]):
            raise ValueError("Original 30k per-locus counts contradict independent full125 summary")
        taxa[row["accepted_key"]].append(row)
    names={}
    for key,rr in taxa.items():
        source={norm(x["S1_source_taxon"]) for x in rr}
        if len(source)!=1:
            raise ValueError("WCVP matched key has multiple original source names")
        names[key]=next(iter(source))
    matching={name:key for key,name in names.items() if name in TARGET_NAMES}
    if set(matching)!=TARGET_NAMES:
        raise ValueError("Seven original 120k deep sample taxa do not all have per-locus Q20 source evidence")
    chosen={}
    for name,key in sorted(matching.items()):
        data=taxa[key]
        successes=[x for x in data if bysource[x["run"]]["source_30k_to_120k_valid"]=="True"]
        if not successes:
            raise ValueError("Seven deep species must have source verified 120k run")
        best=max(successes,key=lambda r:(
            int(bysource[r["run"]]["first120k_n_loci_dp3_ge50pct"]),r["run"]))
        count=int(bysource[best["run"]]["first120k_n_loci_dp3_ge50pct"])
        if count<100:
            raise ValueError("Seven high-coverage species membership changed")
        chosen[name]=best
    expected_loci={z["locus"] for z in next(iter(runs.values()))["loci"]}
    if len(expected_loci)!=N_LOCI or any({z["locus"] for z in x["loci"]}!=expected_loci
         for x in runs.values()):
        raise ValueError("Mixed different original author low-paralog 241-locus sets")
    lookup={run:{l["locus"]:l for l in row["loci"]} for run,row in runs.items()}
    matrix=[]
    for locus in sorted(expected_loci):
        best_for_taxon={}
        for key,rr in taxa.items():
            rmax=max((lookup[x["run"]][locus] for x in rr),
                     key=lambda z:(z["dp3_breadth"],z["dp1_breadth"]))
            best_for_taxon[key]=rmax
        seven={name:lookup[r["run"]][locus] for name,r in chosen.items()}
        for rec in seven.values():
            if not rec["reference_choice_maximizes_coverage_not_orthology"]:
                raise ValueError("Single original reference per locus safeguard not maintained")
        refs={v["best_original_reference"] for v in seven.values()}
        matrix.append({
          "locus":locus,
          "n_original_source_taxa":len(best_for_taxon),
          "n_original_taxa_30k_dp3_ge25pct":sum(x["dp3_breadth"]>=.25 for x in best_for_taxon.values()),
          "n_original_taxa_30k_dp3_ge50pct":sum(x["dp3_breadth"]>=.5 for x in best_for_taxon.values()),
          "n_deep7_30k_dp3_ge25pct":sum(x["dp3_breadth"]>=.25 for x in seven.values()),
          "n_deep7_30k_dp3_ge50pct":sum(x["dp3_breadth"]>=.5 for x in seven.values()),
          "deep7_all_same_single_original_reference":len(refs)==1,
          "deep7_seven_refs":"|".join(sorted(refs)),
          "deep7_30k_four_fifths_taxa_qualify_25pct":sum(x["dp3_breadth"]>=.25 for x in seven.values())>=6,
          "deep7_30k_all7_qualify_25pct":all(x["dp3_breadth"]>=.25 for x in seven.values()),
          "not_120k_common_locus_coverage":True,
          "locus_is_original_author_241_candidate_not_ortholog_confirmation":True,
        })
    selected=[{"source_taxon":k,"accepted_key":r["accepted_key"],
               "run":r["run"],"biosample":r["biosample"],
               "same_run_120k_dp3_ge50pct":int(
                   bysource[r["run"]]["first120k_n_loci_dp3_ge50pct"]),
               "source_name_candidate_voucher_unverified":True}
              for k,r in sorted(chosen.items())]
    summary={
       "schema":"azami_v4_all113_source_taxa_and_seven_120k_deep_Q20_original241_shared_locus_screen_v1",
       "original_source_taxa":N_TAXA,
       "original_public_SRA_runs":N_RUNS,
       "all_original_candidate_loci":N_LOCI,
       "source_S1_photo_WCVP_keys_retained_in_master":244,
       "source_S1_original_name_units_retained_in_master":251,
       "n_seven_deep_source_taxa":TARGET_N,
       "seven_real_120k_selected_source_specimens":selected,
       "original_30k_loci_all7_ge25pct":sum(x["deep7_30k_all7_qualify_25pct"] for x in matrix),
       "original_30k_loci_at_least6of7_ge25pct":sum(x["deep7_30k_four_fifths_taxa_qualify_25pct"] for x in matrix),
       "original_30k_loci_all7_ge25pct_same_ref":sum(
           x["deep7_30k_all7_qualify_25pct"] and x["deep7_all_same_single_original_reference"] for x in matrix),
       "original_30k_loci_113taxa_at_least90_ge25pct":sum(
           x["n_original_taxa_30k_dp3_ge25pct"]>=90 for x in matrix),
       "all7_real120k_locus_overlap_not_calculated":True,
       "reference_copy_consistency_does_not_resolve_orthology":True,
       "nuclear_tree_fit_allowed_from_breadth_only":False,
       "true_biological_posture_of_7_unscored_taxa_not_inferred":True,
       "n_loci_reused_from_author_350_not_assumed":True,
    }
    return matrix,selected,summary

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-shards",type=Path,required=True)
    p.add_argument("--paired120k",type=Path,required=True)
    p.add_argument("--out-dir",type=Path,required=True)
    a=p.parse_args()
    source=read_perrun(a.source_shards)
    pair=read_csv(a.paired120k)
    data,taxa,s=compute(source,pair)
    a.out_dir.mkdir(parents=True,exist_ok=True)
    for name,rows in (
        ("AUTHOR241_SHARED_DEPTH_30K_ACROSS113_AND_DEEP7.csv",data),
        ("REAL120K_DEEP7_SELECTED_SPECIMENS.csv",taxa)):
        with (a.out_dir/name).open("w",newline="",encoding="utf-8") as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),extrasaction="raise")
            w.writeheader();w.writerows(rows)
    (a.out_dir/"AUTHOR241_SHARED_DEPTH_READYNESS.json").write_text(
        json.dumps(s,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    print("T2_SHARED241_ACTUAL_30K_113TAXA_AND_SEVEN120K_NUCLEAR_TARGETS="+json.dumps(s,sort_keys=True,ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()
