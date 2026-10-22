"""Source-verified aggregation of all 125 real 241-locus site-depth audits.

Do not confuse >=3 high-quality mapped bases on a best original reference
copy with orthology-certified Cirsium gene reconstructions or full-species
phylogenetic results. All 244 S1 WCVP keys remain in the frozen census.
"""
from __future__ import annotations
import argparse
from collections import defaultdict,Counter
import csv
import json
from pathlib import Path
import statistics

RUNS=125
N_WCVP=113
SHARD_IDS=set(range(10))
LOCI=241
KEYS=("n_loci_dp1_any","n_loci_dp3_any","n_loci_dp5_any",
      "n_loci_dp3_100bp","n_loci_dp3_ge10pct","n_loci_dp3_ge25pct",
      "n_loci_dp3_ge50pct","n_loci_dp3_ge75pct",
      "total_best_reference_sites_dp1","total_best_reference_sites_dp3",
      "total_best_reference_sites_dp5")
CRITICAL=("Cirsium_conspicuum","Cirsium_scopulorum","Cirsium_arvense","Cirsium_vulgare")


def load_receipts(root):
    found=sorted(Path(root).rglob("shard_*.json"))
    if len(found)!=len(SHARD_IDS):
        raise ValueError("Missing one or more of ten original 125-read shards")
    reports=[json.loads(p.read_text(encoding="utf-8")) for p in found]
    if {r["shard_id"] for r in reports}!=SHARD_IDS:
        raise ValueError("Repeated or missing 125-run shard ID")
    if any(r["n_shards"]!=10 or r["read_data_in_artifact"] or
           r["BAM_in_artifact"] for r in reports):
        raise ValueError("Invalid real-data source proof or raw read disclosure")
    return reports


def validate(source,reports):
    if len(source)!=RUNS or len({x["run"] for x in source})!=RUNS:
        raise ValueError("Original S1 Moreyra 125-source run IDs changed")
    if len({x["accepted_key"] for x in source})!=N_WCVP:
        raise ValueError("Original S1 Moreyra 113 taxon keys changed")
    expected={x["run"]:x for x in source}
    observed={}
    failures=[]
    detail=[]
    for report in reports:
        if report["requested"]!=len(report["outcomes"]):
            raise ValueError("One source SRA run vanished within a result shard")
        for x in report["outcomes"]:
            run=x["run"]
            if run not in expected or run in observed:
                raise ValueError("Wrong or duplicated SRA accession in depth pilot")
            observed[run]=x
            e=expected[run]
            if (x["biosample"]!=e["biosample"] or
                x["accepted_key"]!=e["accepted_key"]):
                raise ValueError("Run/BioSample/S1 accepted key drift")
            if x["status"]=="failed":
                failures.append({"run":run,"reason":x.get("error","")})
                continue
            if (x["status"]!="real_ENA_reference_site_depth_q20_Q20"
                    or not 0<x["real_ena_read_pairs"]<=30000
                    or len(x["loci"])!=LOCI):
                raise ValueError("Real 241 locus depth receipt missing/unverified")
            if len({l["locus"] for l in x["loci"]})!=LOCI:
                raise ValueError("Duplicate/missing public no-warning locus")
            if (x["n_loci_dp3_ge75pct"]>x["n_loci_dp3_ge50pct"] or
                x["n_loci_dp3_ge50pct"]>x["n_loci_dp3_ge25pct"] or
                x["n_loci_dp3_ge25pct"]>x["n_loci_dp3_ge10pct"]):
                raise ValueError("Impossible depth breadth threshold order")
            for key in KEYS:
                if not isinstance(x[key],int) or x[key]<0:
                    raise ValueError("Invalid source real numeric depth count "+key)
            for l in x["loci"]:
                if (not l["reference_choice_maximizes_coverage_not_orthology"] or
                    not 0<=l["dp5_breadth"]<=l["dp3_breadth"]<=l["dp1_breadth"]<=1):
                    raise ValueError("Invalid best-reference high-quality locus breadth")
            detail.append(x)
    if len(observed)!=RUNS:
        raise ValueError("Not all 125 source run IDs represented in depth scan")
    if failures:
        raise ValueError(f"{len(failures)} original public run depth scans failed: "+
                         json.dumps(failures[:6],sort_keys=True))
    first=set(detail[0]["loci"][i]["locus"] for i in range(LOCI))
    if any({l["locus"] for l in x["loci"]}!=first for x in detail):
        raise ValueError("Mixed different 241 locus reference sets across runs")
    return detail


def summarize(source,results):
    if len(results)!=RUNS:
        raise ValueError("All 125 real read outputs required")
    tax=defaultdict(list)
    by_locus=defaultdict(list)
    for r in results:
        tax[r["accepted_key"]].append(r)
        for z in r["loci"]:
            by_locus[z["locus"]].append(z)
    if len(tax)!=N_WCVP or len(by_locus)!=LOCI:
        raise ValueError("Original taxon and author locus universe incomplete")
    # Taxa can have >1 independent biological sample. Report max as optimistic
    # rather than pooling different specimens or treating runs as distinct taxa.
    species=[]
    for key,rx in sorted(tax.items()):
        name={x["S1_source_taxon"] for x in rx}
        if len(name)!=1:
            raise ValueError("Multiple source taxon names on one S1 accepted key")
        species.append({
            "accepted_key":key,"source_taxon":next(iter(name)),
            "n_independent_public_runs":len(rx),
            "median_per_run_loci_dp3_ge25pct":statistics.median(
                 x["n_loci_dp3_ge25pct"] for x in rx),
            "minimum_per_run_loci_dp3_ge25pct":min(
                 x["n_loci_dp3_ge25pct"] for x in rx),
            "maximum_per_run_loci_dp3_ge25pct":max(
                 x["n_loci_dp3_ge25pct"] for x in rx),
            "maximum_per_run_loci_dp3_ge50pct":max(
                 x["n_loci_dp3_ge50pct"] for x in rx),
            "maximum_per_run_loci_dp3_100bp":max(
                 x["n_loci_dp3_100bp"] for x in rx),
            "not_same_sample_sequence_concatenation":True,
        })
    loci=[]
    for locus,rr in sorted(by_locus.items()):
        if len(rr)!=RUNS:
            raise ValueError("Depth summary lacks public runs at a frozen locus")
        loci.append({
            "locus":locus,"n_public_runs":len(rr),
            "n_runs_dp3_ge25pct":sum(x["dp3_breadth"]>=.25 for x in rr),
            "n_runs_dp3_ge50pct":sum(x["dp3_breadth"]>=.5 for x in rr),
            "n_runs_dp3_100bp":sum(x["dp3_at_least100bp"] for x in rr),
            "median_best_ref_dp3_breadth":statistics.median(x["dp3_breadth"] for x in rr),
            "not_verified_single_copy":True,
        })
    per_run=[{k:r[k] for k in ("run","biosample","accepted_key",
                "S1_source_taxon","real_ena_read_pairs",*KEYS)} for r in results]
    def quantiles(rows,name):
        x=sorted(z[name] for z in rows)
        return {"min":min(x),"q25":x[int((len(x)-1)*.25)],
                "median":statistics.median(x),"q75":x[int((len(x)-1)*.75)],
                "max":max(x)}
    name_index={x["S1_source_taxon"].replace(" ","_"):x for x in species}
    critical={name: {k:name_index[name][k] for k in (
        "n_independent_public_runs","maximum_per_run_loci_dp3_ge25pct",
        "maximum_per_run_loci_dp3_ge50pct","maximum_per_run_loci_dp3_100bp")}
        for name in CRITICAL if name in name_index}
    summary={
      "schema":"azami_v4_all125_Moreyra_241_real_site_depth_20261022_v1",
      "status":"ACTUAL_REAL_PUBLIC_SITE_DEPTH_CENSUS_NOT_HYBPIPER_OR_GENUS_TREE",
      "source_S1_original_WCVP_accepted_keys":244,
      "source_S1_original_photos":170066,
      "source_SRA_name_candidate_taxon_keys":N_WCVP,
      "actual_source_BioSample_run_count":RUNS,
      "original_author_no_warning_high_occupancy_locus_universe":LOCI,
      "real_prefix_pairs_total":sum(x["real_ena_read_pairs"] for x in results),
      "quality_filter":"samtools depth -q 20 (base Q) -Q 20 (mapping Q); absent sites are not imputed",
      "per_original_reference_locus_choose_max_single_ref_dp3_coverage":True,
      "dp3_ge25pct_loci_per_run":quantiles(results,"n_loci_dp3_ge25pct"),
      "dp3_ge50pct_loci_per_run":quantiles(results,"n_loci_dp3_ge50pct"),
      "dp3_atleast100bp_loci_per_run":quantiles(results,"n_loci_dp3_100bp"),
      "s2_and_nodding_source_taxa_if_exact_name":critical,
      "run_level_metrics_not_independent_species":True,
      "strict_241_ortholog_alignments_completed":False,
      "eazami294_concatenated_tree_completed":False,
      "eazami294_astral_species_tree_completed":False,
      "real_read_prefix_is_not_random_sample":True,
      "taxon_rank_and_voucher_matching_not_completed":True,
      "possible_next_gate":"If >=100 source 241 loci retain high depth and >=50% breadth for multiple vetted taxa, recover allele-aware consensus, inspect copy-number and orthology, and run paired locus phylogenetic QC.",
    }
    return per_run,species,loci,summary


def write_csv(path,rows):
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),extrasaction="raise")
        w.writeheader()
        w.writerows(rows)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--shards",type=Path,required=True)
    p.add_argument("--out-dir",type=Path,required=True)
    a=p.parse_args()
    with a.source.open(newline="",encoding="utf-8-sig") as f:
        source=list(csv.DictReader(f))
    receipts=load_receipts(a.shards)
    full=validate(source,receipts)
    runs,sp,loci,report=summarize(source,full)
    a.out_dir.mkdir(parents=True,exist_ok=True)
    write_csv(a.out_dir/"ALL125_SOURCE_241_LOCUS_DEPTH_RUNS.csv",runs)
    write_csv(a.out_dir/"ALL113_TAXON_BEST_RUN_241_LOCUS_DEPTH.csv",sp)
    write_csv(a.out_dir/"AUTHOR241_LOCI_ACROSS125_REAL_RUNS.csv",loci)
    (a.out_dir/"ALL125_AUTHOR241_REAL_DEPTH_SUMMARY.json").write_text(
       json.dumps(report,ensure_ascii=False,sort_keys=True,indent=2)+"\n",
       encoding="utf-8")
    print("T2_REAL_ALL125_241_LOCUS_DP3_METRICS="+json.dumps(report,
          ensure_ascii=False,sort_keys=True),flush=True)


if __name__=="__main__":
    main()
