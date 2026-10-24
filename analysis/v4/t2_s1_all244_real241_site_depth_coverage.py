"""Keep all 244 S1 source taxa while attaching actual 241-locus real read depth.

Unobserved public SRA sequence coverage is NA, NOT 0 loci or no genes.
Actual 30k read-prefix site depth is NOT orthologous sequence recovery or
biological head nodding; this is a cross-sectional availability/precision
ledger for eventual global independent nuclear placement.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import statistics

TOTAL=244
OBS=170066
SOURCE_MATCHED=113
COLLISIONS=7
MINIMUMS=(5,10)


def read_csv(path):
    with Path(path).open(newline="",encoding="utf-8-sig") as handle:
        rows=list(csv.DictReader(handle))
    if not rows:
        raise ValueError("Missing original source data rows")
    return rows


def build(census,depth):
    if len(census)!=TOTAL or len({x["accepted_key"] for x in census})!=TOTAL:
        raise ValueError("Original full-S1-244 taxon source lost or duplicated")
    if len(depth)!=SOURCE_MATCHED or len({x["accepted_key"] for x in depth})!=SOURCE_MATCHED:
        raise ValueError("Incomplete actual callable 241-gene depth source")
    if sum(int(x["eligible_photo_observations"]) for x in census)!=OBS:
        raise ValueError("Original source 170066 photo-level observations changed")
    bykey={x["accepted_key"]:x for x in depth}
    table=[]
    for s in census:
        key=s["accepted_key"]
        collision=s["source_concept_collision"].lower()=="true"
        candidate=s["Moreyra2025_SRA_taxon_name_candidate"].lower()=="true"
        if s["source_concept_collision"].lower() not in ("true","false") or \
           s["Moreyra2025_SRA_taxon_name_candidate"].lower() not in ("true","false"):
            raise ValueError("S1 explicit source taxon/Moreyra matching status missing")
        hit=bykey.get(key)
        if candidate!=(hit is not None):
            raise ValueError("S1 Moreyra matching and actual 125-read depth identity disagree")
        if candidate and collision:
            raise ValueError("WCVP source-name collision cannot become a nuclear sample")
        if hit is not None:
            norm=lambda x:" ".join(str(x).replace("_"," ").split())
            if norm(hit["source_taxon"])!=norm(s["normalized_exact_source_name_candidate"]):
                raise ValueError("Wrong nuclear source taxon attached to S1 accepted key")
        v=lambda k:int(hit[k]) if hit is not None else None
        values={
            "accepted_key":key,"original_S1_taxon":s["source_names"],
            "source_original_photo_count":int(s["eligible_photo_observations"]),
            "native_S1_mean_image_angle_deg":float(s["S1_mean_EXIF_angle_deg"]),
            "S1_BIO12_mean_mm":float(s["chelsa_bio12_photo_eligible_mean"]),
            "source_WCVP_names_conflict":collision,
            "nuclear_original_125_SRA_name_candidate":candidate,
            "n_source_public_SRA_runs":v("n_independent_public_runs"),
            "best_single_run_241_loci_dp3_ge25pct":v("maximum_per_run_loci_dp3_ge25pct"),
            "best_single_run_241_loci_dp3_ge50pct":v("maximum_per_run_loci_dp3_ge50pct"),
            "best_single_run_241_loci_dp3_100bp":v("maximum_per_run_loci_dp3_100bp"),
            "nuclear_read_coverage_observed":bool(hit),
            "nuclear_missing_reason":(
                "" if hit else "WCVP_source_name_collision" if collision
                else "no_exact_name_matched_Moreyra_public_SRA_run"),
            "preliminary_next_stage":(
                "source_name_candidate_requires_orthology_and_voucher_validation"
                if hit and int(hit["maximum_per_run_loci_dp3_ge50pct"])>=100
                else "low_prefix_depth_not_nuclear_gene_absence" if hit
                else "no_qualifying_public_SRA_in_frozen_125_run_cohort"),
            "nuclear_taxon_identity_verified":False,
            "species_nuclear_gene_tree_completed":False,
        }
        if hit and any(v(k)>241 for k in ("maximum_per_run_loci_dp3_ge25pct",
                                         "maximum_per_run_loci_dp3_ge50pct",
                                         "maximum_per_run_loci_dp3_100bp")):
            raise ValueError("More than 241 independent original author loci counted")
        table.append(values)
    if len(table)!=TOTAL or sum(x["source_WCVP_names_conflict"] for x in table)!=COLLISIONS:
        raise ValueError("Original full S1 source alias distribution changed")
    main=[x for x in table if not x["source_WCVP_names_conflict"]]
    high=[x for x in main if x["nuclear_read_coverage_observed"]]
    size_groups={}
    for n in MINIMUMS:
        retained=[x for x in main if x["source_original_photo_count"]>=n]
        matches=[x for x in retained if x["nuclear_read_coverage_observed"]]
        size_groups[str(n)]={
            "n_source_photo_taxa":len(retained),"n_real_nuclear_run_sample_taxa":len(matches),
            "n_without_public_sra":len(retained)-len(matches),
            "n_ge100_loci_dp3_ge50pct":sum(x["best_single_run_241_loci_dp3_ge50pct"]>=100 for x in matches),
        }
    summary={
       "status":"ALL_ORIGINAL_244_SOURCE_TAXA_REAL_AUTHOR241_SITE_DEPTH_COVERAGE_LEDGER",
       "original_source_photo_taxa":TOTAL,
       "original_source_photo_observations":OBS,
       "WCVP_accepted_alias_collisions":COLLISIONS,
       "original_clean_accepted_taxa":len(main),
       "actual_public_read_name_candidate_taxa":len(high),
       "actual_read_supported_taxa_ge100_loci_dp3_ge50pct":sum(
          x["best_single_run_241_loci_dp3_ge50pct"]>=100 for x in high),
       "actual_read_supported_taxa_below100_loci_dp3_ge50pct":sum(
          x["best_single_run_241_loci_dp3_ge50pct"]<100 for x in high),
       "photo_depth_sensitivity":size_groups,
       "unmatched_taxa_count_not_gene_absence":True,
       "author241_read_prefix_qualifying_is_not_one_to_one_orthology":True,
       "multi_voucher_runs_not_pooled":True,
       "all244_on_independent_nuclear_tree":False,
       "no_generated_fasta_or_tree":True,
    }
    return table,summary


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--s1-full244",type=Path,required=True)
    ap.add_argument("--author241-real-taxon-depth",type=Path,required=True)
    ap.add_argument("--out-dir",type=Path,required=True)
    a=ap.parse_args()
    data,summary=build(read_csv(a.s1_full244),read_csv(a.author241_real_taxon_depth))
    a.out_dir.mkdir(parents=True,exist_ok=True)
    with (a.out_dir/"S1_ALL244_REAL241_BASE_DEPTH_COVERAGE.csv").open(
        "w",newline="",encoding="utf-8") as h:
        w=csv.DictWriter(h,fieldnames=list(data[0]),extrasaction="raise")
        w.writeheader();w.writerows(data)
    (a.out_dir/"S1_ALL244_REAL241_SITE_DEPTH_COVERAGE_SUMMARY.json").write_text(
       json.dumps(summary,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print("T2_S1_ALL244_REAL241_SITE_DEPTH_COVERAGE="+json.dumps(summary,
          ensure_ascii=False,sort_keys=True),flush=True)


if __name__=="__main__":
    main()
