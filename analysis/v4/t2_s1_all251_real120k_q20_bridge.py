"""Integrate the REAL 125-run nested 30k/120k Q20 depth into ALL 251 S1 units.

All 170066 original S1 photo measurements (244 accepted keys, 251 original-name
units) survive. Each 120k sample is the SAME run/BioSample compared to its
own 30k observation. A shorter archival FASTQ is a technical run-size ceiling,
not zero alleles/genes. Original EAzami294 tips are sample sources, not an
inferred species tree. Post-hoc 100-locus/50%-breadth is a sequencing gate,
NOT proven single-copy ortholog recovery or botanical posture evolution.
"""
from __future__ import annotations

import argparse
from collections import Counter,defaultdict
import csv
import json
from pathlib import Path

N_UNITS=251
N_ACCEPTED=244
N_COLLISION_KEYS=7
N_IMAGES=170066
N_NUCLEAR_KEYS=113
N_RUNS=125
N_120K_RUNS=123
N_SHORT=2
N_120K_KEYS=111
MINIMUM_DP3_HALF=100
S2={"Cirsium_arvense","Cirsium_vulgare"}
TROPICAL={"Cirsium_brevicaule","Cirsium_irumtiense"}
SHARED_WCVP="206745753"

def norm(x):
    return "_".join(str(x or "").strip().split()).replace(" ","_")

def load_csv(path):
    with Path(path).open(newline="",encoding="utf-8-sig") as h:
        rows=list(csv.DictReader(h))
    if not rows:
        raise ValueError("Empty original S1/EAzami/real nested nuclear source file")
    return rows

def present_bool(value,field):
    x=str(value).lower()
    if x not in {"true","false"}:
        raise ValueError("Unspecified source flag "+field)
    return x=="true"

def validate_run_comparisons(runs):
    if len(runs)!=N_RUNS or len({x["run"] for x in runs})!=N_RUNS or \
       len({x["accepted_key"] for x in runs})!=N_NUCLEAR_KEYS:
        raise ValueError("All 125 original source SRA run/113 matched taxon keys required")
    grouped=defaultdict(list)
    ok=short=0
    for r in runs:
        status=r["full_120k_comparison_status"]
        success=present_bool(r["source_30k_to_120k_valid"],"same source run success")
        if success!=(status=="real_ENA_120k_quality20_depth"):
            raise ValueError("Nested exact-run 120k success status is inconsistent")
        before=int(r["first30k_n_loci_dp3_ge50pct"])
        if before<0 or before>241 or int(r["original_first30k_read_pairs"])>30000:
            raise ValueError("Original 30k source quality evidence invalid")
        if success:
            now=int(r["first120k_n_loci_dp3_ge50pct"])
            if not before<=now<=241:
                raise ValueError("Nested mapped Q20 depth cannot decrease within a run")
            ok+=1
        else:
            if r["first120k_n_loci_dp3_ge50pct"] not in {"","None"}:
                raise ValueError("Unavailable full 120k run must have NULL site depth")
            if r["source_run_unavailability_reason"]!="ENA source full FASTQ has fewer than 120k paired reads":
                raise ValueError("Unexplained source-run 120k coverage failure")
            short+=1
        grouped[r["accepted_key"]].append(r)
    if ok!=N_120K_RUNS or short!=N_SHORT:
        raise ValueError("123 genuine depth-120k and 2 short-archive source contract drift")
    return grouped

def merge(units,depth113,runs,summary):
    if (len(units)!=N_UNITS or
        len({(r["accepted_key"],r["S1_original_source_name"]) for r in units})!=N_UNITS or
        len({r["accepted_key"] for r in units})!=N_ACCEPTED or
        sum(int(r["original_S1_photo_count"]) for r in units)!=N_IMAGES):
        raise ValueError("Original 170066-photo/244-key/251-source-name cohort changed")
    if (len(depth113)!=N_NUCLEAR_KEYS or
        len({x["accepted_key"] for x in depth113})!=N_NUCLEAR_KEYS):
        raise ValueError("One or more original 113 source-matched nuclear taxon rows lost")
    bykey={x["accepted_key"]:x for x in depth113}
    runby=validate_run_comparisons(runs)
    if set(bykey)!=set(runby):
        raise ValueError("113 source-matched taxa do not match actual same-run evidence")
    if (summary.get("n_runs_with_recovered_120k_read_depth")!=N_120K_RUNS or
        summary.get("n_runs_without_complete_120k_prefix")!=N_SHORT or
        summary.get("n_source_taxa_with_any_120k_test")!=N_120K_KEYS):
        raise ValueError("Original all125 Q20/nested 120k source summary changed")
    if sum(len(v) for v in runby.values())!=N_RUNS:
        raise ValueError("Original run mapping incompletely grouped")
    overlaps=Counter(r["accepted_key"] for r in units)
    if sum(n==2 for n in overlaps.values())!=N_COLLISION_KEYS or \
       any(n>2 for n in overlaps.values()):
        raise ValueError("Seven original WCVP name collisions not conserved")
    outputs=[]
    for row in units:
        key=row["accepted_key"]
        full=bykey.get(key)
        ownruns=runby.get(key,[])
        present=present_bool(row["real_241_site_depth_available"],"original real source depth")
        if present!=(full is not None):
            raise ValueError("S1 251 original source cohort and real241 identity drift")
        original_name=norm(row["S1_original_source_name"])
        original_depth=row["Q20_dp3_ge50pct_loci_best_single_run"]
        if full:
            if overlaps[key]!=1 or norm(full["source_name"])!=original_name:
                raise ValueError("Cannot match nuclear run onto collided or different S1 name")
            old=int(full["best_30k_single_run_dp3_ge50pct"])
            if old!=int(original_depth):
                raise ValueError("Baseline 30k real site-depth did not match prior 251 bridge")
            real_name={norm(r["source_taxon"]) for r in ownruns}
            if real_name!={original_name}:
                raise ValueError("Source run-to-photo name mismatch")
            usable=[r for r in ownruns
                    if present_bool(r["source_30k_to_120k_valid"],"120k coverage")]
            if len(usable)!=int(full["successfully_depth120k_mapped_runs"]):
                raise ValueError("Taxon 120k successes and source-run evidence differ")
            observed=max((int(r["first120k_n_loci_dp3_ge50pct"]) for r in usable),
                         default=None)
            declared=full["best_120k_single_run_dp3_ge50pct"]
            if (None if declared in ("","None") else int(declared))!=observed:
                raise ValueError("Taxon best 120k-depth metric contradicts exact same-run results")
            # Do not compare best 30k from an unavailable 120k accession to
            # best 120k of another accession as if a nested sample changed.
            paired_gain=[int(r["first120k_n_loci_dp3_ge50pct"])-
                         int(r["first30k_n_loci_dp3_ge50pct"]) for r in usable]
            status=("measured_public_prefix120k_Q20_unverified_orthologs" if usable
                    else "no120k_due_to_short_public_fastq")
            n_run=len(ownruns)
        else:
            observed=None
            paired_gain=[]
            status=("no_Moreyra125_SRA_in_source_name_candidate_panel")
            n_run=0
        record={
          **row,
          "source_125_public_sra_run_count":n_run,
          "source_same_biosample_run_120k_Q20_successes":len(paired_gain),
          "source_same_biosample_run_120k_Q20_short_archive":n_run-len(paired_gain),
          "Q20_dp3_ge50pct_loci_best_120k_single_source_run":observed,
          "original_30k_best_single_run_loci_Q20":int(original_depth) if full else None,
          "min_paired_run_30k_to_120k_gain_loci":min(paired_gain) if paired_gain else None,
          "max_paired_run_30k_to_120k_gain_loci":max(paired_gain) if paired_gain else None,
          "Q20_120k_max_ge100_loci":bool(observed is not None and observed>=MINIMUM_DP3_HALF),
          "Q20_120k_read_state":status,
          "no_full_241_gene_consensus_or_validated_ortholog_alignment":True,
          "no_nuclear_species_tree_estimated":True,
        }
        outputs.append(record)
    if sum(x["source_same_biosample_run_120k_Q20_successes"] for x in outputs)!=N_120K_RUNS:
        raise ValueError("Source 123 120k run outcomes lost in original251 bridge")
    if sum(x["Q20_120k_max_ge100_loci"] for x in outputs)!=summary.get(
            "n_source_taxa_with_best_single_120k_dp3_ge50pct_ge100"):
        raise ValueError("100 covered loci candidate count drift from real nested source")
    pair=[r for r in outputs if r["accepted_key"]==SHARED_WCVP]
    if len(pair)!=2 or {x["S1_original_source_name"] for x in pair}!=TROPICAL:
        raise ValueError("C. brevicaule/C. irumtiense still merged by WCVP")
    if any(x["Q20_120k_max_ge100_loci"] or x["real_241_site_depth_available"]=="True"
           for x in pair):
        raise ValueError("Do not impute Moreyra125 SRA loci to separately sampled Chang taxa")
    meets=[x for x in outputs if x["Q20_120k_max_ge100_loci"]]
    distribution={}
    for threshold in (1,5,10):
        good=[x for x in outputs if int(x["original_S1_photo_count"])>=threshold]
        distribution[str(threshold)]={
          "n_S1_original_source_name_units":len(good),
          "n_with_125_Moreyra_SRA_source_evidence":sum(
              present_bool(x["real_241_site_depth_available"],"depth") for x in good),
          "n_source_names_with_real_120k_depth":sum(x["source_same_biosample_run_120k_Q20_successes"]>0 for x in good),
          "n_with_at_least100_loci_dp3_ge50pct":sum(x["Q20_120k_max_ge100_loci"] for x in good),
        }
    result={
       "schema":"azami_v4_S1_all251_real_Moreyra125_source_Q20_30k120k_integrated_v1",
       "S1_photo_observations":N_IMAGES,
       "S1_WCVP_accepted_keys":N_ACCEPTED,
       "S1_original_source_name_units":N_UNITS,
       "S1_original_taxon_key_collisions":N_COLLISION_KEYS,
       "EAzami294_is_sample_panel_not_inferred_global_nuclear_tree":True,
       "source_Moreyra125_real_baseline30k_taxon_keys":N_NUCLEAR_KEYS,
       "source_complete_exact_run_depth120k":N_120K_RUNS,
       "short_public_SRA_read_archives":N_SHORT,
       "S1_source_taxa_with_any_120k":sum(x["source_same_biosample_run_120k_Q20_successes"]>0 for x in outputs),
       "S1_source_taxa_with_at_least100_loci_dp3_half_length":len(meets),
       "qualified_source_units":[{
           "accepted_key":x["accepted_key"],
           "original_name":x["S1_original_source_name"],
           "original_S1_photo_count":int(x["original_S1_photo_count"]),
           "best_single_run_120k_Q20_loci":x["Q20_dp3_ge50pct_loci_best_120k_single_source_run"],
           "source_BIO12_mm":float(x["S1_CHELSA_BIO12_mm"]),
           "original_WCVP_collapsed_pair":present_bool(x["original_WCVP_collapsed_pair"],"WCVP conflict"),
       } for x in meets],
       "original_photo_depth_sensitivity":distribution,
       "unmeasured_138_source_name_units_are_not_zero_nuclear_genes":True,
       "same_run_before_after_compared_on_identical_biosample":True,
       "observed_120k_not_single_copy_validated":True,
       "homolog_consensus_sequences_created":False,
       "biological_U_D_orientation_not_inferred_from_EXIF":True,
       "directional_Mk_BM_OU_on_global_true_nuclear_tree_not_executed":True,
       "claims_boundary":"Read depth at one reference copy is not orthology or an accepted nuclear tip; source-name matching not voucher/accepted-species verification; no recurrent rainfall-driven nodding adaptation demonstrated",
    }
    return outputs,result

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source251",type=Path,required=True)
    p.add_argument("--taxon120k",type=Path,required=True)
    p.add_argument("--run125",type=Path,required=True)
    p.add_argument("--summary120k",type=Path,required=True)
    p.add_argument("--out-dir",type=Path,required=True)
    args=p.parse_args()
    outputs,summary=merge(load_csv(args.source251),load_csv(args.taxon120k),
                          load_csv(args.run125),
                          json.loads(args.summary120k.read_text(encoding="utf-8")))
    args.out_dir.mkdir(parents=True,exist_ok=True)
    with (args.out_dir/"S1_ALL251_ORIGINAL_NAMES_REAL_120K_Q20_SOURCE_COVERAGE.csv").open(
         "w",newline="",encoding="utf-8") as h:
        w=csv.DictWriter(h,fieldnames=list(outputs[0]),extrasaction="raise")
        w.writeheader();w.writerows(outputs)
    (args.out_dir/"S1_ALL251_REAL_120K_Q20_NUCLEAR_SOURCE_RESULT.json").write_text(
        json.dumps(summary,indent=2,sort_keys=True,ensure_ascii=False,allow_nan=False)+"\n",
        encoding="utf-8")
    print("T2_FULL251_REAL120K_Q20_NUCLEAR_ALIGNMENT_READINESS="+json.dumps({
        "n_source_name_units":summary["S1_original_source_name_units"],
        "n_public_SRA_name_matched":summary["source_Moreyra125_real_baseline30k_taxon_keys"],
        "n_120k_source_units":summary["S1_source_taxa_with_any_120k"],
        "n_at_least100_dp3_half":summary["S1_source_taxa_with_at_least100_loci_dp3_half_length"],
        "qualifying_taxa":summary["qualified_source_units"],
    },ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":
    main()
