"""Attach ACTUAL all-111 original241 reference-copy coverage to ALL S1 source names.

Read-only join of source-verified public taxon-level artifacts:
251 original (WCVP-key, source-name) photo/environment units and the
111-taxon x 241 original same-run Q20 site-depth/copy matrix. All missing
nuclear sources stay explicit NA, never zero loci or non-nodding. The 118
near-common ref-copy loci are potential alignment inputs, NOT orthologs or
an accepted global species tree.
"""
from __future__ import annotations
import argparse
import csv
import json
from collections import Counter,defaultdict
from pathlib import Path
import statistics

N_SOURCE=251
N_ACCEPTED=244
N_PHOTOS=170066
N_COLLISIONS=7
N_REAL_NAME_CANDIDATES=113
N_MEASURED=111
N_LOCI=241
N_MATCHED_CORE=118
N_MISSING_SHORT=2


def normalize(text):
    return " ".join(str(text).replace("_"," ").split())


def read_csv(path):
    with Path(path).open("r",encoding="utf-8-sig",newline="") as handle:
        rows=list(csv.DictReader(handle))
    if not rows:
        raise ValueError("Missing original census / actual Q20 nuclear matrix")
    return rows


def compile_coverage(original,chosen,matrix,source_loci):
    keys=Counter(r["accepted_key"] for r in original)
    if (len(original)!=N_SOURCE or len(keys)!=N_ACCEPTED
        or len({(r["accepted_key"],normalize(r["S1_original_source_name"])) for r in original})!=N_SOURCE
        or sum(int(r["original_S1_photo_count"]) for r in original)!=N_PHOTOS
        or sum(c==2 for c in keys.values())!=N_COLLISIONS
        or any(c>2 for c in keys.values())):
        raise ValueError("251 source-name / 244 accepted WCVP / 170066 native observations drift")
    sample={r["accepted_key"]:r for r in chosen}
    if len(chosen)!=N_MEASURED or len(sample)!=N_MEASURED:
        raise ValueError("111 original source Q20 public BioSamples not unique/complete")
    if len(source_loci)!=N_LOCI or len({x["locus"] for x in source_loci})!=N_LOCI:
        raise ValueError("Frozen original author241 locus copy panel missing")
    minimum=89  # ceil(0.8*111) source 120k/Q20 copy breadth gate
    if N_MEASURED!=111:
        raise ValueError("80-percent common copy criterion refers to original 111 taxa")
    qualifying={}
    for r in source_loci:
        if int(r["taxa_120k_dp3_ge25pct_same_ref"])>=minimum:
            name=r["modal_120k_dp3_ge25pct_ref"]
            if not name or not name.endswith("-"+r["locus"]):
                raise ValueError("Invalid single original reference label in 111-common locus panel")
            qualifying[r["locus"]]=name
    if len(qualifying)!=N_MATCHED_CORE:
        raise ValueError("Real observed same-copy 111-taxon 118 source loci changed")

    run_index={}
    for r in matrix:
        key=r["accepted_key"]
        locus=r["author241_locus"]
        if key not in sample:
            raise ValueError("Depth matrix contains public run not in 111 measured taxon panel")
        if r["public_run"]!=sample[key]["run"] or r["biosample"]!=sample[key]["biosample"]:
            raise ValueError("Not the selected same verified SRA run/BioSample")
        tag=(key,locus)
        if tag in run_index:
            raise ValueError("Duplicate taxon-locus high-quality row")
        run_index[tag]=r
    if len(run_index)!=N_MEASURED*N_LOCI:
        raise ValueError("Incomplete source 111x241 nuclear site depth panel")
    taxon_occupancy={}
    for key in sample:
        retained=[run_index[(key,locus)] for locus in sorted(qualifying)]
        c25=sum(
            float(r["120k_dp3_breadth"])>=.25 and
            r["120k_best_reference"]==qualifying[r["author241_locus"]]
            for r in retained)
        c50=sum(
            float(r["120k_dp3_breadth"])>=.50 and
            r["120k_best_reference"]==qualifying[r["author241_locus"]]
            for r in retained)
        taxon_occupancy[key]=(c25,c50)
        if not (0<=c50<=c25<=N_MATCHED_CORE):
            raise ValueError("Source depth gate count exceeds 118 original loci")
    outputs=[]
    for original_row in original:
        row=original_row.copy()
        key=row["accepted_key"]
        measured=sample.get(key)
        old=row["Q20_120k_read_state"]
        if not old:
            raise ValueError("Missing original S1 source name 120k state")
        n_same=int(row["source_same_biosample_run_120k_Q20_successes"])
        if n_same>0 and measured is None:
            raise ValueError("Original 251 source has measured 120k but no authentic 111 sample")
        if measured is not None:
            if keys[key]!=1:
                raise ValueError("One source Moreyra run assigned to pooled WCVP names")
            if normalize(measured["name"])!=normalize(row["S1_original_source_name"]):
                raise ValueError("Cirsium source name differs from measured 120k same-run tip")
            if n_same<1:
                raise ValueError("Conflicting observed source 120k status")
            row["selected_real120k_SRA_run"]=measured["run"]
            row["selected_real120k_BioSample"]=measured["biosample"]
            row["shared_original_reference_118_loci_dp3_ge25pct"]=taxon_occupancy[key][0]
            row["shared_original_reference_118_loci_dp3_ge50pct"]=taxon_occupancy[key][1]
            row["nuclear120k_source_observation_status"]="measured_real_Q20_source_read_prefix"
        else:
            row["selected_real120k_SRA_run"]=""
            row["selected_real120k_BioSample"]=""
            row["shared_original_reference_118_loci_dp3_ge25pct"]=None
            row["shared_original_reference_118_loci_dp3_ge50pct"]=None
            if old=="no120k_due_to_short_public_fastq":
                row["nuclear120k_source_observation_status"]="public_run_shorter_than_120k_not_gene_loss"
            else:
                row["nuclear120k_source_observation_status"]="no_source_verified_120k_Moreyra_sample"
        row["original_same_copy_reference_not_certified_ortholog"]=True
        row["origin_of_botanical_nodding_not_inferred_from_EXIF_image"]=True
        row["all294_independent_nuclear_tree_not_yet_inferred"]=True
        outputs.append(row)
    states=Counter(r["nuclear120k_source_observation_status"] for r in outputs)
    if states!={
        "measured_real_Q20_source_read_prefix":N_MEASURED,
        "public_run_shorter_than_120k_not_gene_loss":N_MISSING_SHORT,
        "no_source_verified_120k_Moreyra_sample":N_SOURCE-N_REAL_NAME_CANDIDATES
    }:
        raise ValueError("Original full 251 source identity or real Q20 completeness drift")

    focal={}
    for name in ("Cirsium arvense","Cirsium vulgare",
                 "Cirsium conspicuum","Cirsium scopulorum",
                 "Cirsium brevicaule","Cirsium irumtiense"):
        matches=[r for r in outputs if normalize(r["S1_original_source_name"])==name]
        if len(matches)!=1:
            raise ValueError("Original focal Cirsium source name lost or duplicated")
        rec=matches[0]
        focal[name]={
            "original_S1_photo_count":int(rec["original_S1_photo_count"]),
            "real120k_nuclear_source_status":rec["nuclear120k_source_observation_status"],
            "n_original118_same_ref_loci_dp3_ge25pct":
                rec["shared_original_reference_118_loci_dp3_ge25pct"],
        }
    q=sorted(x[0] for x in taxon_occupancy.values())
    summary={
        "status":"ALL_251_ORIGINAL_S1_NAMES_WITH_REAL111_Q20_118_SHARED_REF_LOCUS_READINESS",
        "original_photo_observations":N_PHOTOS,
        "original_WCVP_accepted_keys":N_ACCEPTED,
        "original_source_name_units":N_SOURCE,
        "original_collided_accepted_keys":N_COLLISIONS,
        "original_SRA_source_name_candidates":N_REAL_NAME_CANDIDATES,
        "S1_units_with_real_Q20_120k_NA_source":states["measured_real_Q20_source_read_prefix"],
        "S1_units_short_original_SRA_archives":states["public_run_shorter_than_120k_not_gene_loss"],
        "S1_units_without_original_Moreyra125_source":states["no_source_verified_120k_Moreyra_sample"],
        "original_author_no_warning_candidates":N_LOCI,
        "core_reference_copy_Q20_shared_atleast80pct_111_taxa":N_MATCHED_CORE,
        "taxon_Q20_dp3_ge25pct_loci_within_core118":{
           "min":min(q),"median":statistics.median(q),"max":max(q),
           "n_taxa_ge100":sum(x>=100 for x in q),
           "n_taxa_ge80":sum(x>=80 for x in q),
        },
        "source_focus":focal,
        "n_original_source_photo_units_dropped":0,
        "no_missing_sequences_mislabeled_as_zero":True,
        "reference_compatibility_not_ortholog_certificate":True,
        "species_phylogeny_constructed":False,
        "causal_rainfall_evolution_estimated":False,
    }
    return outputs,summary


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source251",type=Path,required=True)
    ap.add_argument("--source111",type=Path,required=True)
    ap.add_argument("--matrix111",type=Path,required=True)
    ap.add_argument("--loci241",type=Path,required=True)
    ap.add_argument("--out-dir",type=Path,required=True)
    a=ap.parse_args()
    result,receipt=compile_coverage(
        read_csv(a.source251),read_csv(a.source111),
        read_csv(a.matrix111),read_csv(a.loci241))
    a.out_dir.mkdir(parents=True,exist_ok=True)
    with (a.out_dir/"FULL251_ORIGINAL_CIRSIUM_REAL120K_118_REF_COPY_READINESS.csv").open(
        "w",encoding="utf-8",newline="") as out:
        writer=csv.DictWriter(out,fieldnames=list(result[0]),extrasaction="raise")
        writer.writeheader();writer.writerows(result)
    (a.out_dir/"FULL251_ORIGINAL_Q20_118_SHARED_REF_READINESS_SUMMARY.json").write_text(
        json.dumps(receipt,sort_keys=True,ensure_ascii=False,indent=2)+"\n",
        encoding="utf-8")
    print("T2_REAL111_CORE118_COVERAGE_JOINED_TO_FULL251="+json.dumps(
        receipt,sort_keys=True,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
