"""Compare SAME 111 S1 taxa and SAME public SRA runs at 30k vs 120k reads.

The actual 120k/241 per-run depth matrices are acquired from original ENA
short-read files in an independently source-validated Actions workflow.
This analysis:
  1. keeps all 125 run identities / 113 original WCVP keys;
  2. selects exactly ONE unchanged 120k public run per measured taxon;
  3. compares 30k vs 120k Q20 per-locus depth on identical public specimens;
  4. quantifies common *original reference-copy* occupancy for all 111;
  5. reports phylogenetic input feasibility without inventing aligned
     orthologous sequences or a 294-tip nuclear tree.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import json
import math
from pathlib import Path

from analysis.v4 import t2_moreyra125_120k_shared_locus_occupancy as real120
from analysis.v4 import t2_113taxa_author241_shared_locus_screen as real30

N_SOURCE_RUNS=125
N_SOURCE_KEYS=113
N_MEASURED=111
N_LOCI=241
SHORT={"SRR30887185","SRR30887248"}
LEVELS=(.25,.50)
TAXON_PCT=(.50,.65,.80,.90)


def choose_fixed_sample(source,real):
    """Same specimen at both depths; high-depth best-run optimization is QC only."""
    by_tax=defaultdict(list)
    if len(source)!=N_SOURCE_RUNS or len(real)!=N_SOURCE_RUNS:
        raise ValueError("Complete source run count changed")
    for row in source:
        run=row["run"]
        if run not in real:
            raise ValueError("Missing original public run")
        entry=real[run]
        if (row["biosample"],row["accepted_key"])!=(entry["biosample"],entry["accepted_key"]):
            raise ValueError("Original run / specimen / accepted-key identity disagreement")
        by_tax[row["accepted_key"]].append((row,entry))
    if len(by_tax)!=N_SOURCE_KEYS:
        raise ValueError("Original 113 source accepted keys changed")
    chosen=[]
    missing=[]
    for key,items in sorted(by_tax.items()):
        present=[(row,r) for row,r in items if r["status"].startswith("real_")]
        names={str(x[0]["S1_source_taxon"]).replace("_"," ").strip() for x in items}
        if len(names)!=1:
            raise ValueError("Source accepted key has competing Cirsium names")
        if not present:
            missing.append({"accepted_key":key,
                            "name":next(iter(names)),
                            "source_run_ids":"|".join(sorted(x[0]["run"] for x in items)),
                            "reason":"NO_120K_SOURCE_PREFIX_NOT_GENETIC_ABSENCE"})
            continue
        best=sorted(present,key=lambda pair:(
            -int(pair[1]["n_loci_dp3_ge50pct"]),
            -int(pair[1]["n_loci_dp3_ge25pct"]),
            pair[0]["run"]))[0]
        row,r=best
        chosen.append({
           "accepted_key":key,"name":next(iter(names)),
           "run":row["run"],"biosample":row["biosample"],
           "n_source_runs":len(items),
           "real_n_loci_Q20_120k_dp3_ge50pct":r["n_loci_dp3_ge50pct"],
           "real_n_loci_Q20_120k_dp3_ge25pct":r["n_loci_dp3_ge25pct"],
           "same_sample_for_both_depths":True,
           "taxon_concept_verified":False,
        })
    if len(chosen)!=N_MEASURED or len(missing)!=N_SOURCE_KEYS-N_MEASURED:
        raise ValueError("Original 111 measured / 2 short source taxon coverage changed")
    if {x for d in missing for x in d["source_run_ids"].split("|")}!=SHORT:
        raise ValueError("Two original short source runs changed")
    return chosen,missing


def build_locus_matrix(source, old_runs, new_runs):
    selected,missing=choose_fixed_sample(source,new_runs)
    if len(old_runs)!=N_SOURCE_RUNS or set(old_runs)!=set(new_runs):
        raise ValueError("Original same 125 run identity 30k-vs-120k source missing")
    matrix=[]
    combined=defaultdict(list)
    for sample in selected:
        run=sample["run"]
        old=old_runs[run]
        new=new_runs[run]
        if old["biosample"]!=new["biosample"] or old["accepted_key"]!=new["accepted_key"]:
            raise ValueError("30k versus 120k original biological voucher disagreement")
        if (not isinstance(old.get("loci"),list) or
            not isinstance(new.get("per_locus_241"),list) or
            len(old["loci"])!=N_LOCI or len(new["per_locus_241"])!=N_LOCI):
            raise ValueError("Both real-read depth levels require 241 named loci")
        old_by={d["locus"]:d for d in old["loci"]}
        new_by={d["locus"]:d for d in new["per_locus_241"]}
        if len(old_by)!=N_LOCI or set(old_by)!=set(new_by):
            raise ValueError("Original author241 named source-locus set disagrees")
        for name in sorted(old_by):
            a,b=old_by[name],new_by[name]
            if not all(x.get("reference_choice_maximizes_coverage_not_orthology") is True for x in (a,b)):
                raise ValueError("Unverified reference-copy selection")
            for z in (a,b):
                if not 0<=float(z["dp5_breadth"])<=float(z["dp3_breadth"])<=float(z["dp1_breadth"])<=1:
                    raise ValueError("Bad Q20 single-reference depth fraction")
            row={
                "accepted_key":sample["accepted_key"],
                "source_taxon":sample["name"],
                "public_run":run,
                "biosample":sample["biosample"],
                "author241_locus":name,
                "30k_best_reference":a["best_original_reference"],
                "120k_best_reference":b["best_original_reference"],
                "same_reference_copy_at_both_depths":
                    a["best_original_reference"]==b["best_original_reference"],
                "30k_dp3_breadth":float(a["dp3_breadth"]),
                "120k_dp3_breadth":float(b["dp3_breadth"]),
                "30k_dp3_sites":int(a["dp3_callable_sites"]),
                "120k_dp3_sites":int(b["dp3_callable_sites"]),
                "taxon_identity_verified":False,
                "Cirsium_orthology_certified":False,
            }
            matrix.append(row)
            combined[name].append(row)
    if len(matrix)!=N_MEASURED*N_LOCI or len(combined)!=N_LOCI:
        raise ValueError("Incomplete actual 111x241 locus-by-taxon matrix")
    return selected,missing,matrix,combined


def summarize(source,old,new):
    selected,missing,matrix,by_locus=build_locus_matrix(source,old,new)
    locus_summary=[]
    for name,rows in sorted(by_locus.items()):
        if len(rows)!=N_MEASURED:
            raise ValueError("SRA source locus not represented in every measured taxon")
        metrics={"locus":name}
        for r in (30000,120000):
            breadth=f"{r//1000}k_dp3_breadth"
            ref=f"{r//1000}k_best_reference"
            for threshold in LEVELS:
                eligible=[x for x in rows if x[breadth]>=threshold]
                counts=Counter(x[ref] for x in eligible)
                t=int(threshold*100)
                metrics[f"taxa_{r//1000}k_dp3_ge{t}pct"]=len(eligible)
                metrics[f"taxa_{r//1000}k_dp3_ge{t}pct_same_ref"]=max(counts.values(),default=0)
                metrics[f"modal_{r//1000}k_dp3_ge{t}pct_ref"]=(
                    sorted(counts,key=lambda x:(-counts[x],x))[0] if counts else "")
        metrics["same_ref_copy_30k_vs120k_fraction"]=sum(
           x["same_reference_copy_at_both_depths"] for x in rows)/N_MEASURED
        metrics["no_orthology_certified"]=True
        locus_summary.append(metrics)
    common={}
    for r in (30,120):
        for pct in (25,50):
            for n in TAXON_PCT:
                minimum=math.ceil(N_MEASURED*n)
                series=lambda col:sum(x[col]>=minimum for x in locus_summary)
                tag=f"{r}k_dp3_{pct}pct_atleast{round(n*100)}pct_taxa"
                common[tag]={
                    "n_loci_any_reference":series(f"taxa_{r}k_dp3_ge{pct}pct"),
                    "n_loci_one_modal_original_reference":series(
                         f"taxa_{r}k_dp3_ge{pct}pct_same_ref"),
                    "n_taxa_required":minimum,
                }
    gain={key: {
            "n_loci_any_ref_gain_120k_minus30k":
                v["n_loci_any_reference"]-common["30"+key[3:]]["n_loci_any_reference"],
            "n_loci_same_ref_gain_120k_minus30k":
                v["n_loci_one_modal_original_reference"]-
                common["30"+key[3:]]["n_loci_one_modal_original_reference"],
        } for key,v in common.items() if key.startswith("120k_")}
    focal={}
    for name in ("Cirsium conspicuum","Cirsium scopulorum","Cirsium arvense","Cirsium vulgare"):
        taxa=[row for row in selected if row["name"]==name]
        if len(taxa)>1:
            raise ValueError("More than one sample for strict photo source name")
        if taxa:
            accession=taxa[0]["run"]
            rows=[x for x in matrix if x["public_run"]==accession]
            focal[name]={
                "run":accession,"n_loci_30k_dp3_ge25pct":
                   sum(x["30k_dp3_breadth"]>=.25 for x in rows),
                "n_loci_120k_dp3_ge25pct":
                   sum(x["120k_dp3_breadth"]>=.25 for x in rows),
                "n_loci_30k_dp3_ge50pct":
                   sum(x["30k_dp3_breadth"]>=.5 for x in rows),
                "n_loci_120k_dp3_ge50pct":
                   sum(x["120k_dp3_breadth"]>=.5 for x in rows),
                "published_botanical_U_D_independent_of_S1_photo":(
                    "nodding" if name in ("Cirsium conspicuum","Cirsium scopulorum")
                    else "not_assigned_here"),
            }
    summary={
        "status":"REAL_PUBLIC_111_S1_TAXON_SAME_RUN_30K_120K_AUTHOR241_REFERENCE_COPY_OCCUPANCY",
        "source_original_S1_accepted_keys":244,
        "source_original_S1_named_units":251,
        "source_original_S1_photos":170066,
        "source_original_SRA_public_runs":N_SOURCE_RUNS,
        "source_original_SRA_name_matched_keys":N_SOURCE_KEYS,
        "measured_same_original_run_120k_taxa":N_MEASURED,
        "short_public_runs_not_genetically_absent":missing,
        "real_30k_vs120k_Q20_3x_241_loci":len(by_locus),
        "real_111x241_taxon_locus_comparisons":len(matrix),
        "occupation_same_reference_does_not_ensure_Cirsium_orthology":True,
        "across_taxa_occupancy":common,
        "increase_on_same_111_taxa_and_241_targets":gain,
        "focal_real_public_run_30k_vs120k":focal,
        "phylogenetic_tree_fitted_from_this_matrix":False,
        "de_novo_moreyra_350_recreated":False,
        "biological_orientation_caused_by_BIO12":False,
        "S1_photo_axis_not_gravitational_nodding":True,
    }
    return selected,missing,matrix,locus_summary,summary


def write_csv(path,rows):
    with Path(path).open("w",encoding="utf-8",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]),extrasaction="raise")
        writer.writeheader();writer.writerows(rows)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--original30k-shards",type=Path,required=True)
    p.add_argument("--real120k-shards",type=Path,required=True)
    p.add_argument("--out-dir",type=Path,required=True)
    a=p.parse_args()
    with a.source.open(encoding="utf-8-sig",newline="") as handle:
        source=list(csv.DictReader(handle))
    old=real30.read_perrun(a.original30k_shards)
    new=real120.reconcile(source,real120.read_shards(a.real120k_shards))
    chosen,missing,matrix,loci,summary=summarize(source,old,new)
    a.out_dir.mkdir(parents=True,exist_ok=True)
    write_csv(a.out_dir/"111_SOURCE_TAXA_BIOSAMPLE_SELECTION.csv",chosen)
    write_csv(a.out_dir/"111_TAXA_X_241_REAL_SITE_DEPTH_COPY_MATRIX.csv",matrix)
    write_csv(a.out_dir/"241_AUTHOR_LOCI_SAME_REF_OCCUPANCY_30K_120K.csv",loci)
    (a.out_dir/"111_TAXON_REAL_COMMON241_30K120K_SUMMARY.json").write_text(
        json.dumps(summary,sort_keys=True,ensure_ascii=False,indent=2)+"\n",
        encoding="utf-8")
    print("T2_ORIGINAL111_TAXA_SAME_RUN_30K120K_REFERENCE_COPY_GAIN="+json.dumps({
       "run_count":summary["source_original_SRA_public_runs"],
       "measured_taxa":summary["measured_same_original_run_120k_taxa"],
       "loci":summary["real_30k_vs120k_Q20_3x_241_loci"],
       "common":summary["across_taxa_occupancy"],
       "focal":summary["focal_real_public_run_30k_vs120k"],
       "short":summary["short_public_runs_not_genetically_absent"],
    },sort_keys=True,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
