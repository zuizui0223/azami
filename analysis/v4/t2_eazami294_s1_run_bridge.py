"""Link every original S1 taxon to EAzami's already audited 294-sample panel.

Both records carry PUBLIC source accessions. This prevents silently downloading
S1 125 Moreyra public runs as a second duplicate panel instead of reusing
EAzami's 256 Moreyra plus 38 Chang biological samples.

Sample/run identity, published tree_code and taxon-name candidates are separate.
The 294 panel is NOT a completed 294-tip phylogeny or admitted S1 accepted
concept phylogenetic result.
"""
from __future__ import annotations

import argparse
from collections import defaultdict, Counter
import csv
import hashlib
import json
from pathlib import Path

EA_PANEL_FIELDS = {"biosample","run_accessions","source_studies",
                   "source_taxon_label","analysis_taxon_label","panel_id",
                   "assay","name_review_required"}
S1_FIELDS = {"accepted_key","source_names","source_concept_collision",
             "normalized_exact_source_name_candidate",
             "eligible_photo_observations"}
SRA_FIELDS = {"accepted_key","S1_source_taxon","original_tree_code",
              "biosample","run","voucher_and_herbarium"}
EXPECTED = {"s1_taxa":244,"s1_photos":170066,
            "eazami_samples":294,"eazami_runs":295,
            "eazami_moreyra":256,"eazami_chang":38,
            "s1_sra_rows":125,"s1_sra_taxa":113}

def read_csv(path, columns):
    with Path(path).open(newline="",encoding="utf-8-sig") as handle:
        rows=list(csv.DictReader(handle))
    if not rows or not columns.issubset(rows[0]):
        raise ValueError(f"{path}: missing fields or empty table")
    return rows

def normalize_name(name):
    return " ".join(name.replace("_"," ").split()).strip()

def parse_eazami_panel(rows):
    by_run={}
    by_bio={}
    by_name=defaultdict(list)
    counts=Counter()
    for row in rows:
        bio=row["biosample"].strip()
        runs=row["run_accessions"].split("|")
        if not bio or not runs or any(not acc for acc in runs):
            raise ValueError("EAzami sample/run inventory contains empty IDs")
        if bio in by_bio:
            raise ValueError("EAzami duplicates biological BioSample")
        if int(row.get("run_count",len(runs)))!=len(runs):
            raise ValueError("EAzami multiple-run sample run-count mismatch")
        by_bio[bio]=row
        for run in runs:
            if run in by_run:
                raise ValueError("EAzami duplicates source run")
            by_run[run]=(bio,row)
        counts["moreyra" if row["source_studies"]=="Moreyra2025" else "chang"]+=1
        name=normalize_name(row["analysis_taxon_label"])
        if name:
            by_name[name].append(row)
    if (len(rows)!=EXPECTED["eazami_samples"] or
        len(by_run)!=EXPECTED["eazami_runs"] or
        counts!={"moreyra":256,"chang":38}):
        raise ValueError("EAzami frozen 294/295, 256+38 inventory changed")
    return by_run,by_bio,by_name

def link(source_taxa, sra, panel):
    if (len(source_taxa)!=EXPECTED["s1_taxa"] or
        sum(int(x["eligible_photo_observations"]) for x in source_taxa)!=EXPECTED["s1_photos"] or
        len({x["accepted_key"] for x in source_taxa})!=len(source_taxa)):
        raise ValueError("Original S1 full 244/170066 inventory not authenticated")
    if (len(sra)!=EXPECTED["s1_sra_rows"] or
        len({x["accepted_key"] for x in sra})!=EXPECTED["s1_sra_taxa"]):
        raise ValueError("Original 125-run / 113-source candidate SRA table changed")
    by_run,by_bio,by_name=parse_eazami_panel(panel)
    # Keep sample-level and run-level match separate: impossible mappings fail.
    per_run=[]
    index=defaultdict(list)
    for x in sra:
        acc=x["run"]
        hit=by_run.get(acc)
        status="absent_from_EAzami_294"
        panel_id=""
        assay=""
        panel_bio=""
        sample_taxon=""
        if hit:
            panel_bio,row=hit
            if panel_bio!=x["biosample"]:
                raise ValueError("S1 Moreyra / EAzami 294 same run different BioSample")
            if row["source_studies"]!="Moreyra2025":
                raise ValueError("Original Moreyra source run mapped to Chang-only sample")
            if row["source_taxon_label"]!=x["original_tree_code"]:
                raise ValueError("EAzami Moreyra original tree_code differs from source run")
            status="exact_run_BioSample_tree_code"
            panel_id=row["panel_id"]
            assay=row["assay"]
            sample_taxon=row["analysis_taxon_label"]
        indexed={
            "accepted_key":x["accepted_key"],
            "S1_source_taxon":x["S1_source_taxon"],
            "run":acc,"biosample":x["biosample"],
            "EAzami_294_panel_id":panel_id,
            "EAzami_294_assay":assay,
            "EAzami_294_taxon_label":sample_taxon,
            "sample_run_identity_status":status,
            "accepted_concept_voucher_identity_verified":False,
            "inferred_global_species_tree_tip":False,
        }
        per_run.append(indexed)
        index[x["accepted_key"]].append(indexed)
    if len({x["run"] for x in per_run})!=125:
        raise ValueError("S1 Moreyra original source run unexpectedly repeated")
    taxon_rows=[]
    for row in source_taxa:
        key=row["accepted_key"]
        original=normalize_name(row["normalized_exact_source_name_candidate"])
        conflicted=row["source_concept_collision"].lower()=="true"
        if row["source_concept_collision"].lower() not in {"true","false"}:
            raise ValueError("Original S1 collision field must be explicit")
        linked=index.get(key,[])
        hits=[x for x in linked if x["sample_run_identity_status"]=="exact_run_BioSample_tree_code"]
        namehits=[] if conflicted or not original else by_name.get(original,[])
        moreyra_namehits=[x for x in namehits if x["source_studies"]=="Moreyra2025"]
        chang_namehits=[x for x in namehits if x["source_studies"]!="Moreyra2025"]
        if conflicted and linked:
            raise ValueError("Source-collided WCVP key not allowed inside strict 113 name-match cohort")
        taxon_rows.append({
            "accepted_key":key,
            "source_names":row["source_names"],
            "n_original_photo":int(row["eligible_photo_observations"]),
            "WCVP_source_alias_collision":conflicted,
            "n_Moreyra_SRA_source_runs":len(linked),
            "n_exact_run_BioSample_linked_to_EAzami294":len(hits),
            "n_original_SRA_not_in_EAzami294":len(linked)-len(hits),
            "n_EAzami294_Moreyra_label_only_samples":len(moreyra_namehits),
            "n_EAzami294_Chang_label_only_samples":len(chang_namehits),
            "label_only_EAzami294_sample_IDs":"|".join(sorted(
                {x["panel_id"] for x in namehits})),
            "original_SRA_run_BioSample_verified":len(linked)>0 and len(linked)==len(hits),
            "accepted_taxonomic_concept_reviewed":False,
            "completed_global_nuclear_phylogeny":False,
        })
    stats=Counter(x["sample_run_identity_status"] for x in per_run)
    status=Counter()
    for r in taxon_rows:
        if r["n_Moreyra_SRA_source_runs"]>0:
            status["Moreyra_SRA_candidate_taxa"]+=1
        if r["n_exact_run_BioSample_linked_to_EAzami294"]>0:
            status["confirmed_EAzami294_same_run_taxa"]+=1
        if r["n_EAzami294_Chang_label_only_samples"]>0:
            status["Chang_label_only_taxa"]+=1
        if r["n_EAzami294_Moreyra_label_only_samples"]>0:
            status["Moreyra_label_only_taxa"]+=1
    return taxon_rows,per_run,{
        "schema":"azami_T2_S1_244_EAzami_294_exact_sample_run_bridge_v1",
        "status":"SOURCE_PANEL_OVERLAP_AUDIT_NO_294_TIP_TREE",
        "S1_original_photo_taxa":len(source_taxa),
        "S1_original_photo_observations":sum(r["n_original_photo"] for r in taxon_rows),
        "EAzami_original_biological_samples":len(panel),
        "EAzami_original_run_accessions":len(by_run),
        "EAzami_original_Moreyra_samples":256,
        "EAzami_original_Chang_samples":38,
        "S1_source_candidate_Moreyra_run_rows":len(per_run),
        "source_run_identity_status":dict(stats),
        "taxon_level_status":dict(status),
        "S1_and_EAzami_sample_runs_reacquired_or_resequenced":False,
        "EAzami_294_global_nuclear_tree_completed":False,
        "Moreyra_published_original_final350_tree_recovered":False,
        "accepted_taxonomic_concept_and_voucher_verified":False,
        "no_auto_combining_Moreyra_Chang_assays":True,
        "requires_original_HPC_294_tree_input_and_orthology_validation":True,
        "main_next_action":"Use EAzami frozen 294-tip Compositae1061 HPC assembly once; join verified S1 taxon concepts onto resulting independently inferred branch-length phylogeny.",
    }

def output_csv(path,rows):
    with path.open("w",newline="",encoding="utf-8") as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]),extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--s1-census",type=Path,required=True)
    parser.add_argument("--s1-sra",type=Path,required=True)
    parser.add_argument("--eazami-panel",type=Path,required=True)
    parser.add_argument("--out-dir",type=Path,required=True)
    a=parser.parse_args()
    f=[a.s1_census,a.s1_sra,a.eazami_panel]
    provenance={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in f}
    source=read_csv(a.s1_census,S1_FIELDS)
    sra=read_csv(a.s1_sra,SRA_FIELDS)
    panel=read_csv(a.eazami_panel,EA_PANEL_FIELDS)
    all_taxa,all_runs,summary=link(source,sra,panel)
    summary["source_file_sha256"]=provenance
    a.out_dir.mkdir(parents=True,exist_ok=True)
    output_csv(a.out_dir/"T2_ALL244_TO_EAzami294_TAXON_COVERAGE.csv",all_taxa)
    output_csv(a.out_dir/"T2_ALL125_MOREYRA_TO_EAzami294_EXACT_RUNS.csv",all_runs)
    (a.out_dir/"T2_EAzami294_S1_244_BRIDGE_SUMMARY.json").write_text(
        json.dumps(summary,indent=2,sort_keys=True,ensure_ascii=False)+"\n",
        encoding="utf-8")
    print("T2_EAzami294_S1_244_REAL_SAMPLE_OVERLAP="+json.dumps(summary,
                                               sort_keys=True,ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()
