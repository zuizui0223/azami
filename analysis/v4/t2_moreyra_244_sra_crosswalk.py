"""Reconcile every original S1 photographic taxon to public Moreyra2025 SRA vouchers.

The input is the exact, existing owner-authenticated 244-key census; no private
images/observer data. The raw SRA ledger is EAzami's frozen, SHA-verified 258
reconciled tree-code records. Names alone do NOT establish accepted taxon
concept identity or independently estimated 350-locus branch lengths.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import hashlib
import io
import json
from pathlib import Path
import re
import urllib.request

from analysis.v4.t2_native_s1_all_taxa_registry import (
    EA_ROOT, MOREYRA_DIR, generic_taxon_label,
)

CENSUS_N=244
PHOTO_N=170066
SRA_ROWS_N=258
EXPECTED_MATCH_N=113
S2_NAMES={"Cirsium arvense","Cirsium vulgare"}
SOURCE_MANIFEST="manifest.json"
SRA_FIELDS=("accepted_key","S1_source_taxon","original_tree_code",
            "biosample","run","experiment","voucher_and_herbarium",
            "scope_class","name_reconciliation_priority",
            "tree_code_vs_sra_name","accepted_concept_verified")
OUTPUT_FIELDS=("accepted_key","S1_source_taxon","n_S1_photo",
               "S1_source_name_collision","sra_name_candidate",
               "n_public_SRA_samples","n_unique_public_BioSamples",
               "source_sra_biosamples","source_sra_runs",
               "sra_taxon_name_disagreement_count",
               "voucher_concept_verification",
               "reconstruction_status","is_S2_focal")


def download_bytes(url):
    with urllib.request.urlopen(url,timeout=90) as handle:
        return handle.read()


def read_manifest_and_samples(fetch=download_bytes):
    manifest_raw=fetch(EA_ROOT+MOREYRA_DIR+SOURCE_MANIFEST)
    manifest=json.loads(manifest_raw.decode("utf-8"))
    if (manifest.get("rows")!=SRA_ROWS_N or
            len(manifest.get("parts",[]))!=8 or
            manifest.get("contract_version")!="moreyra2025_cirsium_reconciliation_v1"):
        raise ValueError("EAzami original Moreyra manifest changed")
    all_samples=[]
    part_digests=[]
    for part in manifest["parts"]:
        raw=fetch(EA_ROOT+MOREYRA_DIR+part["file"])
        digest=hashlib.sha256(raw).hexdigest()
        if (digest!=part["sha256"] or len(raw)!=part["bytes"]):
            raise ValueError("SHA/size mismatch in EAzami public Moreyra sample ledger")
        rows=list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
        if len(rows)!=part["rows"] or set(rows[0])!=set(manifest["columns"]):
            raise ValueError("Author SRA part changed shape")
        all_samples.extend(rows)
        part_digests.append({"name":part["file"],"sha256":digest,"n":len(rows)})
    if len(all_samples)!=SRA_ROWS_N:
        raise ValueError("Frozen Moreyra reconciliation row count changed")
    labelled=defaultdict(list)
    for row in all_samples:
        if row["sra_link_status"]!="linked_runinfo":
            continue
        stem=re.sub(r" \d+$","",row["tree_code"].strip())
        label=generic_taxon_label(stem)
        if label:
            labelled[label].append(row)
    return manifest,labelled,part_digests


def source_census(path):
    with Path(path).open(newline="",encoding="utf-8-sig") as f:
        rows=list(csv.DictReader(f))
    if (len(rows)!=CENSUS_N or
        sum(int(row["eligible_photo_observations"]) for row in rows)!=PHOTO_N or
        len({row["accepted_key"] for row in rows})!=CENSUS_N):
        raise ValueError("Source-authenticated S1 all-244 roster changed")
    required={"accepted_key","normalized_exact_source_name_candidate",
              "source_names","source_concept_collision",
              "eligible_photo_observations",
              "Moreyra2025_SRA_taxon_name_candidate"}
    if not required<=set(rows[0]):
        raise ValueError("S1 census missing original taxonomic evidence columns")
    return rows


def reconcile(census, samples):
    output=[]
    per_run=[]
    matched=0
    for row in census:
        name=row["normalized_exact_source_name_candidate"]
        collided=row["source_concept_collision"]=="True"
        if row["source_concept_collision"] not in {"True","False"}:
            raise ValueError("Original WCVP collision status missing")
        actual=[] if collided or not name else samples.get(name,[])
        flagged=row["Moreyra2025_SRA_taxon_name_candidate"]=="True"
        if row["Moreyra2025_SRA_taxon_name_candidate"] not in {"True","False"}:
            raise ValueError("Original Moreyra SRA candidate flag missing")
        if flagged!=bool(actual):
            raise ValueError("S1/EAzami Moreyra source name matching drifted")
        bios=sorted({r["biosample"] for r in actual if r["biosample"]})
        runs=sorted({r["run"] for r in actual if r["run"]})
        if flagged:
            matched+=1
        mismatches=sum(r["tree_code_vs_sra_name"]!="exact" for r in actual)
        status=("WCVP_source_name_collision" if collided else
                "SRA_source_name_candidate_voucher_concept_unverified"
                if actual else "no_linked_public_Moreyra_SRA_name")
        output.append({
            "accepted_key":row["accepted_key"],
            "S1_source_taxon":row["source_names"],
            "n_S1_photo":int(row["eligible_photo_observations"]),
            "S1_source_name_collision":collided,
            "sra_name_candidate":bool(actual),
            "n_public_SRA_samples":len(actual),
            "n_unique_public_BioSamples":len(bios),
            "source_sra_biosamples":"|".join(bios),
            "source_sra_runs":"|".join(runs),
            "sra_taxon_name_disagreement_count":mismatches,
            "voucher_concept_verification":"pending_exact_WCVP_name_rank_voucher",
            "reconstruction_status":status,
            "is_S2_focal":name in S2_NAMES,
        })
        for z in actual:
            per_run.append({
                "accepted_key":row["accepted_key"],
                "S1_source_taxon":row["source_names"],
                "original_tree_code":z["tree_code"],
                "biosample":z["biosample"],
                "run":z["run"],
                "experiment":z["experiment"],
                "voucher_and_herbarium":z["voucher_and_herbarium"],
                "scope_class":z["scope_class"],
                "name_reconciliation_priority":z["name_reconciliation_priority"],
                "tree_code_vs_sra_name":z["tree_code_vs_sra_name"],
                "accepted_concept_verified":"no",
            })
    if matched!=EXPECTED_MATCH_N:
        raise ValueError("113 source-eligible name-matched S1 taxa changed")
    s2={row["S1_source_taxon"]:row for row in output if row["is_S2_focal"]}
    if set(s2)!=S2_NAMES or any(not x["sra_name_candidate"] for x in s2.values()):
        raise ValueError("Both S2 comparison species must retain source-vouched SRA candidates")
    if len({(r["accepted_key"],r["biosample"],r["run"]) for r in per_run})!=len(per_run):
        raise ValueError("Duplicate exact S1×public SRA run linkage")
    summary={
        "status":"AUTHENTIC_PUBLIC_SRA_NAME_AND_VOUCHER_LEDGER_NOT_NUCLEAR_TREE",
        "all_S1_taxon_keys":len(output),
        "total_S1_source_photo_observations":sum(x["n_S1_photo"] for x in output),
        "S1_WCVP_source_name_collisions":sum(x["S1_source_name_collision"] for x in output),
        "S1_with_Moreyra_public_SRA_name_candidate":matched,
        "matching_public_SRA_sample_records":len(per_run),
        "matching_unique_public_SRA_runs":len({x["run"] for x in per_run}),
        "matching_unique_BioSamples":len({x["biosample"] for x in per_run}),
        "candidate_records_where_submitted_name_differs":sum(
            x["tree_code_vs_sra_name"]!="exact" for x in per_run),
        "S2_focal_both_have_SRA_sample_candidates":True,
        "reconstructed_113_tip_nuclear_tree":False,
        "original_350_locus_Moreyra_tree_recovered":False,
        "accepted_taxonomic_concept_and_same_voucher_verified":False,
        "sra_metadata_only_no_sequence_read_download":True,
        "source": "EAzami frozen 258-row public Moreyra 2025 SRA reconciliation",
    }
    return output,per_run,summary


def write_csv(path,header,rows):
    with path.open("w",newline="",encoding="utf-8") as handle:
        writer=csv.DictWriter(handle,fieldnames=header,extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--all244",type=Path,required=True)
    p.add_argument("--out-dir",type=Path,required=True)
    args=p.parse_args()
    manifest,sra,parts=read_manifest_and_samples()
    roster=source_census(args.all244)
    output,per_run,summary=reconcile(roster,sra)
    summary["EAzami_source_run"]=manifest["source_workflow_run"]
    summary["source_parts"]=parts
    summary["source_manifest_sha256"]=hashlib.sha256(
        json.dumps(manifest,sort_keys=True).encode("utf-8")).hexdigest()
    args.out_dir.mkdir(parents=True,exist_ok=True)
    write_csv(args.out_dir/"S1_244_MOREYRA_SRA_NAME_CROSSWALK.csv",
              OUTPUT_FIELDS,output)
    write_csv(args.out_dir/"S1_MOREYRA_PUBLIC_SRA_RUN_LINKS.csv",
              SRA_FIELDS,per_run)
    (args.out_dir/"S1_244_MOREYRA_SRA_CROSSWALK_SUMMARY.json").write_text(
        json.dumps(summary,indent=2,sort_keys=True,ensure_ascii=False)+"\n",
        encoding="utf-8")
    print("T2_MOREYRA_244_TO_PUBLIC_SRA_SOURCE_VOUCHERS="+json.dumps({
        k:v for k,v in summary.items()
        if k not in {"source_parts","source_manifest_sha256"}},
        sort_keys=True,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
