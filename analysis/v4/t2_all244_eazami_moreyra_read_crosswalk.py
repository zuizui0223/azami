"""Exact-name S1 all244 vs EAzami pinned Moreyra-2025 258 public nuclear read ledger.

Never invent a 248-tip Newick or treat raw SRA read access as validated
independent nuclear-phylogeny placement. Do not collapse varieties into a
binomial. Only terminal integer replicate codes are normalized.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import urllib.request
from pathlib import Path

import pandas as pd

EA_REPO="zuizui0223/EAzami"
EA_COMMIT="c9870e91604ed1514344ddc8f765eae5237dbf27"
EA_DIR="data/evidence/moreyra2025_cirsium_reconciliation_v1"
S2_NAMES={"Cirsium arvense","Cirsium vulgare"}
EXPECTED_ORIGINAL_244=244
EXPECTED_RECONCILIATION_ROWS=258


def read_pinned_eazami_bytes(name):
    url=f"https://raw.githubusercontent.com/{EA_REPO}/{EA_COMMIT}/{EA_DIR}/{name}"
    with urllib.request.urlopen(url,timeout=90) as r:
        return r.read()


def load_moreyra():
    manifest_raw=read_pinned_eazami_bytes("manifest.json")
    manifest=json.loads(manifest_raw)
    if (manifest["contract_version"]!="moreyra2025_cirsium_reconciliation_v1"
            or manifest["rows"]!=EXPECTED_RECONCILIATION_ROWS
            or len(manifest["parts"])!=8):
        raise ValueError("EAzami Moreyra public SRA reconciliation manifest changed")
    combined=[]
    part_digests=[]
    for meta in manifest["parts"]:
        raw=read_pinned_eazami_bytes(meta["file"])
        if (len(raw)!=meta["bytes"]
            or hashlib.sha256(raw).hexdigest()!=meta["sha256"]):
            raise ValueError("EAzami Moreyra original-SRA part SHA mismatch")
        rows=list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
        if len(rows)!=meta["rows"] or set(rows[0])!=set(manifest["columns"]):
            raise ValueError("EAzami Moreyra public sample schema drift")
        combined.extend(rows)
        part_digests.append(meta["sha256"])
    if len(combined)!=EXPECTED_RECONCILIATION_ROWS:
        raise ValueError("258 public SRA-linked Cirsium sample scope drift")
    if len(set(row["run"] for row in combined))!=len(combined):
        raise ValueError("Duplicate public raw run ID in Moreyra sample ledger")
    if any(row["sra_link_status"]!="linked_runinfo" for row in combined):
        raise ValueError("Moreyra sample row without an original SRA link")
    return manifest,combined,hashlib.sha256(manifest_raw).hexdigest()


def name_representation(raw):
    """Allow a final sample number; never equate species with varieties."""
    raw=re.sub(r"\s+"," ",str(raw).strip().replace("_"," "))
    stripped=re.sub(r"\s+[0-9]+$","",raw)
    if not stripped.startswith("Cirsium "):
        return None,None
    return raw,stripped


def analyse(coverage:pd.DataFrame,moreyra):
    needed={"accepted_key","source_taxon_names","n_source_orientation_photos",
            "tree194_covered","strict_phylogenetic_eligible","s2_focal"}
    if not needed<=set(coverage) or len(coverage)!=EXPECTED_ORIGINAL_244:
        raise ValueError("S1 original 244 taxon roster missing or incomplete")
    if (coverage.accepted_key.nunique()!=244
            or int(coverage.n_source_orientation_photos.sum())!=170066
            or int(coverage.tree194_covered.sum())!=194
            or int(coverage.strict_phylogenetic_eligible.sum())!=191):
        raise ValueError("S1 original source cohort checksum/coverage contract changed")
    grouped={}
    for r in moreyra:
        raw,standard=name_representation(r["tree_code"])
        if raw is None:
            continue
        grouped.setdefault(standard,[]).append(r)
    out=[]
    for row in coverage.to_dict("records"):
        names=[name_representation(x)[0] for x in row["source_taxon_names"].split("|")]
        if any(not name for name in names):
            raise ValueError("Unexpected non-Cirsium S1 source taxon name")
        collision=bool(row["source_name_collision"])
        if collision and len(names)<2:
            raise ValueError("Source-name collision flag without multiple original names")
        matching=grouped.get(names[0],[]) if len(names)==1 and not collision else []
        # One exact species name to a public run is evidence for raw sequence
        # availability only; not congruence of vouchers or taxonomy treatment.
        exact=[r for r in matching if name_representation(r["tree_code"])[0]==names[0]]
        replicate=[r for r in matching if name_representation(r["tree_code"])[0]!=names[0]]
        if collision:
            status="blocked_source_taxonomic_collision"
        elif exact:
            status="public_original_tree_code_exact_name"
        elif replicate:
            status="public_original_tree_code_numbered_replicate_only"
        else:
            status="no_public_original_tree_code_name_match"
        out.append({
            "accepted_key":str(row["accepted_key"]),
            "original_source_taxon_names":"|".join(names),
            "S1_eligible_photos":int(row["n_source_orientation_photos"]),
            "original_grafted_194_tree_covered":bool(row["tree194_covered"]),
            "strict_191_tree_eligible":bool(row["strict_phylogenetic_eligible"]),
            "S2_same_source_species":bool(row["s2_focal"]),
            "moreyra_read_name_match_status":status,
            "n_matched_public_moreyra_reads":len(matching),
            "moreyra_runs":"|".join(sorted({r["run"] for r in matching})),
            "moreyra_biosamples":"|".join(sorted({r["biosample"] for r in matching})),
            "original_moreyra_tree_codes":"|".join(sorted({r["tree_code"] for r in matching})),
            "sample_name_match_is_validated_tree_tip":False,
            "source_taxonomic_concept_equivalence_validated":False,
        })
    frame=pd.DataFrame(out)
    status={k:int(v) for k,v in frame.moreyra_read_name_match_status.value_counts().items()}
    hasreads=frame.loc[frame.n_matched_public_moreyra_reads>0]
    s2=frame.loc[frame.S2_same_source_species]
    if len(s2)!=2 or set(s2.original_source_taxon_names)!=S2_NAMES:
        raise ValueError("Both native S1/S2 focus species must remain separately named")
    return frame,{
        "schema":"azami_v4_s1_all244_vs_pinned_EAzami_Moreyra2025_public_read_inventory_v1",
        "source_original_S1_native_taxa":244,
        "source_original_S1_orientation_photos":170066,
        "source_S1_current_tree_covered":194,
        "source_S1_strict_clean_tree":191,
        "source_Moreyra_public_Cirsium_SRA_sample_rows":258,
        "n_S1_taxa_with_public_raw_read_name_match":int(len(hasreads)),
        "n_S1_taxa_with_grafted_tree_and_Moreyra_public_reads":
           int((hasreads.original_grafted_194_tree_covered).sum()),
        "n_S1_taxa_without_grafted_tree_but_with_Moreyra_public_reads":
           int((~hasreads.original_grafted_194_tree_covered).sum()),
        "n_S1_taxa_without_Moreyra_public_raw_name_match":
           int((frame.n_matched_public_moreyra_reads==0).sum()),
        "taxon_coverage_status_counts":status,
        "S2_two_species":[{
            "taxon":r["original_source_taxon_names"],
            "S1_in_191":r["strict_191_tree_eligible"],
            "n_public_nuclear_reads":r["n_matched_public_moreyra_reads"],
            "name_match_status":r["moreyra_read_name_match_status"]
        } for r in s2.to_dict("records")],
        "disallowed_claims":[
            "All 244 original S1 concepts already have genuine nuclear sequences",
            "SRA original run-name match means a validated orthologous nuclear species tree",
            "Moreyra final exact 350 loci tree or dated tree is available in these parts",
            "Numbered replicate labels are new independently evolving taxon concepts",
            "Unmatched original tree code proves nuclear sequence is globally unavailable",
        ],
        "next_scientific_gate":"Acquire/derive an independently rooted nuclear supertree with audited tip/voucher/accepted-key mapping, branch lengths, topology uncertainty; separately adjudicate botanical orientation.",
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--all244",type=Path,required=True)
    p.add_argument("--out-dir",type=Path,required=True)
    args=p.parse_args()
    manifest,moreyra,manifest_sha=load_moreyra()
    frame,summary=analyse(pd.read_csv(args.all244),moreyra)
    summary["EAzami_commit"]=EA_COMMIT
    summary["EAzami_manifest_sha256"]=manifest_sha
    summary["EAzami_original_reconciliation_sha256"]=manifest["source_full_reconciliation_sha256"]
    args.out_dir.mkdir(parents=True,exist_ok=True)
    frame.to_csv(args.out_dir/"S1_ALL244_MOREYRA_PUBLIC_READ_CROSSWALK.csv",index=False)
    (args.out_dir/"S1_ALL244_MOREYRA_PUBLIC_READ_COVERAGE.json").write_text(
        json.dumps(summary,sort_keys=True,ensure_ascii=False,indent=2)+"\n")
    print("ALL244_MOREYRA_PUBLIC_READS="+json.dumps(summary,sort_keys=True),flush=True)


if __name__=="__main__":
    main()
