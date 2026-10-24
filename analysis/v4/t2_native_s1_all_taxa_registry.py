"""Export ALL 244 S1 image-eligible accepted keys and public nuclear-source coverage.

Reads owner-authenticated native source *only* in Actions. Exports taxon-level
summaries, never records, coordinates, photos or observers. Phylogenetic source
name matches are candidates, NOT accepted-concept or final-tree validation.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import io
import json
from pathlib import Path
import re
import urllib.request

import pandas as pd

from analysis.v4.joint_s1_t2_phylogeny import (
    ENDPOINT, EXPOSURES, MODEL_PATH, COVERAGE_PATH, MAPPING_PATH,
    actual_native_roster, source_concept_index,
)
from analysis.v4.t2_pgls_placement_screen import download_author_states

EA_COMMIT="c9870e91604ed1514344ddc8f765eae5237dbf27"
EA_ROOT=("https://raw.githubusercontent.com/zuizui0223/EAzami/"
         +EA_COMMIT+"/data/evidence/")
EA_JAPAN="japan38_comp1061_concept_map_v1.csv"
EA_EASTASIA="orientation_comp1061_20tip_source_crosswalk_v1.csv"
MOREYRA_DIR="moreyra2025_cirsium_reconciliation_v1/"
EA_JAPAN_N=38
EA_EASTASIA_N=20
MOREYRA_ROWS=258


def fetch_bytes(url):
    with urllib.request.urlopen(url,timeout=100) as r:
        return r.read()


def generic_taxon_label(value):
    """Remove a known author suffix but never collapse infraspecific ranks."""
    x=value.replace("_"," ").strip().split()
    if len(x)<2 or x[0]!="Cirsium" or not x[1].islower():
        return ""
    short=" ".join(x[:2])
    for i,token in enumerate(x[2:],start=2):
        if token in ("var.","subsp.","ssp.","f.") and i+1<len(x):
            short+=" "+("subsp." if token=="ssp." else token)+" "+x[i+1]
            break
    return short


def load_nuclear_candidate_names():
    states,sha=download_author_states()
    na={generic_taxon_label(t) for t in states}
    ja=list(csv.DictReader(io.StringIO(fetch_bytes(EA_ROOT+EA_JAPAN).decode("utf-8-sig"))))
    ea=list(csv.DictReader(io.StringIO(fetch_bytes(EA_ROOT+EA_EASTASIA).decode("utf-8-sig"))))
    if len(ja)!=EA_JAPAN_N or len(ea)!=EA_EASTASIA_N:
        raise ValueError("EAzami Japan38 or EastAsia20 source roster changed")
    japan={generic_taxon_label(x["paper_taxon_concept"]) for x in ja}
    east={generic_taxon_label(x["accepted_taxon"]) for x in ea}
    manifest=json.loads(fetch_bytes(EA_ROOT+MOREYRA_DIR+"manifest.json"))
    if manifest["rows"]!=MOREYRA_ROWS or len(manifest["parts"])!=8:
        raise ValueError("Moreyra SRA reconciliation provenance drifted")
    moreyra=defaultdict(set)
    seen_records=0
    for item in manifest["parts"]:
        raw=fetch_bytes(EA_ROOT+MOREYRA_DIR+item["file"])
        if hashlib.sha256(raw).hexdigest()!=item["sha256"]:
            raise ValueError("Moreyra source panel SHA256 mismatch")
        rows=list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
        if len(rows)!=item["rows"]:
            raise ValueError("Moreyra source panel row count drift")
        for row in rows:
            if row["sra_link_status"]!="linked_runinfo":
                continue
            stem=re.sub(r" \d+$","",row["tree_code"].strip())
            label=generic_taxon_label(stem)
            if label:
                moreyra[label].add(row["biosample"])
        seen_records+=len(rows)
    if seen_records!=MOREYRA_ROWS:
        raise ValueError("Frozen original 258-run metadata audit changed")
    return {
        "north_america_nuclear":na,
        "EAzami_japan38":japan,
        "EAzami_eastasia20":east,
        "Moreyra2025_SRA_candidate":moreyra,
        "provenance":{
           "EAzami_commit":EA_COMMIT,
           "NA_published_trait_sha256":sha,
           "Moreyra_reconciled_run_metadata_rows":seen_records,
           "Moreyra_original_final_phylogeny_present":False,
           "candidate_name_match_is_verified_taxon_concept":False,
           "japan38_nuclear_tip_is_a_sample_not_identical_to_taxon_concept":True
        }
    }


def census(plain):
    models=json.loads(MODEL_PATH.read_text(encoding="utf-8"))
    coverage=json.loads(COVERAGE_PATH.read_text(encoding="utf-8"))[ENDPOINT]
    if hashlib.sha256(plain.read_bytes()).hexdigest()!=models["input_sha256"]:
        raise ValueError("Original S1 170066 selected observation SHA mismatch")
    if tuple(models["selection"]["selected"])!=EXPOSURES:
        raise ValueError("Source S1 fixed climatic axes changed")
    table=pd.read_csv(plain,low_memory=False)
    roster=actual_native_roster(table,models,coverage)
    mapping=source_concept_index(json.loads(MAPPING_PATH.read_text(encoding="utf-8")))
    nuclear=load_nuclear_candidate_names()
    if len(roster)!=244 or sum(r["n"] for r in roster.values())!=170066:
        raise ValueError("Full 244 source concepts not recovered")
    out=[]
    for key,info in sorted(roster.items(),key=lambda kv:int(kv[0])):
        original=info["source_names"]
        if len(original)==1:
            name=generic_taxon_label(original[0])
        else:
            name=""
        record=mapping.get(key,{"status":"unmapped","tree_matches":set()})
        candidate=bool(name)
        moreyra=nuclear["Moreyra2025_SRA_candidate"].get(name,set()) if candidate else set()
        item={
            "accepted_key":key,
            "source_names":"|".join(original),
            "source_concept_collision":len(original)!=1,
            "normalized_exact_source_name_candidate":name,
            "eligible_photo_observations":info["n"],
            "S1_mean_EXIF_angle_deg":info["mean_image_angle_deg"],
            "S1_median_EXIF_angle_deg":info["median_image_angle_deg"],
            "grafted_placement_status":record["status"],
            "grafted_placement_tip_candidates":"|".join(sorted(record["tree_matches"])),
            "NA63_nuclear_taxon_name_candidate":candidate and name in nuclear["north_america_nuclear"],
            "EAzami_JP38_nuclear_taxon_name_candidate":candidate and name in nuclear["EAzami_japan38"],
            "EAzami_EA20_nuclear_taxon_name_candidate":candidate and name in nuclear["EAzami_eastasia20"],
            "Moreyra2025_SRA_taxon_name_candidate":bool(moreyra),
            "Moreyra2025_linked_source_BioSample_candidates":len(moreyra),
        }
        for e in EXPOSURES:
            item[e+"_photo_eligible_mean"]=info["mean_frozen_five_exposures"][e]
        out.append(item)
    tree_matched=sum(x["grafted_placement_status"]=="matched" for x in out)
    if tree_matched!=194:
        raise ValueError("S1 matched tree taxonomy coverage changed")
    return out, {
        "total_original_S1_image_eligible_WCVP_keys":len(out),
        "total_original_S1_image_observations":sum(r["eligible_photo_observations"] for r in out),
        "n_grafted_tree_matched":tree_matched,
        "n_not_matched_to_52_grafted_placements":len(out)-tree_matched,
        "source_names_with_intra_key_collision":sum(r["source_concept_collision"] for r in out),
        "source_name_candidate_coverage":{
            k:sum(r[k] for r in out)
            for k in ("NA63_nuclear_taxon_name_candidate",
                      "EAzami_JP38_nuclear_taxon_name_candidate",
                      "EAzami_EA20_nuclear_taxon_name_candidate",
                      "Moreyra2025_SRA_taxon_name_candidate")
        },
        "T2_RAIN_ANGLE_UNIVERSE":"ALL native-range S1 trait-eligible concepts; no after-outcome species selection",
        "source_identity":"owner-authenticated original S1 decoded source SHA + public taxon-level aggregate",
        "phylogeny_rules":"52 grafted placement scenarios are not molecular; named SRA candidate not final Moreyra tree",
        "missing_nuclear_species_not_imputed":True,
        "separately_verified_S2_native_introduced":{"Cirsium arvense":True,"Cirsium vulgare":True},
        "external_phylo_sources":nuclear["provenance"],
        "full_species_adaptive_radiation_confirmed":False,
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--decoded-native-table",required=True,type=Path)
    p.add_argument("--out-dir",required=True,type=Path)
    args=p.parse_args()
    rows,summary=census(args.decoded_native_table)
    args.out_dir.mkdir(parents=True,exist_ok=True)
    with (args.out_dir/"S1_ALL_244_ACCEPTED_TAXON_CENSUS.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)
    (args.out_dir/"S1_ALL_244_SOURCE_CONCEPT_PHYLO_COVERAGE.json").write_text(
        json.dumps(summary,sort_keys=True,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("T2_FULL_S1_CENSUS="+json.dumps(summary,ensure_ascii=False,sort_keys=True))


if __name__=="__main__":
    main()
