"""Original-sourced botanical head-state sampling gate after 125 REAL SRA mappings.

Only exact author names (terminal integer specimen label stripped) are matched.
No synonym/rank collapse. This is NOT a test of rainfall-driven evolution or
nuclear tree placement, and image EXIF angles never become U/D states.
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

TRAIT_REPO="carol-siniscalchi/cirsium"
TRAIT_SHA="c925d4af86160064448f794a9715f70b917a2786"
TRAIT_PATH="traits/traits.csv"
TRAIT_BLOB_SHA="fa127f0587b676e6e004ea121424d26069586bcb"
EA_REPO="zuizui0223/EAzami"
EA_SHA="c9870e91604ed1514344ddc8f765eae5237dbf27"
EA_DIR="data/evidence/moreyra2025_cirsium_reconciliation_v1"
N_TAXA=63
N_NODDING=5
MOREYRA_ROWS=258
S1_RUNS=125


def get_bytes(url):
    with urllib.request.urlopen(url,timeout=95) as r:
        return r.read()


def blob_sha(content):
    return hashlib.sha1(b"blob "+str(len(content)).encode()+b"\x00"+content).hexdigest()


def download_published_states(fetch=get_bytes):
    raw=fetch(f"https://raw.githubusercontent.com/{TRAIT_REPO}/{TRAIT_SHA}/{TRAIT_PATH}")
    if blob_sha(raw)!=TRAIT_BLOB_SHA:
        raise ValueError("Author botanical HeadPosition original Git blob drift")
    rows=list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
    if len(rows)!=N_TAXA or len({z["Species"] for z in rows})!=N_TAXA:
        raise ValueError("63 author-scored North American taxa missing/duplicated")
    states={}
    for row in rows:
        if row["HeadPosition"] not in {"0","1"}:
            raise ValueError("Published author HeadPosition no longer binary")
        states[row["Species"]]=int(row["HeadPosition"])
    if sum(states.values())!=N_NODDING:
        raise ValueError("Published five independently scored nodding species changed")
    return states


def normalize_author_name(source):
    """Remove terminal specimen NUMBER ONLY, never 'var.' or subsp/rank."""
    name=re.sub(r"\s+[0-9]+$","",source.strip())
    return "_".join(name.split())


def download_original_sra_index(fetch=get_bytes):
    url=f"https://raw.githubusercontent.com/{EA_REPO}/{EA_SHA}/{EA_DIR}/"
    raw=fetch(url+"manifest.json")
    manifest=json.loads(raw.decode("utf-8-sig"))
    if manifest["rows"]!=MOREYRA_ROWS or len(manifest["parts"])!=8:
        raise ValueError("Original EAzami 258 run reconciliation changed")
    by_run=defaultdict(list)
    n=0
    for part in manifest["parts"]:
        content=fetch(url+part["file"])
        if hashlib.sha256(content).hexdigest()!=part["sha256"] or len(content)!=part["bytes"]:
            raise ValueError("EAzami author source SRA manifest/digest drift")
        rows=list(csv.DictReader(io.StringIO(content.decode("utf-8-sig"))))
        if len(rows)!=part["rows"]:
            raise ValueError("EAzami original source SRA part count drift")
        for row in rows:
            if row["sra_link_status"]=="linked_runinfo":
                by_run[row["run"]].append(row)
        n+=len(rows)
    if n!=MOREYRA_ROWS:
        raise ValueError("Missing EAzami source SRA rows")
    return by_run


def evaluate(receipt,original,states):
    if (receipt.get("read_results",{}).get("mapped_run_count")!=S1_RUNS or
        receipt.get("read_results",{}).get("failed_run_count")!=0 or
        receipt.get("workflow_conclusion")!="success" or
        len(receipt.get("per_original_public_run",[]))!=S1_RUNS or
        receipt.get("source_moreyra_name_matched_taxon_keys")!=113):
        raise ValueError("125-run empirical source mapping proof missing")
    join=[]
    source_conflicts=[]
    for row in receipt["per_original_public_run"]:
        run=row["run"]
        orig=original.get(run)
        if not orig:
            raise ValueError("S1 real mapped run missing from Moreyra original SRA")
        names={normalize_author_name(z["tree_code"]) for z in orig}
        if len(names)!=1:
            if any(name in states for name in names):
                raise ValueError("One SRA run has ambiguous original tree-code botanical HeadPosition")
            source_conflicts.append({
                "run":run,"original_tree_code_names":sorted(names),
                "botanical_state_attribution":"blocked_ambiguous_source_labels",
                "read_mapping_success_preserved":True
            })
            continue
        name=next(iter(names))
        if name not in states:
            continue
        join.append({
            "run":run,"original_source_taxon":name,
            "author_coded_botanical_nodding":bool(states[name]),
            "read_pairs":int(row["pairs"]),
            "reference_loci_with_any_read":int(row["reference_loci_at_least_one_alignment"]),
            "source_taxon_concept_verified":False,
            "original_Moreyra_350_locus_tree_built":False,
        })
    taxa=defaultdict(list)
    for r in join:taxa[r["original_source_taxon"]].append(r)
    down=[x for x,v in taxa.items() if states[x]==1]
    up=[x for x,v in taxa.items() if states[x]==0]
    if (len(join)!=32 or len(taxa)!=26 or len(down)!=2 or
        len(up)!=24 or sum(x["author_coded_botanical_nodding"] for x in join)!=3):
        raise ValueError("Independent botanical head-state × original 125 run intersection drift")
    detailed=[]
    for name in sorted(taxa):
        a=sorted(taxa[name],key=lambda x:x["run"])
        detailed.append({
            "taxon":name,
            "botanical_state":"nodding" if states[name] else "erect",
            "n_public_SRA_runs":len(a),
            "source_SRA_run_ids":[x["run"] for x in a],
            "real_reference_locus_counts":[x["reference_loci_with_any_read"] for x in a],
            "reference_locus_min":min(x["reference_loci_with_any_read"] for x in a),
            "reference_locus_max":max(x["reference_loci_with_any_read"] for x in a),
            "source_concept_taxonomic_voucher_verification_pending":True,
        })
    scop=next(x for x in detailed if x["taxon"]=="Cirsium_scopulorum")
    if (scop["reference_locus_min"],scop["reference_locus_max"])!=(589,965):
        raise ValueError("Published nodding C. scopulorum source replicate evidence drift")
    return {
        "schema":"azami_v4_t2_na63_original_head_position_Moreyra125_real_read_evidence_gate_v1",
        "status":"SOURCE_AUTHENTIC_NA_BOTANICAL_ORIENTATION_READ_ACCESS_OVERLAP_ONLY",
        "original_author_scored_north_american_taxa":N_TAXA,
        "original_botanical_nodding_scored_taxa":N_NODDING,
        "all_125_original_sra_runs_have_real_mappings":True,
        "multiple_source_tree_code_conflict_runs":source_conflicts,
        "n_unassignable_source_conflict_runs":len(source_conflicts),
        "n_original_na_taxa_with_exact_named_S1_Moreyra_mapped_runs":len(taxa),
        "n_original_na_nodding_taxa_with_exact_named_runs":len(down),
        "n_original_na_erect_taxa_with_exact_named_runs":len(up),
        "n_exact_mapped_runs_with_author_botanical_state":len(join),
        "n_runs_nodding":sum(x["author_coded_botanical_nodding"] for x in join),
        "n_runs_erect":sum(not x["author_coded_botanical_nodding"] for x in join),
        "nodding_taxa_exact_name":sorted(down),
        "unmatched_to_any_source_name_is_not_gene_absence":True,
        "missing_author_scored_nodding_taxa_are_not_reclassified_as_erect":True,
        "S1_image_angle_is_not_used_as_botanical_U_D":True,
        "one_read_to_reference_not_an_ortholog_or_phylogenetic_reconstruction":True,
        "cannot_infer_four_independent_nodding_origins_from_two_nodding_taxa":True,
        "source_sra_run_to_taxon_original_voucher_unverified":True,
        "species_details":detailed,
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-receipt",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    receipt=json.loads(a.source_receipt.read_text(encoding="utf-8"))
    original=download_original_sra_index()
    states=download_published_states()
    out=evaluate(receipt,original,states)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("T2_REAL_MAPPED_NA63_AUTHOR_NODDING_SCOPE="+json.dumps({
       k:out[k] for k in ("n_original_na_taxa_with_exact_named_S1_Moreyra_mapped_runs",
          "n_original_na_nodding_taxa_with_exact_named_runs",
          "n_exact_mapped_runs_with_author_botanical_state",
          "n_runs_nodding")},sort_keys=True),flush=True)

if __name__=="__main__":
    main()
