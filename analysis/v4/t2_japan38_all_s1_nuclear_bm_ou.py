"""All exact-name-candidate S1 photographic taxa on EAzami Japan38 genuine nuclear tree.

EAzami Japan38 has 39 ingroup biological specimens representing 38 nominal
taxa; the original JPN20 replicated non-monophyly and JPN31 source identity
conflict must not be hidden by taxon averages. Name matches are candidates
until exact WCVP concept/voucher reconciliation, not a completed global tree.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from copy import deepcopy
import hashlib
import io
import json
from pathlib import Path
import urllib.request

from Bio import Phylo
import numpy as np

from analysis.v4.t2_full_species_bm_ou import fit_three_models, FIXED_S1
from analysis.v4.t2_native_s1_all_taxa_registry import (
    EA_ROOT, EA_JAPAN, generic_taxon_label,
)
from analysis.v4.s1_continuous_phylogenetic_history import (
    audited_pruned_tree, shared_path_matrix,
)

ORIGINAL_TREE="japan38_comp1061_primary_tree_v1.nwk"
ACCEPTANCE="japan38_comp1061_primary_tree_acceptance_v1.json"
MIN_TREE_TIPS=7


def load_source():
    with urllib.request.urlopen(EA_ROOT+EA_JAPAN,timeout=100) as url:
        cmap_raw=url.read()
    with urllib.request.urlopen(EA_ROOT+ORIGINAL_TREE,timeout=100) as url:
        tree_raw=url.read()
    with urllib.request.urlopen(EA_ROOT+ACCEPTANCE,timeout=100) as url:
        manifest=json.load(url)
    if (manifest["tree_sha256"]!=hashlib.sha256(tree_raw).hexdigest()
            or manifest["tree_tip_count"]!=40
            or manifest["root_tip"]!="OUTGROUP_saff"):
        raise ValueError("EAzami recovered nuclear ML tree SHA/source audit drift")
    cmap=list(csv.DictReader(io.StringIO(cmap_raw.decode("utf-8-sig"))))
    if len(cmap)!=38 or len({r["paper_japan_member_id"] for r in cmap})!=38:
        raise ValueError("EAzami 38 taxonomy source drift")
    tree=Phylo.read(io.StringIO(tree_raw.decode("utf-8-sig")),"newick")
    tips=[t.name for t in tree.get_terminals()]
    if len(tips)!=40 or len(set(tips))!=40 or "OUTGROUP_saff" not in tips:
        raise ValueError("Wrong Japan38 independent nuclear topology")
    tree.root_with_outgroup({"name":"OUTGROUP_saff"})
    return tree,cmap, {
        "tree_sha256":hashlib.sha256(tree_raw).hexdigest(),
        "concept_source_sha256":hashlib.sha256(cmap_raw).hexdigest(),
        "tree_loci":manifest["rootable_loci"],
        "tree_alignment_bp":manifest["alignment_length_bp"],
        "branch_unit":manifest["branch_length_semantics"],
        "multiple_samples_nonmonophyly_resolved_by_name_collapse":False,
    }


def join_original_244(full244,cmap):
    with Path(full244).open(newline="",encoding="utf-8") as f:
        rows=list(csv.DictReader(f))
    if len(rows)!=244 or sum(int(x["eligible_photo_observations"]) for x in rows)!=170066:
        raise ValueError("Authentic all-244 census input changed")
    by_name=defaultdict(list)
    for row in rows:
        if row["source_concept_collision"]=="True":
            continue
        label=row["normalized_exact_source_name_candidate"]
        if label:
            by_name[label].append(row)
    used=[]
    status=[]
    for mapping in cmap:
        label=generic_taxon_label(mapping["paper_taxon_concept"])
        tips=mapping["tip_ids"].split("|")
        reason=None
        if mapping["trait_asr_primary_allowed"]!="true":
            reason="EAzami_disallowed_primary_identity_or_concept"
        elif len(tips)!=1 or mapping["n_biological_samples"]!="1":
            reason="EAzami_nonmonophyletic_replications_not_collapsed"
        elif not label or len(by_name.get(label,[]))!=1:
            reason=("S1_name_unmatched" if not by_name.get(label) else
                    "ambiguous_multi_WCVP_concepts_for_same_label")
        else:
            sample=by_name[label][0]
            if sample["source_concept_collision"]!="False":
                reason="source_WCVP_name_collision"
            else:
                used.append({
                    "tip":tips[0],
                    "EAzami_JPN_concept_id":mapping["paper_japan_member_id"],
                    "EAzami_source_concept":mapping["paper_taxon_concept"],
                    "S1_accepted_key":sample["accepted_key"],
                    "n_S1_photo_eligible":int(sample["eligible_photo_observations"]),
                    "image_angle_deg":float(sample["S1_mean_EXIF_angle_deg"]),
                    "image_median_deg":float(sample["S1_median_EXIF_angle_deg"]),
                    "five_climate":{e:float(sample[e+"_photo_eligible_mean"]) for e in FIXED_S1},
                })
        status.append({
            "EAzami_JPN_concept_id":mapping["paper_japan_member_id"],
            "EAzami_source_concept":mapping["paper_taxon_concept"],
            "source_label_candidate":label,
            "status":"included_source_name_candidate_not_verified_accepted_concept" if reason is None else reason
        })
    if len({r["tip"] for r in used})!=len(used) or len({r["S1_accepted_key"] for r in used})!=len(used):
        raise ValueError("Japanese one-to-one tip/WCVP key violated")
    return used,status


def run(census_file):
    tree,cmap,provenance=load_source()
    rows,ledger=join_original_244(census_file,cmap)
    n=len(rows)
    if n < MIN_TREE_TIPS:
        return {"status":"BLOCKED_TOO_FEW_NUCLEAR_AND_S1_MATCHING_TAXA",
                "n_candidate_exact_names":n,"provenance":provenance,
                "full_japan38_ledger":ledger}
    tips=sorted(r["tip"] for r in rows)
    phy,audit=audited_pruned_tree(tree,tips)
    C=shared_path_matrix(phy,tips)[0]
    lookup={r["tip"]:r for r in rows}
    X0=np.ones((n,1))
    angle=np.array([lookup[t]["image_angle_deg"] for t in tips])
    models={"all_available_S1_image_angle_intercept_only":
            fit_three_models(angle,X0,np.diag(C),C)}
    # Even 1-image taxa stay in the main all-available panel, and are flagged.
    # The 5 covariates need enough taxa for an identifiable n-p-1 AICc.
    if n >= 11:
        five=np.array([[lookup[t]["five_climate"][v] for v in FIXED_S1] for t in tips])
        sd=five.std(axis=0)
        if np.any(sd<1e-12):
            models["five_climate_adjusted"]={"status":"blocked_constant_climate"}
        else:
            XX=np.column_stack([np.ones(n),(five-five.mean(axis=0))/sd])
            try:
                models["five_climate_adjusted"]=fit_three_models(
                    angle,XX,np.diag(C),C)
            except ValueError as e:
                models["five_climate_adjusted"]={"status":"blocked_identifiability","detail":str(e)}
    else:
        models["five_climate_adjusted"]={"status":"blocked_n_fewer_than_11"}
    deep=[r for r in rows if r["n_S1_photo_eligible"]>=10]
    models["n10_eligible_taxa_count"]=len(deep)
    return {
        "status":"SOURCE_NAME_MATCHED_NOT_YET_EXACT_ACCEPTED_CONCEPT_VERIFIED",
        "n_all_38_nominal_EAzami":38,
        "n_one_to_one_noncollided_name_matched":n,
        "n_one_image_species":sum(r["n_S1_photo_eligible"]==1 for r in rows),
        "n10_photo_taxa":len(deep),
        "source_matched_taxa":rows,
        "model_results":models,
        "full_japan38_ledger":ledger,
        "model_assumption":"original EAzami independent 176-locus nuclear branch lengths substitutions/site",
        "physiological_orientation_validity":"S1 photo-axis proxy, NOT botanical U/D nor gravity",
        "EAzami_discrete_Mk_literature":"19 resolved Japan38 taxonomy concepts, ML minimum changes 6; UFBoot 4-6; separate original analyses",
        "not_an_independent_global_origin_count":True,
        "not_a_causal_rain_evolutionary_model":True,
        "root_outgroup_is_Safflower":True,
        "zero_branches_regularized":audit["zero_length_edges_regularised"],
        "provenance":provenance,
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--s1-full244",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    out=run(a.s1_full244)
    path=Path(a.output)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(out,indent=2,ensure_ascii=False,sort_keys=True,
                               allow_nan=False)+"\n",encoding="utf-8")
    print("T2_JAPAN38_FULL_AVAILABLE_SOURCE_RESULT="+json.dumps({
        "status":out["status"],
        "n_photo_matched":out["n_one_to_one_noncollided_name_matched"] if "n_one_to_one_noncollided_name_matched" in out else out["n_candidate_exact_names"],
        "models":{k:v.get("preferred",v.get("status","count"))
                  for k,v in out.get("model_results",{}).items()
                  if isinstance(v,dict)}
    },ensure_ascii=False,sort_keys=True))


if __name__=="__main__":
    main()
