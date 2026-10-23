"""Independent North American 63-tip nuclear topology concordance for real111 Q20 ML.

The real111 technical source tree may retain 60-111 exact original WCVP keys;
the independently published Siniscalchi2023 nuclear 63-tip topology is
author-authenticated by the exact Git blob SHA in earlier Azami audits.

Compare ONLY exact original source names with accepted 111 run/S1 photo
measurements, with no guessed synonym or infraspecific identity. Prune
both trees to identical tips and report unrooted split RF and quartet
agreement. Quartet samples share tips/branches and are NOT independent
replications or binomial selection tests. This cannot infer direction
of botanical U->D transitions, dated adaptive radiation, or BM/OU.
"""
from __future__ import annotations
import argparse
from collections import Counter
from copy import deepcopy
import csv
from itertools import combinations
import json
from pathlib import Path

from Bio import Phylo
from analysis.v4.t2_real7_NA63_quartet_crosscheck import load_na63

MASTER_NAMES=251
WCVP_KEYS=244
SAMPLE_UNIVERSE=111
NA63=63
MIN_COMMON=4


def normalize(name):
    return "_".join(str(name).replace(" ","_").split("_")).strip("_")


def phylo_names(technical,samples):
    if len(samples)!=SAMPLE_UNIVERSE or len({x["accepted_key"] for x in samples})!=SAMPLE_UNIVERSE:
        raise ValueError("Full111 DNA source-authenticated taxon matching table missing")
    bytip={}
    for row in samples:
        tip=row["phylogenetic_tip"]
        if tip!="WCVP"+row["accepted_key"]:
            raise ValueError("Fictitious/WCVP misassigned sample label")
        name=normalize(row["taxon_name"])
        if name in bytip.values() or not name.startswith("Cirsium_"):
            raise ValueError("Duplicate original Cirsium source name for two tree tips")
        bytip[tip]=name
    copy=deepcopy(technical)
    original={z.name for z in copy.get_terminals()}
    if not original<=set(bytip):
        raise ValueError("Technical DNA tree has missing or invented source accepted-key tips")
    if len(copy.get_terminals())!=len(original):
        raise ValueError("Technical DNA tree duplicate named sample")
    for tip in copy.get_terminals():
        tip.name=bytip[tip.name]
    return copy


def restrict_tree(tree,names):
    clone=deepcopy(tree)
    available={r.name for r in clone.get_terminals()}
    if not names<=available:
        raise ValueError("Species not found in exact independent nuclear tree")
    for taxon in sorted(available-names):
        clone.prune(taxon)
    assert len(clone.get_terminals())==len(names)
    return clone


def unrooted_splits(tree,names):
    n=len(names)
    if n<4:
        raise ValueError("Unrooted topology needs at least four overlapping taxa")
    parts=set()
    for c in tree.get_nonterminals():
        desc={z.name for z in c.get_terminals()}
        if 2<=len(desc)<=n-2:
            one=tuple(sorted(desc))
            two=tuple(sorted(names-desc))
            parts.add(min(one,two))
    return parts


def distance_pairs(tree,names):
    return {(a,b):float(tree.distance(a,b)) for a,b in combinations(names,2)}


def quartet_choice(distance,quad):
    a,b,c,d=quad
    def dist(x,y):
        return distance[tuple(sorted((x,y)))]
    pairs=(
        (dist(a,b)+dist(c,d),(a,b)),
        (dist(a,c)+dist(b,d),(a,c)),
        (dist(a,d)+dist(b,c),(a,d)),
    )
    ranked=sorted(pairs,key=lambda x:x[0])
    gap=ranked[1][0]-ranked[0][0]
    if gap<=max(1e-9,1e-6*max(abs(x[0]) for x in ranked)):
        return None
    return tuple(sorted(ranked[0][1]))


def compare(technical,reference,samples):
    names=phylo_names(technical,samples)
    refnames={z.name for z in reference.get_terminals()}
    if len(refnames)!=NA63:
        raise ValueError("Independent published nuclear NA63 original taxonomy changed")
    got={x.name for x in names.get_terminals()}
    shared=got&refnames
    if len(shared)<MIN_COMMON:
        return {
            "status":"NOT_ENOUGH_EXACT_NAMES_TO_COMPARE_INDEPENDENT_TREES",
            "exact_original_taxon_overlap":len(shared),
            "technical_unrooted_NA63_validation_completed":False,
            "exact_taxonomic_voucher_identity_verified":False,
            "original_S1_244_keys_and_251_named_units_preserved_elsewhere":True,
        }
    source=restrict_tree(names,shared)
    gold=restrict_tree(reference,shared)
    ss=unrooted_splits(source,shared)
    gs=unrooted_splits(gold,shared)
    ordered=sorted(shared)
    d1=distance_pairs(source,ordered)
    d2=distance_pairs(gold,ordered)
    tally=Counter()
    for quad in combinations(ordered,4):
        a=quartet_choice(d1,quad)
        b=quartet_choice(d2,quad)
        if a is None or b is None:
            tally["unresolved_distance"]+=1
        elif a==b:
            tally["same_unrooted_quartet"]+=1
        else:
            tally["discordant_unrooted_quartet"]+=1
    total=sum(tally.values())
    if total!=len(list(combinations(shared,4))):
        raise ValueError("Missing one source phylogenetic quartet")
    resolved=tally["same_unrooted_quartet"]+tally["discordant_unrooted_quartet"]
    return {
        "schema":"azami_v4_111_realQ20_118copy_vs_63_published_nuclear_source_topology_v1",
        "status":"INDEPENDENT_PUBLISHED_NUCLEAR_TOPOLOGY_CONCORDANCE_TECHNICAL_ONLY",
        "technical_tree_original_tip_count":len(got),
        "published_independent_tree_original_tips":NA63,
        "exact_original_taxon_overlap":len(shared),
        "overlap_names":ordered,
        "technical_unrooted_splits":len(ss),
        "published_na63_unrooted_splits":len(gs),
        "shared_unrooted_splits":len(ss&gs),
        "normalized_unrooted_RF":len(ss^gs)/max(1,len(ss)+len(gs)),
        "n_possible_quartets":total,
        "n_resolved_same_quartets":tally["same_unrooted_quartet"],
        "n_resolved_discordant_quartets":tally["discordant_unrooted_quartet"],
        "n_unresolved_quartets":tally["unresolved_distance"],
        "resolved_quartet_concordance_fraction":
            tally["same_unrooted_quartet"]/resolved if resolved else None,
        "quartets_are_heavily_nonnindependent_no_binomial_test":True,
        "two_author_coded_botanical_nodding_in_global_panel":
            sorted(shared&{"Cirsium_conspicuum","Cirsium_scopulorum"}),
        "independent_true_nuclear_tree_is_NA63_only":True,
        "real111_technical_tree_is_not_single_copy_orthology_species_tree":True,
        "technical_tree_unrooted_not_chronogram":True,
        "original_244_accepted_and_251_names_preserved":True,
        "photo_EXIF_angle_not_botanical_gravity_nodding":True,
        "adaptive_radiation_or_environmental_causality_tested":False,
        "taxonomic_voucher_same_concept_across_studies_verified":False,
    }


def load_csv(path):
    with Path(path).open(newline="",encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("--technical-tree",type=Path,required=True)
    a.add_argument("--source-samples",type=Path,required=True)
    a.add_argument("--output",type=Path,required=True)
    args=a.parse_args()
    tech=Phylo.read(str(args.technical_tree),"newick")
    gold=load_na63()
    result=compare(tech,gold,load_csv(args.source_samples))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,sort_keys=True,ensure_ascii=False,indent=2)+"\n",
                           encoding="utf-8")
    print("T2_REAL111_VS_NA63_INDEPENDENT_UNROOTED="+json.dumps({
        "overlap":result["exact_original_taxon_overlap"],
        "status":result["status"],
        "quartet_concordance":result.get("resolved_quartet_concordance_fraction"),
        "normalized_RF":result.get("normalized_unrooted_RF")},sort_keys=True),flush=True)


if __name__=="__main__":
    main()
