"""Actual 131-tip Senecioneae linked-rDNA genealogy and direct nodding evidence.

Reads newly computed original GenBank ITS+ETS IQ-TREE from GitHub Actions
39961377067 (NOT a published whole-nuclear coalescent species tree). Compares
genuinely observed D6 and 125 UNKNOWN, exact source original taxa, nearest
molecular neighbor, root-independent genus exclusivity and node confidence.

Crucially no independently scored flowering-stage upright heads were matched
on these 131 DNA tips. Thus "different genera" is NOT "independent origins",
an ITS/ETS pattern is NOT a reticulation-proof species phylogeny, and this
genetic tree alone says NOTHING about precipitation-adaptive fitness.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
import re
from collections import Counter
from pathlib import Path

from Bio import Phylo

TREE_SHA="fb1c2019b06e3691fac13e19d2f250ad7a34287fbd29abd209b9405d30748310"
TREE_RUN=39961377067
N_GENUINE_GB_ITS_ETS=131
D_N=6
UNKNOWN_N=125
STRONG_ALRT=80.
STRONG_UFBOOT=95.

def read_source(treefile,crosswalk):
    if hashlib.sha256(Path(treefile).read_bytes()).hexdigest()!=TREE_SHA:
        raise ValueError("True newly inferred original 131 Senecioneae GenBank ML gene tree SHA drift")
    tree=Phylo.read(str(treefile),"newick")
    tips=tree.get_terminals()
    names=[x.name.replace("_"," ") for x in tips]
    if len(tips)!=N_GENUINE_GB_ITS_ETS or len(set(names))!=N_GENUINE_GB_ITS_ETS:
        raise ValueError("Original true 131 linked-rDNA taxon tree changed")
    with Path(crosswalk).open(encoding="utf-8-sig",newline="") as f:
        table=list(csv.DictReader(f))
    if len(table)!=N_GENUINE_GB_ITS_ETS or set(names)!={r["GenBank_original_published_species"] for r in table}:
        raise ValueError("Phylogenetic original source and independently scored morphology are not the same specimen-name list")
    if any(r["linked_multicopy_rDNA_not_single_copy_independent_genes"]!="True" or
           r["DNA_voucher_equals_flora_specimen"]!="False" for r in table):
        raise ValueError("ITS ETS nuclear linkage or DNA vs flora vouchers were falsely called independent")
    state={r["GenBank_original_published_species"]:r["original_species_morphology_state"] for r in table}
    if Counter(state.values())!={"D":D_N,"?":UNKNOWN_N}:
        raise ValueError("Direct head posture must remain only six D, 125 unmeasured")
    return tree,state

def root_independent_split_is_exclusive(tree,original_genus):
    names={t.name for t in tree.get_terminals()}
    group={x for x in names if x.split("_")[0]==original_genus}
    if not 2<=len(group)<=len(names)-2:return False
    for node in tree.find_clades():
        subset={x.name for x in node.get_terminals()}
        if subset==group or names-subset==group:
            return True
    return False

def audit(tree,state):
    tips={x.name:x for x in tree.get_terminals()}
    if len(tips)!=N_GENUINE_GB_ITS_ETS:return None
    D=[name.replace(" ","_") for name,v in state.items() if v=="D"]
    D.sort()
    if len(D)!=D_N or not set(D)<=set(tips):
        raise ValueError("Six truly documented Senecioneae head postures not mapped to 131 nuclear sequence labels")
    genera=Counter(x.split("_")[0] for x in D)
    if genera!={"Parasenecio":4,"Cremanthodium":2}:
        raise ValueError("Original published direct flora-matched two genera marker coverage drift")
    scores=[]
    for internal in tree.get_nonterminals():
        if not internal.name:continue
        match=re.fullmatch(r"([0-9]+(?:\.[0-9]+)?)/([0-9]+(?:\.[0-9]+)?)",internal.name)
        if match:scores.append(tuple(float(a) for a in match.groups()))
    if len(scores)!=128:
        raise ValueError("Real 131-tip published-source ML did not preserve standard IQTREE support")
    neighbors=[]
    for original in D:
        other=sorted(
            ((float(tree.distance(original,x)),x)
             for x in tips if x!=original),
            key=lambda x:(x[0],x[1]))
        close_D=sorted(((float(tree.distance(original,x)),x)
                        for x in D if x!=original))
        if not other or not close_D:
            raise ValueError("No observed nuclear genealogy nearest-tip pattern")
        neighbors.append({
            "source_D_species":original.replace("_"," "),
            "nearest_full_nuclear_ITS_ETS_GenBank_tip":other[0][1].replace("_"," "),
            "nearest_full_tip_is_floristically_scored_D":state[other[0][1].replace("_"," ")]=="D",
            "nearest_original_gene_tree_patristic_substitutions_per_site":other[0][0],
            "nearest_another_D_species_name":close_D[0][1].replace("_"," "),
            "nearest_another_D_gene_tree_patristic_distance":close_D[0][0],
            "floral_phenotype_not_inferred_for_unscored_nearest_tip":True,
            "this_is_not_independent_adaptive_gain_event":True,
        })
    return {
        "schema":"azami_T3_2025_real_genbank_linked_ITS_ETS_genealogy_senecioneae_D_posture_structure_v1",
        "status":"REAL_131_SOURCE_GENBANK_LINKED_RDNA_TOPOLOGY_6_D_125_UNKNOWN_NOT_ORIGIN_TEST",
        "original_2025_genbank_linked_rDNA_sequence_taxa":N_GENUINE_GB_ITS_ETS,
        "original_D_sequenced_taxa":D_N,
        "head_orientation_unknown_taxa":UNKNOWN_N,
        "documented_downward_ITS_ETS_taxa_by_genus":dict(genera),
        "genus_scoped_marker_tip_count":{
            g:sum(x.startswith(g+"_") for x in tips)
            for g in ("Cremanthodium","Parasenecio","Ligularia","Sinosenecio")
        },
        "exact_genus_tip_sets_exclusive_under_root_invariant_gene_tree_splits":{
            g:root_independent_split_is_exclusive(tree,g)
            for g in ("Cremanthodium","Parasenecio","Ligularia","Sinosenecio")
        },
        "actual_IQTREE_internal_nodes_with_support_labels":len(scores),
        "actual_node_with_SHaLRT_atleast80_and_UFBoot_atleast95":
            sum(a>=STRONG_ALRT and b>=STRONG_UFBOOT for a,b in scores),
        "actual_node_with_UFBoot_atleast95":sum(b>=STRONG_UFBOOT for a,b in scores),
        "nodding_specimen_nearest_real_genealogical_tip":neighbors,
        "genuinely_sequenced_but_head_orientation_unknown_not_counted_as_erect":True,
        "genus_phylogenetic_nonexclusivity_in_linked_rDNA_not_synonymy_or_accepted_species_history":True,
        "root_direction_species_tree_and_ancestral_D_origins_not_identified":True,
        "linked_multicopy_ITS_ETS_not_independent_coalescent_genes":True,
        "rainfall_at_original_head_transition_or_fitness_not_measured":True,
        "tree_source_sha256":TREE_SHA,
        "tree_original_public_2025_genbank_ml_action":TREE_RUN,
    }

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tree",type=Path,required=True)
    p.add_argument("--head-crosswalk",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    t,states=read_source(a.tree,a.head_crosswalk)
    result=audit(t,states)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,sort_keys=True,
                                ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    print("T3_TRUE_131_D_GENBANK_GENEALOGY="+json.dumps({
        "D_genus":result["documented_downward_ITS_ETS_taxa_by_genus"],
        "n_strong_node":result["actual_node_with_SHaLRT_atleast80_and_UFBoot_atleast95"],
        "n_erect_source_tips":0,
        "n_real_taxa":result["original_2025_genbank_linked_rDNA_sequence_taxa"],
        "nearest_genetically_to_C_campanulatum":next(z["nearest_full_nuclear_ITS_ETS_GenBank_tip"] for z in result["nodding_specimen_nearest_real_genealogical_tip"] if z["source_D_species"]=="Cremanthodium campanulatum"),
    },ensure_ascii=False),flush=True)

if __name__=="__main__":main()
