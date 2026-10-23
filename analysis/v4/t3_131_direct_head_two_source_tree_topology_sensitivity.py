"""T3: exact DNA topology sensitivity of Senecioneae directly described heads.

Uses original same 131 GenBank ITS+ETS nuclear ribosomal tips under:
(1) IQ-TREE maximum likelihood .treefile,
(2) IQ-TREE UFBoot consensus .contree (same gene/data),
(3) BioNJ initial distance .bionj (same gene/data, NOT independent marker).
The original Chinese Flora independently documented down head of Ligularia
purdomii (D9 among 3 genera), 2 other true D in Cremanthodium (one
functionally manipulated to show rain/UV shelter), 6 in Parasenecio,
1 documented ONTOGENETIC V in P. koualapensis; 121 otherwise unknown.
No studied taxon has an independently confirmed fixed U on this 131-tip DNA
set. LCA size from displayed Newick root is NOT independent gain count.
Use unrooted split membership when comparing topologies.

BioNJ retaining/dropping the 4-tip split is a topology-sensitivity probe,
not an independent sequence or study replication. ITS/ETS are linked copies
of multicopy nuclear rDNA; no polyploid homeolog, gene flow or specimen-level
flower/voucher verification. Rain-shelter ancestral transitions NOT inferred.
"""
from __future__ import annotations
import argparse
import csv
from collections import Counter
import hashlib
import io
import json
from pathlib import Path
from Bio import Phylo

from analysis.v4.t3_131_chinese_flora_new_Ligularia_nodding import (
    reconcile, SHELTER, FLORA, REQUIRED_LCA_TIPS,
)
from analysis.v4.t3_2025_real_nuclear_rDNA_phylogeny_posture_topology import TREE_SHA

PINNED={
    "max_likelihood":{
        "suffix":"treefile",
        "sha256":"fb1c2019b06e3691fac13e19d2f250ad7a34287fbd29abd209b9405d30748310",
    },
    "ultrafast_bootstrap_consensus":{
        "suffix":"contree",
        "sha256":"6f674d4a256f3eefc3ec840458fed328adda0896d71dc76b47a2cad206b00761",
    },
    "BioNJ_initial_distance_tree":{
        "suffix":"bionj",
        "sha256":"719239a205533ebc1f18b07bc9f84548761c56da3f24c7a3a73b4959b0b26a27",
    }
}
PAIR_STAGED={"Parasenecio_cyclotus","Parasenecio_koualapensis"}
PAIR_OTHER={"Cremanthodium_campanulatum","Ligularia_purdomii"}
TIPS=131


def read_tree(path,sha):
    payload=Path(path).read_bytes()
    if hashlib.sha256(payload).hexdigest()!=sha:
        raise ValueError("Original 131 nuclear rDNA source topology SHA changed")
    t=Phylo.read(io.StringIO(payload.decode("utf-8-sig")),"newick")
    terminals=[x.name for x in t.get_terminals()]
    if len(terminals)!=TIPS or len(set(terminals))!=TIPS:
        raise ValueError("Gene tree with 131 biological ITS+ETS original tips lost")
    return t


def unrooted_split_exists(t,expected):
    alltips={x.name for x in t.get_terminals()}
    if not expected<=alltips or not expected or len(expected)>=len(alltips):
        raise ValueError("Observed real source taxon split not a strict molecular subset")
    # Pair splits on displayed topology are true undirected branches, not roots.
    for node in t.find_clades():
        names={x.name for x in node.get_terminals()}
        if names==expected or alltips-names==expected:
            return True
    return False


def nearest_tip(t,source):
    name=source.replace(" ","_")
    candidates=sorted(
        [(float(t.distance(name,z.name)),z.name) for z in t.get_terminals()
         if z.name!=name],
        key=lambda x:(x[0],x[1]))
    return {"species":candidates[0][1].replace("_"," "),
            "distance_substitutions_per_site":candidates[0][0]}


def actual_topology(treepath,conspath,bionjpath,original131,genbank163,florafile,stagefile):
    files={"max_likelihood":treepath,
           "ultrafast_bootstrap_consensus":conspath,
           "BioNJ_initial_distance_tree":bionjpath}
    species_source=reconcile(treepath,original131,florafile,genbank163)
    with Path(stagefile).open(encoding="utf-8-sig",newline="") as f:
        stage=list(csv.DictReader(f))
    if (len(stage)!=1 or
        stage[0]["2025_original_species"]!="Parasenecio koualapensis" or
        stage[0]["source_posture"]!="V" or
        stage[0]["developmental_scope"]!="early_spreading_later_pendulous" or
        stage[0]["is_true_fixed_flowering_erect_control"]!="no"):
        raise ValueError("Original independent floral ontogenetic observation changed")
    results={}
    reference_tips=None
    for key,meta in PINNED.items():
        if Path(files[key]).suffix!="."+meta["suffix"]:
            raise ValueError("Unknown genuine IQTREE source tree format")
        tree=read_tree(files[key],meta["sha256"])
        tips={x.name for x in tree.get_terminals()}
        if reference_tips is not None and tips!=reference_tips:
            raise ValueError("Cannot compare same 131-tip nuclear gene trees with disjoint taxa")
        reference_tips=tips
        lca=tree.common_ancestor(*sorted(PAIR_OTHER))
        lca2=tree.common_ancestor(*sorted(PAIR_STAGED))
        results[key]={
            "source_GenBank_true_original_ITS_ETS_tip_set":TIPS,
            "source_tree_sha256":meta["sha256"],
            "unrooted_source_four_tip_D_D_clade_split_exists":unrooted_split_exists(
                tree,REQUIRED_LCA_TIPS),
            "observed_four_tip_clade_members":sorted(
                x.replace("_"," ") for x in REQUIRED_LCA_TIPS),
            "displayed_root_LCA_tip_count_for_D_D":len(lca.get_terminals()),
            "displayed_root_LCA_tip_count_for_D_and_developmental_V":
                len(lca2.get_terminals()),
            "original_source_sister_Parasenecio_D_V_two_tip_split_exists":
                unrooted_split_exists(tree,PAIR_STAGED),
            "experiment_species_closest_GenBank_tip":nearest_tip(tree,SHELTER),
            "direct_flora_D_Ligularia_closest_GenBank_tip":nearest_tip(tree,FLORA),
            "absolute_historical_head_ancestor_identifiable_from_source_root":False,
        }
    # Predeclared original source scenarios have different actual results;
    # never turn the starting distance topology into "replicate" evolution.
    if ([results[x]["unrooted_source_four_tip_D_D_clade_split_exists"] for x in PINNED]
           !=[True,True,False] or
        not all(x["original_source_sister_Parasenecio_D_V_two_tip_split_exists"]
                for x in results.values()) or
        [results[x]["displayed_root_LCA_tip_count_for_D_D"] for x in PINNED]
           !=[4,4,6] or
        results["BioNJ_initial_distance_tree"][
            "experiment_species_closest_GenBank_tip"]["species"]==FLORA):
        raise ValueError("Original source 2025 gene tree reconstruction sensitivity not reproducible")
    return {
        "schema":"azami_T3_source_authentic_linked_rDNA_131_direct_D9_V1_topology_modes_v1",
        "status":"EMPIRICAL_ORIGINAL_T3_NUCLEAR_RDNA_TOPOLOGY_ECOLOGY_SENSITIVITY",
        "original_real_sequence_phylogeny_taxa":TIPS,
        "source_genbank_ITS_and_ETS_linked_not_independent_coalescent_loci":True,
        "documented_published_direct_species_head_down":species_source["new_authority_D"],
        "documented_direct_D_genera":species_source["new_D_across_genera"],
        "documented_developmental_V":1,
        "unreviewed_true_head_states":121,
        "botanically_erect_species_U":0,
        "shelter_experiment_in_direct_D_species":SHELTER,
        "previously_unknown_English_FOC_but_documented_Chinese_FRPS_D_species":FLORA,
        "all_real_tree_topologies_from_exact_same_ITS_ETS_alignment_not_independent_replicates":True,
        "UFBoot_consensus_shares_same_IQTREE_source_alignment":True,
        "BioNJ_is_initial_same_alignment_distance_start_not_independent_evolution":True,
        "unrooted_between_genus_D_D_clade_split_in_ML_and_consensus_but_not_BioNJ":True,
        "D_V_Parasenecio_two_tip_split_in_all_three_true_rDNA_topologies":True,
        "2021_target_capture_lowcopy_tree_not_observed_here":True,
        "no_independent_erect_U_for_cross_tribe_transition_test":True,
        "n_nodding_origin_events_estimated":False,
        "phenotypic_rain_water_protection_in_Ligularia_experimentally_measured":False,
        "source_original_source_citation_and_genbank_gate":species_source["sources"],
        "source_authentic_tree_sensitivity":results,
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ("ML-tree","consensus-tree","BioNJ-tree","original131",
                "genbank163","chinese-flora","stage-flora","out"):
        p.add_argument("--"+arg,required=True,type=Path)
    a=p.parse_args()
    out=actual_topology(a.ML_tree,a.consensus_tree,a.BioNJ_tree,
                        a.original131,a.genbank163,a.chinese_flora,a.stage_flora)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,ensure_ascii=False,sort_keys=True,
                                indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print("T3_REAL131_TRUE_TOPOLOGY_SENSITIVITY="+json.dumps({
        "D":out["documented_published_direct_species_head_down"],
        "D_V_cherry_all_topologies":out["D_V_Parasenecio_two_tip_split_in_all_three_true_rDNA_topologies"],
        "cross_genus_D_D_four_tip_split_only_ML_and_consensus":
            out["unrooted_between_genus_D_D_clade_split_in_ML_and_consensus_but_not_BioNJ"],
        "ML_closest_to_shelter":out["source_authentic_tree_sensitivity"][
            "max_likelihood"]["experiment_species_closest_GenBank_tip"]["species"],
        "BioNJ_closest_to_shelter":out["source_authentic_tree_sensitivity"][
            "BioNJ_initial_distance_tree"]["experiment_species_closest_GenBank_tip"]["species"]
    },ensure_ascii=False,sort_keys=True),flush=True)


if __name__=="__main__":
    main()
