"""T3: independently published Chinese Flora adds a true head D to REAL 131-tip tree.

The OWN GenBank nuclear linked ITS+ETS gene tree already has 131 authentic
Senecioneae taxa. Original T3 author-flora D8, ?123. The native-language
Flora Reipublicae Popularis Sinicae ORIGINAL SPECIES text for *Ligularia
purdomii* explicitly says "头状花序多数，盘状，下垂" (capitula numerous,
discoid, NODDING). ENGLISH Flora of China describes erect STEM ONLY, not
an erect head. It is therefore an INDEPENDENT direct D species entry,
not a taxon genus projection, photo-camera-frame posture, or tree tip ASR.

In the authentic 131-tip ITS+ETS gene tree the nearest tip to the
experimentally sheltered *Cremanthodium campanulatum* is L. purdomii. Yet
their immediate LCA has FOUR tips including TWO unscored Ligularia
and is not a cherry; this unrooted nuclear-rDNA relationship is NOT a
two-species independent E->D gain, nor verified evolutionary homology.
L. purdomii has a genuine original ITS AY723257 but no ETS. The
2021 48-sample Ren low-copy target-capture set did not sample it.

No DNA voucher has been matched to the particular flowering head described
in Chinese Flora; no species tree/homoeolog correction or ancestral rain
selective fitness proof is derived.
"""
from __future__ import annotations
import argparse
import csv
from collections import Counter
import hashlib
import json
from pathlib import Path
import re

from analysis.v4.t3_2025_real_nuclear_rDNA_phylogeny_posture_topology import (
    read_source, TREE_SHA, N_GENUINE_GB_ITS_ETS
)

FLORA="Ligularia purdomii"
SHELTER="Cremanthodium campanulatum"
ITS_ACC="AY723257"
REQUIRED_LCA_TIPS={
    "Cremanthodium_campanulatum",
    "Ligularia_yunnanensis",
    "Ligularia_curvisquama",
    "Ligularia_purdomii",
}
EXPECTED_FOC_HEAD="头状花序多数，盘状，下垂"
FRPS_CHINESE_URL="https://www.iflora.cn/search_result.aspx?p=AA797647D57B0EB4C1FABEBC584460E8"
FOC_ENGLISH_URL="https://www.efloras.org/florataxon.aspx?flora_id=2&taxon_id=200024274"


def csv_rows(path):
    with Path(path).open(newline="",encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def reconcile(treefile,original131crosswalk,florafile,genbank163):
    tree,states=read_source(treefile,original131crosswalk)
    if (len(states)!=N_GENUINE_GB_ITS_ETS or
        Counter(states.values())!={"D":8,"?":123}):
        raise ValueError("Previously independently verified 131-tip D8 morphology changed")
    flora=csv_rows(florafile)
    if len(flora)!=1:
        raise ValueError("Single independently traceable Chinese Flora source expected")
    fl=flora[0]
    if (
        fl["taxon"]!=FLORA or fl["head_posture"]!="D" or
        fl["chinese_species_head_phrase"]!=EXPECTED_FOC_HEAD or
        fl["species_description_source_URL"]!=FRPS_CHINESE_URL or
        fl["source_plant_stage"]!="stage_unspecified" or
        fl["2025_original_ITS_accession"]!=ITS_ACC or
        fl["2025_original_ETS_accession"] or
        fl["does_english_FOC_assert_nodding"]!="no" or
        fl["english_FOC_describes_erect_stem_not_head"]!="yes" or
        any(fl[n]!="no" for n in (
            "DNA_and_flora_same_flowering_specimen_verified",
            "double_blind_head_review_completed",
            "independent_U_to_D_origin_or_rain_utility_demonstrated"
        ))
    ):
        raise ValueError("Original Chinese flora nodding CAPITULUM or scientific stops not satisfied")
    archived=csv_rows(genbank163)
    if len(archived)!=163:
        raise ValueError("Actual Taiwania 2025 original 163 public GenBank accession rows changed")
    exact=[x for x in archived if x["original_published_scientific_name"]==FLORA]
    if len(exact)!=1 or exact[0]["ITS_GenBank_accession"]!=ITS_ACC or exact[0]["ETS_GenBank_accession"]:
        raise ValueError("Cannot verify genuine published SAME species GenBank ITS-only accession")
    if states[FLORA]!="?":
        raise ValueError("New independently discovered Chinese flora cannot overwrite already coded head")
    states={**states,FLORA:"D"}
    counts=Counter(states.values())
    if counts!={"D":9,"?":122}:
        raise ValueError("Original 131 D8 source not properly extended by only ONE direct Chinese flora")
    terminal={t.name:t for t in tree.get_terminals()}
    stem1=SHELTER.replace(" ","_")
    stem2=FLORA.replace(" ","_")
    if stem1 not in terminal or stem2 not in terminal:
        raise ValueError("Original real ITS ETS nuclear tree lacks exact protected and Chinese flora tips")
    clade=tree.common_ancestor(stem1,stem2)
    members={z.name for z in clade.get_terminals()}
    if members!=REQUIRED_LCA_TIPS:
        raise ValueError("Original D-D near tip relationship is NOT conserved authentic four-tip clade")
    match=re.fullmatch(r"([0-9]+(?:\.[0-9]+)?)/([0-9]+(?:\.[0-9]+)?)",
                       clade.name or "")
    if not match:
        raise ValueError("Actual 4-tip nuclear gene-tree node source IQ-TREE support unrecorded")
    sh,uf=map(float,match.groups())
    distances=sorted(
        ((float(tree.distance(stem1,x)),x)
         for x in terminal if x!=stem1),
        key=lambda row:(row[0],row[1]),
    )
    if distances[0][1]!=stem2:
        raise ValueError("Genuinely nearest ML-tip to experimentally anchored capitulum changed")
    if not len(members)>2:
        raise ValueError("Source root-invariant clade falsely promoted to a true two-tip cherry")
    desc={
        "schema":"azami_T3_original_FRPS_Ligularia_purdomii_Senecioneae_131_nuclear_linked_rDNA_v1",
        "status":"ORIGINAL_CHINESE_SPECIES_FLORA_CONFIRMS_NODDING_THIRD_GENUS_EXACT_ON_REAL_131_NUCLEAR_TIPS",
        "original_2025_genbank_linked_rDNA_taxa":131,
        "previous_authority_documented_D":8,
        "additional_source_direct_genbank_same_species_D":1,
        "new_authority_D":9,
        "new_D_across_genera":dict(Counter(x.split()[0] for x in states if states[x]=="D")),
        "source_unknown_heads_retained":122,
        "source_verified_U":0,
        "original_FRPS_morphology_taxon":FLORA,
        "original_species_capitula_nodding_Chinese_phrase":EXPECTED_FOC_HEAD,
        "English_flora_narrates_ERECT_STEM_NOT_ERECT_CAPITULUM":True,
        "original_same_species_ITS_only_accession":ITS_ACC,
        "L_purdomii_rDNA_ETS_gene_available_in_2025_source":False,
        "identical_morphological_voucher_and_GenBank_nuclear_tip_verified":False,
        "direct_rain_protection_experiment_species":SHELTER,
        "experiment_species_nearest_nuclear_genealogy_taxon":FLORA,
        "experiment_species_nearest_original_gene_tree_distance_substitutions_per_site":distances[0][0],
        "experiment_species_nearest_real_two_genus_heads_are_both_source_D":True,
        "original_4tip_clade_member_names":[z.replace("_"," ") for z in sorted(members)],
        "clade_IQTREE_SHaLRT":sh,
        "clade_IQTREE_UFBoot":uf,
        "within_real_clade_other_two_members_have_UNKNOWN_morphology":True,
        "experiment_and_new_D_species_are_A_CHERRY_two_tip_sisters":False,
        "unrooted_nuclear_rDNA_genealogy_not_homoeolog_corrected_species_tree":True,
        "one_Chinese_species_source_without_independent_double_blind_flora_review":True,
        "no_imputed_erect_controls_in_senescent_or_unobserved_stage":True,
        "no_gravity_angle_or_stage_from_S1_image_orientation_inferred":True,
        "number_of_independent_down_head_origins_inferred":False,
        "rainfall_at_head_origin_or_effect_on_pollen_fitness_in_L_purdomii_measured":False,
        "sources":{
            "original_FRPS_direct_species_head_URL":FRPS_CHINESE_URL,
            "English_Flora_of_China_stage_unspecified_URL":FOC_ENGLISH_URL,
            "original_Taiwania2025_163_accessions_SHA256":hashlib.sha256(Path(genbank163).read_bytes()).hexdigest(),
            "original_2025_nuclear_131_tree_SHA256":TREE_SHA,
            "original_Chinese_flora_admission_sha256":hashlib.sha256(Path(florafile).read_bytes()).hexdigest(),
        }
    }
    if desc["new_D_across_genera"]!={"Cremanthodium":2,"Parasenecio":6,"Ligularia":1}:
        raise ValueError("An extra source taxon would falsely imply extra independent evolutionary event")
    return desc


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("tree","original131","chinese-flora","genbank163","out"):
        p.add_argument("--"+key,type=Path,required=True)
    a=p.parse_args()
    result=reconcile(a.tree,a.original131,a.chinese_flora,a.genbank163)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,sort_keys=True,ensure_ascii=False,indent=2,allow_nan=False)+"\n",
                     encoding="utf-8")
    print("T3_REAL_FRPS_HEAD_D_NUCLEAR131="+json.dumps({
        "n_D":result["new_authority_D"],
        "genera":result["new_D_across_genera"],
        "near_experimental":result["experiment_species_nearest_nuclear_genealogy_taxon"],
        "LCA_size":len(result["original_4tip_clade_member_names"]),
        "SHaLRT":result["clade_IQTREE_SHaLRT"],
        "UFBoot":result["clade_IQTREE_UFBoot"],
    },sort_keys=True,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
