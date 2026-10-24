"""Real molecular sister contrast: two Parasenecio original source DNA samples.

2025 new 131 GenBank ITS+ETS source tree contains a strongly supported
P. cyclotus – P. koualapensis unrooted cherry (actual SHaLRT/UFBoot).
Both species independently occur in 2021 Ren PRJNA672812 capture SRA
source sample roster. Independent FOC species descriptions:
P. cyclotus: capitula pendulous, stage not reported;
P. koualapensis: peduncles spreading, pendulous LATER.
This does NOT say P. cyclotus never moves; both may have analogous
ontogenetic sequences. Time-of-anthesis of later-down posture unknown.
2025 ITS+ETS is linked multicopy rDNA, neither a gene-wise independent
two-locus species tree nor a rainfall-adaptation/independent-gain count.
The 2021 genuine low-copy capture tips are known but its readable tree
has not been recovered; no 2021 resolved topology is inferred here.
"""
from __future__ import annotations
import argparse
import csv
import json
import hashlib
from pathlib import Path
from collections import Counter
import re

from analysis.v4.t3_2025_real_nuclear_rDNA_phylogeny_posture_topology import (
    read_source, TREE_SHA, N_GENUINE_GB_ITS_ETS
)

D_SPECIES="Parasenecio cyclotus"
V_SPECIES="Parasenecio koualapensis"
SRA={"Parasenecio cyclotus":("SRR12967040","SAMN16588858"),
     "Parasenecio koualapensis":("SRR12967039","SAMN16588859")}
SOURCE_GENBANK={"Parasenecio cyclotus":("KY979028","KY978948"),
                "Parasenecio koualapensis":("KY979045","KY978965")}


def load_csv(path):
    with Path(path).open(newline="",encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def audit(tree_file,original131_crosswalk,original48_sra,stage_source,GenBank163):
    tree,original_states=read_source(tree_file,original131_crosswalk)
    if Counter(original_states.values())!={"D":8,"?":123}:
        raise ValueError("True original 131 genotype/morphology D8 + unscored123 changed")
    staged=load_csv(stage_source)
    if len(staged)!=1:
        raise ValueError("Only one independently checked ontogenetic source")
    v=staged[0]
    if (v["2025_original_species"]!=V_SPECIES or
        v["source_posture"]!="V" or
        v["developmental_scope"]!="early_spreading_later_pendulous" or
        "pendulous later" not in v["capitulum_bearing_peduncle_source_phrase"] or
        not v["original_direct_species_flora"].startswith("https://efloras.org/") or
        any(v[x]!="no" for x in (
            "DNA_voucher_equals_morphological_flora_specimen",
            "is_true_fixed_flowering_erect_control",
            "ancestral_state_or_rain_shelter_fitness_proven"))):
        raise ValueError("Direct SRA species head-stage label not uniquely sourced as variable")
    if original_states.get(V_SPECIES)!="?" or original_states.get(D_SPECIES)!="D":
        raise ValueError("Staged V cannot overwrite fixed genuine D or known U in source genealogy")
    original2021=load_csv(original48_sra)
    if len(original2021)!=48:
        raise ValueError("Original Ren2021 true target capture must retain all 48 independent real runs")
    genes=load_csv(GenBank163)
    if len(genes)!=163:
        raise ValueError("Original Taiwania 2025 public 163 GenBank accession rows missing")
    anchor={}
    for name,(run,biosample) in SRA.items():
        c=[x for x in original2021 if x["original_species_name_exact_from_ENA"]==name]
        if (len(c)!=1 or c[0]["original_nuclear_SRA_run"]!=run or
            c[0]["original_nuclear_BioSample"]!=biosample or
            c[0]["DNA_sequence_reads_or_original_phylo_Newick_downloaded"]!="False"):
            raise ValueError("The original 2021 low-copy DNA BioSample species identity differs")
        g=[x for x in genes if x["original_published_scientific_name"]==name
           and x["ITS_GenBank_accession"]==SOURCE_GENBANK[name][0]
           and x["ETS_GenBank_accession"]==SOURCE_GENBANK[name][1]]
        if len(g)!=1 or g[0]["ITS_and_ETS_independent_unlinked_gene_loci"]!="no":
            raise ValueError("Original Taiwania 2025 SAME linked nuclear rDNA ITS-ETS source incorrectly relabeled")
        anchor[name]={"Ren2021_original_target_capture_run":run,
                      "Ren2021_original_BioSample":biosample,
                      "Taiwania2025_nuclear_ITS_accession":g[0]["ITS_GenBank_accession"],
                      "Taiwania2025_nuclear_ETS_accession":g[0]["ETS_GenBank_accession"],
                      "original_taxa_matched_at_species_name_not_same_voucher":True}
    n=[x.name for x in tree.get_terminals()]
    if len(n)!=N_GENUINE_GB_ITS_ETS:
        raise ValueError("Original 2025 novel ITS-ETS rDNA nuclear genealogy tips missing")
    a=D_SPECIES.replace(" ","_")
    b=V_SPECIES.replace(" ","_")
    node=tree.common_ancestor(a,b)
    descendants={x.name for x in node.get_terminals()}
    if descendants!={a,b}:
        raise ValueError("Original documented fixed D–developmental V are NOT a two-species nuclear tip cherry")
    s=re.fullmatch(r"(\d+(?:\.\d+)?)/(\d+(?:\.\d+)?)",node.name or "")
    if not s:
        raise ValueError("Original direct ITS/ETS ML internal sister node has no support label")
    alrt,uf=map(float,s.groups())
    branch=float(tree.distance(a,b))
    if not 0<branch<.02:
        raise ValueError("Authentic original GenBank nuclear rDNA sister cherry branch scale changed")
    staged_states={**original_states,V_SPECIES:"V"}
    if Counter(staged_states.values())!={"D":8,"V":1,"?":122}:
        raise ValueError("Developmental floral posture was silently made fixed or missing species changed")
    genus={k.split()[0] for k,v in staged_states.items() if v=="D"}
    # 131 other nuclear genealogy tips are real, but are not automatically erect.
    return {
        "schema":"azami_T3_Ren2021_Taiwania2025_real_nuclear_Pcyclotus_Pkoualapensis_posture_developmental_sister_v1",
        "status":"EMPIRICAL_TRUE_2025_RDNA_SISTER_AND_2021_LOWCOPY_SOURCE_BIOSAMPLES_DIRECT_STAGE_EVIDENCE",
        "real_2025_ITS_ETS_molecular_original_taxa":N_GENUINE_GB_ITS_ETS,
        "reconciled_flowering_morphology_states":{"D":8,"V":1,"unknown":122,"U":0},
        "same_species_in_2021_public_lowcopy_DNA_and_2025_linked_nuclear_rDNA":anchor,
        "DNA_original_ITS_ETS_SRA_voucher_pairs_are_not_identical_phenotyped_individuals":True,
        "documented_D_species":D_SPECIES,
        "documented_variable_peduncle_posture_species":V_SPECIES,
        "original_V_documentation":"peduncles spreading, pendulous later; anthesis stage not explicit",
        "P_cyclotus_documented_head_pendulous_ontogeny_unmeasured":True,
        "gene_tree_cherry_pair_phylogenetic_distance_substitutions_per_site":branch,
        "actual_original_IQTREE_sister_node_SHaLRT":alrt,
        "actual_original_IQTREE_sister_node_UFBoot":uf,
        "original_genealogy_is_rooted_2021_coalescent_species_tree":False,
        "n_true_D_genera_in_131_linked_rDNA":len(genus),
        "no_true_erect_species_tips_described_in_131":True,
        "unscored_source_taxa_are_NOT_erect":True,
        "no_inferred_U_to_D_independent_gains":True,
        "no_rainfall_fitness_or_historical_precipitation_inferred":True,
        "technical_same_sibling_pair_2025_only_not_confirmed_in_Ren2021_tree":True,
        "next_discriminating_biological_measurements":[
            "Observe same capitulum P cyclotus and P koualapensis through bud-anthesis-fruiting",
            "Measure gravitational head inclination and peduncle attachment under controlled viewpoints",
            "Test exposure under natural/artificial rainfall on matched flower stages",
            "Check 2021 multi-gene capture clade topology and voucher identity independently",
        ],
        "source":{
            "original_2025_131_Linked_nuclear_rDNA_SHA256":TREE_SHA,
            "original_2021_48_SRA_source_SHA256":hashlib.sha256(Path(original48_sra).read_bytes()).hexdigest(),
            "original_2025_163_genbank_accession_table_sha256":hashlib.sha256(Path(GenBank163).read_bytes()).hexdigest(),
            "original_direct_species_flora_URL":v["original_direct_species_flora"],
        },
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tree",type=Path,required=True)
    p.add_argument("--source131",type=Path,required=True)
    p.add_argument("--source2021",type=Path,required=True)
    p.add_argument("--GenBank163",type=Path,required=True)
    p.add_argument("--stage-source",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    result=audit(a.tree,a.source131,a.source2021,a.stage_source,a.GenBank163)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True,
                                   allow_nan=False)+"\n",encoding="utf-8")
    print("T3_REAL_NUCLEAR_SISTER_DOCUMENTED_DEVELOPMENT="+json.dumps({
        "sister_pair":[D_SPECIES,V_SPECIES],
        "SHaLRT":result["actual_original_IQTREE_sister_node_SHaLRT"],
        "UFBoot":result["actual_original_IQTREE_sister_node_UFBoot"],
        "branch_length":result["gene_tree_cherry_pair_phylogenetic_distance_substitutions_per_site"],
        "morphology_states":result["reconciled_flowering_morphology_states"],
    },sort_keys=True,ensure_ascii=False),flush=True)


if __name__=="__main__":main()
