"""Compare source-identical real110 public nuclear GTR+G4 ML and Q20 NJ.

Both methods used exactly the same original111 fixed 108 candidate nuclear
loci, original >=89-of-111 29310 Q20 positions and 110 source BioSamples.
One source-referenced low-coverage C. tioganum was excluded explicitly only
from the posthoc110 sensitivity; original111 failure remains in the record.

An unrooted RF split comparison is purely TECHNICAL estimator sensitivity:
these 108 public reference-copy candidates have unverified paralogy/homeology
and phylogenetic reticulation. Neither tree is accepted as an ancestral
Cirsium species tree nor grounds BM/OU or nodding transition events.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from Bio import Phylo

TIPS=110
ORIGINAL=111
NUCLEAR_LOCI=108
Q20_POS=29310
INDEPENDENT_NA63_OVERLAP=26


def split_set(tree):
    tips=[z.name for z in tree.get_terminals()]
    if len(tips)!=len(set(tips)) or not all(tips):
        raise ValueError("Unrecognized duplicated/unnamed original source tips")
    names=set(tips)
    splits=set()
    for node in tree.get_nonterminals():
        part={x.name for x in node.get_terminals()}
        if 2<=len(part)<=len(names)-2:
            a=tuple(sorted(part))
            b=tuple(sorted(names-part))
            splits.add(min(a,b))
    return names,splits


def source_checks(ml_tree,nj_tree,original_receipt,ml_na63,nj_na63,ml_verified,nj_verified):
    data=json.loads(Path(original_receipt).read_text(encoding="utf-8"))
    q_ml=json.loads(Path(ml_verified).read_text(encoding="utf-8"))
    q_nj=json.loads(Path(nj_verified).read_text(encoding="utf-8"))
    na_ml=json.loads(Path(ml_na63).read_text(encoding="utf-8"))
    na_nj=json.loads(Path(nj_na63).read_text(encoding="utf-8"))
    if (data.get("posthoc_sensitivity_taxa")!=TIPS or
        data.get("original_all111_Q20_successful_real_SRA_BioSamples")!=ORIGINAL or
        data.get("original111_technical_gate_passed") is not False or
        data.get("original_111_source_108_retained_loci")!=NUCLEAR_LOCI or
        data.get("observed_Q20_reference_positions_reused_unchanged")!=Q20_POS or
        data.get("original_111_predeclared_89_of111_quorum_unmodified")!=89 or
        data.get("original_111_reference_position_selection_recomputed_on110") is not False or
        data.get("original_S1_WCVP_master_accepted_keys")!=244 or
        data.get("original_S1_source_name_operational_units")!=251 or
        data.get("source111_single_species_Q20_outlier",{}).get("tip")!="WCVP206740970"):
        raise ValueError("Original111 failed Q20 locus/site mask modified")
    if (q_ml.get("recovered_technical_ML_tree_tip_count")!=TIPS or
        q_ml.get("original_observed_Q20_nucleotide_sites")!=Q20_POS or
        q_ml.get("precipitation_driven_adaptive_radiation_tested") is not False or
        q_nj.get("source_posthoc_Q20_BioSamples")!=TIPS or
        q_nj.get("original_Q20_reference_positions")!=Q20_POS or
        q_nj.get("original_failed111_quality_gate_relabelled_as_success") is not False or
        not q_nj.get("quartets_are_not_independent_replicates") or
        not q_ml.get("unrooted_substitutions_per_site_not_clock")):
        raise ValueError("Source-aligned real Q20 ML and NJ method receipt missing")
    if q_ml["treefile_sha256"]!=hashlib.sha256(Path(ml_tree).read_bytes()).hexdigest():
        raise ValueError("IQ-TREE tree SHA differs from verified genuine real110")
    if (na_ml.get("exact_original_taxon_overlap")!=INDEPENDENT_NA63_OVERLAP or
        na_nj.get("exact_original_taxon_overlap")!=INDEPENDENT_NA63_OVERLAP):
        raise ValueError("Published independent NA63 comparison scope drift")
    if (na_ml.get("real111_technical_tree_is_not_single_copy_orthology_species_tree") is not True or
        na_nj.get("real111_technical_tree_is_not_single_copy_orthology_species_tree") is not True):
        raise ValueError("Unwarranted acceptance of unvetted technical phylogeny")
    mt=Phylo.read(str(ml_tree),"newick")
    nt=Phylo.read(str(nj_tree),"newick")
    names_m,splits_m=split_set(mt)
    names_n,splits_n=split_set(nt)
    if len(names_m)!=TIPS or names_m!=names_n or "WCVP206740970" in names_m:
        raise ValueError("Actual ML/NJ source BioSample/WCVP tip universes differ")
    if len(splits_m)<TIPS//2 or len(splits_n)<TIPS//2:
        raise ValueError("Insufficient observed internal splits for real110 topology screen")
    shared=splits_m&splits_n
    union=splits_m^splits_n
    return {
        "schema":"azami_T2_real110_originalQ20_same_DNA_ML_NJ_unrooted_sensitivity_v1",
        "status":"SOURCE_AUTHENTIC_110_ML_VS_NJ_UNROOTED_ROBUSTNESS_NOT_EVOLUTIONARY_SPECIES_TREE",
        "original111_quality_source_passed":False,
        "one_sparse_C_tioganum_excluded_posthoc":True,
        "original_WCVP_photo_keys":244,
        "original_source_name_units":251,
        "posthoc_public_BioSamples":TIPS,
        "same_original_108_Q20_reference_loci":True,
        "same_original_29310_site_89of111_mask":True,
        "ML_unrooted_internal_splits":len(splits_m),
        "NJ_unrooted_internal_splits":len(splits_n),
        "ML_and_NJ_shared_unrooted_internal_splits":len(shared),
        "ML_vs_NJ_symmetric_Robinson_Foulds":len(union),
        "ML_vs_NJ_normalized_unrooted_RF":len(union)/(len(splits_m)+len(splits_n)),
        "independent_published_NA63_exact_source_taxa":INDEPENDENT_NA63_OVERLAP,
        "ML_vs_NA63_quartet_fraction":na_ml["resolved_quartet_concordance_fraction"],
        "NJ_vs_NA63_quartet_fraction":na_nj["resolved_quartet_concordance_fraction"],
        "ML_vs_NA63_normalized_RF":na_ml["normalized_unrooted_RF"],
        "NJ_vs_NA63_normalized_RF":na_nj["normalized_unrooted_RF"],
        "NA63_quartets_heavily_nonindependent":True,
        "both_technical_not_true_orthology_or_homeolog_certified":True,
        "full_294_tip_EAzami_nuclear_tree_completed":False,
        "botanical_U_to_D_change_number_or_rain_causality_inferred":False,
        "source_NJ_tree_SHA256":hashlib.sha256(Path(nj_tree).read_bytes()).hexdigest(),
        "source_ML_tree_SHA256":q_ml["treefile_sha256"],
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("ml-tree","nj-tree","original-receipt","ml-na63","nj-na63","ml-verified","nj-verified","output"):
        p.add_argument("--"+name,required=True,type=Path)
    a=p.parse_args()
    report=source_checks(a.ml_tree,a.nj_tree,a.original_receipt,a.ml_na63,
                         a.nj_na63,a.ml_verified,a.nj_verified)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(report,sort_keys=True,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    print("T2_SOURCE110_GTRG4_ML_VS_NJ_IDENTICAL_29310_Q20="+json.dumps({
        k:report[k] for k in (
            "ML_unrooted_internal_splits","NJ_unrooted_internal_splits",
            "ML_and_NJ_shared_unrooted_internal_splits",
            "ML_vs_NJ_normalized_unrooted_RF",
            "ML_vs_NA63_quartet_fraction","NJ_vs_NA63_quartet_fraction"
        )},sort_keys=True),flush=True)


if __name__=="__main__":
    main()
