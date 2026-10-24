"""Fail-closed verification of post-result 110-tip REAL Q20 technical ML.

This is a source-conditioned technical IQ-TREE GTR+G4 phylogram. It preserves
the original 111-tip failure: one original C. tioganum specimen had only
119/29310 observed Q20 positions. The other 110 samples inherit the ORIGINAL
108-locus/29310-site 89/111 mask. This is not ortholog/homeolog-verified,
not a dated/outgroup-rooted Cirsium species tree, and supports neither U->D
event number nor precipitation adaptation, even if compared with NA63.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path

from Bio import Phylo

N_TAXA=110
N_LOCI=108
N_SITES=29310
EXCLUDED="WCVP206740970"
REQUIRED_NODDING={"Cirsium conspicuum","Cirsium scopulorum"}


def read_rows(path):
    with Path(path).open(newline="",encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def verify(tree_file,source_receipt,source_samples,iqtree_report):
    report=json.loads(Path(source_receipt).read_text(encoding="utf-8"))
    samples=read_rows(source_samples)
    if (report.get("posthoc_sensitivity_taxa")!=N_TAXA or
        report.get("original_all111_Q20_successful_real_SRA_BioSamples")!=111 or
        report.get("original111_technical_gate_passed") is not False or
        report.get("original_111_source_108_retained_loci")!=N_LOCI or
        report.get("observed_Q20_reference_positions_reused_unchanged")!=N_SITES or
        report.get("original_111_predeclared_89_of111_quorum_unmodified")!=89 or
        report.get("original_111_reference_position_selection_recomputed_on110") is not False or
        report.get("original_S1_WCVP_master_accepted_keys")!=244 or
        report.get("original_S1_source_name_operational_units")!=251 or
        report.get("source111_single_species_Q20_outlier",{}).get("tip")!=EXCLUDED or
        report.get("source111_single_species_Q20_outlier",{}).get("n_real_Q20_positions")!=119 or
        report.get("copy_number_paralogy_polyploidy_review_completed") is not False or
        report.get("global_EAzami294_true_species_tree_accepted") is not False or
        report.get("true_Mk_nodding_transition_or_rainfall_adaptation_supported") is not False or
        report.get("technical_GTRG4_ML_source_sensitivity_permitted") is not True):
        raise ValueError("Original 111 failed gate or strict posthoc110 real Q20 source mask changed")
    if len(samples)!=N_TAXA:
        raise ValueError("Posthoc110 source sample inventory incomplete")
    ids=[r["source_tip"] for r in samples]
    bios=[r["source_BioSample"] for r in samples]
    runs=[r["original_run"] for r in samples]
    if (len(set(ids))!=N_TAXA or EXCLUDED in ids or
        len(set(bios))!=N_TAXA or len(set(runs))!=N_TAXA or
        any(not tip.startswith("WCVP") for tip in ids) or
        any(int(r["n_actual_Q20_called"])<1000 for r in samples) or
        not REQUIRED_NODDING<={r["taxon_name"] for r in samples}):
        raise ValueError("Source BioSample, origin direction or retained original Q20 sites changed")
    tree=Phylo.read(str(tree_file),"newick")
    tips=tree.get_terminals()
    tip_names=[t.name for t in tips]
    if len(tip_names)!=N_TAXA or set(tip_names)!=set(ids):
        raise ValueError("Inferred technical IQ-TREE tip IDs do not equal authentic 110 BioSamples")
    lengths=[clade.branch_length for clade in tree.find_clades()
             if clade.branch_length is not None]
    if len(lengths)<N_TAXA-1 or any(v<0 for v in lengths) or not any(v>0 for v in lengths):
        raise ValueError("Observed technical tree missing valid substitution-distance branches")
    model=Path(iqtree_report).read_text(encoding="utf-8",errors="replace")
    if "GTR" not in model:
        raise ValueError("A substituted or absent IQ-TREE GTR model report")
    return {
        "schema":"azami_T2_real111_to_posthoc110_observedQ20_IQTREE_GTRG4_v1",
        "status":"ACTUAL_110_SOURCE_Q20_POSTHOC_TECHNICAL_ML_PHYLOGRAM_NOT_ACCEPTED_CIRSIUM_SPECIES_TREE",
        "original_real_source_independent_BioSamples":111,
        "posthoc_real_source_BioSamples":N_TAXA,
        "source_WCVP_keys_retained_elsewhere":244,
        "source_original_S1_name_units_retained_elsewhere":251,
        "single_excluded_lowQ20_specimen":EXCLUDED,
        "excluded_original_Q20_sites":119,
        "original_locus_selection_108_no_posthoc_refit":True,
        "original_111_site_selection_89_of111_no_posthoc_refit":True,
        "original_observed_Q20_nucleotide_sites":N_SITES,
        "recovered_technical_ML_tree_tip_count":len(tips),
        "ML_tree_branch_length_count":len(lengths),
        "ML_tree_positive_branch_count":sum(x>0 for x in lengths),
        "treefile_sha256":hashlib.sha256(Path(tree_file).read_bytes()).hexdigest(),
        "IQTREE_report_sha256":hashlib.sha256(Path(iqtree_report).read_bytes()).hexdigest(),
        "original_two_author_scored_nodding_taxa_kept":sorted(REQUIRED_NODDING),
        "unrooted_substitutions_per_site_not_clock":True,
        "taxon_concept_and_voucher_cross_study_verified":False,
        "per_locus_single_copy_ortholog_and_homeolog_verified":False,
        "NA63_author_independent_tree_is_separate_test":True,
        "no_reclassification_of_unsampled_Cirsium_as_erect":True,
        "no_evolutionary_transition_count_from_this_technical_tree":True,
        "precipitation_driven_adaptive_radiation_tested":False,
        "EAzami294_full_nuclear_species_tree_finished":False,
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("tree","original-receipt","source-samples","iqtree-report","output"):
        p.add_argument("--"+name,type=Path,required=True)
    a=p.parse_args()
    result=verify(a.tree,a.original_receipt,a.source_samples,a.iqtree_report)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,sort_keys=True,ensure_ascii=False,indent=2)+"\n",
                        encoding="utf-8")
    print("T2_REAL110_ORIGINAL111_POSTHOC_GENUINE_Q20_IQTREE="+json.dumps({
        k:result[k] for k in (
            "original_real_source_independent_BioSamples",
            "posthoc_real_source_BioSamples","original_observed_Q20_nucleotide_sites",
            "ML_tree_positive_branch_count","per_locus_single_copy_ortholog_and_homeolog_verified"
        )},sort_keys=True),flush=True)


if __name__=="__main__":
    main()
