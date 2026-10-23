"""Validate a source-backed, non-destructive S1 posture extension, 66 -> 72.

The 6 added literature head-orientation claims are single-pass candidates only.
Two variable/ontogenetic concepts stay excluded from binary ecology.
This audits coverage and descriptive precipitation means, NOT a phylogenetic
rate effect, independent gains, or rainfall adaptation.
"""
from __future__ import annotations
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
from statistics import mean

BASE = Path(__file__).resolve().parent / "results"
CENSUS = BASE / "20261023_T2_S1_ALL251_66_TO72_PROVISIONAL_POSTURE_OVERLAY.csv"
CANDIDATES = BASE / "20261023_T2_S1_8_EXTERNAL_BOTANICAL_POSTURE_CANDIDATES.csv"
QUEUE = BASE / "20261023_T2_S1_177_UNKNOWN_FLORA_ANNOTATION_PRIORITY.csv"
SOURCE_COUNTS = {"U":54, "D":12, "?":185}
EXTRA_COUNTS = {"U":3, "D":3, "V":2}


def load(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError("Missing original species-level evidence")
    return rows


def checked_summary(all_rows, external, queue):
    if len(all_rows) != 251 or len({r["accepted_key"] for r in all_rows}) != 244:
        raise ValueError("Original 251/244 source labels were lost or relabelled")
    if len({(r["accepted_key"], r["original_source_name"]) for r in all_rows}) != 251:
        raise ValueError("One WCVP source-name evidence record duplicated")
    if sum(int(r["original_photo_count"]) for r in all_rows) != 170066:
        raise ValueError("170066 photographic observations missing")
    if Counter(r["botanical_U_D_state"] for r in all_rows) != SOURCE_COUNTS:
        raise ValueError("Original independent author U/D baseline changed")
    if Counter(r["source_posture"] for r in external) != EXTRA_COUNTS:
        raise ValueError("Only 3 U, 3 D, 2 variable published source records allowed")
    if len({r["original_source_name"] for r in external}) != 8:
        raise ValueError("Duplicated additional botanical literature taxon")
    evidence = {r["original_source_name"]:r for r in external}
    seen = set()
    for record in all_rows:
        taxon = record["original_source_name"]
        original = record["botanical_U_D_state"]
        follow = record["expanded_posture_sensitivity_U_D_or_unknown"]
        new = evidence.get(taxon)
        if (record["frozen_20261023_botanical_U_D_original"] != original or
            record["botanical_state_not_exif_angle_threshold"] != "True" or
            record["nuclear_independent_tree_tip_and_taxon_concept_verified"] != "False" or
            record["newly_evolutionary_gain_or_rain_adaptation_established"] != "False"):
            raise ValueError("Original botanical source/phylogeny stop flags changed")
        if not new:
            if follow != original or record["new_single_source_literature_posture"]:
                raise ValueError("Unreviewed source name was silently assigned U/D")
            if record["developmental_or_taxon_polymorphism_flag"] != "False":
                raise ValueError("Fictitious variable-stage source label")
            continue
        seen.add(taxon)
        if original != "?" or new["source_review_level"] not in (
            "single_authority_direct_capitulum", "single_authority_variable_state",
            "single_authority_developmental_state"):
            raise ValueError("New literature must not replace already scored traits")
        if (new["strict_dual_blind_taxon_concept_checked"] != "no"
            or not new["primary_source_url"].startswith("https://")
            or record["new_single_source_literature_URL"] != new["primary_source_url"]
            or record["new_single_source_literature_posture"] != new["source_posture"]
            or record["source_review_stage"] != "additional_single_authority_not_dual_blinded"):
            raise ValueError("Posthoc single-authority results cannot be called validated")
        expected = new["source_posture"] if new["source_posture"] in ("U","D") else "?"
        if follow != expected or (
            record["developmental_or_taxon_polymorphism_flag"] == "True"
        ) != (new["source_posture"] == "V"):
            raise ValueError("Stage-varying posture cannot be coded as fixed U/D")
    if seen != set(evidence):
        raise ValueError("Published sources not matched to EXACT S1 original species names")
    strict = [r for r in all_rows if r["botanical_U_D_state"] in ("U","D")]
    proposed = [r for r in all_rows if r["expanded_posture_sensitivity_U_D_or_unknown"] in ("U","D")]
    variables = [r for r in all_rows if r["developmental_or_taxon_polymorphism_flag"] == "True"]
    remaining = [r for r in all_rows if
        r["expanded_posture_sensitivity_U_D_or_unknown"] == "?" and
        r["developmental_or_taxon_polymorphism_flag"] == "False"]
    if len(strict) != 66 or len(proposed) != 72 or len(variables) != 2 or len(remaining) != 177:
        raise ValueError("Full251 source cohort conservation violated")
    if set(r["original_source_name"] for r in remaining) != set(r["original_source_name"] for r in queue):
        raise ValueError("177 still-unscored source taxa lost from next-annotation queue")
    for r in queue:
        if not r["mandatory_evidence"].startswith("original_species_flora") or r["current_annotation"] != "unknown":
            raise ValueError("Missing appropriate future taxonomic source review")
    n = lambda r: r["nuclear_real120k_BioSample"] == "True"
    groups = {
        "original_independent_authorities_66": (strict,"botanical_U_D_state"),
        "expanded_single_authority_sensitivity_72": (proposed,"expanded_posture_sensitivity_U_D_or_unknown"),
    }
    reports = {}
    for label,(rows,key) in groups.items():
        down=[r for r in rows if r[key]=="D"]
        up=[r for r in rows if r[key]=="U"]
        d12 = mean(float(r["chelsa_bio12"]) for r in down)
        u12 = mean(float(r["chelsa_bio12"]) for r in up)
        reports[label] = {
            "n_U":len(up),"n_D":len(down),
            "n_original_source_images":sum(int(r["original_photo_count"]) for r in rows),
            "n_original_real111_source_nuclear_samples":sum(n(r) for r in rows),
            "n_nodding_nuclear_source_samples":sum(n(r) for r in down),
            "descriptive_mean_BIO12_D_mm":d12,
            "descriptive_mean_BIO12_U_mm":u12,
            "descriptive_difference_D_minus_U_mm":d12-u12,
            "phylogenetic_independent_tip_effect_fitted":False,
        }
    if (reports["original_independent_authorities_66"]["n_original_real111_source_nuclear_samples"]!=38 or
        reports["expanded_single_authority_sensitivity_72"]["n_original_real111_source_nuclear_samples"]!=42 or
        sum(n(r) for r in all_rows)!=111):
        raise ValueError("Real111 DNA coverage must remain candidate-only and sample-conserved")
    return {
        "status":"SOURCE_REVIEW_EXPANDED_NOT_GLOBAL_PHYLOGENETIC_RAIN_TEST",
        "original_total_taxon_source_name_units":251,
        "original_photographic_observations":170066,
        "original_accepted_WCVP_keys":244,
        "source_authority_cohorts":reports,
        "new_binary_candidates":6,
        "source_documented_variable_or_ontogenetic_taxa":[r["original_source_name"] for r in variables],
        "remaining_original_source_name_units_unscored":177,
        "remaining_unscored_original_photo_observations":sum(int(r["original_photo_count"]) for r in remaining),
        "remaining_unscored_real111_source_nuclear_candidates":sum(n(r) for r in remaining),
        "priority_one_real111_nuclear_with_10photos":sum(r["priority"]=="P1_nuclear_and_photo10plus" for r in queue),
        "all_source_name_units_kept_even_when_unknown":True,
        "new_six_original_S1_photo_EXIF_values_not_used_to_assign_U_D":True,
        "external_authorities_not_double_blind_or_voucher_matched":True,
        "global_independent_molecular_species_tree_not_constructed":True,
        "rain_induced_evolutionary_transition_not_inferred":True,
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--census",type=Path,default=CENSUS)
    p.add_argument("--external",type=Path,default=CANDIDATES)
    p.add_argument("--queue",type=Path,default=QUEUE)
    p.add_argument("--out",type=Path)
    a=p.parse_args()
    summary=checked_summary(load(a.census),load(a.external),load(a.queue))
    contents=json.dumps(summary,ensure_ascii=False,sort_keys=True,indent=2,allow_nan=False)+"\n"
    if a.out:
        a.out.parent.mkdir(parents=True,exist_ok=True)
        a.out.write_text(contents,encoding="utf-8")
    print(contents,end="")


if __name__=="__main__":
    main()
