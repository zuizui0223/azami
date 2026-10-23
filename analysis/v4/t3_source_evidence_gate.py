"""T3 publication-source head posture candidates, not inferred convergence."""
from __future__ import annotations
import argparse
from collections import Counter
import csv
import json
from pathlib import Path

SCHEMA = (
    "record_id", "taxon", "tribe", "source_posture", "life_stage_scope",
    "source_kind", "brief_source_phrase", "source_url", "evidence_basis",
    "taxon_concept_independently_reviewed", "blind_double_annotation_completed",
    "same_tip_rooted_nuclear_tree_verified", "strict_T3_admissible",
    "coding_limit", "taxonomic_status"
)
STATES = frozenset(("U", "D", "?", "V"))

def read_evidence(path):
    with Path(path).open(newline="", encoding="utf-8-sig") as h:
        reader = csv.DictReader(h)
        if tuple(reader.fieldnames or ()) != SCHEMA:
            raise ValueError("T3 source evidence schema drift")
        rows = list(reader)
    if not rows:
        raise ValueError("No species-level posture evidence")
    seen_ids, seen_taxa = set(), set()
    for r in rows:
        if r["record_id"] in seen_ids or r["taxon"] in seen_taxa:
            raise ValueError("Duplicated source species or ID would inflate T3")
        seen_ids.add(r["record_id"])
        seen_taxa.add(r["taxon"])
        if (r["source_posture"] not in STATES
                or r["tribe"] != "Senecioneae"
                or len(r["taxon"].split()) < 2):
            raise ValueError("Unexpected real Senecioneae taxon or head state")
        if not r["source_url"].startswith("https://"):
            raise ValueError("Missing independently traceable source URL")
        for k in (
            "taxon_concept_independently_reviewed",
            "blind_double_annotation_completed",
            "same_tip_rooted_nuclear_tree_verified",
            "strict_T3_admissible",
        ):
            if r[k] not in ("yes", "no"):
                raise ValueError("Invalid independent evidence flag")
        if r["source_posture"] == "V" and r["life_stage_scope"] != "developmentally_changes":
            raise ValueError("Developmental posture switch must not be fixed-tip binary")
        if r["life_stage_scope"] == "age_dependent_tendency" and r["source_posture"] != "?":
            raise ValueError("Age-related tendency cannot become a fixed-tip U/D state")
        if r["taxonomic_status"] not in ("no_conflict_identified", "contested_2018_2022", "contested_2016_2026"):
            raise ValueError("Unreviewed or malformed botanical taxonomy status")
        if r["taxonomic_status"] != "no_conflict_identified" and r["strict_T3_admissible"] == "yes":
            raise ValueError("Contested author taxon concept cannot enter independent T3 phylogeny")
        if r["strict_T3_admissible"] == "yes":
            valid = (
                r["source_posture"] in ("U", "D")
                and r["life_stage_scope"] == "anthesis_independently_verified"
                and r["taxon_concept_independently_reviewed"] == "yes"
                and r["blind_double_annotation_completed"] == "yes"
                and r["same_tip_rooted_nuclear_tree_verified"] == "yes"
            )
            if not valid:
                raise ValueError("Prematurely promoted posture to strict T3 inference")
    return rows

def summarize(rows):
    states = Counter(r["source_posture"] for r in rows)
    genus = {}
    for r in rows:
        genus.setdefault(r["taxon"].split()[0], Counter())[r["source_posture"]] += 1
    strict = [r for r in rows if r["strict_T3_admissible"] == "yes"]
    return {
        "schema": "azami_T3_Senecioneae_public_posture_sources_20261023",
        "source_taxa": len(rows),
        "nodding_species_descriptions_not_independent_gains": states["D"],
        "erect_species_descriptions_not_ancestral_states": states["U"],
        "qualified_uncertain_species": states["?"],
        "ontogenetically_variable_species": states["V"],
        "qualified_age_dependent_nonbinary_species": sum(
            r["life_stage_scope"] == "age_dependent_tendency" for r in rows
        ),
        "source_taxon_concept_contests": sum(
            r["taxonomic_status"] != "no_conflict_identified" for r in rows
        ),
        "source_genus_counts": {k: dict(sorted(v.items())) for k,v in sorted(genus.items())},
        "strict_reviewed_same_tip_phenotypes": len(strict),
        "ready_tribes": sorted({r["tribe"] for r in strict}),
        "genuine_rooted_cross_tribe_nuclear_tree_join_verified": False,
        "historical_repeated_nodding_origins_estimated": False,
        "rain_adaptation_supported_by_this_screen": False,
        "status": "REAL_SPECIES_DESCRIPTIONS_ADDED_T3_EVOLUTIONARY_TEST_BLOCKED",
        "reason": "No jointly independent-reviewed anthesis-only two-tribe nuclear tree phenotype panel"
    }

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, default=Path(__file__).resolve().parent / "results" / "20261023_T3_SENECIONEAE_POSTURE_SOURCE_CANDIDATES.csv")
    p.add_argument("--output", type=Path)
    args = p.parse_args()
    result = summarize(read_evidence(args.source))
    text = json.dumps(result, sort_keys=True, ensure_ascii=False, allow_nan=False, indent=2)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    print(text, end="")

if __name__ == "__main__":
    main()
