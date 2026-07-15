#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Iterable


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def as_float(value: str | None) -> float:
    if value is None or value == "":
        return math.nan
    return float(value)


def bh_adjust(values: Iterable[float]) -> list[float]:
    values = list(values)
    finite = [(i, value) for i, value in enumerate(values) if math.isfinite(value)]
    result = [math.nan] * len(values)
    if not finite:
        return result
    ordered = sorted(finite, key=lambda item: item[1])
    m = len(ordered)
    running = 1.0
    for rank_from_end, (index, value) in enumerate(reversed(ordered), start=1):
        rank = m - rank_from_end + 1
        adjusted = min(running, value * m / rank, 1.0)
        result[index] = adjusted
        running = adjusted
    return result


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--molecular-dir", type=Path, required=True)
    parser.add_argument("--signal-dir", type=Path, required=True)
    args = parser.parse_args()

    coverage = read_csv(args.molecular_dir / "cirsium_molecular_database_coverage.csv")
    deterministic = read_csv(args.signal_dir / "phylogenetic_signal_deterministic_and_direct.csv")
    random_rows = read_csv(args.signal_dir / "phylogenetic_signal_random_tree_summary.csv")

    coverage_columns = [
        ("nuccore_all", "Any nucleotide record"),
        ("nuccore_its", "ITS-related nucleotide record"),
        ("nuccore_plastid", "Plastid/chloroplast record"),
        ("nuccore_complete_plastome", "Complete plastome"),
        ("sra_all", "SRA record"),
    ]
    n_taxa = len(coverage)
    coverage_summary: list[dict[str, object]] = []
    for column, label in coverage_columns:
        counts = [int(float(row[column])) for row in coverage]
        present = sum(value > 0 for value in counts)
        coverage_summary.append({
            "database_category": column,
            "label": label,
            "n_taxa_with_records": present,
            "n_taxa_total": n_taxa,
            "coverage_percent": round(100 * present / n_taxa, 3),
            "total_records_returned": sum(counts),
            "median_records_per_taxon": sorted(counts)[n_taxa // 2] if n_taxa else 0,
            "maximum_records_for_one_taxon": max(counts, default=0),
        })
    write_csv(
        args.molecular_dir / "cirsium_molecular_database_coverage_summary.csv",
        coverage_summary,
        list(coverage_summary[0]),
    )

    s1 = [row for row in deterministic if row["scenario"] == "S1"]
    groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in s1:
        groups[row["information_tier"]].append(row)
    for rows in groups.values():
        k_q = bh_adjust(as_float(row["K_p_value"]) for row in rows)
        lambda_q = bh_adjust(as_float(row["lambda_p_value"]) for row in rows)
        for row, k_value, lambda_value in zip(rows, k_q, lambda_q):
            row["K_q_value_BH"] = "" if not math.isfinite(k_value) else f"{k_value:.12g}"
            row["lambda_q_value_BH"] = "" if not math.isfinite(lambda_value) else f"{lambda_value:.12g}"

    random_by_endpoint = {row["endpoint"]: row for row in random_rows}
    s1_by_key = {(row["endpoint"], row["information_tier"]): row for row in s1}
    evidence_rows: list[dict[str, object]] = []
    for endpoint, random_row in sorted(random_by_endpoint.items()):
        all_tip = s1_by_key[(endpoint, "all_tips")]
        direct = s1_by_key[(endpoint, "direct_backbone_only")]
        all_k_q = as_float(all_tip.get("K_q_value_BH"))
        all_lambda_q = as_float(all_tip.get("lambda_q_value_BH"))
        direct_k_q = as_float(direct.get("K_q_value_BH"))
        direct_lambda_q = as_float(direct.get("lambda_q_value_BH"))
        all_supported = (math.isfinite(all_k_q) and all_k_q < 0.05) or (
            math.isfinite(all_lambda_q) and all_lambda_q < 0.05
        )
        direct_supported = (math.isfinite(direct_k_q) and direct_k_q < 0.05) or (
            math.isfinite(direct_lambda_q) and direct_lambda_q < 0.05
        )
        if direct_supported:
            tier = "direct_backbone_supported"
            manuscript_role = "candidate_main_or_supplement"
        elif all_supported:
            tier = "grafting_sensitive"
            manuscript_role = "supplement_only"
        else:
            tier = "no_FDR_robust_signal"
            manuscript_role = "supplement_only"
        if endpoint.startswith("corolla_hue_"):
            endpoint_note = "Hue sine/cosine component; interpret jointly, not as an independent biological trait."
            tier = "requires_joint_circular_test"
            manuscript_role = "supplement_only"
        else:
            endpoint_note = ""
        evidence_rows.append({
            "analysis_tier": random_row["analysis_tier"],
            "trait_group": random_row["trait_group"],
            "endpoint": endpoint,
            "interpretation": random_row["interpretation"],
            "n_species_all_tips": all_tip["n_species"],
            "n_species_direct_backbone": direct["n_species"],
            "S1_all_K": all_tip["K"],
            "S1_all_K_p": all_tip["K_p_value"],
            "S1_all_K_q_BH": all_tip.get("K_q_value_BH", ""),
            "S1_all_lambda": all_tip["lambda"],
            "S1_all_lambda_p": all_tip["lambda_p_value"],
            "S1_all_lambda_q_BH": all_tip.get("lambda_q_value_BH", ""),
            "direct_K": direct["K"],
            "direct_K_p": direct["K_p_value"],
            "direct_K_q_BH": direct.get("K_q_value_BH", ""),
            "direct_lambda": direct["lambda"],
            "direct_lambda_p": direct["lambda_p_value"],
            "direct_lambda_q_BH": direct.get("lambda_q_value_BH", ""),
            "random_tree_K_median": random_row["K_median"],
            "random_tree_K_q025": random_row["K_q025"],
            "random_tree_K_q975": random_row["K_q975"],
            "random_tree_lambda_median": random_row["lambda_median"],
            "random_tree_lambda_q025": random_row["lambda_q025"],
            "random_tree_lambda_q975": random_row["lambda_q975"],
            "evidence_tier": tier,
            "recommended_manuscript_role": manuscript_role,
            "endpoint_note": endpoint_note,
        })
    write_csv(
        args.signal_dir / "phylogenetic_signal_evidence_summary.csv",
        evidence_rows,
        list(evidence_rows[0]),
    )

    tier_counts: dict[str, int] = defaultdict(int)
    for row in evidence_rows:
        tier_counts[str(row["evidence_tier"])] += 1
    decision = {
        "n_endpoints": len(evidence_rows),
        "evidence_tier_counts": dict(sorted(tier_counts.items())),
        "primary_decision_rule": (
            "A signal is considered directly supported only when Blomberg K or Pagel lambda "
            "remains significant after BH correction in the 54 direct-backbone tips."
        ),
        "manuscript_recommendation": (
            "Keep phylogenetic signal as a sensitivity analysis unless at least one non-circular "
            "endpoint is direct_backbone_supported. All-tip-only results are grafting-sensitive."
        ),
        "limitations": [
            "The backbone is not a resolved nuclear Cirsium species tree.",
            "Random grafting does not model reticulation, chloroplast capture, or allopolyploid parentage.",
            "Hue sine and cosine are components of one circular trait and require joint interpretation.",
        ],
    }
    (args.signal_dir / "phylogenetic_signal_decision.json").write_text(
        json.dumps(decision, indent=2), encoding="utf-8"
    )
    print(json.dumps({"coverage": coverage_summary, "decision": decision}, indent=2))


if __name__ == "__main__":
    main()
