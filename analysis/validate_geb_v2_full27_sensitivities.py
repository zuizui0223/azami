#!/usr/bin/env python3
"""Fail-closed validation for v2 full-27 spatial and historical gates."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo


KEYS = ["unit_id", "predictor"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--atlas-dir", required=True, type=Path)
    parser.add_argument("--spatial-dir", required=True, type=Path)
    parser.add_argument("--historical-dir", required=True, type=Path)
    parser.add_argument("--tree-dir", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    return parser.parse_args()


def truth(values: pd.Series) -> pd.Series:
    return values.astype(str).str.lower().isin({"true", "1", "yes"})


def pair_set(frame: pd.DataFrame) -> set[tuple[str, str]]:
    return set(frame[KEYS].astype(str).itertuples(index=False, name=None))


def add_check(
    checks: list[dict[str, object]], name: str, condition: bool, detail: object
) -> None:
    checks.append({"check": name, "passed": bool(condition), "detail": detail})


def recompute_spatial_pass(frame: pd.DataFrame) -> pd.Series:
    direction = frame["inferential_unit"].ne("linear_endpoint") | truth(
        frame["same_linear_direction_as_base"]
    )
    return (
        frame["status"].eq("ok")
        & pd.to_numeric(frame["spatial_permutation_p_value"], errors="coerce").lt(0.05)
        & direction
        & pd.to_numeric(frame["residual_morans_p_value"], errors="coerce").ge(0.05)
    )


def count_placement_trees(tree_dir: Path) -> int:
    total = 0
    for path in sorted(tree_dir.glob("*.tre")):
        total += sum(1 for _ in Phylo.parse(str(path), "newick"))
    for path in sorted(tree_dir.glob("*.trees")):
        total += sum(1 for _ in Phylo.parse(str(path), "newick"))
    return total


def main() -> int:
    args = parse_args()
    within = pd.read_csv(args.atlas_dir / "v2_full27_environment_within.csv")
    among = pd.read_csv(args.atlas_dir / "v2_full27_environment_among.csv")
    spatial_within = pd.read_csv(args.spatial_dir / "v2_full27_spatial_within.csv")
    spatial_among = pd.read_csv(args.spatial_dir / "v2_full27_spatial_among.csv")
    models = pd.read_csv(args.historical_dir / "v2_full27_historical_placement_models.csv")
    summary = pd.read_csv(args.historical_dir / "v2_full27_historical_placement_summary.csv")
    checks: list[dict[str, object]] = []

    expected_within = within[truth(within["fdr_significant_0_05"])]
    expected_among = among[
        among["scope"].eq("among_taxon_min5") & truth(among["fdr_significant_0_05"])
    ]
    add_check(
        checks,
        "spatial_within_entry_set_matches_atlas_q_gate",
        pair_set(expected_within) == pair_set(spatial_within),
        {"expected": len(expected_within), "observed": len(spatial_within)},
    )
    add_check(
        checks,
        "spatial_among_entry_set_matches_primary_atlas_q_gate",
        pair_set(expected_among) == pair_set(spatial_among),
        {"expected": len(expected_among), "observed": len(spatial_among)},
    )
    for label, frame in (("within", spatial_within), ("among", spatial_among)):
        expected_pass = recompute_spatial_pass(frame)
        actual_pass = truth(frame["broad_spatial_sensitivity_pass"])
        add_check(
            checks,
            f"spatial_{label}_pass_rule_recomputed",
            expected_pass.equals(actual_pass),
            {
                "rows": len(frame),
                "passes": int(actual_pass.sum()),
                "mismatches": int((expected_pass != actual_pass).sum()),
            },
        )
        add_check(
            checks,
            f"spatial_{label}_all_rows_executed",
            frame["status"].eq("ok").all(),
            frame["status"].value_counts().to_dict(),
        )

    tree_count = count_placement_trees(args.tree_dir)
    add_check(checks, "historical_resource_contains_52_trees", tree_count == 52, tree_count)
    spatial_pass_pairs = pair_set(
        spatial_among[truth(spatial_among["broad_spatial_sensitivity_pass"])]
    )
    add_check(
        checks,
        "historical_pair_set_matches_spatial_among_passes",
        pair_set(summary) == spatial_pass_pairs and pair_set(models) == spatial_pass_pairs,
        {
            "spatial_pairs": len(spatial_pass_pairs),
            "summary_pairs": len(pair_set(summary)),
            "model_pairs": len(pair_set(models)),
        },
    )
    add_check(
        checks,
        "historical_model_grid_complete_and_successful",
        len(models) == len(spatial_pass_pairs) * tree_count and models["status"].eq("ok").all(),
        {"rows": len(models), "status_counts": models["status"].value_counts().to_dict()},
    )

    summary_mismatches = 0
    for row in summary.itertuples(index=False):
        subset = models[
            models["unit_id"].astype(str).eq(str(row.unit_id))
            & models["predictor"].astype(str).eq(str(row.predictor))
        ]
        successful = subset[subset["status"].eq("ok")]
        expected_pass = (
            len(successful) == tree_count
            and pd.to_numeric(successful["p_value"], errors="coerce").lt(0.05).all()
            and truth(successful["same_linear_direction_as_spatial"]).all()
        )
        consistent = (
            int(row.n_successful_placement_trees) == len(successful)
            and int(row.n_placement_trees_p_lt_0_05)
            == int(pd.to_numeric(successful["p_value"], errors="coerce").lt(0.05).sum())
            and (
                str(row.historical_placement_sensitivity_pass).lower()
                in {"true", "1", "yes"}
            )
            == expected_pass
        )
        summary_mismatches += int(not consistent)
    add_check(
        checks,
        "historical_summary_recomputed",
        summary_mismatches == 0,
        {"rows": len(summary), "mismatches": summary_mismatches},
    )
    add_check(
        checks,
        "candidate_labels_only_after_both_gates",
        truth(summary["historical_placement_sensitivity_pass"]).all()
        and summary["candidate_class"].eq(
            "adaptive_pattern_candidate_under_current_controls"
        ).all(),
        summary["candidate_class"].value_counts().to_dict(),
    )

    passed = all(item["passed"] for item in checks)
    output = {
        "status": "PASS" if passed else "FAIL",
        "n_checks": len(checks),
        "n_passed": sum(bool(item["passed"]) for item in checks),
        "checks": checks,
        "claim_ceiling": (
            "adaptive-pattern candidates under broad-space and historical-placement "
            "sensitivities; not proof of adaptation, selection, convergence or mechanism"
        ),
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(output, ensure_ascii=False, indent=2))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
