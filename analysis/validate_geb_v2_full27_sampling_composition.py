#!/usr/bin/env python3
"""Fail-closed validation for the full-27 sampling-composition sensitivity."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


KEYS = ["scale", "unit_id", "predictor"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--atlas-within", required=True, type=Path)
    parser.add_argument("--atlas-among", required=True, type=Path)
    parser.add_argument("--sampling-dir", required=True, type=Path)
    parser.add_argument("--spatial-within", required=True, type=Path)
    parser.add_argument("--spatial-among", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    return parser.parse_args()


def truth(values: pd.Series) -> pd.Series:
    return values.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def pair_set(frame: pd.DataFrame) -> set[tuple[str, str, str]]:
    return set(frame[KEYS].astype(str).itertuples(index=False, name=None))


def add_check(
    checks: list[dict[str, object]], name: str, condition: bool, detail: object
) -> None:
    checks.append({"check": name, "passed": bool(condition), "detail": detail})


def recompute_direction_stability(frame: pd.DataFrame) -> pd.Series:
    result = pd.Series(False, index=frame.index, dtype=bool)
    ok = frame["status"].eq("ok")
    linear = ok & frame["inferential_unit"].eq("linear_endpoint")
    baseline = pd.to_numeric(frame["baseline_beta_std"], errors="coerce")
    estimate = pd.to_numeric(frame["beta_std_sensitivity"], errors="coerce")
    result.loc[linear] = (
        baseline.loc[linear].ne(0)
        & estimate.loc[linear].ne(0)
        & np.sign(baseline.loc[linear]).eq(np.sign(estimate.loc[linear]))
    )
    circular = ok & ~frame["inferential_unit"].eq("linear_endpoint")
    alignment = pd.to_numeric(frame["circular_vector_alignment_cosine"], errors="coerce")
    result.loc[circular] = alignment.loc[circular].gt(0)
    return result


def recompute_class(part: pd.DataFrame) -> str:
    all_evaluable = True
    all_stable = True
    any_unstable = False
    for _, family in part.groupby("sensitivity_family", sort=True):
        ok = family["status"].eq("ok")
        stable = truth(family.loc[ok, "direction_stable"])
        all_evaluable = all_evaluable and bool(ok.all())
        all_stable = all_stable and bool(len(stable) > 0 and stable.all())
        any_unstable = any_unstable or bool((~stable).any())
    if any_unstable:
        return "direction_unstable"
    if all_evaluable and all_stable:
        return "stable_all_declared_scenarios"
    return "stable_where_evaluable_incomplete"


def main() -> int:
    args = parse_args()
    within = pd.read_csv(args.atlas_within, low_memory=False)
    among = pd.read_csv(args.atlas_among, low_memory=False)
    scenarios = pd.read_csv(
        args.sampling_dir / "v2_full27_sampling_composition_scenarios.csv",
        low_memory=False,
    )
    summary = pd.read_csv(
        args.sampling_dir / "v2_full27_sampling_composition_summary.csv",
        low_memory=False,
    )
    sampling_report = json.loads(
        (args.sampling_dir / "v2_full27_sampling_composition_report.json").read_text(
            encoding="utf-8"
        )
    )
    spatial_within = pd.read_csv(args.spatial_within, low_memory=False)
    spatial_among = pd.read_csv(args.spatial_among, low_memory=False)
    checks: list[dict[str, object]] = []

    expected_within = within[
        within["status"].eq("ok")
        & pd.to_numeric(within["q_fdr_bh_global_family"], errors="coerce").lt(0.05)
    ].assign(scale="within_taxon")
    expected_among = among[
        among["scope"].eq("among_taxon_min5")
        & among["status"].eq("ok")
        & pd.to_numeric(among["q_fdr_bh_global_family"], errors="coerce").lt(0.05)
    ].assign(scale="among_taxon")
    expected_pairs = pair_set(pd.concat([expected_within, expected_among], ignore_index=True))
    add_check(
        checks,
        "entry_pair_set_matches_scale_specific_global_q_gate",
        pair_set(scenarios) == expected_pairs and pair_set(summary) == expected_pairs,
        {
            "expected": len(expected_pairs),
            "scenario_pairs": len(pair_set(scenarios)),
            "summary_pairs": len(pair_set(summary)),
        },
    )

    n_top = len(sampling_report["top_taxa"])
    n_regions = len(sampling_report["broad_regions"])
    expected_scenarios = {
        "within_taxon": 1 + n_top + 1 + n_regions + 1,
        "among_taxon": n_top + 1 + n_regions + 1,
    }
    counts = scenarios.groupby(KEYS).size()
    scenario_count_ok = all(
        int(count) == expected_scenarios[str(key[0])]
        for key, count in counts.items()
    )
    add_check(
        checks,
        "complete_declared_scenario_grid",
        scenario_count_ok and len(counts) == len(expected_pairs),
        {
            "expected_per_scale": expected_scenarios,
            "observed_min": int(counts.min()),
            "observed_max": int(counts.max()),
            "rows": len(scenarios),
        },
    )
    add_check(
        checks,
        "no_post_selection_p_or_q_values",
        not any(
            column in scenarios.columns
            for column in ("p_value", "q_value", "q_fdr_bh_global_family")
        ),
        sorted(scenarios.columns),
    )
    recomputed_direction = recompute_direction_stability(scenarios)
    observed_direction = truth(scenarios["direction_stable"])
    add_check(
        checks,
        "direction_rule_recomputed",
        recomputed_direction.equals(observed_direction),
        {
            "rows": len(scenarios),
            "unstable": int((~observed_direction).sum()),
            "mismatches": int((recomputed_direction != observed_direction).sum()),
        },
    )

    recomputed_classes = {
        key: recompute_class(part)
        for key, part in scenarios.groupby(KEYS, sort=True)
    }
    class_mismatches = 0
    for row in summary.itertuples(index=False):
        key = (str(row.scale), str(row.unit_id), str(row.predictor))
        class_mismatches += int(
            recomputed_classes.get(key) != str(row.sampling_composition_stability_class)
        )
    add_check(
        checks,
        "pair_summary_classes_recomputed",
        class_mismatches == 0,
        {"rows": len(summary), "mismatches": class_mismatches},
    )
    add_check(
        checks,
        "all_scenarios_executed_or_explicitly_unevaluable",
        scenarios["status"].astype(str).str.len().gt(0).all(),
        scenarios["status"].value_counts().to_dict(),
    )

    spatial_within_pass = spatial_within[
        truth(spatial_within["broad_spatial_sensitivity_pass"])
    ].assign(scale="within_taxon")
    spatial_among_pass = spatial_among[
        truth(spatial_among["broad_spatial_sensitivity_pass"])
    ].assign(scale="among_taxon")
    spatial_pass = pd.concat([spatial_within_pass, spatial_among_pass], ignore_index=True)[
        KEYS
    ].merge(summary, on=KEYS, how="left", validate="one_to_one")
    add_check(
        checks,
        "every_broad_space_pass_has_sampling_composition_result",
        spatial_pass["sampling_composition_stability_class"].notna().all(),
        {
            "broad_space_passes": len(spatial_pass),
            "stable_all": int(
                spatial_pass["sampling_composition_stability_class"]
                .eq("stable_all_declared_scenarios")
                .sum()
            ),
            "direction_unstable": int(
                spatial_pass["sampling_composition_stability_class"]
                .eq("direction_unstable")
                .sum()
            ),
        },
    )

    passed = all(item["passed"] for item in checks)
    output = {
        "status": "PASS" if passed else "FAIL",
        "n_checks": len(checks),
        "n_passed": sum(bool(item["passed"]) for item in checks),
        "checks": checks,
        "scientific_summary": {
            "pair_stability_classes": summary[
                "sampling_composition_stability_class"
            ].value_counts().to_dict(),
            "broad_space_pass_sampling_classes": spatial_pass[
                "sampling_composition_stability_class"
            ].value_counts().to_dict(),
        },
        "claim_ceiling": (
            "directional stability in a retrospective selected-row sensitivity; "
            "not correction to global representativeness and not a new discovery family"
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
