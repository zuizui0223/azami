#!/usr/bin/env python3
"""Summarize direct among-vs-within integration contrast from the frozen common-cohort bootstrap.

This is a secondary v3 synthesis. It does not refit the frozen v2 endpoint atlas,
change any multiplicity family, or redefine the two manuscript headline candidates.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--pairwise", type=Path, required=True)
    p.add_argument("--bootstrap", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    return p.parse_args()


def q(series: pd.Series, p: float) -> float:
    return float(series.quantile(p))


def main() -> int:
    args = parse_args()
    pairs = pd.read_csv(args.pairwise)
    boot = pd.read_csv(args.bootstrap)

    required_pairs = {
        "within_taxon_rv", "among_taxon_rv", "delta_among_minus_within"
    }
    required_boot = {
        "median_within_rv", "median_among_rv", "relations_stronger_among"
    }
    if not required_pairs.issubset(pairs.columns):
        raise SystemExit(f"pairwise columns missing: {required_pairs - set(pairs.columns)}")
    if not required_boot.issubset(boot.columns):
        raise SystemExit(f"bootstrap columns missing: {required_boot - set(boot.columns)}")
    if len(pairs) != 36:
        raise SystemExit(f"expected 36 construct relations, got {len(pairs)}")
    if len(boot) < 100:
        raise SystemExit(f"bootstrap too small: {len(boot)}")

    boot = boot.copy()
    boot["median_rv_difference_among_minus_within"] = (
        boot["median_among_rv"] - boot["median_within_rv"]
    )
    delta = boot["median_rv_difference_among_minus_within"]
    stronger = boot["relations_stronger_among"]

    report = {
        "analysis_id": "ch1_v3_construct_scale_contrast_summary_20260910",
        "claim_boundary": (
            "secondary direct scale contrast from the existing complete-18 taxon bootstrap; "
            "frozen v2 endpoint conclusions, multiplicity families and headline candidates unchanged"
        ),
        "observed": {
            "relations": int(len(pairs)),
            "median_within_rv": float(pairs["within_taxon_rv"].median()),
            "median_among_rv": float(pairs["among_taxon_rv"].median()),
            "difference_of_medians_among_minus_within": float(
                pairs["among_taxon_rv"].median() - pairs["within_taxon_rv"].median()
            ),
            "median_pairwise_delta_among_minus_within": float(
                pairs["delta_among_minus_within"].median()
            ),
            "relations_stronger_among": int(
                (pairs["among_taxon_rv"] > pairs["within_taxon_rv"]).sum()
            ),
        },
        "taxon_bootstrap": {
            "replicates": int(len(boot)),
            "difference_of_median_rv_median": float(delta.median()),
            "difference_of_median_rv_low95": q(delta, 0.025),
            "difference_of_median_rv_high95": q(delta, 0.975),
            "probability_median_among_exceeds_within": float((delta > 0).mean()),
            "relations_stronger_among_median": float(stronger.median()),
            "relations_stronger_among_low95": q(stronger, 0.025),
            "relations_stronger_among_high95": q(stronger, 0.975),
            "probability_majority_of_relations_stronger_among": float(
                (stronger > (len(pairs) / 2)).mean()
            ),
        },
        "interpretation": (
            "This quantifies the scale contrast already visible in the common-cohort matrices. "
            "It supports stronger visible-phenotype integration among taxa if the bootstrap "
            "difference remains positive, but does not identify evolutionary, developmental, "
            "genetic or causal mechanisms for that difference."
        ),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
