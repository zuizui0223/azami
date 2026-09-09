#!/usr/bin/env python3
"""Run the frozen v2 full-27 atlas with a supplied subset of frozen predictors.

Only the predictor test family is filtered.  The statistical implementation is
imported directly from analysis/run_geb_v2_full27_environment_atlas.py and its
within/among estimators, permutation counts, endpoint handling and BH procedure
are unchanged.
"""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from analysis import run_geb_v2_full27_environment_atlas as atlas  # noqa: E402

FROZEN_NINE = [
    "chelsa_bio01",
    "chelsa_bio04",
    "chelsa_bio12",
    "chelsa_bio15",
    "chelsa_rsds_mean",
    "chelsa_vpd_mean",
    "chelsa_sfcwind_mean",
    "chelsa_gsp",
    "chelsa_npp",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--endpoint-contract", required=True, type=Path)
    parser.add_argument("--analysis-contract", required=True, type=Path)
    parser.add_argument("--traits-long", required=True, type=Path)
    parser.add_argument("--environment", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--predictors", required=True, nargs="+")
    parser.add_argument("--lane", required=True, choices=["within_taxon", "among_taxon"])
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    selected = [name for name in FROZEN_NINE if name in args.predictors]
    if selected != args.predictors:
        raise SystemExit(
            f"Predictors must be an ordered subset of the frozen nine: requested={args.predictors}, normalized={selected}"
        )
    if not selected:
        raise SystemExit("At least one predictor is required")

    contract = json.loads(args.analysis_contract.read_text(encoding="utf-8"))
    contract = copy.deepcopy(contract)
    contract["analysis_id"] = f"geb_v2_full27_native_vifstep10_{args.lane}_v1"
    contract["status"] = "user_requested_vifstep10_comparative_reanalysis"
    contract["purpose"] = (
        "Repeat the frozen v2 full-27 scale-specific atlas after a phenotype-blind VIFstep<=10 filter; "
        "the downstream statistical model is unchanged."
    )
    contract["environment"]["predictors"] = selected
    filtered_blocks = []
    for block in contract["environment"]["blocks"]:
        keep = [name for name in block["predictors"] if name in selected]
        if keep:
            updated = dict(block)
            updated["predictors"] = keep
            filtered_blocks.append(updated)
    contract["environment"]["blocks"] = filtered_blocks
    contract["multiplicity"]["families"] = (
        "one BH family across all successful trait units and the VIFstep<=10-retained predictors, "
        "separately for within taxon, among min5 and among min2"
    )
    contract["claim_boundary"] = list(contract.get("claim_boundary", [])) + [
        "This comparison uses VIFstep as a predictor-family filter; the frozen v2 primary atlas did not use VIF for model selection.",
        "VIF selection is scale-specific and phenotype-blind; no trait outcome enters predictor retention.",
    ]

    args.out_dir.mkdir(parents=True, exist_ok=True)
    filtered_contract = args.out_dir / "selected_predictor_analysis_contract.json"
    filtered_contract.write_text(json.dumps(contract, indent=2) + "\n", encoding="utf-8")

    atlas.EXPECTED_PREDICTORS = selected
    saved_argv = sys.argv
    try:
        sys.argv = [
            "run_geb_v2_full27_environment_atlas.py",
            "--contract", str(args.endpoint_contract),
            "--analysis-contract", str(filtered_contract),
            "--traits-long", str(args.traits_long),
            "--environment", str(args.environment),
            "--out-dir", str(args.out_dir),
        ]
        result = atlas.main()
    finally:
        sys.argv = saved_argv

    metadata = {
        "lane": args.lane,
        "selected_predictors": selected,
        "n_selected_predictors": len(selected),
        "implementation": "direct import of analysis/run_geb_v2_full27_environment_atlas.py with only EXPECTED_PREDICTORS and matching contract subset changed",
    }
    (args.out_dir / "selected_predictor_runner_report.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(metadata, indent=2))
    return int(result or 0)


if __name__ == "__main__":
    raise SystemExit(main())
