"""Audit the ALREADY-EXECUTED native Gate B aggregate, without private rows.

Input: native_integration.json from azami Actions run 39428065666,
artifact 12266812186, archive member SHA256 below. This reinterprets a
previously fitted predictor-complete 157654/174 cohort; it does NOT fit
new species/observer-balanced covariance, KDE, environmental or genetic models.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path

SOURCE_RUN = 39428065666
ARTIFACT_ID = 12266812186
MEMBER_SHA256 = "1545037041b62b2494934d022c7b21e7a7109f54f3582579c6557dfadf998752"
CORRECTED_INPUT_SHA256 = "342d78e9a8a960e8d9bdcf0bc2cdacd13bee48214a8b962419d5a71d5aef46ed"
EXPECTED_N = 157654
EXPECTED_TAXA = 174
CONSTRUCT_GROUP = {
    "angle": "presentation",
    "lightness": "colour",
    "chroma": "colour",
    "hue": "colour",
    "aspect": "outline",
    "circularity": "outline",
    "solidity": "outline",
    "width_profile": "outline",
}
EXPECTED_PAIR_NAMES = {
    frozenset(pair) for pair in itertools.combinations(CONSTRUCT_GROUP, 2)
}


def verify_aggregate_file(path: Path) -> dict:
    contents = path.read_bytes()
    if hashlib.sha256(contents).hexdigest() != MEMBER_SHA256:
        raise ValueError("Wrong native integration aggregate bytes (fail-closed)")
    data = json.loads(contents)
    if str(data["run_id"]) != str(SOURCE_RUN):
        raise ValueError("Wrong archived source workflow")
    if data["input_sha256"] != CORRECTED_INPUT_SHA256:
        raise ValueError("Uncorrected or substituted native environmental input")
    return data


def summarize_existing_aggregate(report: dict) -> dict:
    """Re-group matched constructed trait pairs; NO new ecological estimation."""
    common = report["common_cohort"]
    if common.get("status") != "estimated":
        raise ValueError("Gate B aggregate is not an executed result")
    if (common.get("n"), common.get("species")) != (EXPECTED_N, EXPECTED_TAXA):
        raise ValueError("Gate B cohort changed")
    pairs = common["pairs"]
    if len(pairs) != 28:
        raise ValueError("Unexpected number of registered construct pairs")
    if {frozenset((p["left"], p["right"])) for p in pairs} != EXPECTED_PAIR_NAMES:
        raise ValueError("Missing, duplicated or unregistered construct pairs")

    counts = {}
    notable = []
    def category(a: str, b: str) -> str:
        return "/".join(sorted((CONSTRUCT_GROUP[a], CONSTRUCT_GROUP[b])))
    for p in pairs:
        a, b = p["left"], p["right"]
        group = category(a, b)
        row = counts.setdefault(group, {
            "pairs": 0,
            "raw_positive_delta": 0,
            "raw_pointwise_interval_entirely_positive": 0,
            "adjusted_positive_delta": 0,
            "adjusted_pointwise_interval_entirely_positive": 0,
            "adjusted_pointwise_interval_entirely_negative": 0,
        })
        row["pairs"] += 1
        for mode in ("raw", "adjusted"):
            est = p["estimates"]
            within = float(est[mode + "/within"]["rv"])
            between = float(est[mode + "/between"]["rv"])
            record = est[mode + "/between_minus_within_rv"]
            delta = float(record["estimate"])
            bounds = record["interval"]
            lo, hi = float(bounds["lower"]), float(bounds["upper"])
            if not all(map(math.isfinite, (within, between, delta, lo, hi))):
                raise ValueError("Nonfinite effect or interval")
            if not 0 <= within <= 1 or not 0 <= between <= 1:
                raise ValueError("Invalid RV coefficient")
            if abs((between - within) - delta) > 1e-8:
                raise ValueError("Between minus within identity broken")
            if lo > hi or not -1 <= lo <= 1 or not -1 <= hi <= 1:
                raise ValueError("Invalid bootstrap interval")
            if mode == "raw":
                row["raw_positive_delta"] += int(delta > 0)
                row["raw_pointwise_interval_entirely_positive"] += int(lo > 0)
            else:
                row["adjusted_positive_delta"] += int(delta > 0)
                row["adjusted_pointwise_interval_entirely_positive"] += int(lo > 0)
                row["adjusted_pointwise_interval_entirely_negative"] += int(hi < 0)
                if lo > 0 or hi < 0:
                    notable.append({
                        "pair": a + "/" + b,
                        "module_pair": group,
                        "adjusted_among_minus_within_rv": delta,
                        "pointwise_interval_95": [lo, hi],
                    })
    total = {field: sum(entry[field] for entry in counts.values()) for field in next(iter(counts.values()))}
    if total["pairs"] != 28:
        raise AssertionError("Pair census changed")
    cross = {key: value for key, value in counts.items() if key.split("/")[0] != key.split("/")[1]}
    cross_total = {field: sum(entry[field] for entry in cross.values()) for field in total}
    return {
        "status": "VERIFIED_PREVIOUS_NATIVE_COMMON_COHORT_NOT_NEW_FIT",
        "source_run": SOURCE_RUN,
        "source_artifact_id": ARTIFACT_ID,
        "source_member_sha256": MEMBER_SHA256,
        "corrected_native_input_sha256": CORRECTED_INPUT_SHA256,
        "cohort": {"observations": EXPECTED_N, "accepted_taxa": EXPECTED_TAXA,
                   "selection": "predictor complete and minimum five accepted taxon observations"},
        "constructs": len(CONSTRUCT_GROUP),
        "image_outputs": 9,
        "pointwise_bootstrap_95_intervals_not_multiplicity_adjusted": True,
        "by_module_pair": dict(sorted(counts.items())),
        "all_pairs": total,
        "cross_module_only": cross_total,
        "pointwise_adjusted_interval_nonzero_pairs": sorted(notable, key=lambda x: x["pair"]),
        "already_executed_equal_n_sampling_range": common["equal_n_diagnostic"],
        "limitations": [
            "Original sample contains repeated source observations, not proven independent plants.",
            "Original observer IDs and matched taxon-by-cell-by-observer counts absent.",
            "Selected original 174 taxa are not all 242 taxa in nine-output measurement cohort.",
            "Colour measurements share pixels and hue uses joint sine/cosine.",
            "RV is unsigned and n-sensitive; within/among RV are not genetic integration.",
            "Interval exclusions are pointwise exploratory and cannot imply FDR-controlled support.",
            "Equal-n resampling range holds between means fixed; it is not an ecological confidence interval.",
            "The adjusted fit holds environmental coefficients fixed in bootstrap.",
            "Full observer/region-balanced Gate B on original corrected rows remains unexecuted.",
        ],
    }


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--integration-json", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    args = p.parse_args()
    result = summarize_existing_aggregate(verify_aggregate_file(args.integration_json))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
                        encoding="utf-8")
    print(json.dumps({"status": result["status"], "all_pairs": result["all_pairs"],
                      "cross_module_only": result["cross_module_only"]}))


if __name__ == "__main__":
    main()
