"""Outcome-blind native-source attrition audit, before ecological fitting.

Inputs are source membership/exposures and terminal endpoint processing states,
never phenotype values. Missing ledger rows mean UNREPORTED, not failed images.
No overlap statistic or retention threshold is an ecological authorization.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .workflow import digest

EXPOSURES = ("pr_month", "rsds_month", "vpd_month", "tasmax_month", "sfcWind_month")
STATES = {
    "not_scheduled", "rights_restricted", "metadata_unresolved", "scheduled",
    "download_failed", "downloaded", "no_detection", "detected",
    "measurement_failed", "qc_unusable", "usable",
}
MILESTONES = {
    "attempted": {"download_failed", "downloaded", "no_detection", "detected", "measurement_failed", "qc_unusable", "usable"},
    "downloaded": {"downloaded", "no_detection", "detected", "measurement_failed", "qc_unusable", "usable"},
    "detected": {"detected", "measurement_failed", "qc_unusable", "usable"},
    "usable": {"usable"},
}


def _identity(frame: pd.DataFrame, keys: list[str], label: str) -> None:
    if not set(keys) <= set(frame):
        raise ValueError(f"{label}: required identity columns absent")
    if frame[keys].isna().any().any():
        raise ValueError(f"{label}: missing identity")
    if frame[keys].astype(str).apply(lambda col: col.str.strip().eq("")).any().any():
        raise ValueError(f"{label}: blank identity")
    if frame[keys].astype(str).duplicated().any():
        raise ValueError(f"{label}: duplicate identity")


def audit(source: pd.DataFrame, states: pd.DataFrame, modules: dict[str, list[str]]) -> dict:
    """Each state is observation x endpoint, aggregated from all known photos.

    A producer must not write ``no_detection`` while another linked photo remains
    unprocessed. This reader checks identity and states, not image-level history.
    Weights stay fixed at 1/source observations per taxon throughout attrition.
    """
    _identity(source, ["obs_id"], "source")
    if source.empty or not {"accepted_key", "native_range_status"} <= set(source):
        raise ValueError("Nonempty native source with accepted taxon keys required")
    if source["accepted_key"].isna().any() or source["accepted_key"].astype(str).str.strip().eq("").any():
        raise ValueError("Missing accepted taxon key")
    if not source["native_range_status"].eq("native").all():
        raise ValueError("Audit baseline must be the eligible native ecological source")
    if not modules or any(not name or not endpoints or len(set(endpoints)) != len(endpoints)
                          for name, endpoints in modules.items()):
        raise ValueError("Nonempty modules with unique endpoint lists required")
    endpoints = sorted({e for values in modules.values() for e in values})
    _identity(states, ["obs_id", "endpoint"], "states")
    if "state" not in states or not states["state"].isin(STATES).all():
        raise ValueError("Missing or unrecognized processing state")
    if not states["endpoint"].isin(endpoints).all():
        raise ValueError("Unknown endpoint in state ledger")
    x = source.copy()
    x["obs_id"] = x["obs_id"].astype(str)
    x["accepted_key"] = x["accepted_key"].astype(str)
    states = states.copy()
    states["obs_id"] = states["obs_id"].astype(str)
    if not states["obs_id"].isin(x["obs_id"]).all():
        raise ValueError("State ledger contains observations outside the frozen source")
    x = x.set_index("obs_id")
    weights = 1 / x.groupby("accepted_key")["accepted_key"].transform("size")
    state_matrix = states.pivot(index="obs_id", columns="endpoint", values="state")
    state_matrix = state_matrix.reindex(index=x.index, columns=endpoints).fillna("unreported")
    denominator = len(x)
    taxon_mass = float(weights.sum())

    def retention(mask: pd.Series) -> dict:
        return {"observations": int(mask.sum()), "fraction_of_source": float(mask.mean()),
                "source_taxon_standardized_fraction": float(weights[mask].sum() / taxon_mass),
                "represented_taxa": int(x.loc[mask, "accepted_key"].nunique())}

    endpoint_reports = {}
    for endpoint in endpoints:
        col = state_matrix[endpoint]
        endpoint_reports[endpoint] = {
            "terminal_states": {str(k): int(v) for k, v in col.value_counts().sort_index().items()},
            "milestones": {name: retention(col.isin(eligible)) for name, eligible in MILESTONES.items()},
        }
    joint_masks = {name: state_matrix[values].eq("usable").all(axis=1) for name, values in modules.items()}
    environmental = {}
    for variable in EXPOSURES:
        if variable not in x:
            environmental[variable] = {"status": "not_evaluable_exposure_absent"}
            continue
        values = pd.to_numeric(x[variable], errors="coerce")
        finite = pd.Series(np.isfinite(values.to_numpy(dtype=float)), index=x.index)
        if not finite.any():
            environmental[variable] = {"status": "not_evaluable_no_finite_source_exposure", "missing": denominator}
            continue
        # Source-defined bins; collapse tied boundaries without consulting outcomes.
        edges = np.unique(np.quantile(values[finite], [0, .2, .4, .6, .8, 1]))
        if len(edges) == 1:
            codes = pd.Series(0, index=values[finite].index)
        else:
            codes = pd.cut(values[finite], bins=edges, labels=False, include_lowest=True)
        strata = []
        masks = [(f"source_bin_{int(code)}", x.index.isin(codes.index[codes.eq(code)]))
                 for code in sorted(codes.unique())]
        masks.append(("missing_exposure", ~finite.to_numpy()))
        for label, array in masks:
            mask = pd.Series(array, index=x.index)
            count = int(mask.sum())
            mass = float(weights[mask].sum())
            strata.append({
                "stratum": label, "source_observations": count,
                "module_retention": {name: {
                    "usable_observations": int((mask & usable).sum()),
                    "fraction": float((mask & usable).sum() / count) if count else None,
                    "source_taxon_standardized_fraction": float(weights[mask & usable].sum() / mass) if mass else None,
                } for name, usable in joint_masks.items()},
            })
        environmental[variable] = {"status": "source_binned_retention_reported", "source_bin_edges": edges.tolist(), "strata": strata}
    unreported = int(state_matrix.eq("unreported").to_numpy().sum())
    return {
        "schema_version": 1,
        "status": "COVERAGE_AUDIT_WITH_UNREPORTED_STATES" if unreported else "COVERAGE_AUDIT_REPORTED",
        "source_observations": denominator, "source_taxa": int(x["accepted_key"].nunique()),
        "expected_observation_endpoint_slots": denominator * len(endpoints),
        "unreported_slots": unreported,
        "endpoints": endpoint_reports,
        "joint_modules": {name: retention(mask) for name, mask in joint_masks.items()},
        "environmental_support": environmental,
        "trait_values_read": 0, "ecological_models_executed": 0,
        "ecological_fitting_authorized": False,
        "limits": [
            "Stages summarize a supplied observation-level processing ledger; photo-level completeness requires its separate provenance check.",
            "Missing state rows are unreported, not unavailable images, detector negatives or QC failures.",
            "Taxon-standardized fractions use fixed source weights, not weights recomputed among survivors.",
            "Marginal source-bin retention is not a multivariate hypervolume or proof of unbiased sampling.",
            "Selection relative to recorded native source support does not identify unrecorded plants or recording effort.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--states", type=Path, required=True)
    parser.add_argument("--modules", type=Path, required=True, help="JSON mapping module names to endpoint IDs")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise ValueError("Preserve prior audit; output already exists")
    allowed = {"obs_id", "accepted_key", "native_range_status", *EXPOSURES}
    source = pd.read_csv(args.source, usecols=lambda c: c in allowed,
                         dtype={"obs_id": str, "accepted_key": str})
    states = pd.read_csv(args.states, usecols=["obs_id", "endpoint", "state"], dtype=str)
    report = audit(source, states, json.loads(args.modules.read_text(encoding="utf-8")))
    report["input_sha256"] = {"source": digest(args.source), "states": digest(args.states), "modules": digest(args.modules)}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({k: report[k] for k in ("status", "source_observations", "source_taxa", "unreported_slots")}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
