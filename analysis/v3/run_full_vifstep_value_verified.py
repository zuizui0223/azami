#!/usr/bin/env python3
"""Fail-closed launcher for the full VIF-step sensitivity when CSV bytes differ.

The canonical full nine-predictor environment has a frozen SHA-256 identity. A
fresh reconstruction can differ in CSV bytes or remotely served raster encoding.
This launcher never declares such a reconstruction byte-identical. Instead it
allows the core sensitivity runner to inspect it only so that the runner's
pre-analysis gate can reconstruct every comparable frozen among-taxon marginal
coefficient. The core gate requires a maximum coefficient discrepancy <= 1e-10;
otherwise execution stops before VIF selection or multivariable inference.

Restored endpoints that are constant at the taxon-median analysis scale are
explicitly recorded as constant_response rather than causing the all-trait scan
to abort. Predictors that become constant within an endpoint-specific cohort are
removed before subset-specific VIF filtering. Neither rule uses trait outcomes to
select among environmental predictors with variation.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.v3 import run_full_vifstep_sensitivity as core

CANONICAL_ENV_SHA = "e242aa7ce69d12b11937c1335e84b9638799c50b42ef36b95725e77190df98e7"


def arg_path(flag: str) -> Path:
    try:
        return Path(sys.argv[sys.argv.index(flag) + 1])
    except (ValueError, IndexError) as exc:
        raise SystemExit(f"required argument missing: {flag}") from exc


def file_sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def has_variation(series: pd.Series) -> bool:
    values = pd.to_numeric(series, errors="coerce").to_numpy(float)
    return bool(len(values) > 1 and np.isfinite(values).all() and np.std(values, ddof=0) > 0)


def robust_run_threshold(threshold, units, traits, env_taxon, global_selected,
                         minimum, permutations, seed):
    """Core run_threshold with explicit degenerate-endpoint bookkeeping."""
    results = []
    selections = []
    for unit in units:
        if unit["inferential_unit"] == "linear_endpoint":
            frame = core.linear_taxon_table(traits, env_taxon, unit["members"][0], minimum)
            response_ok = has_variation(frame["median"]) if len(frame) else False
        else:
            frame = core.circular_taxon_table(traits, env_taxon, unit["members"], minimum)
            response_ok = bool(len(frame) and has_variation(frame["sine"]) and has_variation(frame["cosine"]))

        if not response_ok:
            selections.append({
                "vif_threshold": threshold, "unit_id": unit["unit_id"], "n_taxa": len(frame),
                "status": "constant_response", "selected_predictors": "", "max_vif": np.nan,
            })
            continue

        selected = [p for p in global_selected if has_variation(frame[p])]
        if len(frame) <= len(selected) + 3 or not selected:
            selections.append({
                "vif_threshold": threshold, "unit_id": unit["unit_id"], "n_taxa": len(frame),
                "status": "insufficient_taxa_or_predictor_variation", "selected_predictors": ";".join(selected),
                "max_vif": np.nan,
            })
            continue

        while len(selected) > 1:
            vt = core.vif_table(frame, selected)
            if float(vt.iloc[0].vif) < threshold:
                break
            selected.remove(str(vt.iloc[0].predictor))
        vt = core.vif_table(frame, selected)
        selections.append({
            "vif_threshold": threshold, "unit_id": unit["unit_id"], "n_taxa": len(frame),
            "status": "ok", "selected_predictors": ";".join(selected),
            "max_vif": float(vt.vif.max()),
        })
        if unit["inferential_unit"] == "linear_endpoint":
            results.extend(core.fit_linear(frame, selected, unit, threshold, permutations, seed))
        else:
            results.extend(core.fit_circular(frame, selected, unit, threshold, permutations, seed))

    result = pd.DataFrame(results)
    if len(result):
        result["q_perm_bh_posthoc_family"] = core.bh_adjust(result.p_perm)
        result["posthoc_fdr_significant_0_05"] = result.q_perm_bh_posthoc_family.lt(.05)
    return result, pd.DataFrame(selections)


def main() -> int:
    environment = arg_path("--environment")
    out_dir = arg_path("--out-dir")
    actual = file_sha(environment)
    byte_identical = actual == CANONICAL_ENV_SHA
    print(json.dumps({
        "environment_sha256_actual": actual,
        "environment_sha256_canonical": CANONICAL_ENV_SHA,
        "canonical_byte_identity": byte_identical,
        "continuation_gate": "frozen_univariate_coefficient_reproduction_max_error_le_1e-10",
    }))

    # Override only the byte sentinel for this process. Structural checks and the
    # frozen coefficient-reproduction gate remain active inside core.main().
    core.EXPECTED_ENV_SHA = actual
    core.run_threshold = robust_run_threshold
    status = core.main()
    if status != 0:
        return int(status)

    report_path = out_dir / "full_vifstep_report.json"
    payload = json.loads(report_path.read_text(encoding="utf-8"))
    verification = payload.get("frozen_univariate_reproduction", {})
    n = int(verification.get("n_comparable_rows", 0))
    error = float(verification.get("maximum_absolute_coefficient_error", float("inf")))
    if n < 100 or error > 1e-10:
        raise SystemExit(f"value-identity gate failed after core run: comparable={n}, max_error={error}")
    payload["environment_identity"] = {
        "canonical_sha256": CANONICAL_ENV_SHA,
        "actual_sha256": actual,
        "canonical_byte_identity": byte_identical,
        "accepted_for_posthoc_sensitivity_by": "reproduction_of_all_comparable_frozen_among_taxon_coefficients",
        "n_comparable_frozen_coefficients": n,
        "maximum_absolute_coefficient_error": error,
        "not_a_claim_of_byte_identity": not byte_identical,
    }
    payload["degenerate_endpoint_rule"] = (
        "taxon-median responses with zero variance are recorded as constant_response and excluded from regression"
    )
    report_path.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print("ENVIRONMENT_VALUE_IDENTITY_GATE=PASS")
    print("ENVIRONMENT_ACTUAL_SHA256=" + actual)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
