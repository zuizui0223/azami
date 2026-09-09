#!/usr/bin/env python3
"""Fail-closed launcher for the full VIF-step sensitivity when CSV bytes differ.

The canonical full nine-predictor environment has a frozen SHA-256 identity. A
fresh reconstruction can differ in CSV bytes or remotely served raster encoding.
This launcher never declares such a reconstruction byte-identical. Instead it
allows the core sensitivity runner to inspect it only so that the runner's
pre-analysis gate can reconstruct every comparable frozen among-taxon marginal
coefficient. The core gate requires a maximum coefficient discrepancy <= 1e-10;
otherwise execution stops before VIF selection or multivariable inference.

Because standardized coefficients and VIF are invariant to positive affine unit
rescaling, this gate is suitable for diagnosing transport/scale differences but
is not a replacement for the canonical archived input. The resulting analysis
remains a post-hoc sensitivity and does not modify frozen v2 outputs.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

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

    # The core runner otherwise rejects before reaching its stronger scientific
    # value gate. Override only the byte sentinel for this process; all structural
    # checks and the frozen coefficient reproduction remain active.
    core.EXPECTED_ENV_SHA = actual
    status = core.main()
    if status != 0:
        return int(status)

    report_path = out_dir / "full_vifstep_report.json"
    payload = json.loads(report_path.read_text(encoding="utf-8"))
    verification = payload.get("frozen_univariate_reproduction", {})
    n = int(verification.get("n_comparable_rows", 0))
    error = float(verification.get("maximum_absolute_coefficient_error", float("inf")))
    if n < 100 or error > 1e-10:
        raise SystemExit(
            f"value-identity gate failed after core run: comparable={n}, max_error={error}"
        )
    payload["environment_identity"] = {
        "canonical_sha256": CANONICAL_ENV_SHA,
        "actual_sha256": actual,
        "canonical_byte_identity": byte_identical,
        "accepted_for_posthoc_sensitivity_by": "reproduction_of_all_comparable_frozen_among_taxon_coefficients",
        "n_comparable_frozen_coefficients": n,
        "maximum_absolute_coefficient_error": error,
        "not_a_claim_of_byte_identity": not byte_identical,
    }
    report_path.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print("ENVIRONMENT_VALUE_IDENTITY_GATE=PASS")
    print("ENVIRONMENT_ACTUAL_SHA256=" + actual)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
