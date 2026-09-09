"""Verify the minimal v3 upstream patch and v2 reuse path.

No traits are fit and no images are downloaded.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RECEIPT = ROOT / "reproducibility" / "v3_minimal_reuse_receipt_20260909.json"
CONTRACT = ROOT / "analysis" / "v3" / "workflow_contract.json"

V2_REUSE = [
    "analysis/run_geb_v2_full27_environment_atlas.py",
    "analysis/run_geb_v2_full27_spatial_sensitivity.py",
    "analysis/run_geb_v2_full27_sampling_composition_sensitivity.py",
    "analysis/audit_v2_environment_collinearity.py",
]
V3_PATCH = [
    "analysis/v3/reconcile_sources.py",
    "analysis/v3/audit_processing_history.py",
    "analysis/v3/recover_display_composition.py",
]


def validate() -> dict:
    receipt = json.loads(RECEIPT.read_text(encoding="utf-8"))
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    if receipt["status"] != "MINIMAL_V3_UPSTREAM_PATCH_VERIFIED_FROM_PRIOR_EXECUTIONS":
        raise ValueError("Minimal reuse receipt is not verified")
    if contract["canonical_flow"] != ["source_provenance", "historical_measurement_reuse", "v2_ecology", "output"]:
        raise ValueError("Minimal canonical flow changed")
    if receipt["source_reconciliation"]["retained_unique_observations"] != 665139:
        raise ValueError("Recovered source observation count changed")
    if receipt["historical_processing_reuse"]["detected_heads"] != 1255791:
        raise ValueError("Historical head denominator changed")
    v2 = receipt["frozen_v2_primary_view"]
    if (v2["observations"], v2["taxa"]) != (46276, 259):
        raise ValueError("Frozen v2 primary cohort changed")
    patch = receipt["historical_measurement_patch"]
    if patch["source_heads"] != 1255791 or patch["recovered_observation_values_per_field"] != 347608:
        raise ValueError("Recovered historical measurement counts changed")
    if receipt["canonical_reuse"]["new_ecological_model"] or receipt["canonical_reuse"]["new_primary_cohort"] or receipt["canonical_reuse"]["new_image_remeasurement"]:
        raise ValueError("Minimal revision cannot promote a new analysis universe")
    missing = [path for path in [*V2_REUSE, *V3_PATCH] if not (ROOT / path).is_file()]
    if missing:
        raise ValueError("Required reused file missing: " + ", ".join(missing))
    return {
        "status": "MINIMAL_V3_TO_V2_REUSE_PATH_READY",
        "canonical_flow": contract["canonical_flow"],
        "source_observations_traced": 665139,
        "historical_detected_heads_reused": 1255791,
        "v2_primary_observations": 46276,
        "v2_primary_taxa": 259,
        "recovered_omitted_fields": 5,
        "new_image_remeasurement_required": False,
        "new_ecological_model_required": False,
        "next_step": "reuse frozen v2 inputs/code; apply only claim corrections and optional endpoint-structure QC",
    }


if __name__ == "__main__":
    print(json.dumps(validate(), indent=2))
