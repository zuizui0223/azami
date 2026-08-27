#!/usr/bin/env python3
"""Validate the machine-readable Azami Chapter 1 hypothesis-recovery ledger.

This audit protects the final pattern/mechanism boundary. It does not recompute
scientific results; it checks that the submission-facing status resolves all
registered hypotheses, preserves immutable provenance and keeps external
validation gates open rather than silently promoting observational patterns.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STATUS_PATH = ROOT / "analysis" / "ch1" / "hypothesis_recovery_status.json"

EXPECTED_IDS = {
    "FOUNDATION_REPEATED_CONTINUOUS_PHENOTYPES",
    "FOUNDATION_BELOW_TAXON_DIVERSITY",
    "FOUNDATION_ONE_UNIVERSAL_SYNDROME",
    "H1_MEASUREMENT_MODULE_ORGANIZATION",
    "H2_PARTIAL_CROSS_SCALE_CORRESPONDENCE",
    "H3_PROCESS_BLOCK_ENVIRONMENT_REPRESENTATION",
    "H4_STABLE_CROSS_SCALE_COEFFICIENT_ROTATION",
    "H5_PHENOTYPE_SPACE_HANDOFF",
    "CANDIDATE_INVOLUCRE_ENVIRONMENT_STRUCTURE",
    "LEGACY_NEGATIVE_LABILITY_RESPONSIVENESS",
    "FUNCTIONAL_CAUSAL_MECHANISMS",
}

ALLOWED_STATUSES = {
    "operationally_complete_externally_unvalidated",
    "supported_primary_endpoints",
    "rejected",
    "supported",
    "supported_scale_specific",
    "inconclusive",
    "completed_workflow_boundary",
    "partially_supported_within_only",
    "rejected_withdrawn",
    "not_tested_in_azami",
}

EXPECTED_PROVENANCE = {
    "continuous_trait_run": 32975451732,
    "continuous_trait_head": "f4a6fd5e01a2befd4f49174984a99e53856c2330",
    "continuous_trait_artifact": 9612943217,
    "continuous_trait_digest": "sha256:101e996b638996a0c5ae79d358bf51293c3585f0e84c4a961b91dcbedf96211e",
    "multilevel_run": 33035785120,
    "multilevel_head": "227c0e7b8c338894806785b8545c7c77c8724de1",
    "multilevel_artifact": 9632715852,
    "multilevel_digest": "sha256:51e7a26b5bd09e030b67b9342586699abaaf46e630f45b6bb4ee7bfc9152ced6",
}

EXPECTED_SYNTHESIS = (
    "The thistle capitulum is a multidimensional, partially organized phenotype "
    "whose environmental alignment changes across biological scales."
)


def fail(message: str) -> None:
    raise SystemExit(message)


def main() -> int:
    data = json.loads(STATUS_PATH.read_text(encoding="utf-8"))

    if data.get("schema_version") != 1:
        fail("hypothesis recovery schema_version must equal 1")
    if data.get("overall_status") != "computational_pattern_estimand_complete_external_validation_open":
        fail("overall status must retain the computational-complete/external-validation-open boundary")
    if data.get("chapter_synthesis") != EXPECTED_SYNTHESIS:
        fail("Chapter 1 synthesis changed without a versioned result update")
    if data.get("frozen_sources") != EXPECTED_PROVENANCE:
        fail("immutable Azami result provenance changed")

    hypotheses = data.get("hypotheses", [])
    ids = [row.get("id") for row in hypotheses]
    if len(ids) != len(set(ids)):
        fail("duplicate hypothesis IDs")
    if set(ids) != EXPECTED_IDS:
        fail(f"unexpected hypothesis registry: {sorted(set(ids) ^ EXPECTED_IDS)}")

    by_id = {row["id"]: row for row in hypotheses}
    for row in hypotheses:
        if row.get("status") not in ALLOWED_STATUSES:
            fail(f"unsupported status for {row.get('id')}: {row.get('status')}")
        for field in ("scope", "key_result", "allowed_claim", "not_allowed"):
            if not str(row.get(field, "")).strip():
                fail(f"{row.get('id')} is missing {field}")

    required_fixed = {
        "FOUNDATION_ONE_UNIVERSAL_SYNDROME": "rejected",
        "H1_MEASUREMENT_MODULE_ORGANIZATION": "supported",
        "H2_PARTIAL_CROSS_SCALE_CORRESPONDENCE": "supported",
        "H3_PROCESS_BLOCK_ENVIRONMENT_REPRESENTATION": "supported_scale_specific",
        "H4_STABLE_CROSS_SCALE_COEFFICIENT_ROTATION": "inconclusive",
        "H5_PHENOTYPE_SPACE_HANDOFF": "completed_workflow_boundary",
        "LEGACY_NEGATIVE_LABILITY_RESPONSIVENESS": "rejected_withdrawn",
        "FUNCTIONAL_CAUSAL_MECHANISMS": "not_tested_in_azami",
    }
    for hypothesis_id, expected in required_fixed.items():
        if by_id[hypothesis_id]["status"] != expected:
            fail(f"{hypothesis_id} must remain {expected}")

    completion = data.get("completion", {})
    if completion.get("computational_pattern_layers", {}).get("status") != "complete":
        fail("computational pattern layer must be marked complete")
    if completion.get("independent_scientific_validation", {}).get("status") != "open_submission_blocking":
        fail("independent validation must remain submission blocking")
    if completion.get("independent_scientific_validation", {}).get("open_items") != 6:
        fail("expected six open independent validation gates")
    gates = data.get("open_validation_gates", [])
    if len(gates) != 6 or len(set(gates)) != 6:
        fail("open_validation_gates must contain six unique entries")

    forbidden_promotions = (
        "functional modularity is supported",
        "genetic modularity is supported",
        "plasticity is supported",
        "adaptation is supported",
        "pollinator causation is supported",
        "defensive efficacy is supported",
    )
    serialized = json.dumps(data, ensure_ascii=False).lower()
    for phrase in forbidden_promotions:
        if phrase in serialized:
            fail(f"forbidden causal promotion found: {phrase}")

    if "do not add an uncontracted correlational raster screen" not in data.get("next_step_rule", "").lower():
        fail("next-step rule must prohibit uncontracted raster fishing")

    print(json.dumps({
        "status": "ok",
        "hypotheses": len(hypotheses),
        "supported_or_partially_supported": sum(
            row["status"] in {
                "supported_primary_endpoints", "supported", "supported_scale_specific",
                "completed_workflow_boundary", "partially_supported_within_only",
            }
            for row in hypotheses
        ),
        "rejected_or_withdrawn": sum(
            row["status"] in {"rejected", "rejected_withdrawn"}
            for row in hypotheses
        ),
        "inconclusive": sum(row["status"] == "inconclusive" for row in hypotheses),
        "not_tested_in_azami": sum(row["status"] == "not_tested_in_azami" for row in hypotheses),
        "external_validation_gates_open": len(gates),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
