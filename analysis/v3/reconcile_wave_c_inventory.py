"""Reconcile Wave C completion from protected numerical asset metadata only.

A GitHub Actions job failure is not a measurement failure: a completed chunk may
already have been uploaded to the unpublished numerical draft before a return
verification timed out.  This audit therefore uses the protected final-asset
inventory, not Actions job conclusions, to decide which chunks need resume.

No protected bundle is downloaded and no observation, coordinate, trait or
environment value is inspected.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

from .protected_artifacts import DraftStore, new_json, require

FINAL_DIGEST = re.compile(r"sha256:[0-9a-f]{64}")


def reconcile(batch: dict, assets: list[dict]) -> dict:
    require(batch.get("schema_version") == 1, "Unknown Wave C batch schema")
    require(batch.get("status") == "bounded_native_raw_measurement_batch_no_ecology", "Wrong Wave C batch scope")
    chunks = sorted(batch["chunks"])
    require(chunks == [f"c{i:06d}" for i in range(18, 146)], "Wave C chunk universe changed")
    plan_id = batch["plan_id"]
    require(re.fullmatch(r"[0-9a-f]{64}", plan_id) is not None, "Invalid Wave C plan id")
    prefix = f"v3-raw-{plan_id[:16]}-"

    finals: dict[str, list[dict]] = {chunk: [] for chunk in chunks}
    partials: dict[str, int] = {chunk: 0 for chunk in chunks}
    for asset in assets:
        name = asset.get("name", "")
        match = re.fullmatch(re.escape(prefix) + r"(c[0-9]{6})\.zip", name)
        if match and match.group(1) in finals:
            chunk = match.group(1)
            require(asset.get("state") == "uploaded", f"Protected final asset is not uploaded: {chunk}")
            require(type(asset.get("size")) is int and asset["size"] > 0, f"Protected final asset has invalid size: {chunk}")
            require(FINAL_DIGEST.fullmatch(asset.get("digest", "")) is not None, f"Protected final asset lacks server SHA-256: {chunk}")
            finals[chunk].append(asset)
            continue
        partial = re.fullmatch(re.escape(prefix) + r"(c[0-9]{6})-incomplete-[0-9]+-[0-9]+\.zip", name)
        if partial and partial.group(1) in partials:
            require(asset.get("state") == "uploaded", f"Protected checkpoint asset is not uploaded: {partial.group(1)}")
            partials[partial.group(1)] += 1

    duplicates = [chunk for chunk, rows in finals.items() if len(rows) != len({row.get("id") for row in rows}) or len(rows) > 1]
    require(not duplicates, "Duplicate protected final assets for Wave C chunks: " + ",".join(duplicates))

    completed = [chunk for chunk in chunks if len(finals[chunk]) == 1]
    missing = [chunk for chunk in chunks if not finals[chunk]]
    missing_with_checkpoint = [chunk for chunk in missing if partials[chunk] > 0]
    missing_without_checkpoint = [chunk for chunk in missing if partials[chunk] == 0]

    return {
        "schema_version": 1,
        "status": (
            "WAVE_C_PROTECTED_FINAL_INVENTORY_COMPLETE_128_OF_128_NO_RESUME"
            if not missing
            else "WAVE_C_PROTECTED_FINAL_INVENTORY_INCOMPLETE_RESUME_ONLY_MISSING"
        ),
        "plan_id_prefix": plan_id[:16],
        "planned_chunks": len(chunks),
        "protected_final_chunks": len(completed),
        "missing_final_chunks": len(missing),
        "missing_chunks": missing,
        "missing_with_protected_checkpoint": missing_with_checkpoint,
        "missing_without_protected_checkpoint": missing_without_checkpoint,
        "resume_chunks": missing,
        "aggregate_source_observations": batch["aggregate"]["observations"],
        "aggregate_source_photo_jobs": batch["aggregate"]["photo_jobs"],
        "aggregate_request_slots": batch["aggregate"]["requests"],
        "inventory_basis": "authenticated unpublished-draft asset metadata; final asset name + uploaded state + positive byte size + GitHub SHA-256",
        "actions_job_conclusion_used_for_completion": False,
        "protected_numerical_payload_downloaded": False,
        "source_images_downloaded": 0,
        "source_images_persisted": 0,
        "observation_identifiers_read": 0,
        "coordinates_read": 0,
        "trait_values_read": 0,
        "environment_values_read": 0,
        "empirical_trait_environment_values_read": 0,
        "ecological_models_executed": 0,
        "ecological_fitting_authorized": False,
        "limits": [
            "A protected final asset proves a completed numerical chunk transaction was committed, not biological measurement accuracy.",
            "This inventory audit does not restore or inspect numerical rows and does not authorize ecological inference.",
            "Only chunks lacking a unique protected final asset may be resumed; completed chunks must not be remeasured."
        ],
    }


def run(batch_path: Path, out: Path) -> dict:
    require(not out.exists(), "Preserve previous Wave C reconciliation output")
    batch = json.loads(batch_path.read_text(encoding="utf-8"))
    contract = {
        "repository": batch["repository"],
        "release_id": batch["release_id"],
        "tag_name": batch["tag_name"],
    }
    store = DraftStore(contract)
    try:
        release = store.check()
        report = reconcile(batch, release["assets"])
    finally:
        store.close()
    out.mkdir(parents=True)
    new_json(out / "public_report.json", report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch", type=Path, default=Path("analysis/v3/native_measurement_wave_20260908_c.json"))
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.batch, args.out), sort_keys=True))


if __name__ == "__main__":
    main()
