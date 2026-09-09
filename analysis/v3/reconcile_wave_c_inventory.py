"""Reconcile Wave C completion from protected numerical asset metadata only.

A GitHub Actions job failure is not a measurement failure: a completed chunk may
already have been uploaded to the unpublished numerical draft before a return
verification timed out. This audit therefore uses the protected final-asset
inventory, not Actions job conclusions, to decide which chunks need resume.

No protected bundle is downloaded and no observation, coordinate, trait or
environment value is inspected.
"""
from __future__ import annotations

import argparse
import json
import os
import re
from pathlib import Path

import requests

FINAL_DIGEST = re.compile(r"sha256:[0-9a-f]{64}")
REPOSITORY = "zuizui0223/azami"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def new_json(path: Path, value: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as handle:
        json.dump(value, handle, sort_keys=True, indent=2, allow_nan=False)
        handle.write("\n")


def protected_asset_inventory(batch: dict) -> list[dict]:
    """Read only authenticated draft-release metadata with anonymous denial."""
    require(batch["repository"] == REPOSITORY, "Wave C repository differs")
    release_id = batch["release_id"]
    require(type(release_id) is int and release_id > 0, "Invalid Wave C release id")
    require(batch["tag_name"].startswith("private-v3-numerical-"), "Unexpected protected tag")
    token = os.environ.get("GH_TOKEN")
    require(bool(token), "GitHub token required for protected inventory")

    base = f"https://api.github.com/repos/{REPOSITORY}"
    url = f"{base}/releases/{release_id}"
    session = requests.Session()
    anonymous = requests.Session()
    anonymous.trust_env = False
    session.headers.update({"Authorization": "Bearer " + token, "X-GitHub-Api-Version": "2022-11-28"})
    try:
        response = session.get(url, timeout=(15, 60))
        response.raise_for_status()
        release = response.json()
        require(release.get("draft") is True and release.get("published_at") is None, "Protected release is published or not draft; STOP")
        require(release.get("tag_name") == batch["tag_name"] and release.get("id") == release_id, "Protected draft identity differs")
        public = anonymous.get(url, timeout=(15, 60), allow_redirects=False)
        require(public.status_code == 404, "Protected draft is not anonymously denied; STOP")

        assets, seen = [], set()
        for page in range(1, 1001):
            reply = session.get(url + "/assets", params={"per_page": 100, "page": page}, timeout=(15, 60))
            reply.raise_for_status()
            rows = reply.json()
            require(isinstance(rows, list) and len(rows) <= 100, "Invalid protected asset page")
            for row in rows:
                require(type(row.get("id")) is int and row["id"] not in seen, "Duplicate or invalid protected asset identity")
                seen.add(row["id"])
                assets.append(row)
            if len(rows) < 100:
                break
        else:
            raise ValueError("Protected inventory exceeds bounded pagination; STOP")
        return assets
    finally:
        session.close()
        anonymous.close()


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
    report = reconcile(batch, protected_asset_inventory(batch))
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
