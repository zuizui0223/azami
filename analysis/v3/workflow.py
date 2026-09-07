#!/usr/bin/env python3
"""Loss-accounted v3 entry point from the full recovered photo metadata.

Standard library only. Inventory does not download images, thin observations,
select native records, run a detector, or fit ecological models.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import platform
import sqlite3

ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "analysis/v3/workflow_contract.json"
NATIVE_SHA = "c01eeb9ff245d7f73da1a12fa4eede904dd9770467655f20e3d85de2ac8dd84a"
FIELDS = ["source_row", "obs_id", "photo_id", "photo_index", "taxon_id", "taxon_name", "taxon_rank",
          "coordinate_status", "captive_status", "photo_license_code", "legacy_v2_observation",
          "native_status_v2", "metadata_chunk_source", "identity_status"]


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def canonical_digest(data: dict) -> str:
    return hashlib.sha256(json.dumps(data, sort_keys=True, separators=(",", ":"),
                                     ensure_ascii=False, allow_nan=False).encode()).hexdigest()


def text_digest(path: Path) -> str:
    return hashlib.sha256(path.read_text(encoding="utf-8").encode()).hexdigest()


def validate_contract(contract: dict) -> None:
    if contract["schema_version"] != 2:
        raise ValueError("Unsupported full-source schema")
    if contract["study_status"] != "retrospective_redesign_after_v2_results_seen":
        raise ValueError("Prior outcome exposure must remain explicit")
    if contract["legacy_v2_native_reference_sha256"] != NATIVE_SHA:
        raise ValueError("Legacy reference identity changed")
    preservation = contract["preservation"]
    for field in ("retain_every_source_row", "retain_all_photos_per_observation"):
        if preservation[field] is not True:
            raise ValueError("Full-source preservation required: " + field)
    for field in ("drop_rows_for_qc", "thin_at_acquisition", "first_photo_only", "native_only_at_acquisition",
                  "require_coordinates_for_image_measurement", "delete_detector_negative_photos"):
        if preservation[field] is not False:
            raise ValueError("A downstream filter cannot define the acquisition universe: " + field)
    if not all(contract["reporting"][field] is True for field in
               ("no_required_result_direction", "no_requirement_that_candidates_survive")):
        raise ValueError("Evaluation must not predetermine findings")
    if [s["id"] for s in contract["stages"]] != ["inventory", "recover", "measure", "analyse", "report"]:
        raise ValueError("Full-source inventory must precede measurement and analysis")


def boolean_state(value: str) -> str:
    value = value.strip().lower()
    if value in ("true", "1"):
        return "true"
    if value in ("false", "0"):
        return "false"
    return "unknown"  # The original value remains in the immutable source CSV.


def load_legacy_native(path: Path | None) -> dict:
    if path is None:
        return {}
    if digest(path) != NATIVE_SHA:
        raise ValueError("Legacy v2 native table hash mismatch")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len({r["obs_id"] for r in rows}) != len(rows):
        raise ValueError("Duplicate legacy observation ID")
    return {r["obs_id"]: (r["taxon_name"], r["native_range_status"]) for r in rows}


def inventory_records(metadata: Path, out: Path, contract: dict,
                      legacy: dict | None = None) -> dict:
    """Stream every source row into a local relational ledger without filtering.

    `legacy` is an optional identity-linked comparison only, never an eligibility
    condition. Unknown source identities and duplicate links remain visible.
    """
    validate_contract(contract)
    out = out.resolve()
    if out == ROOT or (ROOT in out.parents and not any(
            out.is_relative_to(ROOT / allowed) for allowed in ("local_data", "outputs"))):
        raise ValueError("Use a new external directory or local_data/outputs inside the repository")
    if out.exists():
        raise ValueError("Output directory already exists; preserve earlier runs")
    source_hash = digest(metadata)
    if source_hash != contract["source"]["member_sha256"]:
        raise ValueError("Source metadata SHA-256 mismatch")
    legacy_supplied = legacy is not None
    legacy = legacy or {}
    out.mkdir(parents=True, exist_ok=False)
    database = out / "source_ledger.sqlite"
    try:
        with sqlite3.connect(database) as db, metadata.open(encoding="utf-8-sig", newline="") as handle:
            db.execute("CREATE TABLE photo_records (" + ",".join(
                f"{name} {'INTEGER PRIMARY KEY' if name == 'source_row' else 'INTEGER' if name == 'legacy_v2_observation' else 'TEXT'}"
                for name in FIELDS) + ")")
            source = csv.DictReader(handle)
            required = {"obs_id", "photo_id", "photo_index", "taxon_name", "taxon_rank",
                        "coordinate_usable_for_environment", "captive", "photo_license_code"}
            if not required.issubset(source.fieldnames or []):
                raise ValueError("Missing source columns: " + str(sorted(required - set(source.fieldnames or []))))
            batch = []
            count = 0
            for count, row in enumerate(source, start=1):
                if None in row or any(v is None for v in row.values()):
                    raise ValueError(f"Malformed source CSV row {count}; retain source and investigate")
                obs, photo = row["obs_id"], row["photo_id"]
                old = legacy.get(obs)
                native = (old[1] if old and old[0] == row["taxon_name"] else
                          "not_assessed_taxon_changed" if old else "not_assessed_outside_v2")
                batch.append((count, obs, photo, row["photo_index"], row.get("taxon_id", ""),
                              row["taxon_name"], row["taxon_rank"],
                              boolean_state(row["coordinate_usable_for_environment"]), boolean_state(row["captive"]),
                              row["photo_license_code"], int(old is not None), native,
                              row.get("metadata_chunk_source", ""), "linked" if obs and photo else "missing_identifier"))
                if len(batch) >= 20000:
                    db.executemany("INSERT INTO photo_records VALUES (" + ",".join(["?"] * len(FIELDS)) + ")", batch)
                    batch.clear()
            if batch:
                db.executemany("INSERT INTO photo_records VALUES (" + ",".join(["?"] * len(FIELDS)) + ")", batch)
            db.execute("CREATE INDEX observation_lookup ON photo_records(obs_id)")
            db.execute("CREATE INDEX photo_lookup ON photo_records(photo_id)")
            db.execute("""CREATE TABLE observation_records AS SELECT obs_id, COUNT(*) AS n_source_photo_rows,
                COUNT(DISTINCT NULLIF(photo_id,'')) AS n_unique_photo_ids,
                COUNT(DISTINCT NULLIF(taxon_name,'')) AS n_source_taxon_labels,
                MAX(legacy_v2_observation) AS legacy_v2_observation
                FROM photo_records WHERE obs_id <> '' GROUP BY obs_id""")
            scalar = lambda sql: db.execute(sql).fetchone()[0]
            n_observations = scalar("SELECT COUNT(*) FROM observation_records")
            if count != contract["source"]["expected_photo_rows"] or n_observations != contract["source"]["expected_observations_with_photos"]:
                raise ValueError("Source photo/observation counts differ from the pinned snapshot")
            counts = {
                "source_photo_rows": count,
                "retained_photo_record_rows": scalar("SELECT COUNT(*) FROM photo_records"),
                "unique_photo_ids": scalar("SELECT COUNT(DISTINCT photo_id) FROM photo_records WHERE photo_id <> ''"),
                "observations_with_photos": n_observations,
                "observations_with_multiple_photos": scalar("SELECT COUNT(*) FROM observation_records WHERE n_unique_photo_ids > 1"),
                "additional_photos_beyond_one_per_observation": scalar("SELECT COALESCE(SUM(n_unique_photo_ids-1),0) FROM observation_records WHERE n_unique_photo_ids > 1"),
                "species_rank_taxa": scalar("SELECT COUNT(DISTINCT taxon_name) FROM photo_records WHERE taxon_rank = 'species' AND taxon_name <> ''"),
                "rows_with_missing_identifiers": scalar("SELECT COUNT(*) FROM photo_records WHERE identity_status <> 'linked'"),
                "duplicate_observation_photo_rows": scalar("SELECT COALESCE(SUM(n-1),0) FROM (SELECT COUNT(*) n FROM photo_records GROUP BY obs_id,photo_id)"),
                "photo_ids_linked_to_multiple_observations": scalar("SELECT COUNT(*) FROM (SELECT photo_id FROM photo_records WHERE photo_id <> '' GROUP BY photo_id HAVING COUNT(DISTINCT obs_id)>1)"),
                "observations_with_conflicting_taxon_labels": scalar("SELECT COUNT(*) FROM observation_records WHERE n_source_taxon_labels>1"),
                "rows_removed_in_v3_inventory": 0,
            }
            if counts["retained_photo_record_rows"] != count:
                raise ValueError("Source-row conservation failed")
            groups = {}
            for column in ("coordinate_status", "captive_status", "taxon_rank", "photo_license_code", "native_status_v2"):
                groups[column] = dict(db.execute(f"SELECT {column}, COUNT(*) FROM photo_records GROUP BY {column} ORDER BY {column}"))
            overlap = scalar("SELECT COUNT(*) FROM observation_records WHERE legacy_v2_observation=1")
            legacy_summary = {"status": "LINKED_COMPARISON_ONLY" if legacy_supplied else "NOT_SUPPLIED",
                              "source_v2_observations": len(legacy) if legacy_supplied else None,
                              "overlapping_observations": overlap if legacy_supplied else None,
                              "missing_from_photo_snapshot": len(legacy) - overlap if legacy_supplied else None}
            with (out / "observation_ledger.csv").open("w", encoding="utf-8", newline="") as handle:
                cursor = db.execute("SELECT * FROM observation_records ORDER BY obs_id")
                writer = csv.writer(handle, lineterminator="\n")
                writer.writerow([c[0] for c in cursor.description])
                writer.writerows(cursor)
        report = {
            "workflow_id": contract["workflow_id"], "status": "FULL_PHOTO_METADATA_INVENTORIED_NOT_IMAGE_VALIDATED",
            "contract_sha256_canonical_json": canonical_digest(contract), "metadata_sha256": source_hash,
            "implementation_sha256_text_lf": text_digest(Path(__file__)),
            "counts": counts, "photo_record_counts_by_status": groups, "legacy_v2_comparison": legacy_summary,
            "legacy_mapping_sha256_canonical_json": canonical_digest(legacy) if legacy_supplied else None,
            "stages": {"inventory": "COMPLETED_FOR_RECOVERED_PHOTO_SNAPSHOT", "recover": "NOT_EXECUTED",
                       "measure": "NOT_EXECUTED", "analyse": "NOT_EXECUTED", "report": "NOT_EXECUTED"},
            "image_download_completeness": "NOT_ASSESSED_FROM_METADATA", "detector_negative_count": None,
            "image_operations_performed": False, "ecological_models_executed": False,
            "upstream_merge_loss": contract["source"]["earlier_source_limit"],
            "raw_observations_without_photos": "NOT_ENUMERATED_BY_PHOTO_METADATA",
            "software": {"python": platform.python_version(), "sqlite": sqlite3.sqlite_version},
        }
        (out / "workflow_contract.json").write_text(json.dumps(contract, indent=2) + "\n", encoding="utf-8", newline="\n")
        report["output_sha256"] = {name: digest(out / name) for name in
                                   ("source_ledger.sqlite", "observation_ledger.csv", "workflow_contract.json")}
        (out / "source_inventory_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8", newline="\n")
        return report
    except Exception as error:
        (out / "incomplete_run.json").write_text(json.dumps({"status": "INCOMPLETE_DO_NOT_USE",
            "error_type": type(error).__name__, "reason": str(error)}, indent=2) + "\n", encoding="utf-8")
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["plan", "inventory"])
    parser.add_argument("--contract", type=Path, default=CONTRACT)
    parser.add_argument("--metadata", type=Path)
    parser.add_argument("--legacy-v2-native", type=Path)
    parser.add_argument("--out-dir", type=Path)
    args = parser.parse_args()
    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    validate_contract(contract)
    if args.command == "plan":
        result = {"status": "DESIGN_CHECKED_NOT_EXECUTED", "workflow_id": contract["workflow_id"],
                  "source": contract["source"], "stages": contract["stages"]}
    else:
        if args.metadata is None or args.out_dir is None:
            parser.error("inventory requires --metadata and --out-dir")
        legacy = load_legacy_native(args.legacy_v2_native) if args.legacy_v2_native else None
        result = inventory_records(args.metadata, args.out_dir, contract, legacy)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
