"""Versioned, outcome-blind enrichment of the unchanged native source cohort.

All three inputs require caller-pinned exact hashes. Available sibling receipts
are checked too. The old cohort and source databases are opened read-only; this
does not rebuild taxonomy, change membership, or authorize image/trait execution.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import date
import hashlib
import json
import math
from pathlib import Path
import re
import sqlite3

from .prepare_observation_annotations import annotate
from .workflow import ROOT, digest, text_digest

VERSION = "v3_enriched_native_source_cohort_v1"
STATUS = "ENRICHED_SOURCE_COHORT_VERIFIED_NO_ECOLOGICAL_AUTHORIZATION"
COHORT_FIELDS = ["obs_id", "accepted_key", "accepted_name", "source_taxon_name", "source_taxon_rank",
                 "analysis_latitude", "analysis_longitude", "observation_month", "position_accuracy_status",
                 "position_accuracy_m", "quarter_degree_cell", "native_range_status"]
CALENDAR_FIELDS = ["observed_year", "doy", "days_in_year", "sin_doy", "cos_doy", "hemisphere",
                   "south_indicator", "south_sin", "south_cos"]
EXTRA_FIELDS = ["observed_on", *CALENDAR_FIELDS, "date_status", "coordinate_status", "captive_state",
                "dependence_component_id", "dependence_component_source_observations", "prior_annotation_component_id",
                "reconciled_source_photo_count", "source_user_id", "source_quality_grade", "source_created_at",
                "source_updated_at", "source_annotation_records_json", "source_preference", "source_conflicted_fields_json",
                "image_quality_covariates_status"]
ROUNDTRIP_ABSOLUTE_TOLERANCE = 1e-10


def pinned_hash(path: Path, expected: str, role: str) -> str:
    if not re.fullmatch(r"[0-9a-fA-F]{64}", expected or ""):
        raise ValueError(f"{role} requires a caller-pinned SHA-256")
    actual = digest(path)
    if actual != expected.lower():
        raise ValueError(f"{role} input identity differs from pinned SHA-256")
    return actual


def checked_receipts(paths: dict[str, Path], hashes: dict[str, str]) -> dict:
    specs = {"cohort": ("ecological_source_cohort_report.json", "cohort_sha256"),
             "annotations": ("observation_annotations_report.json", "output_database_sha256"),
             "reconciliation": ("source_reconciliation_report.json", "ledger_sha256")}
    receipts = {}
    for role, (name, key) in specs.items():
        if role not in paths:
            continue
        path = paths[role].with_name(name)
        if not path.exists():
            receipts[role] = {"present": False, "verification": "caller_pinned_input_sha256"}
            continue
        report = json.loads(path.read_text(encoding="utf-8"))
        if report.get(key) != hashes[role]:
            raise ValueError(f"{role} sibling receipt input identity mismatch")
        if role == "annotations" and report.get("execution_contract", {}).get("source_reconciliation_sha256") != hashes["reconciliation"]:
            raise ValueError("Annotation receipt references a different reconciliation")
        receipts[role] = {"present": True, "receipt_sha256": digest(path), "verification": "receipt_and_caller_pin"}
    return receipts


def probe(annotations: Path, reconciliation: Path, *, expected_annotations_sha256: str,
          expected_reconciliation_sha256: str) -> dict:
    """Verify available source/dependence inputs without claiming a native cohort."""
    paths = {"annotations": annotations.resolve(), "reconciliation": reconciliation.resolve()}
    hashes = {role: pinned_hash(paths[role], pin, role) for role, pin in (
        ("annotations", expected_annotations_sha256), ("reconciliation", expected_reconciliation_sha256))}
    receipts = checked_receipts(paths, hashes)
    with sqlite3.connect(paths["annotations"].as_uri() + "?mode=ro", uri=True) as ann, sqlite3.connect(paths["reconciliation"].as_uri() + "?mode=ro", uri=True) as source:
        _, _, _, counts = source_components(source, ann, set())
    return {"schema_version": 1, "status": "SOURCE_INTEGRITY_VERIFIED_NATIVE_COHORT_UNCHECKED",
            "input_sha256": hashes, "checked_upstream_receipts": receipts, "dependence": counts,
            "implementation_sha256_text_lf": text_digest(Path(__file__)),
            "native_cohort_verified": False, "ecological_fitting_authorized": False,
            "production_image_execution_authorized": False, "trait_files_read": 0,
            "limits": ["This probe verifies source identities, full membership and known dependence only; it does not verify the missing native cohort or native authority chain.",
                       "Known dependence components do not prove independence of all other records."]}


def read_cohort(path: Path) -> tuple[dict[str, dict], list[str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)) or set(fields) != set(COHORT_FIELDS):
            raise ValueError("Old cohort columns differ from the source-only schema; do not admit trait columns")
        rows = {}
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError("Malformed cohort source row")
            obs = row["obs_id"].strip()
            if not obs or obs in rows:
                raise ValueError("Old cohort observation identities are missing or duplicated")
            if row["native_range_status"] != "native" or row["source_taxon_rank"] not in {"species", "subspecies", "variety"} or not row["accepted_key"]:
                raise ValueError("Old cohort does not meet its unchanged native taxonomic scope")
            rows[obs] = row
    if not rows:
        raise ValueError("Old cohort is empty")
    return rows, fields


class SourceComponents:
    def __init__(self, observations):
        self.parent = {obs: obs for obs in observations}

    def find(self, obs):
        if obs not in self.parent:
            raise ValueError("Dependence link has an observation outside the reconciled source")
        while self.parent[obs] != obs:
            self.parent[obs] = self.parent[self.parent[obs]]
            obs = self.parent[obs]
        return obs

    def union(self, left, right):
        left, right = self.find(left), self.find(right)
        if left != right:
            self.parent[max(left, right)] = min(left, right)

    def identifier(self, obs):
        # The least source observation is deterministic; exact input hashes bind
        # the full component membership, including bridges outside this cohort.
        return "source_component_" + hashlib.sha256((VERSION + "|" + self.find(obs)).encode()).hexdigest()


def source_components(source, annotations, cohort_ids):
    raw_source_rows = [row[0] for row in source.execute("SELECT obs_id FROM source_observations ORDER BY obs_id")]
    if any(obs is None for obs in raw_source_rows):
        raise ValueError("Reconciled source observation identities are invalid")
    source_rows = [str(obs) for obs in raw_source_rows]
    if any(not obs for obs in source_rows) or len(source_rows) != len(set(source_rows)):
        raise ValueError("Reconciled source observation identities are invalid")
    groups = SourceComponents(source_rows)
    annotation_rows = [str(row[0]) for row in annotations.execute("SELECT obs_id FROM annotations")]
    if len(annotation_rows) != len(set(annotation_rows)) or set(annotation_rows) != set(source_rows):
        raise ValueError("Full annotation and reconciled source observation membership differ")
    if not set(cohort_ids).issubset(groups.parent):
        raise ValueError("A cohort observation is absent from the reconciled source")
    photo_counts = Counter()
    previous_photo = previous_obs = first_obs = None
    link_count = 0
    for photo, obs in source.execute("SELECT photo_id,obs_id FROM source_links ORDER BY photo_id,obs_id"):
        photo, obs = str(photo or ""), str(obs or "")
        if not photo or not obs:
            raise ValueError("Reconciled source contains an empty photo link")
        if photo == previous_photo and obs == previous_obs:
            raise ValueError("Reconciled source contains a duplicate photo link")
        groups.find(obs)
        if photo == previous_photo:
            groups.union(first_obs, obs)
        else:
            previous_photo, first_obs = photo, obs
        previous_obs = obs
        link_count += 1
        if obs in cohort_ids:
            photo_counts[obs] += 1
    columns = {row[1] for row in annotations.execute("PRAGMA table_info(annotations)")}
    prior_present = "component_id" in columns
    if prior_present:
        previous_component = first_obs = None
        for component, obs in annotations.execute("SELECT component_id,obs_id FROM annotations ORDER BY component_id,obs_id"):
            if component is None or not str(component):
                raise ValueError("Existing annotation dependence component is missing")
            if component == previous_component:
                groups.union(first_obs, str(obs))
            else:
                previous_component, first_obs = component, str(obs)
    component_sizes = Counter(groups.find(obs) for obs in source_rows)
    return groups, component_sizes, photo_counts, {
        "full_source_observations": len(source_rows), "reconciled_source_photo_links": link_count,
        "known_dependence_components": len(component_sizes), "prior_annotation_components_preserved": prior_present,
        "full_source_components_with_multiple_observations": sum(size > 1 for size in component_sizes.values()),
    }


def same_number(left, right) -> bool:
    try:
        return math.isfinite(float(left)) and math.isfinite(float(right)) and math.isclose(float(left), float(right), rel_tol=0, abs_tol=ROUNDTRIP_ABSOLUTE_TOLERANCE)
    except (TypeError, ValueError):
        return False


def enrich_row(row: dict, annotation: dict, component: str, component_size: int, photo_count: int) -> dict:
    fields = json.loads(annotation["selected_fields_json"])
    conflicts = json.loads(annotation["conflicted_fields_json"])
    if not isinstance(fields, dict) or not isinstance(conflicts, list):
        raise ValueError("Invalid annotation source provenance")
    derived = annotate(fields, conflicts, annotation["preferred_source_kind"])
    if derived["date_status"] != "exact_day" or derived["coordinate_status"] != "public_location_present_precision_not_gated" or derived["captive_state"] != "false":
        raise ValueError("Annotation date, privacy or wild-status conflicts with eligible old cohort")
    for name in ("source_taxon_name", "source_taxon_rank", "date_status", "coordinate_status", "captive_state", "hemisphere", "position_accuracy_status"):
        if str(annotation.get(name)) != str(derived[name]):
            raise ValueError("Saved annotation state differs from source provenance")
    for name in ("source_taxon_name", "source_taxon_rank"):
        if row[name] != derived[name]:
            raise ValueError("Cohort taxon assignment differs from source annotation")
    for name in ("analysis_latitude", "analysis_longitude"):
        if not same_number(row[name], derived[name]) or not same_number(annotation.get(name), derived[name]):
            raise ValueError("Cohort and source annotation coordinates differ")
    for name in CALENDAR_FIELDS:
        if name != "hemisphere" and not same_number(annotation.get(name), derived[name]):
            raise ValueError("Saved calendar terms differ from source exact date")
    if derived["position_accuracy_m"] is not None or annotation.get("position_accuracy_m") is not None:
        if not same_number(derived["position_accuracy_m"], annotation.get("position_accuracy_m")):
            raise ValueError("Saved positional accuracy differs from source provenance")
    observed_on = fields["observed_on"]
    if not same_number(row["observation_month"], date.fromisoformat(observed_on).month):
        raise ValueError("Old cohort month differs from source exact date")
    if row["position_accuracy_status"] != annotation.get("position_accuracy_status"):
        raise ValueError("Cohort positional-accuracy status differs from annotation")
    if row["position_accuracy_m"] or annotation.get("position_accuracy_m") is not None:
        if not same_number(row["position_accuracy_m"], annotation.get("position_accuracy_m")):
            raise ValueError("Cohort positional accuracy differs from annotation")
    return {**row, "observed_on": observed_on, **{name: derived[name] for name in CALENDAR_FIELDS},
            **{name: derived[name] for name in ("date_status", "coordinate_status", "captive_state")},
            "dependence_component_id": component, "dependence_component_source_observations": component_size,
            "prior_annotation_component_id": annotation.get("component_id") or "",
            "reconciled_source_photo_count": photo_count,
            "source_user_id": fields.get("user_id", ""), "source_quality_grade": fields.get("quality_grade", ""),
            "source_created_at": fields.get("created_at", ""), "source_updated_at": fields.get("updated_at", ""),
            "source_annotation_records_json": fields.get("annotation_records", ""),
            "source_preference": annotation["preferred_source_kind"],
            "source_conflicted_fields_json": annotation["conflicted_fields_json"],
            "image_quality_covariates_status": "not_joined_requires_endpoint_matched_measurements"}


def build(cohort: Path, annotations: Path, reconciliation: Path, out: Path, *,
          expected_cohort_sha256: str, expected_annotations_sha256: str, expected_reconciliation_sha256: str) -> dict:
    paths = {"cohort": cohort.resolve(), "annotations": annotations.resolve(), "reconciliation": reconciliation.resolve()}
    out = out.resolve()
    if out == ROOT or (ROOT in out.parents and not any(out.is_relative_to(ROOT / name) for name in ("local_data", "outputs"))):
        raise ValueError("Private cohort outputs require an external or ignored output directory")
    if out.exists():
        raise ValueError("Output exists; preserve every previous cohort version")
    hashes = {role: pinned_hash(paths[role], pin, role) for role, pin in
              (("cohort", expected_cohort_sha256), ("annotations", expected_annotations_sha256), ("reconciliation", expected_reconciliation_sha256))}
    receipts = checked_receipts(paths, hashes)
    rows, columns = read_cohort(paths["cohort"])
    with sqlite3.connect(paths["annotations"].as_uri() + "?mode=ro", uri=True) as ann, sqlite3.connect(paths["reconciliation"].as_uri() + "?mode=ro", uri=True) as source:
        ann.row_factory = sqlite3.Row
        annotation_columns = {r[1] for r in ann.execute("PRAGMA table_info(annotations)")}
        needed = {"obs_id", "selected_fields_json", "conflicted_fields_json", "preferred_source_kind", "source_taxon_name", "source_taxon_rank",
                  "date_status", "coordinate_status", "captive_state", "analysis_latitude", "analysis_longitude", "position_accuracy_status", "position_accuracy_m", *CALENDAR_FIELDS}
        if not needed.issubset(annotation_columns):
            raise ValueError("Annotation schema lacks the exact-date/source provenance fields")
        groups, sizes, photo_counts, counts = source_components(source, ann, rows)
        selected_columns = sorted(needed | ({"component_id"} if "component_id" in annotation_columns else set()))
        out.mkdir(parents=True)
        try:
            hemisphere_counts = Counter()
            components = set()
            written = 0
            csv_path = out / "enriched_ecological_source_cohort_private.csv"
            with csv_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=[*columns, *EXTRA_FIELDS], lineterminator="\n")
                writer.writeheader()
                for record in ann.execute("SELECT " + ",".join(selected_columns) + " FROM annotations ORDER BY obs_id"):
                    obs = str(record["obs_id"])
                    if obs not in rows:
                        continue
                    root = groups.find(obs)
                    enriched = enrich_row(rows[obs], dict(record), groups.identifier(obs), sizes[root], photo_counts[obs])
                    writer.writerow(enriched)
                    hemisphere_counts[enriched["hemisphere"]] += 1
                    components.add(root)
                    written += 1
            if written != len(rows):
                raise ValueError("Enriched cohort membership differs from the unchanged source cohort")
            report = {"schema_version": 1, "version": VERSION, "status": STATUS,
                      "input_sha256": hashes, "checked_upstream_receipts": receipts,
                      "implementation_sha256_text_lf": text_digest(Path(__file__)),
                      "annotation_helper_sha256_text_lf": text_digest(Path(__file__).with_name("prepare_observation_annotations.py")),
                      "cohort_rows": written, "accepted_taxa": len({r["accepted_key"] for r in rows.values()}),
                      "cohort_membership_changed": False, "source_rows_deleted": 0,
                      "dependence": {**counts, "cohort_known_components": len(components),
                                     "cohort_observations_with_source_photo_links": sum(photo_counts[obs] > 0 for obs in rows)},
                      "calendar": {"exact_date_rows": written, "hemisphere_counts": dict(hemisphere_counts),
                                   "recomputed_and_compared_to_source_annotation": True,
                                   "numeric_roundtrip_absolute_tolerance": ROUNDTRIP_ABSOLUTE_TOLERANCE},
                      "output_csv_sha256": digest(csv_path), "output_csv_bytes": csv_path.stat().st_size,
                      "trait_files_read": 0, "environment_values_read": 0, "ecological_models_executed": 0,
                      "ecological_fitting_authorized": False, "production_image_execution_authorized": False,
                      "limits": ["Known shared-photo and existing annotation components are retained across the full source; this does not prove all remaining observations independent.",
                                 "Exact dates and hemisphere harmonics describe calendar timing, not observed developmental stage.",
                                 "Endpoint-matched resolution, sharpness and colour context are not present until joined from verified measurement products.",
                                 "Native classification is inherited from the pinned old cohort; no taxonomy or distribution decision was changed."]}
            (out / "enriched_source_cohort_report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
            return report
        except Exception as error:
            (out / "incomplete_run.json").write_text(json.dumps({"status": "INCOMPLETE_DO_NOT_USE", "error_type": type(error).__name__}) + "\n", encoding="utf-8")
            raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-only-probe", action="store_true")
    parser.add_argument("--cohort", type=Path)
    parser.add_argument("--expected-cohort-sha256")
    for name in ("annotations", "reconciliation", "out-dir"):
        parser.add_argument("--" + name, type=Path, required=True)
    for role in ("annotations", "reconciliation"):
        parser.add_argument("--expected-" + role + "-sha256", required=True)
    args = parser.parse_args()
    if args.source_only_probe:
        if args.cohort or args.expected_cohort_sha256:
            parser.error("Source-only probe cannot verify or accept a cohort")
        if args.out_dir.exists():
            raise ValueError("Preserve prior source probe")
        report = probe(args.annotations, args.reconciliation,
                       expected_annotations_sha256=args.expected_annotations_sha256,
                       expected_reconciliation_sha256=args.expected_reconciliation_sha256)
        args.out_dir.mkdir(parents=True)
        (args.out_dir / "source_integrity_probe.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        print(json.dumps(report, indent=2))
        return 0
    if not args.cohort or not args.expected_cohort_sha256:
        parser.error("Enrichment requires the original cohort and its pinned SHA-256")
    report = build(args.cohort, args.annotations, args.reconciliation, args.out_dir,
                   expected_cohort_sha256=args.expected_cohort_sha256, expected_annotations_sha256=args.expected_annotations_sha256,
                   expected_reconciliation_sha256=args.expected_reconciliation_sha256)
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
