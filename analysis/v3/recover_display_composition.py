"""Recover five measured head fields omitted from the historical atlas input.

The source measured these fields; the historical nine-field aggregator did not
carry them forward. Preserve all head records and derive explicitly versioned,
equal-photo observation means. This does not change frozen v2 atlas results.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import sqlite3

from .workflow import ROOT, digest, text_digest

HEAD_SHA = "b3d13ba7014100d4e44bf63d7732e21361c863db949f304e513be140a5011aa8"
FIELDS = ["corolla_visible_fraction", "corolla_white_fraction", "corolla_redmagenta_fraction", "corolla_purple_fraction", "corolla_yellow_fraction"]
ENDPOINTS = ["visible_floret_fraction", "corolla_white_pixel_fraction", "corolla_redmagenta_pixel_fraction", "corolla_purple_pixel_fraction", "corolla_yellow_pixel_fraction"]
COMPOSITION_TOLERANCE = 1e-6


def finite(value):
    try:
        number = float(value)
        return number if math.isfinite(number) else None
    except (ValueError, TypeError):
        return None


def classify(row):
    values = [finite(row[field]) for field in FIELDS]
    colour_ok = row["colour_status"] == "usable"
    visible_ok = colour_ok and values[0] is not None and 0 <= values[0] <= 1
    complete = all(v is not None for v in values[1:])
    bounded = complete and all(0 <= v <= 1 for v in values[1:])
    closed = bounded and abs(sum(values[1:]) - 1) <= COMPOSITION_TOLERANCE
    return values, int(visible_ok), int(colour_ok and closed), int(colour_ok and not closed)


def recover(heads, out, expected_sha=HEAD_SHA, expected_rows=1255791):
    out = out.resolve()
    if out == ROOT or (ROOT in out.parents and not any(out.is_relative_to(ROOT / p) for p in ("local_data", "outputs"))):
        raise ValueError("Use an external or ignored local output directory")
    if out.exists():
        raise ValueError("Output exists; preserve previous executions")
    if digest(heads) != expected_sha:
        raise ValueError("Historical head source identity mismatch")
    out.mkdir(parents=True)
    try:
        with sqlite3.connect(out / "display_composition.sqlite") as db, heads.open(encoding="utf-8-sig", newline="") as handle:
            db.execute("CREATE TABLE heads (source_row INTEGER PRIMARY KEY, head_id TEXT UNIQUE, photo_id TEXT, obs_id TEXT, colour_status TEXT, " +
                       ",".join(f"{field} REAL" for field in FIELDS) + ",visible_eligible INTEGER,composition_eligible INTEGER,composition_qc_conflict INTEGER,raw_fields_json TEXT)")
            reader = csv.DictReader(handle)
            required = set(FIELDS + ["annotation_unit_id", "photo_id", "obs_id", "colour_status"])
            if not required.issubset(reader.fieldnames or []):
                raise ValueError("Required original head fields are absent")
            batch = []
            count = 0
            for count, row in enumerate(reader, 1):
                if None in row or any(v is None for v in row.values()):
                    raise ValueError("Malformed historical head record")
                if not all(row[k] for k in ("annotation_unit_id", "photo_id", "obs_id")):
                    raise ValueError("Missing head/photo/observation identifier")
                values, visible, composition, conflict = classify(row)
                batch.append((count, row["annotation_unit_id"], row["photo_id"], row["obs_id"], row["colour_status"],
                              *values, visible, composition, conflict, json.dumps({f: row[f] for f in FIELDS}, sort_keys=True)))
                if len(batch) == 20000:
                    db.executemany("INSERT INTO heads VALUES (" + ",".join(["?"] * 14) + ")", batch)
                    batch.clear()
            db.executemany("INSERT INTO heads VALUES (" + ",".join(["?"] * 14) + ")", batch)
            if count != expected_rows:
                raise ValueError("Historical head denominator mismatch")
            db.execute("CREATE INDEX head_photo ON heads(obs_id,photo_id)")
            means = [f"AVG(CASE WHEN {'visible' if i == 0 else 'composition'}_eligible=1 THEN {f} END) AS {f}" for i,f in enumerate(FIELDS)]
            db.execute("CREATE TABLE photo_values AS SELECT obs_id,photo_id,COUNT(*) AS n_heads_total,SUM(visible_eligible) AS n_visible_heads,SUM(composition_eligible) AS n_composition_heads," +
                       ",".join(means) + " FROM heads GROUP BY obs_id,photo_id")
            db.execute("CREATE UNIQUE INDEX photo_value ON photo_values(obs_id,photo_id)")
            db.execute("CREATE TABLE observation_values AS SELECT obs_id,COUNT(*) AS n_photos_total,SUM(n_heads_total) AS n_heads_total," +
                       "SUM(n_visible_heads>0) AS n_visible_photos,SUM(n_composition_heads>0) AS n_composition_photos," +
                       ",".join(f"AVG({f}) AS {f}" for f in FIELDS) + " FROM photo_values GROUP BY obs_id")
            scalar = lambda sql: db.execute(sql).fetchone()[0]
            fields = {}
            for endpoint, field in zip(ENDPOINTS, FIELDS):
                fields[endpoint] = {"source_field": field, "finite_head_values": scalar(f"SELECT COUNT(*) FROM heads WHERE {field} IS NOT NULL"),
                                    "finite_colour_qc_usable_heads": scalar(f"SELECT COUNT(*) FROM heads WHERE colour_status='usable' AND {field} IS NOT NULL"),
                                    "recovered_photo_values": scalar(f"SELECT COUNT(*) FROM photo_values WHERE {field} IS NOT NULL"),
                                    "recovered_observation_values": scalar(f"SELECT COUNT(*) FROM observation_values WHERE {field} IS NOT NULL")}
            closure = "+".join(FIELDS[1:])
            max_error = scalar(f"SELECT MAX(ABS(({closure})-1)) FROM observation_values WHERE n_composition_photos>0")
            if max_error is not None and max_error > COMPOSITION_TOLERANCE:
                raise ValueError("Observation composition does not sum to one")
            counts = {"source_heads": count, "retained_head_records": scalar("SELECT COUNT(*) FROM heads"),
                      "source_photos": scalar("SELECT COUNT(*) FROM photo_values"), "source_observations": scalar("SELECT COUNT(*) FROM observation_values"),
                      "colour_qc_usable_heads": scalar("SELECT COUNT(*) FROM heads WHERE colour_status='usable'"),
                      "eligible_visible_heads": scalar("SELECT SUM(visible_eligible) FROM heads"),
                      "eligible_composition_heads": scalar("SELECT SUM(composition_eligible) FROM heads"),
                      "colour_qc_usable_with_invalid_composition": scalar("SELECT SUM(composition_qc_conflict) FROM heads"),
                      "head_records_deleted": 0}
            for table in ("photo_values", "observation_values"):
                with (out / f"{table}.csv").open("w", encoding="utf-8", newline="") as target:
                    cursor = db.execute("SELECT * FROM " + table + " ORDER BY obs_id" + (",photo_id" if table == "photo_values" else ""))
                    writer = csv.writer(target, lineterminator="\n")
                    writer.writerow([d[0] for d in cursor.description])
                    writer.writerows(cursor)
        report = {"status": "FIVE_PREVIOUSLY_MEASURED_FIELDS_RECOVERED_TO_VERSIONED_OBSERVATION_VIEW",
                  "historical_source_artifact": 8269246732, "source_head_sha256": expected_sha,
                  "counts": counts, "endpoints": fields, "composition_max_abs_sum_error": max_error,
                  "implementation_sha256_text_lf": text_digest(Path(__file__)),
                  "aggregation_version": "v3_equal_photo_means_display_composition_v1",
                  "aggregation": "Keep all source rows. For visibility, use colour-QC-usable finite bounded head values. For colour composition, require colour QC plus all four finite bounded fractions summing to one within 1e-6. Average eligible heads within photo, then average photo means within observation with equal photo weight. Keep missing values and support counts; never replace missingness with zero.",
                  "correction": "These five head-level fields were computed, not absent because of unfinished measurement functions. The historical 73_merge_exhaustive_continuous_shards.py TRAITS mapping omitted them from photo/observation aggregation. Their absence from the frozen v2 atlas does not establish that they were never measured.",
                  "new_image_operations": False, "ecological_models_executed": False, "frozen_v2_results_changed": False,
                  "limits": ["This is a new aggregation view, not retroactive execution of the frozen v2 atlas.",
                             "Four colour fractions are one composition, not four independent biological traits; downstream analysis must respect closure.",
                             "Visible floral-pixel fraction is not a validated flowering-stage label, and pixel fractions are not pigment concentrations.",
                             "Original head/photo observation assignments are preserved; recovered shared-photo links are not fanned out into independent trait replicates.",
                             "The complete 27-endpoint v3 measurement and ecological workflow remains to be executed and evaluated."]}
        report["output_sha256"] = {name: digest(out / name) for name in ("display_composition.sqlite", "photo_values.csv", "observation_values.csv")}
        (out / "display_composition_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8", newline="\n")
        return report
    except (Exception, KeyboardInterrupt) as error:
        (out / "incomplete_run.json").write_text(json.dumps({"status": "INCOMPLETE_DO_NOT_USE", "error_type": type(error).__name__, "reason": str(error)}, indent=2) + "\n", encoding="utf-8")
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--heads", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(recover(args.heads, args.out_dir), indent=2))


if __name__ == "__main__":
    main()
