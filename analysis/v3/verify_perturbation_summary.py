"""Independently recompute scalar changes and QC attrition through SQLite.

This verifies arithmetic, not physical accuracy. Rank correlations, quantiles and
joint circular/composition summaries are explicitly outside this SQL cross-check.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import sqlite3

from .workflow import ROOT,digest,text_digest

COUNTS=["scheduled_heads","scheduled_images","scheduled_components","baseline_usable_heads","paired_usable_heads","lost_usable_heads","gained_usable_heads","neither_usable_heads","paired_images","paired_components"]
FLOATS=["mean_signed_change","mean_absolute_change","component_weighted_loss_fraction"]
FILTERS={"all_cached":"1=1","recorded_development_exposed":"development_exposed=1","no_recorded_development_exposure":"development_exposed=0"}


def compare(row,actual):
    for field in COUNTS:
        if int(row[field])!=int(actual[field] or 0):
            raise ValueError("Independent count mismatch: "+field)
    errors={}
    for field in FLOATS:
        expected=float(row[field]) if row[field] else None
        got=actual[field]
        if (expected is None)!=(got is None):
            raise ValueError("Independent availability mismatch: "+field)
        if expected is not None:
            if not math.isclose(expected,got,rel_tol=1e-10,abs_tol=1e-10):
                raise ValueError("Independent value mismatch: "+field)
            errors[field]=abs(expected-got)
    return errors


def verify(perturbation,measurement,dependence,summary,out):
    perturbation,measurement,dependence,summary,out=[p.resolve() for p in (perturbation,measurement,dependence,summary,out)]
    if out==ROOT or (ROOT in out.parents and not any(out.is_relative_to(ROOT/p) for p in ("local_data","outputs"))):
        raise ValueError("Use an external or ignored local output directory")
    if out.exists():
        raise ValueError("Preserve earlier verification versions")
    report=json.loads((summary/"technical_sensitivity_summary_report.json").read_text(encoding="utf-8"))
    if report["status"]!="COMPONENT_WEIGHTED_TECHNICAL_SENSITIVITY_SUMMARIZED":
        raise ValueError("A completed summary is required")
    table=summary/"technical_sensitivity_summary.csv"
    paths=[perturbation/"perturbations.sqlite",measurement/"measurements.sqlite",dependence/"dependence_groups.sqlite"]
    for label,path in zip(("perturbation","measurement","dependence"),paths):
        if digest(path)!=report["input_database_sha256"][label]:
            raise ValueError("Verification input identity changed")
    if digest(table)!=report["summary_sha256"]:
        raise ValueError("Summary table identity changed")
    with table.open(encoding="utf-8",newline="") as handle:
        rows=list(csv.DictReader(handle))
    if len(rows)!=1218 or len({(r["condition"],r["metric_id"],r["exposure_stratum"]) for r in rows})!=1218:
        raise ValueError("Summary grid membership changed")
    scalar_rows=[r for r in rows if r["metric_kind"]=="registered_endpoint"]
    if len(scalar_rows)!=1134:
        raise ValueError("Registered scalar grid differs")
    with sqlite3.connect(paths[0].as_uri()+"?mode=ro",uri=True) as db:
        db.execute("ATTACH DATABASE ? AS base",(paths[1].as_uri()+"?mode=ro",))
        db.execute("ATTACH DATABASE ? AS dep",(paths[2].as_uri()+"?mode=ro",))
        db.executescript("""
            CREATE TEMP TABLE image_checks AS
            WITH linked AS (
                SELECT p.condition,p.endpoint_id,j.sha256,j.component_id,p.value AS changed,b.value AS baseline,
                    (p.status='usable' AND p.value IS NOT NULL AND ABS(p.value)<=1.7976931348623157e308) AS pok,
                    (b.status='usable' AND b.value IS NOT NULL AND ABS(b.value)<=1.7976931348623157e308) AS bok
                FROM endpoints p JOIN base.endpoints b ON p.head_id=b.head_id AND p.endpoint_id=b.endpoint_id JOIN jobs j ON p.head_id=j.head_id
            )
            SELECT condition,endpoint_id,sha256,component_id,COUNT(*) AS n_scheduled,SUM(bok) AS n_base,SUM(bok AND pok) AS n_pair,
                SUM(bok AND NOT pok) AS n_lost,SUM(NOT bok AND pok) AS n_gained,SUM(NOT bok AND NOT pok) AS n_neither,
                AVG(CASE WHEN bok AND pok THEN changed-baseline END) AS signed_change,
                AVG(CASE WHEN bok AND pok THEN ABS(changed-baseline) END) AS absolute_change
            FROM linked GROUP BY condition,endpoint_id,sha256,component_id;
            CREATE TEMP TABLE component_checks AS
            SELECT i.condition,i.endpoint_id,i.component_id,d.development_pool_exposed AS development_exposed,
                SUM(n_scheduled) AS n_scheduled,COUNT(*) AS n_images,SUM(n_base) AS n_base,SUM(n_pair) AS n_pair,SUM(n_lost) AS n_lost,SUM(n_gained) AS n_gained,SUM(n_neither) AS n_neither,
                SUM(n_pair>0) AS n_paired_images,AVG(signed_change) AS signed_change,AVG(absolute_change) AS absolute_change,
                AVG(CASE WHEN n_base>0 THEN n_lost*1.0/n_base END) AS loss_fraction
            FROM image_checks i LEFT JOIN dep.components d ON i.component_id=d.component_id GROUP BY i.condition,i.endpoint_id,i.component_id;
            CREATE INDEX component_condition ON component_checks(condition,endpoint_id);
        """)
        if db.execute("SELECT COUNT(*) FROM component_checks WHERE development_exposed IS NULL").fetchone()[0]:
            raise ValueError("Missing component exposure")
        query="""SELECT SUM(n_scheduled) AS scheduled_heads,SUM(n_images) AS scheduled_images,COUNT(*) AS scheduled_components,
            SUM(n_base) AS baseline_usable_heads,SUM(n_pair) AS paired_usable_heads,SUM(n_lost) AS lost_usable_heads,SUM(n_gained) AS gained_usable_heads,SUM(n_neither) AS neither_usable_heads,
            SUM(n_paired_images) AS paired_images,SUM(n_pair>0) AS paired_components,
            AVG(signed_change) AS mean_signed_change,AVG(absolute_change) AS mean_absolute_change,AVG(loss_fraction) AS component_weighted_loss_fraction
            FROM component_checks WHERE condition=? AND endpoint_id=? AND """
        maximum={key:0. for key in FLOATS}
        for row in scalar_rows:
            cursor=db.execute(query+FILTERS[row["exposure_stratum"]],(row["condition"],row["metric_id"]))
            actual=dict(zip([d[0] for d in cursor.description],cursor.fetchone()))
            for field,error in compare(row,actual).items():
                maximum[field]=max(maximum[field],error)
    result={"status":"SCALAR_PERTURBATION_MEANS_AND_ATTRITION_RECOMPUTED_BY_SQL","scalar_rows_checked":len(scalar_rows),
            "checked_fields":COUNTS+FLOATS,"max_absolute_numerical_differences":maximum,"numerical_comparison_tolerance":{"absolute":1e-10,"relative":1e-10},
            "summary_sha256":report["summary_sha256"],"input_database_sha256":report["input_database_sha256"],
            "implementation_sha256_text_lf":text_digest(Path(__file__)),
            "not_independently_recomputed_here":["rank correlations","component quantiles","joint circular hue changes","joint composition distances"],
            "scientific_limit":"Arithmetic verification is not physical accuracy, real-camera uncertainty calibration or evidence of ecological association."}
    out.mkdir(parents=True)
    (out/"summary_verification_report.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8",newline="\n")
    print(json.dumps(result,indent=2))
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ("perturbation","measurement","dependence","summary","out-dir"):
        parser.add_argument("--"+name,type=Path,required=True)
    args=parser.parse_args()
    verify(args.perturbation,args.measurement,args.dependence,args.summary,args.out_dir)


if __name__=="__main__":
    main()
