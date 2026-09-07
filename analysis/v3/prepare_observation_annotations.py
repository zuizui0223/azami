"""Prepare reversible full-source observation annotations, without trait joins.

Original photo rows and raw API records remain independently traceable. Prefer
available archived API fields over the collector's flattened metadata, whose
boolean defaults lost the distinction between absent and false. Never infer
private coordinates, developmental stage, native status or ecological effects.
"""
from __future__ import annotations

import argparse
import calendar
import csv
from datetime import date
import gzip
import hashlib
import io
from itertools import groupby
import json
import math
from pathlib import Path
import re
import sqlite3
import zipfile

from .workflow import ROOT, canonical_digest, digest, text_digest

FIELDS = ["taxon_id", "taxon_name", "taxon_rank", "observed_on", "created_at", "updated_at",
          "quality_grade", "captive", "user_id", "latitude", "longitude", "positional_accuracy",
          "geoprivacy", "obscured", "coordinate_usable_for_environment", "inat_flowers_annotation",
          "inat_annotation_summary", "annotation_records", "observation_license_code"]
SPECIFICATION = {
    "version":"v3_full_source_observation_annotations_v1",
    "source":"all six pinned original metadata chunks and their available raw API records, not the thinned v2 cohort",
    "retention":"retain every source row locator and version; one annotation row per reconciled observation; no photo, head or observation deletion",
    "preference":"use archived API fields when any API version exists for the observation, otherwise original flattened metadata; never fill API-missing fields silently from flattened defaults",
    "conflicts":"all versions remain saved; a field differing among preferred-source versions becomes unavailable with an explicit conflict flag; no last-row or most-complete-row choice",
    "calendar":"exact YYYY-MM-DD only; phase 2*pi*(DOY-1)/days_in_observed_year; retain year, DOY, sine and cosine; add southern indicator and its sine/cosine interactions only with usable public latitude",
    "hemisphere":"north for latitude>0, south for latitude<0, equatorial for exactly zero; missing or restricted coordinates give unknown; these are calendar terms, not developmental-stage labels",
    "coordinates":"use only public API geojson/location fields or archived flattened public fields; private_* fields are never read; restricted, invalid, unknown-privacy and known-conflicting fields block derived analysis coordinates",
    "precision":"retain reported positional accuracy and distinguish missing, invalid, zero and positive; model-resolution eligibility is not decided here",
    "native_status":"not_assessed_full_source; never inherit native status from a photo, genus or absence of an introduced record",
    "separation":"do not load image measurements, ecological predictors or model results; these annotations do not authorize an ecological cohort or demonstrate confounding removal",
    "design_status":"retrospective v3 preparation after v2 results were seen; no preregistration claim",
}


def text(value):
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value).strip()


def boolean(value):
    token=text(value).lower()
    return "true" if token in ("true","1","yes") else "false" if token in ("false","0","no") else "unknown"


def number(value):
    try:
        result=float(value)
        return result if math.isfinite(result) else None
    except (TypeError,ValueError):
        return None


def csv_fields(row):
    fields={key:text(row.get(key)) for key in FIELDS}
    for key in ("captive","obscured","coordinate_usable_for_environment","inat_flowers_annotation"):
        fields[key]=boolean(row.get(key))
    return fields


def api_fields(record):
    taxon=record.get("taxon") or {}
    user=record.get("user") or {}
    if not isinstance(taxon,dict) or not isinstance(user,dict):
        raise ValueError("Invalid API taxon/user structure")
    fields={key:"" for key in FIELDS}
    for key in ("observed_on","created_at","updated_at","quality_grade","positional_accuracy","geoprivacy"):
        fields[key]=text(record.get(key))
    for key in ("captive","obscured"):
        fields[key]=boolean(record.get(key))
    for key in ("id","name","rank"):
        fields["taxon_"+key]=text(taxon.get(key))
    fields["user_id"]=text(user.get("id"))
    fields["observation_license_code"]=text(record.get("license_code"))
    geo=record.get("geojson") or {}
    coords=geo.get("coordinates") if isinstance(geo,dict) else None
    if isinstance(coords,list) and len(coords)>=2:
        fields["longitude"],fields["latitude"]=text(coords[0]),text(coords[1])
    else:
        location=text(record.get("location"))
        if location.count(",")==1:
            fields["latitude"],fields["longitude"]=[v.strip() for v in location.split(",")]
    fields["coordinate_usable_for_environment"]="not_a_raw_api_field"
    fields["inat_flowers_annotation"]="not_reclassified_from_raw_annotations"
    fields["annotation_records"]=json.dumps(record.get("annotations"),ensure_ascii=False,sort_keys=True,separators=(",",":"))
    return fields


def consensus(versions):
    preferred="api" if any(kind=="api" for kind,_ in versions) else "metadata"
    chosen=[fields for kind,fields in versions if kind==preferred]
    if not chosen:
        return "missing",0,{},list(FIELDS)
    conflicts=[key for key in FIELDS if len({row.get(key,"") for row in chosen})>1]
    values={key:("" if key in conflicts else chosen[0].get(key,"")) for key in FIELDS}
    return preferred,len(chosen),values,conflicts


def annotate(fields,conflicts,kind):
    result={"date_status":"missing","observed_year":None,"doy":None,"days_in_year":None,
            "sin_doy":None,"cos_doy":None,"hemisphere":"unknown","south_indicator":None,
            "south_sin":None,"south_cos":None,"coordinate_status":"missing_or_invalid",
            "analysis_latitude":None,"analysis_longitude":None,"position_accuracy_status":"missing",
            "position_accuracy_m":None,"captive_state":boolean(fields.get("captive")),
            "native_status":"not_assessed_full_source","source_taxon_rank":fields.get("taxon_rank", ""),
            "source_taxon_name":fields.get("taxon_name", "")}
    raw=fields.get("observed_on","")
    if "observed_on" in conflicts:
        result["date_status"]="conflicting_source_versions"
    elif raw:
        result["date_status"]="invalid_or_not_exact_day"
        if re.fullmatch(r"\d{4}-\d{2}-\d{2}",raw):
            try:
                day=date.fromisoformat(raw)
                doy=day.timetuple().tm_yday
                days=366 if calendar.isleap(day.year) else 365
                phase=2*math.pi*(doy-1)/days
                result.update(date_status="exact_day",observed_year=day.year,doy=doy,days_in_year=days,sin_doy=math.sin(phase),cos_doy=math.cos(phase))
            except ValueError:
                pass
    accuracy=number(fields.get("positional_accuracy"))
    if "positional_accuracy" in conflicts:
        result["position_accuracy_status"]="conflicting_source_versions"
    elif fields.get("positional_accuracy"):
        result["position_accuracy_status"]="invalid" if accuracy is None or accuracy<0 else "reported_zero" if accuracy==0 else "reported_positive"
        if accuracy is not None and accuracy>=0:
            result["position_accuracy_m"]=accuracy
    lat,lon=number(fields.get("latitude")),number(fields.get("longitude"))
    privacy=fields.get("geoprivacy","").lower()
    obscured=boolean(fields.get("obscured"))
    sensitive=["latitude","longitude","geoprivacy","obscured","coordinate_usable_for_environment"]
    if any(key in conflicts for key in sensitive):
        result["coordinate_status"]="conflicting_source_versions"
    elif privacy in ("private","obscured") or obscured=="true":
        result["coordinate_status"]="restricted"
    elif privacy not in ("","open") or obscured=="unknown":
        result["coordinate_status"]="privacy_unknown"
    elif kind=="metadata" and boolean(fields.get("coordinate_usable_for_environment"))!="true":
        result["coordinate_status"]="source_flag_not_usable"
    elif lat is not None and lon is not None and -90<=lat<=90 and -180<=lon<=180:
        hemisphere="north" if lat>0 else "south" if lat<0 else "equatorial"
        south=1 if lat<0 else 0
        result.update(coordinate_status="public_location_present_precision_not_gated",analysis_latitude=lat,analysis_longitude=lon,hemisphere=hemisphere,south_indicator=south)
        if result["sin_doy"] is not None:
            result.update(south_sin=south*result["sin_doy"],south_cos=south*result["cos_doy"])
    if "captive" in conflicts:
        result["captive_state"]="conflicting_source_versions"
    return result


def save_batch(db,rows):
    db.executemany("INSERT INTO source_records VALUES (?,?,?,?,?,?)",[(kind,origin,line,obs,canonical_digest(fields),sha) for kind,origin,line,obs,fields,sha in rows])
    db.executemany("INSERT OR IGNORE INTO versions VALUES (?,?,?,?)",[(obs,kind,canonical_digest(fields),json.dumps(fields,ensure_ascii=False,sort_keys=True)) for kind,origin,line,obs,fields,sha in rows])


def ingest(db,archive,spec,chunk):
    if digest(archive)!=spec["archive_sha256"] or archive.stat().st_size!=spec["archive_bytes"]:
        raise ValueError("Original archive identity changed")
    origin=str(spec["artifact_id"])
    counts={"metadata":0,"api":0}
    with zipfile.ZipFile(archive) as zipped:
        if len(zipped.namelist())!=len(set(zipped.namelist())):
            raise ValueError("Duplicate archive members")
        with zipped.open(spec["required_member"]) as handle:
            if hashlib.file_digest(handle,"sha256").hexdigest()!=spec["required_member_sha256"]:
                raise ValueError("Metadata member identity changed")
        with zipped.open(spec["required_member"]) as binary,io.TextIOWrapper(binary,encoding="utf-8-sig",newline="") as handle:
            rows=[]
            for line,row in enumerate(csv.DictReader(handle),1):
                if None in row or any(value is None for value in row.values()) or not row.get("obs_id"):
                    raise ValueError("Malformed metadata record")
                rows.append(("metadata",origin,line,row["obs_id"],csv_fields(row),canonical_digest(row)))
                counts["metadata"]+=1
                if len(rows)==10000:
                    save_batch(db,rows); rows.clear()
            save_batch(db,rows)
        raw=chunk["raw_api_member"]
        if raw:
            fingerprint=hashlib.sha256()
            with zipped.open(raw) as binary:
                handle=gzip.GzipFile(fileobj=binary) if raw.endswith(".gz") else binary
                rows=[]
                for line,payload in enumerate(handle,1):
                    fingerprint.update(payload)
                    record=json.loads(payload)
                    obs=text(record.get("id"))
                    if not obs:
                        raise ValueError("Missing API observation identity")
                    rows.append(("api",origin,line,obs,api_fields(record),hashlib.sha256(payload).hexdigest()))
                    counts["api"]+=1
                    if len(rows)==10000:
                        save_batch(db,rows); rows.clear()
                save_batch(db,rows)
            if fingerprint.hexdigest()!=chunk["raw_api_sha256_uncompressed"]:
                raise ValueError("Raw API stream identity changed")
        if counts["metadata"]!=chunk["metadata_rows"] or counts["api"]!=(chunk["raw_api_observations"] or 0):
            raise ValueError("Source record count changed")
    return {"artifact_id":spec["artifact_id"],**counts}


def build(archives,reconciliation,dependence,out,manifest=None):
    archives,reconciliation,dependence,out=[p.resolve() for p in (archives,reconciliation,dependence,out)]
    if out==ROOT or (ROOT in out.parents and not any(out.is_relative_to(ROOT/p) for p in ("local_data","outputs"))):
        raise ValueError("Use an external or ignored local output directory")
    if out.exists():
        raise ValueError("Preserve earlier annotation versions")
    manifest=manifest or json.loads(Path(__file__).with_name("upstream_sources.json").read_text(encoding="utf-8"))
    source=json.loads((reconciliation/"source_reconciliation_report.json").read_text(encoding="utf-8"))
    groups=json.loads((dependence/"dependence_groups_report.json").read_text(encoding="utf-8"))
    if canonical_digest(manifest)!=source["manifest_sha256_canonical_json"]:
        raise ValueError("Source archive manifest changed")
    if digest(reconciliation/"source_reconciliation.sqlite")!=source["ledger_sha256"] or digest(dependence/"dependence_groups.sqlite")!=groups["output_database_sha256"]:
        raise ValueError("Source or dependence database identity changed")
    out.mkdir(parents=True)
    context={"specification":SPECIFICATION,"implementation_sha256_text_lf":text_digest(Path(__file__)),
             "source_reconciliation_sha256":source["ledger_sha256"],"dependence_groups_sha256":groups["output_database_sha256"],
             "archive_manifest_sha256_canonical_json":canonical_digest(manifest)}
    (out/"execution_contract.json").write_text(json.dumps(context,indent=2)+"\n",encoding="utf-8",newline="\n")
    try:
        with sqlite3.connect(out/"observation_annotations.sqlite",uri=True) as db:
            db.execute("ATTACH DATABASE ? AS source",((reconciliation/"source_reconciliation.sqlite").as_uri()+"?mode=ro",))
            db.execute("ATTACH DATABASE ? AS groups",((dependence/"dependence_groups.sqlite").as_uri()+"?mode=ro",))
            db.executescript("""
                CREATE TABLE observations AS SELECT obs_id,component_id,fold FROM groups.observations;
                CREATE UNIQUE INDEX observation_id ON observations(obs_id);
                CREATE TABLE source_records(kind TEXT,origin TEXT,source_row INTEGER,obs_id TEXT,version_hash TEXT,payload_sha256 TEXT,PRIMARY KEY(kind,origin,source_row));
                CREATE TABLE versions(obs_id TEXT,kind TEXT,version_hash TEXT,fields_json TEXT,PRIMARY KEY(obs_id,kind,version_hash));
            """)
            if db.execute("SELECT COUNT(*) FROM observations o WHERE NOT EXISTS(SELECT 1 FROM source.source_observations s WHERE o.obs_id=s.obs_id)").fetchone()[0] or db.execute("SELECT COUNT(*) FROM source.source_observations s WHERE NOT EXISTS(SELECT 1 FROM observations o WHERE o.obs_id=s.obs_id)").fetchone()[0]:
                raise ValueError("Source/dependence observation membership differs")
            chunks={c["artifact_id"]:c for c in source["chunks"]}
            ingested=[]
            for spec in manifest["archives"]:
                if spec["role"]=="raw_metadata_chunk":
                    ingested.append(ingest(db,archives/str(spec["artifact_id"])/"source.zip",spec,chunks[spec["artifact_id"]]))
                    db.commit()
                    print(json.dumps(ingested[-1]),flush=True)
            if db.execute("SELECT COUNT(*) FROM versions v WHERE NOT EXISTS(SELECT 1 FROM observations o WHERE o.obs_id=v.obs_id)").fetchone()[0]:
                raise ValueError("Unlinked observation metadata")
            checks={
                "metadata_locator_mismatches":db.execute("SELECT COUNT(*) FROM source_records r LEFT JOIN source.original_photo_rows s ON r.origin=s.origin AND r.source_row=s.source_row WHERE r.kind='metadata' AND (s.origin IS NULL OR r.obs_id<>s.obs_id OR r.payload_sha256<>s.payload_sha256)").fetchone()[0],
                "api_locator_mismatches":db.execute("SELECT COUNT(*) FROM source_records r LEFT JOIN source.api_observations s ON r.origin=s.origin AND r.source_row=s.source_line WHERE r.kind='api' AND (s.origin IS NULL OR r.obs_id<>s.obs_id OR r.payload_sha256<>s.raw_line_sha256)").fetchone()[0],
            }
            if any(checks.values()):
                raise ValueError("Source locators or exact record payloads differ")
            derived=list(annotate({},[],"missing"))
            numeric={"observed_year","doy","days_in_year","sin_doy","cos_doy","south_indicator","south_sin","south_cos","analysis_latitude","analysis_longitude","position_accuracy_m"}
            db.execute("CREATE TABLE annotations(obs_id TEXT PRIMARY KEY,component_id TEXT,fold INTEGER,preferred_source_kind TEXT,n_candidate_versions INTEGER,conflicted_fields_json TEXT,selected_fields_json TEXT,"+",".join(key+(" REAL" if key in numeric else " TEXT") for key in derived)+")")
            cursor=db.execute("SELECT o.obs_id,o.component_id,o.fold,v.kind,v.fields_json FROM observations o LEFT JOIN versions v ON o.obs_id=v.obs_id ORDER BY o.obs_id,v.kind,v.version_hash")
            batch=[]
            insert="INSERT INTO annotations VALUES ("+",".join("?" for _ in range(7+len(derived)))+")"
            for obs,records in groupby(cursor,key=lambda row:row[0]):
                records=list(records)
                kind,n,fields,conflicts=consensus([(r[3],json.loads(r[4])) for r in records if r[4] is not None])
                values=annotate(fields,conflicts,kind)
                batch.append((obs,records[0][1],records[0][2],kind,n,json.dumps(conflicts),json.dumps(fields,ensure_ascii=False,sort_keys=True),*[values[key] for key in derived]))
                if len(batch)==10000:
                    db.executemany(insert,batch); batch.clear()
            db.executemany(insert,batch)
            scalar=lambda sql:db.execute(sql).fetchone()[0]
            counts={"observations_retained":scalar("SELECT COUNT(*) FROM observations"),"annotation_rows":scalar("SELECT COUNT(*) FROM annotations"),
                    "source_metadata_rows_retained":scalar("SELECT COUNT(*) FROM source_records WHERE kind='metadata'"),"source_api_rows_retained":scalar("SELECT COUNT(*) FROM source_records WHERE kind='api'"),
                    "source_versions_retained":scalar("SELECT COUNT(*) FROM versions"),"observations_with_conflicted_preferred_fields":scalar("SELECT COUNT(*) FROM annotations WHERE conflicted_fields_json<>'[]'"),
                    "source_rows_deleted":0,"ecological_models_executed":0}
            if counts["annotation_rows"]!=counts["observations_retained"] or counts["source_metadata_rows_retained"]!=source["counts"]["original_metadata_rows"] or counts["source_api_rows_retained"]!=source["counts"]["raw_api_records"]:
                raise ValueError("Annotation denominator changed")
            states={key:dict(db.execute("SELECT "+key+",COUNT(*) FROM annotations GROUP BY "+key)) for key in ("preferred_source_kind","date_status","coordinate_status","position_accuracy_status","hemisphere","captive_state","native_status","source_taxon_rank")}
            db.commit()
        result={"status":"FULL_SOURCE_OBSERVATION_ANNOTATIONS_PREPARED_NOT_MODELLED","execution_contract":context,"counts":counts,"states":states,"integrity_checks":checks,
                "output_database_sha256":digest(out/"observation_annotations.sqlite"),
                "limits":["Public coordinates with reported accuracy are not yet matched to environmental grids or accepted for any spatial-resolution model.",
                          "Calendar harmonics do not identify flowering stage; archived annotations have not been treated as verified stage labels.",
                          "First-chunk raw API records are absent. Its flattened boolean flags may already conflate missing with false; the original record cannot be reconstructed from those flags.",
                          "Source taxon assignment is not a resolved or independently checked botanical identification.",
                          "Native status remains unassessed here, including records previously assessed in a smaller v2 view.",
                          "Separate annotation preparation does not join or inspect image features or ecological associations."]}
        (out/"observation_annotations_report.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8",newline="\n")
        print(json.dumps(result,indent=2))
        return result
    except Exception as error:
        (out/"incomplete_run.json").write_text(json.dumps({"status":"INCOMPLETE_ANNOTATION_PREPARATION","error_type":type(error).__name__,"error":str(error)},indent=2)+"\n",encoding="utf-8",newline="\n")
        raise


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ("archives","reconciliation","dependence","out-dir"):
        parser.add_argument("--"+name,type=Path,required=True)
    args=parser.parse_args()
    build(args.archives,args.reconciliation,args.dependence,args.out_dir)


if __name__=="__main__":
    main()
