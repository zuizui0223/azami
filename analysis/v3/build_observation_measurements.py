"""Link the saved image measurements to reversible full-source observation views.

Do not read coordinates, taxonomy, environments or perturbation outcomes. Keep
the source universe separate from which endpoints have usable image measurements.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sqlite3

from .image_features import registry, REGISTRY
from .workflow import ROOT, canonical_digest, digest, text_digest

COMPOSITION=[f"corolla_{part}_pixel_fraction" for part in ("white","redmagenta","purple","yellow")]
HUE=["corolla_hue_sin","corolla_hue_cos"]
SPECIFICATION={
    "version":"v3_observation_image_means_v1",
    "source_universe":"every reconciled observation, including those without a cached or usable image; upstream immutable ledgers retain all photo links",
    "head_to_image":"arithmetic mean of endpoint-eligible heads within an exact image; separately retain the mean and count of all finite values including QC failures",
    "joint_eligibility":"hue sine/cosine share one finite usable head set; the four colour parts share one usable bounded closed-composition head set; raw eligibility remains recorded",
    "image_identity":"EXIF-oriented decoded-pixel identity; multiple byte encodings must have identical processing summaries or the build stops; photo-ID aliases of the same pixels receive one image weight within an observation",
    "photo_versions":"retain every version link. A photo ID with multiple distinct pixel versions has unresolved version choice and contributes no selected image through that photo ID; no first, best-QC or latest version is silently selected",
    "image_to_observation":"equal weight per distinct eligible image pixel identity; retain support counts and missing means. Different photographs are not demonstrated same-head replicates",
    "shared_observations":"retain shared-image links and dependence-component IDs; shared images are not independent observation evidence or resolved taxonomic assignments",
    "quality_covariates":"endpoint-matched means of head minimum pixel dimension and Laplacian variance; if any eligible head lacks a covariate, that covariate mean is unavailable rather than a different subset mean",
    "colour_pairs":"both registered colour endpoints must be usable, uniform floral pixels present, and context support available with finite lightness/chroma. Compare legacy colour, uniform-floral colour and context using exactly the same heads and images; retain losses separately",
    "colour_contrast":"uniform floral minus context is an uncalibrated paired image contrast, not an illumination correction or physiological measurement",
    "missingness":"an unmeasured endpoint is never zero; unavailable bytes, detector negatives, processing errors and endpoint QC remain distinguishable upstream and in support tables",
    "claim_limit":"descriptive image-feature aggregation; not full-source image execution, a selected ecological cohort, a new independent test sample or a biological variance estimate",
    "design_status":"specified after v2 and cached baseline inspection, before these observation means; no preregistration claim",
}


def finite(value):
    return value is not None and isinstance(value,(int,float)) and math.isfinite(value)


def mean(values):
    values=[float(v) for v in values if finite(v)]
    return math.fsum(values)/len(values) if values else None


def eligible_endpoints(endpoints):
    if len(endpoints)!=27 or set(endpoints)!={r["endpoint_id"] for r in registry()}:
        raise ValueError("A head is missing its registered endpoint slots")
    ok={key:row["status"]=="usable" and finite(row["value"]) for key,row in endpoints.items()}
    if not all(ok[key] for key in HUE):
        for key in HUE:
            ok[key]=False
    complete=all(ok[key] and 0<=endpoints[key]["value"]<=1 for key in COMPOSITION)
    if not complete or abs(math.fsum(endpoints[key]["value"] for key in COMPOSITION)-1)>1e-6:
        for key in COMPOSITION:
            ok[key]=False
    return ok


def aggregate_heads(heads):
    eligibility=[eligible_endpoints(h["endpoints"]) for h in heads]
    endpoints=[]
    for spec in registry():
        key=spec["endpoint_id"]
        all_values=[h["endpoints"][key]["value"] for h in heads]
        selected=[h for h,ok in zip(heads,eligibility) if ok[key]]
        values=[h["endpoints"][key]["value"] for h in selected]
        quality={}
        for name,source in (("head_min_dimension_mean_px","head_min_dimension_px"),("head_sharpness_mean","head_laplacian_variance")):
            numbers=[h["diagnostics"].get(source) for h in selected]
            quality[name]=mean(numbers) if numbers and all(finite(v) for v in numbers) else None
        original_ok=sum(h["endpoints"][key]["status"]=="usable" and finite(h["endpoints"][key]["value"]) for h in heads)
        endpoints.append({"endpoint_id":key,"unit":spec["unit"],"n_heads":len(heads),"n_finite_heads":sum(finite(v) for v in all_values),
                          "n_qc_usable_heads":original_ok,"n_eligible_heads":len(selected),"n_joint_blocked_heads":original_ok-len(selected),
                          "raw_mean":mean(all_values),"eligible_mean":mean(values),**quality})
    colour_heads=[h for h,ok in zip(heads,eligibility) if ok["corolla_lab_chroma"] and ok["corolla_lab_lightness"]]
    pairs=[]
    for context in ("non_head_context","green_non_head_context"):
        selected=[]
        for head in colour_heads:
            paired=head["diagnostics"].get("paired_colour",{})
            flower=paired.get("floral_union",{})
            background=paired.get(context,{})
            if flower.get("n_pixels",0)>0 and background.get("support_status")=="available" and all(finite(part.get(stat)) for part in (flower,background) for stat in ("lab_chroma","lab_lightness")):
                selected.append((head,flower,background))
        for statistic in ("lab_chroma","lab_lightness"):
            pairs.append({"context_kind":context,"statistic":statistic,"n_eligible_colour_heads":len(colour_heads),"n_pairs":len(selected),
                          "legacy_mean":mean([h["endpoints"]["corolla_"+statistic]["value"] for h,f,b in selected]),
                          "floral_mean":mean([f[statistic] for h,f,b in selected]),"context_mean":mean([b[statistic] for h,f,b in selected]),
                          "contrast_mean":mean([f[statistic]-b[statistic] for h,f,b in selected])})
    return endpoints,pairs


def equivalent_pixel_summaries(rows):
    """Never choose a favourable processing result among identical pixels."""
    if len({canonical_digest(row) for row in rows})>1:
        raise ValueError("Identical pixels have inconsistent processing summaries")


def build(workspace,detection,measurement,dependence,out):
    workspace,detection,measurement,dependence,out=[p.resolve() for p in (workspace,detection,measurement,dependence,out)]
    if out==ROOT or (ROOT in out.parents and not any(out.is_relative_to(ROOT/p) for p in ("local_data","outputs"))):
        raise ValueError("Use an external or ignored local output directory")
    if out.exists():
        raise ValueError("Preserve earlier aggregation versions")
    roots=[workspace,detection,measurement,dependence]
    reports=[json.loads((p/name).read_text(encoding="utf-8")) for p,name in zip(roots,["image_workspace_report.json","detection_report.json","measurement_report.json","dependence_groups_report.json"])]
    paths=[p/name for p,name in zip(roots,["image_workspace.sqlite","detection.sqlite","measurements.sqlite","dependence_groups.sqlite"])]
    expected=[reports[i][key] for i,key in enumerate(["workspace_database_sha256","detection_database_sha256","measurement_database_sha256","output_database_sha256"])]
    for path,sha in zip(paths,expected):
        if digest(path)!=sha:
            raise ValueError("Aggregation input identity changed: "+path.name)
    if reports[1]["job_states"].get("pending") or reports[2]["job_states"].get("pending"):
        raise ValueError("Complete the saved detection/measurement schedule first")
    if reports[1]["execution_contract"]["source_workspace_sha256"]!=expected[0] or reports[2]["execution_contract"]["detector_database_sha256"]!=expected[1] or reports[3]["source_workspace_sha256"]!=expected[0]:
        raise ValueError("Aggregation provenance chain differs")
    if reports[2]["execution_contract"]["feature_specification"]["source_sha256_text_lf"][REGISTRY.relative_to(ROOT).as_posix()]!=text_digest(REGISTRY):
        raise ValueError("Endpoint registry changed from measured version")
    out.mkdir(parents=True)
    context={"specification":SPECIFICATION,"input_database_sha256":dict(zip(["workspace","detection","measurement","dependence"],expected)),
             "implementation_sha256_text_lf":text_digest(Path(__file__)),"registry_sha256_text_lf":text_digest(REGISTRY)}
    (out/"execution_contract.json").write_text(json.dumps(context,indent=2)+"\n",encoding="utf-8",newline="\n")
    try:
        with sqlite3.connect(out/"observation_measurements.sqlite",uri=True) as db:
            for alias,path in zip(("source","det","meas","dep"),paths):
                db.execute("ATTACH DATABASE ? AS "+alias,(path.as_uri()+"?mode=ro",))
            db.executescript("""
                CREATE TABLE image_objects AS SELECT o.sha256,o.decoded_rgb_sha256 AS pixel_id,j.status AS detection_status FROM source.objects o LEFT JOIN det.jobs j ON o.sha256=j.sha256;
                CREATE UNIQUE INDEX image_object ON image_objects(sha256);
                CREATE INDEX pixel_objects ON image_objects(pixel_id);
                CREATE TABLE photo_versions AS SELECT v.photo_id,v.sha256,o.pixel_id FROM source.photo_versions v JOIN image_objects o ON v.sha256=o.sha256;
                CREATE TABLE photo_summary AS SELECT photo_id,COUNT(*) AS n_byte_versions,COUNT(DISTINCT pixel_id) AS n_pixel_versions,CASE WHEN COUNT(DISTINCT pixel_id)=1 THEN MIN(pixel_id) END AS selected_pixel_id FROM photo_versions GROUP BY photo_id;
                CREATE UNIQUE INDEX photo_id ON photo_summary(photo_id);
                CREATE TABLE cached_source_links AS SELECT s.obs_id,s.photo_id FROM source.source_links s JOIN photo_summary p ON s.photo_id=p.photo_id;
                CREATE TABLE observation_images AS SELECT DISTINCT s.obs_id,p.selected_pixel_id AS pixel_id FROM cached_source_links s JOIN photo_summary p ON s.photo_id=p.photo_id WHERE p.selected_pixel_id IS NOT NULL;
                CREATE UNIQUE INDEX observation_image ON observation_images(obs_id,pixel_id);
                CREATE INDEX image_observations ON observation_images(pixel_id);
                CREATE TABLE shared_pixels AS SELECT pixel_id,COUNT(*) AS n_observations FROM observation_images GROUP BY pixel_id;
                CREATE UNIQUE INDEX shared_pixel ON shared_pixels(pixel_id);
                CREATE TABLE source_photo_counts AS SELECT obs_id,COUNT(*) AS n_source_photos FROM source.source_links GROUP BY obs_id;
                CREATE UNIQUE INDEX source_photo_count ON source_photo_counts(obs_id);
                CREATE TABLE observations AS SELECT o.obs_id,o.component_id,o.fold,COALESCE(s.n_source_photos,0) AS n_source_photos FROM dep.observations o LEFT JOIN source_photo_counts s ON o.obs_id=s.obs_id;
                CREATE UNIQUE INDEX observation_id ON observations(obs_id);
                CREATE TABLE head_jobs AS SELECT head_id,sha256,status FROM meas.jobs;
                CREATE INDEX head_image ON head_jobs(sha256);
                CREATE TABLE endpoint_registry(endpoint_id TEXT PRIMARY KEY,unit TEXT);
                CREATE TABLE image_endpoints(sha256 TEXT,endpoint_id TEXT,unit TEXT,n_heads INTEGER,n_finite_heads INTEGER,n_qc_usable_heads INTEGER,n_eligible_heads INTEGER,n_joint_blocked_heads INTEGER,raw_mean REAL,eligible_mean REAL,head_min_dimension_mean_px REAL,head_sharpness_mean REAL,PRIMARY KEY(sha256,endpoint_id));
                CREATE TABLE image_pairs(sha256 TEXT,context_kind TEXT,statistic TEXT,n_eligible_colour_heads INTEGER,n_pairs INTEGER,legacy_mean REAL,floral_mean REAL,context_mean REAL,contrast_mean REAL,PRIMARY KEY(sha256,context_kind,statistic));
                CREATE TABLE image_summaries(sha256 TEXT PRIMARY KEY,pixel_id TEXT,processing_fingerprint TEXT,n_heads INTEGER,n_failed_head_jobs INTEGER);
            """)
            db.executemany("INSERT INTO endpoint_registry VALUES(?,?)",[(r["endpoint_id"],r["unit"]) for r in registry()])
            if db.execute("SELECT COUNT(*) FROM image_objects WHERE detection_status IS NULL").fetchone()[0]:
                raise ValueError("Cached image lacks its detection job")
            if db.execute("SELECT COUNT(*) FROM photo_versions").fetchone()[0]!=db.execute("SELECT COUNT(*) FROM source.photo_versions").fetchone()[0]:
                raise ValueError("A cached photo version has no image object")
            if db.execute("SELECT COUNT(*) FROM head_jobs h LEFT JOIN det.detections d ON h.head_id=d.head_id WHERE d.head_id IS NULL OR d.sha256<>h.sha256").fetchone()[0] or db.execute("SELECT COUNT(*) FROM det.detections d LEFT JOIN head_jobs h ON h.head_id=d.head_id WHERE h.head_id IS NULL").fetchone()[0]:
                raise ValueError("Detection/measurement head membership differs")
            for sha,pixel,det_status in db.execute("SELECT sha256,pixel_id,detection_status FROM image_objects ORDER BY sha256"):
                heads=[]
                for head_id,state,payload in db.execute("SELECT h.head_id,h.status,d.raw_and_diagnostics_json FROM head_jobs h LEFT JOIN meas.details d ON h.head_id=d.head_id WHERE h.sha256=? ORDER BY h.head_id",(sha,)):
                    endpoints={key:{"value":value,"status":status} for key,value,status in db.execute("SELECT endpoint_id,value,status FROM meas.endpoints WHERE head_id=?",(head_id,))}
                    heads.append({"endpoints":endpoints,"diagnostics":json.loads(payload).get("diagnostics",{}) if payload else {},"status":state})
                endpoints,pairs=aggregate_heads(heads)
                failures=sum(h["status"]!="measured" for h in heads)
                fingerprint=canonical_digest({"detection_status":det_status,"n_failed_head_jobs":failures,"endpoints":endpoints,"pairs":pairs})
                db.execute("INSERT INTO image_summaries VALUES(?,?,?,?,?)",(sha,pixel,fingerprint,len(heads),failures))
                db.executemany("INSERT INTO image_endpoints VALUES ("+",".join("?" for _ in range(12))+")",[(sha,*row.values()) for row in endpoints])
                db.executemany("INSERT INTO image_pairs VALUES ("+",".join("?" for _ in range(9))+")",[(sha,*row.values()) for row in pairs])
            if db.execute("SELECT COUNT(*) FROM (SELECT pixel_id FROM image_summaries GROUP BY pixel_id HAVING COUNT(DISTINCT processing_fingerprint)>1)").fetchone()[0]:
                raise ValueError("Identical pixels have inconsistent processing summaries")
            db.executescript("""
                CREATE TABLE pixel_representatives AS SELECT pixel_id,MIN(sha256) AS sha256 FROM image_objects GROUP BY pixel_id;
                CREATE UNIQUE INDEX pixel_representative ON pixel_representatives(pixel_id);
                CREATE VIEW pixel_endpoints AS SELECT p.pixel_id,e.* FROM pixel_representatives p JOIN image_endpoints e ON p.sha256=e.sha256;
                CREATE VIEW pixel_pairs AS SELECT p.pixel_id,e.* FROM pixel_representatives p JOIN image_pairs e ON p.sha256=e.sha256;
                CREATE TABLE observation_support AS SELECT o.obs_id,COUNT(i.pixel_id) AS n_selected_pixel_images,COALESCE(SUM(s.n_observations>1),0) AS n_shared_pixel_images,
                    COALESCE(SUM(m.detection_status='no_detection'),0) AS n_detector_negative_images,COALESCE(SUM(m.detection_status='error'),0) AS n_detector_failed_images,
                    COALESCE(SUM(h.n_failed_head_jobs),0) AS n_failed_head_jobs
                    FROM observations o LEFT JOIN observation_images i ON o.obs_id=i.obs_id LEFT JOIN shared_pixels s ON i.pixel_id=s.pixel_id
                    LEFT JOIN pixel_representatives p ON i.pixel_id=p.pixel_id LEFT JOIN image_objects m ON p.sha256=m.sha256 LEFT JOIN image_summaries h ON p.sha256=h.sha256 GROUP BY o.obs_id;
                CREATE UNIQUE INDEX observation_support_id ON observation_support(obs_id);
                CREATE TABLE observation_endpoints AS SELECT i.obs_id,e.endpoint_id,COUNT(*) AS n_selected_images,SUM(e.n_heads) AS n_detected_heads,SUM(e.n_eligible_heads) AS n_eligible_heads,SUM(e.eligible_mean IS NOT NULL) AS n_eligible_images,SUM(e.raw_mean IS NOT NULL) AS n_finite_images,AVG(e.eligible_mean) AS value,AVG(e.raw_mean) AS raw_mean,
                    CASE WHEN SUM(e.eligible_mean IS NOT NULL)=SUM(e.eligible_mean IS NOT NULL AND e.head_min_dimension_mean_px IS NOT NULL) THEN AVG(CASE WHEN e.eligible_mean IS NOT NULL THEN e.head_min_dimension_mean_px END) END AS head_min_dimension_mean_px,
                    CASE WHEN SUM(e.eligible_mean IS NOT NULL)=SUM(e.eligible_mean IS NOT NULL AND e.head_sharpness_mean IS NOT NULL) THEN AVG(CASE WHEN e.eligible_mean IS NOT NULL THEN e.head_sharpness_mean END) END AS head_sharpness_mean
                    FROM observation_images i JOIN pixel_endpoints e ON i.pixel_id=e.pixel_id GROUP BY i.obs_id,e.endpoint_id;
                CREATE UNIQUE INDEX observation_endpoint ON observation_endpoints(obs_id,endpoint_id);
                CREATE TABLE observation_pairs AS SELECT i.obs_id,e.context_kind,e.statistic,COUNT(*) AS n_selected_images,SUM(e.n_eligible_colour_heads) AS n_eligible_colour_heads,SUM(e.n_pairs) AS n_paired_heads,SUM(e.n_pairs>0) AS n_paired_images,AVG(e.legacy_mean) AS legacy_mean,AVG(e.floral_mean) AS floral_mean,AVG(e.context_mean) AS context_mean,AVG(e.contrast_mean) AS contrast_mean
                    FROM observation_images i JOIN pixel_pairs e ON i.pixel_id=e.pixel_id GROUP BY i.obs_id,e.context_kind,e.statistic;
                CREATE VIEW observation_endpoint_inventory AS SELECT o.obs_id,o.component_id,o.fold,o.n_source_photos,r.endpoint_id,r.unit,s.n_selected_pixel_images,s.n_shared_pixel_images,s.n_detector_negative_images,s.n_detector_failed_images,s.n_failed_head_jobs,COALESCE(e.n_detected_heads,0) AS n_detected_heads,COALESCE(e.n_eligible_heads,0) AS n_eligible_heads,COALESCE(e.n_eligible_images,0) AS n_eligible_images,e.value,e.raw_mean,e.head_min_dimension_mean_px,e.head_sharpness_mean,
                    CASE WHEN s.n_selected_pixel_images=0 THEN 'no_selected_cached_image' WHEN s.n_detector_negative_images=s.n_selected_pixel_images THEN 'detector_negative_images_only' WHEN e.n_detected_heads=0 THEN 'no_measured_head_in_selected_images' WHEN e.value IS NULL THEN 'no_eligible_endpoint' ELSE 'image_measurement_available_not_ecologically_admitted' END AS measurement_support
                    FROM observations o CROSS JOIN endpoint_registry r JOIN observation_support s ON o.obs_id=s.obs_id LEFT JOIN observation_endpoints e ON o.obs_id=e.obs_id AND r.endpoint_id=e.endpoint_id;
            """)
            scalar=lambda sql:db.execute(sql).fetchone()[0]
            counts={"source_observations_retained":scalar("SELECT COUNT(*) FROM observations"),"source_photo_links_referenced":scalar("SELECT SUM(n_source_photos) FROM observations"),
                    "cached_objects_retained":scalar("SELECT COUNT(*) FROM image_objects"),"cached_photo_version_links_retained":scalar("SELECT COUNT(*) FROM photo_versions"),
                    "photo_ids_with_unresolved_pixel_version_choice":scalar("SELECT COUNT(*) FROM photo_summary WHERE n_pixel_versions>1"),"distinct_pixel_images":scalar("SELECT COUNT(*) FROM pixel_representatives"),
                    "image_endpoint_rows":scalar("SELECT COUNT(*) FROM image_endpoints"),"observation_endpoint_inventory_slots":scalar("SELECT COUNT(*) FROM observations")*27,
                    "observations_with_selected_cached_images":scalar("SELECT COUNT(*) FROM observation_support WHERE n_selected_pixel_images>0"),"observations_with_shared_selected_images":scalar("SELECT COUNT(*) FROM observation_support WHERE n_shared_pixel_images>0"),
                    "observations_with_any_eligible_endpoint":scalar("SELECT COUNT(DISTINCT obs_id) FROM observation_endpoints WHERE value IS NOT NULL"),
                    "observations_with_paired_green_context":scalar("SELECT COUNT(DISTINCT obs_id) FROM observation_pairs WHERE context_kind='green_non_head_context' AND n_paired_images>0"),
                    "source_observations_deleted":0}
            if counts["source_observations_retained"]!=reports[3]["counts"]["observations_retained"] or counts["source_photo_links_referenced"]!=reports[3]["counts"]["source_links_retained"] or counts["image_endpoint_rows"]!=counts["cached_objects_retained"]*27:
                raise ValueError("Full-source aggregation denominator changed")
            coverage={key:{"observations":n,"observation_image_links":images} for key,n,images in db.execute("SELECT endpoint_id,SUM(value IS NOT NULL),SUM(n_eligible_images) FROM observation_endpoints GROUP BY endpoint_id")}
            for key,n in db.execute("SELECT endpoint_id,SUM(eligible_mean IS NOT NULL) FROM pixel_endpoints GROUP BY endpoint_id"):
                coverage[key]["unique_pixel_images"]=n
            db.commit()
        report={"status":"FULL_SOURCE_LINKED_IMAGE_MEASUREMENT_VIEWS_BUILT_NOT_ECOLOGICALLY_ADMITTED","execution_contract":context,"counts":counts,
                "endpoint_coverage":coverage,"output_database_sha256":digest(out/"observation_measurements.sqlite"),"ecological_models_executed":False,
                "limits":["Cached historical images are not a representative sample of the full source universe.","Shared-image observation values remain linked through dependence components and require explicit treatment in the eventual ecological cohort.",
                          "Repeated photographs do not establish same-head identity or separate photographic from biological variance.","Paired colour contrasts are descriptive and uncalibrated; context may contain non-leaf material or undetected flowers.",
                          "Raw finite means include QC failures and must not be mistaken for eligible ecological measurements.","All 27 columns remain visible; joint hue and the closed composition are not independent biological traits."]}
        (out/"observation_measurements_report.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8",newline="\n")
        print(json.dumps(report,indent=2))
        return report
    except Exception as error:
        (out/"incomplete_run.json").write_text(json.dumps({"status":"INCOMPLETE_OBSERVATION_MEASUREMENT_BUILD","error_type":type(error).__name__,"error":str(error)},indent=2)+"\n",encoding="utf-8",newline="\n")
        raise


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ("workspace","detection","measurement","dependence","out-dir"):
        parser.add_argument("--"+name,type=Path,required=True)
    args=parser.parse_args()
    build(args.workspace,args.detection,args.measurement,args.dependence,args.out_dir)


if __name__=="__main__":
    main()
