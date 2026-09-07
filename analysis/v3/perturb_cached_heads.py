"""Execute the saved technical perturbation grid on every cached detected head."""
from __future__ import annotations

import argparse
from collections import deque
from concurrent.futures import ProcessPoolExecutor
from contextlib import contextmanager
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import sqlite3

import cv2
import numpy as np

from . import image_features as features
from .detect_cached_images import crop_box, decode_object
from .measure_cached_heads import failed_endpoints, software_versions
from .workflow import ROOT, digest, text_digest, canonical_digest

CONTRACT = Path(__file__).with_name("perturbation_contract.json")


def pixel_hash(array):
    payload = array.shape[0].to_bytes(8, "big") + array.shape[1].to_bytes(8, "big")
    return hashlib.sha256(payload + np.ascontiguousarray(array).tobytes()).hexdigest()


def transform(image, box, all_boxes, condition):
    image = image.copy()
    box = list(map(float, box))
    others = [list(map(float, b)) for b in all_boxes]
    kind = condition["kind"]
    sx = sy = 1.
    clipping_fraction = 0.
    if kind == "bbox_shift":
        dx, dy = (box[2]-box[0])*condition["dx"], (box[3]-box[1])*condition["dy"]
        box = [box[0]+dx, box[1]+dy, box[2]+dx, box[3]+dy]
    elif kind == "resolution":
        scale = condition["scale"]
        if not 0 < scale <= 1:
            raise ValueError("Resolution probes cannot upsample")
        old_h, old_w = image.shape[:2]
        width, height = max(1, round(old_w*scale)), max(1, round(old_h*scale))
        image = cv2.resize(image, (width, height), interpolation=cv2.INTER_AREA)
        sx, sy = width/old_w, height/old_h
        box = [box[0]*sx, box[1]*sy, box[2]*sx, box[3]*sy]
        others = [[b[0]*sx, b[1]*sy, b[2]*sx, b[3]*sy] for b in others]
    elif kind in ("gamma", "intensity", "channel_gain"):
        values = image.astype(float)
        if kind == "gamma":
            values = 255*(values/255)**condition["exponent"]
        elif kind == "intensity":
            values *= condition["multiplier"]
        else:
            values *= np.asarray(condition["bgr"], dtype=float)
        clipping_fraction = float(((values < 0) | (values > 255)).mean())
        image = np.clip(np.rint(values), 0, 255).astype(np.uint8)
    elif kind == "blur":
        image = cv2.GaussianBlur(image, (0, 0), condition["sigma"], borderType=cv2.BORDER_REFLECT_101)
    elif kind != "identity":
        raise ValueError("Unknown perturbation kind")
    head_extent, head_clipped = crop_box(box, image.shape[1], image.shape[0], .12)
    context_extent, context_clipped = crop_box(box, image.shape[1], image.shape[0], .8)
    all_head_extents = [head_extent]
    for other in others:
        all_head_extents.append(crop_box(other, image.shape[1], image.shape[0], .12)[0])
    hx1, hy1, hx2, hy2 = head_extent
    cx1, cy1, cx2, cy2 = context_extent
    recipe = {"condition": condition, "actual_scale_x": sx, "actual_scale_y": sy,
              "transformed_image_width": image.shape[1], "transformed_image_height": image.shape[0],
              "transformed_bgr_pixel_sha256": pixel_hash(image), "head_box": head_extent, "context_box": context_extent,
              "all_excluded_head_boxes": all_head_extents, "head_clipped": head_clipped,
              "context_clipped": context_clipped, "out_of_range_channel_fraction_before_clipping": clipping_fraction}
    return image[hy1:hy2, hx1:hx2].copy(), image[cy1:cy2, cx1:cx2].copy(), recipe


def initialize_worker():
    cv2.setNumThreads(1)
    cv2.setRNGSeed(20260907)


def run_head(task):
    workspace, obj, head_id, box, other_boxes, baseline, conditions, source_state = task
    records = []
    if source_state not in ("measured", "measured_with_engine_errors"):
        return {"head_id":head_id,"status":"source_not_evaluable","records":[
            {"condition":c["id"],"status":"source_not_evaluable", "source_state":source_state,
             "result":{"endpoints":baseline if c["id"]=="baseline" else failed_endpoints("not_evaluable_saved_baseline")}} for c in conditions]}
    try:
        image = np.asarray(decode_object(Path(workspace), obj))[:, :, ::-1].copy()
    except Exception as error:
        return {"head_id": head_id, "status": "source_error", "error": {"type":type(error).__name__,"message":str(error)}, "records": [
            {"condition": c["id"], "status": "source_error", "error":{"type":type(error).__name__,"message":str(error)}, "result": {"endpoints": failed_endpoints("source_error")}} for c in conditions]}
    for condition in conditions:
        try:
            head, context, recipe = transform(image, box, other_boxes, condition)
            result, masks = features.measure(head, context, recipe["context_box"], recipe["all_excluded_head_boxes"])
            mask_ids = {name:{"pixel_sha256":pixel_hash(mask.astype(np.uint8)*255), "height":mask.shape[0],"width":mask.shape[1]} for name,mask in masks.items()}
            status = "measured_with_engine_errors" if result["engine_errors"] else "measured"
            record = {"condition":condition["id"],"status":status,"recipe":recipe,"mask_identities":mask_ids,"result":result}
            records.append(record)
            if condition["id"] == "baseline" and result["endpoints"] != baseline:
                records[-1]["status"] = "baseline_replay_mismatch"
                records.extend({"condition": c["id"], "status":"not_run_baseline_mismatch", "result":{"endpoints":failed_endpoints("not_run_baseline_mismatch")}} for c in conditions[1:])
                return {"head_id":head_id,"status":"baseline_replay_mismatch","records":records}
        except Exception as error:
            records.append({"condition":condition["id"],"status":"worker_error", "error":{"type":type(error).__name__,"message":str(error)},"result":{"endpoints":failed_endpoints("worker_error")}})
            if condition["id"] == "baseline":
                records.extend({"condition":c["id"],"status":"not_run_baseline_error","result":{"endpoints":failed_endpoints("not_run_baseline_error")}} for c in conditions[1:])
                return {"head_id":head_id,"status":"baseline_error","records":records}
    return {"head_id":head_id,"status":"completed_with_errors" if any(r["status"]!="measured" for r in records) else "completed","records":records}


@contextmanager
def runner_lock(out):
    with (out/"runner.lock").open("a+b") as handle:
        if handle.tell() == 0:
            handle.write(b"0")
            handle.flush()
        handle.seek(0)
        if os.name == "nt":
            import msvcrt
            msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
        else:
            import fcntl
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        try:
            yield
        finally:
            handle.seek(0)
            if os.name == "nt":
                msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                fcntl.flock(handle, fcntl.LOCK_UN)


def run(workspace, detection, measurement, dependence, out, limit=0, workers=2):
    workspace,detection,measurement,dependence,out = [p.resolve() for p in (workspace,detection,measurement,dependence,out)]
    if out == ROOT or (ROOT in out.parents and not any(out.is_relative_to(ROOT/p) for p in ("local_data","outputs"))):
        raise ValueError("Use an external or ignored local output directory")
    if limit < 0 or workers not in (1,2):
        raise ValueError("Use a nonnegative head limit and one or two CPU workers")
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    conditions = contract["conditions"]
    if len(conditions)!=14 or conditions[0]!={"id":"baseline","kind":"identity"} or len({c["id"] for c in conditions})!=14:
        raise ValueError("The declared perturbation grid changed")
    roots = [workspace,detection,measurement,dependence]
    names = ["image_workspace","detection","measurement","dependence_groups"]
    reports = [json.loads((p/(n+"_report.json")).read_text(encoding="utf-8")) for p,n in zip(roots,names)]
    databases = [workspace/"image_workspace.sqlite",detection/"detection.sqlite",measurement/"measurements.sqlite",dependence/"dependence_groups.sqlite"]
    keys = ["workspace_database_sha256","detection_database_sha256","measurement_database_sha256","output_database_sha256"]
    for path,report,key in zip(databases,reports,keys):
        if digest(path)!=report[key]:
            raise ValueError("Input database identity changed: "+path.name)
    if reports[1]["job_states"].get("pending") or reports[2]["job_states"].get("pending"):
        raise ValueError("Finish the saved baseline schedule before this pass; failed source states remain explicit")
    source_sha=reports[0][keys[0]]
    if reports[1]["execution_contract"]["source_workspace_sha256"]!=source_sha or reports[3]["source_workspace_sha256"]!=source_sha or reports[2]["execution_contract"]["detector_database_sha256"]!=reports[1][keys[1]]:
        raise ValueError("Input provenance chain differs")
    if features.specification()!=reports[2]["execution_contract"]["feature_specification"] or software_versions()!=reports[2]["execution_contract"]["software"]:
        raise ValueError("Measurement implementation or numerical environment differs from baseline")
    context = {"contract":contract,"contract_sha256_canonical_json":canonical_digest(contract),"input_database_sha256":dict(zip(names,[r[k] for r,k in zip(reports,keys)])),
               "worker_sha256_text_lf":text_digest(Path(__file__)),"feature_specification":features.specification(),"python":platform.python_version(),
               "software":software_versions(),"opencv_threads_per_worker":1,"worker_processes":workers,
               "helper_sha256_text_lf":{n:text_digest(Path(__file__).with_name(n)) for n in ("detect_cached_images.py","measure_cached_heads.py","build_image_workspace.py")}}
    out.mkdir(parents=True,exist_ok=True)
    with runner_lock(out):
        lock=out/"execution_contract.json"
        if lock.exists():
            if json.loads(lock.read_text(encoding="utf-8"))!=context:
                raise ValueError("Execution context changed; use a new versioned run")
        else:
            if any(p.name!="runner.lock" for p in out.iterdir()):
                raise ValueError("Unrecognized existing output directory")
            lock.write_text(json.dumps(context,indent=2)+"\n",encoding="utf-8",newline="\n")
        conns=[sqlite3.connect(p.as_uri()+"?mode=ro",uri=True) for p in databases]
        try:
            cache,det,base,group=conns
            objects={r[0]:r for r in cache.execute("SELECT sha256,relative_path,width,height,decoded_rgb_sha256 FROM objects")}
            heads=list(det.execute("SELECT head_id,sha256,box_json,roi_status FROM detections ORDER BY head_id"))
            boxes={}
            for h,s,b,_ in heads:
                if _=="roi_ready":
                    boxes.setdefault(s,[]).append((h,json.loads(b)))
            memberships=dict(group.execute("SELECT sha256,component_id FROM objects"))
            baselines={h:json.loads(payload)["endpoints"] for h,payload in base.execute("SELECT head_id,raw_and_diagnostics_json FROM details")}
            source_states=dict(base.execute("SELECT head_id,status FROM jobs"))
            for h,s,b,roi in heads:
                if h not in baselines:
                    baselines[h]=[dict(zip(("endpoint_id","unit","value","original","mirror","mirror_abs_difference","status"),r)) for r in base.execute("SELECT endpoint_id,unit,value,original,mirror,mirror_abs_difference,status FROM endpoints WHERE head_id=? ORDER BY endpoint_id",(h,))]
                if h not in source_states or source_states[h] not in ("measured","measured_with_engine_errors","invalid_roi","error"):
                    raise ValueError("Unknown baseline job state")
                if len(baselines[h])!=27:
                    raise ValueError("Saved baseline endpoint slots are incomplete")
            if set(baselines)!={r[0] for r in heads} or set(objects)!=set(memberships):
                raise ValueError("Baseline or dependence membership has missing/extra units")
            with sqlite3.connect(out/"perturbations.sqlite") as db:
                db.executescript("""
                    CREATE TABLE IF NOT EXISTS jobs(head_id TEXT PRIMARY KEY,sha256 TEXT,component_id TEXT,status TEXT);
                    CREATE TABLE IF NOT EXISTS results(head_id TEXT,condition TEXT,status TEXT,payload_json TEXT,PRIMARY KEY(head_id,condition));
                    CREATE TABLE IF NOT EXISTS endpoints(head_id TEXT,condition TEXT,endpoint_id TEXT,value REAL,original REAL,mirror REAL,status TEXT,PRIMARY KEY(head_id,condition,endpoint_id));
                """)
                db.executemany("INSERT OR IGNORE INTO jobs VALUES (?,?,?,'pending')",[(h,s,memberships[s]) for h,s,_,_ in heads])
                db.commit()
                if db.execute("SELECT COUNT(*) FROM jobs WHERE status IN ('baseline_replay_mismatch','baseline_error')").fetchone()[0]:
                    raise ValueError("Baseline gate previously failed; do not resume without a new diagnosed version")
                pending={r[0] for r in db.execute("SELECT head_id FROM jobs WHERE status='pending'")}
                scheduled=[r for r in heads if r[0] in pending][:limit or None]
                tasks=((str(workspace),objects[s],h,json.loads(b),[box for other,box in boxes.get(s,[]) if other!=h],baselines[h],conditions,source_states[h]) for h,s,b,_ in scheduled)
                completed=0
                with ProcessPoolExecutor(max_workers=workers,initializer=initialize_worker) as pool:
                    running=deque()
                    for _ in range(workers):
                        task=next(tasks,None)
                        if task is not None:
                            running.append(pool.submit(run_head,task))
                    while running:
                        outcome=running.popleft().result()
                        h=outcome["head_id"]
                        if {r["condition"] for r in outcome["records"]}!={c["id"] for c in conditions} or len(outcome["records"])!=14:
                            raise ValueError("Incomplete condition denominator")
                        for record in outcome["records"]:
                            endpoints=record["result"]["endpoints"]
                            if len(endpoints)!=27 or {e["endpoint_id"] for e in endpoints}!={r["endpoint_id"] for r in features.registry()}:
                                raise ValueError("Incomplete endpoint denominator")
                            db.execute("INSERT INTO results VALUES (?,?,?,?)",(h,record["condition"],record["status"],json.dumps(record,sort_keys=True,allow_nan=False)))
                            db.executemany("INSERT INTO endpoints VALUES (?,?,?,?,?,?,?)",[(h,record["condition"],e["endpoint_id"],e["value"],e["original"],e["mirror"],e["status"]) for e in endpoints])
                        db.execute("UPDATE jobs SET status=? WHERE head_id=?",(outcome["status"],h))
                        db.commit()
                        completed+=1
                        if completed%10==0:
                            print(json.dumps({"heads_completed_this_invocation":completed,"condition_rows_saved_this_invocation":completed*14}),flush=True)
                        if outcome["status"] in ("baseline_replay_mismatch","baseline_error"):
                            for future in running:
                                future.cancel()
                            break
                        task=next(tasks,None)
                        if task is not None:
                            running.append(pool.submit(run_head,task))
                states=dict(db.execute("SELECT status,COUNT(*) FROM jobs GROUP BY status"))
                counts={"scheduled_heads":len(heads),"conditions_per_head":14,"endpoints_per_condition":27,"completed_this_invocation":completed,
                        "condition_rows":db.execute("SELECT COUNT(*) FROM results").fetchone()[0],"endpoint_rows":db.execute("SELECT COUNT(*) FROM endpoints").fetchone()[0]}
                done=len(heads)-states.get("pending",0)
                if counts["condition_rows"]!=done*14 or counts["endpoint_rows"]!=done*14*27:
                    raise ValueError("Saved denominators are incomplete")
        finally:
            for conn in conns:
                conn.close()
        status=("STOP_BASELINE_REPLAY_FAILED" if states.get("baseline_replay_mismatch") or states.get("baseline_error") else "TECHNICAL_PERTURBATION_PARTIAL" if states.get("pending") else "TECHNICAL_PERTURBATION_COMPLETED_WITH_ERRORS" if states.get("source_error") or states.get("completed_with_errors") else "TECHNICAL_PERTURBATION_COMPLETED_WITH_UNEVALUABLE_SOURCES" if states.get("source_not_evaluable") else "TECHNICAL_PERTURBATION_COMPLETED")
        report={"status":status,"execution_contract":context,"counts":counts,"job_states":states,"output_database_sha256":digest(out/"perturbations.sqlite"),
                "summary_and_claim_decision_completed":False,"ecological_models_executed":False,"independent_accuracy_estimated":False,"limits":contract["limits"]}
        (out/"perturbation_report.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8",newline="\n")
        print(json.dumps({"status":status,"counts":counts,"job_states":states},indent=2))
        return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ("workspace","detection","measurement","dependence","out-dir"):
        parser.add_argument("--"+name,type=Path,required=True)
    parser.add_argument("--limit",type=int,default=0)
    parser.add_argument("--workers",type=int,default=2)
    args=parser.parse_args()
    run(args.workspace,args.detection,args.measurement,args.dependence,args.out_dir,args.limit,args.workers)


if __name__=="__main__":
    main()
