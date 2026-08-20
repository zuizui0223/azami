#!/usr/bin/env python3
"""Submission-facing control entry point for the frozen Chapter 1 analysis."""
from __future__ import annotations
import argparse, json, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PIPELINE = Path(__file__).with_name("pipeline.json")

def load_pipeline(): return json.loads(PIPELINE.read_text(encoding="utf-8"))

def check(_):
    p=load_pipeline(); missing=[]
    for stage in p.get("stages",{}).values():
        for key in ("script","selection_script","annotation_app","finalizer","evaluation_script","measurement_script","join_script","residual_export_script","region_script","input_builder","audit_script","summary_exporter","japan38_observation_exporter"):
            value=stage.get(key)
            if value and not (ROOT/value).is_file(): missing.append(value)
        wf=stage.get("workflow")
        if wf and not (ROOT/wf).is_file(): missing.append(wf)
    required=["README.md","manuscript/SUBMISSION_MANUSCRIPT.md","manuscript/final_claims.json","manuscript/COHORT_FLOW_AND_ANALYSIS_LEDGER.md","manuscript/EXTERNAL_COMPLETION_GATES.md"]
    missing.extend(x for x in required if not (ROOT/x).is_file())
    if missing: raise SystemExit("Missing canonical paths:\n- "+"\n- ".join(sorted(set(missing))))
    print(json.dumps({"status":"ok","analysis_status":p.get("status"),"checked_stages":sorted(p.get("stages",{}))},indent=2))

def claims(_):
    script=ROOT/load_pipeline()["stages"]["claims"]["script"]
    subprocess.run([sys.executable,str(script)],cwd=ROOT,check=True)

def summary(_):
    p=load_pipeline(); out={"status":p.get("status"),"submission_story":p.get("submission_story"),"current_submission_gates":p.get("current_submission_gates"),"policy":p.get("policy")}
    print(json.dumps(out,ensure_ascii=False,indent=2))

def main():
    ap=argparse.ArgumentParser(); sp=ap.add_subparsers(dest="command",required=True)
    for name,func in [("check",check),("claims",claims),("summary",summary)]:
        q=sp.add_parser(name); q.set_defaults(func=func)
    a=ap.parse_args(); a.func(a)
if __name__=="__main__": main()
