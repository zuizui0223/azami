"""Summarize all perturbations with explicit attrition and dependence weighting.

Head changes are averaged within image, image summaries within known dependence
component, and components receive equal weight. This is a technical sensitivity
description, not an accuracy estimate or a significance-based admission gate.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sqlite3

import numpy as np
import pandas as pd

from .image_features import clean, registry
from .workflow import ROOT, digest, text_digest

SPECIFICATION = {
    "version":"v3_component_weighted_perturbation_summary_v1",
    "timing":"specified after the perturbation execution began, before inspecting perturbation values; no preregistration claim",
    "eligibility":"both original-baseline and perturbed endpoint values must be finite and QC-usable for paired change estimates; all four QC transitions retain the full scheduled-head denominator",
    "weighting":"arithmetic mean of per-head signed/absolute changes within exact image object; arithmetic mean of image summaries within known dependence component; equal component weight for reported distributions",
    "attrition":"loss is computed among original-baseline-eligible heads within image, then averaged within component and across components; groups with no surviving pairs still contribute to loss",
    "rank":"descriptive Spearman correlation of within-image then within-component paired means, requiring at least three nonconstant components; no p-value or accuracy interpretation",
    "joint_hue":"shortest absolute angular difference in degrees from sine/cosine pairs; no arbitrary linear hue rank or signed arithmetic mean",
    "joint_composition":"total-variation distance, one-half the summed absolute changes of the closed four-part colour composition; no independent-part inference",
    "strata":["all_cached","recorded_development_exposed","no_recorded_development_exposure"],
    "decision":"quantify sensitivity and missingness; no favourable correlation threshold, physical-accuracy promotion or ecological conclusion",
}


def summarize_frame(frame, joint=False):
    """Input has one row per head, matched on identity before any filtering."""
    base_ok=frame["baseline_usable"].astype(bool)
    changed_ok=frame["perturbed_usable"].astype(bool)
    pair=base_ok & changed_ok
    n=int(len(frame))
    result={"scheduled_heads":n,"scheduled_images":int(frame.sha256.nunique()),"scheduled_components":int(frame.component_id.nunique()),
            "baseline_usable_heads":int(base_ok.sum()),"paired_usable_heads":int(pair.sum()),
            "lost_usable_heads":int((base_ok & ~changed_ok).sum()),"gained_usable_heads":int((~base_ok & changed_ok).sum()),
            "neither_usable_heads":int((~base_ok & ~changed_ok).sum()),
            "paired_images":0,"paired_components":0,"mean_signed_change":None,"mean_absolute_change":None,
            "median_component_absolute_change":None,"p95_component_absolute_change":None,"spearman_paired_component_means":None,
            "component_weighted_loss_fraction":None}
    if result["paired_usable_heads"]+result["lost_usable_heads"]+result["gained_usable_heads"]+result["neither_usable_heads"]!=n:
        raise ValueError("QC transition denominator mismatch")
    loss=frame.loc[base_ok,["sha256","component_id"]].copy()
    if len(loss):
        loss["lost"]=(~changed_ok[base_ok]).astype(float).to_numpy()
        per_image=loss.groupby(["component_id","sha256"],sort=True).lost.mean()
        result["component_weighted_loss_fraction"]=float(per_image.groupby(level="component_id").mean().mean())
    pairs=frame.loc[pair].copy()
    if pairs.empty:
        return result
    keys=["signed_change","absolute_change","baseline_value","perturbed_value"]
    image_means=pairs.groupby(["component_id","sha256"],sort=True)[keys].mean()
    component_means=image_means.groupby(level="component_id").mean()
    result.update(paired_images=len(image_means),paired_components=len(component_means),
                  mean_absolute_change=float(component_means.absolute_change.mean()),
                  median_component_absolute_change=float(component_means.absolute_change.median()),
                  p95_component_absolute_change=float(component_means.absolute_change.quantile(.95)))
    if not joint:
        result["mean_signed_change"]=float(component_means.signed_change.mean())
        a,b=component_means.baseline_value,component_means.perturbed_value
        if len(component_means)>=3 and a.nunique()>1 and b.nunique()>1:
            result["spearman_paired_component_means"]=float(a.rank(method="average").corr(b.rank(method="average")))
    return clean(result)


def scalar_frame(table):
    out=table.copy()
    out["baseline_value"]=out.baseline_value.astype(float)
    out["perturbed_value"]=out.perturbed_value.astype(float)
    out["baseline_usable"]=out.baseline_status.eq("usable") & np.isfinite(out.baseline_value)
    out["perturbed_usable"]=out.perturbed_status.eq("usable") & np.isfinite(out.perturbed_value)
    out["signed_change"]=out.perturbed_value-out.baseline_value
    out["absolute_change"]=out.signed_change.abs()
    return out


def joint_frame(table, endpoints, kind):
    base=table[table.endpoint_id.isin(endpoints)].copy()
    identifiers=["head_id","sha256","component_id","development_exposed"]
    units=base[identifiers].drop_duplicates().set_index("head_id")
    b=base.pivot(index="head_id",columns="endpoint_id",values="baseline_value").reindex(columns=endpoints).astype(float)
    p=base.pivot(index="head_id",columns="endpoint_id",values="perturbed_value").reindex(columns=endpoints).astype(float)
    bs=base.pivot(index="head_id",columns="endpoint_id",values="baseline_status").reindex(columns=endpoints)
    ps=base.pivot(index="head_id",columns="endpoint_id",values="perturbed_status").reindex(columns=endpoints)
    units=units.reindex(b.index)
    bok=bs.eq("usable").all(axis=1) & np.isfinite(b).all(axis=1)
    pok=ps.eq("usable").all(axis=1) & np.isfinite(p).all(axis=1)
    if kind=="hue":
        # Endpoint order is sine, cosine.
        ba=np.rad2deg(np.arctan2(b.iloc[:,0],b.iloc[:,1]))
        pa=np.rad2deg(np.arctan2(p.iloc[:,0],p.iloc[:,1]))
        bok &= np.hypot(b.iloc[:,0],b.iloc[:,1])>1e-12
        pok &= np.hypot(p.iloc[:,0],p.iloc[:,1])>1e-12
        difference=((pa-ba+180)%360-180).abs()
    elif kind=="composition":
        bok &= b.ge(0).all(axis=1) & b.le(1).all(axis=1) & b.sum(axis=1).sub(1).abs().le(1e-6)
        pok &= p.ge(0).all(axis=1) & p.le(1).all(axis=1) & p.sum(axis=1).sub(1).abs().le(1e-6)
        difference=(p-b).abs().sum(axis=1)/2
    else:
        raise ValueError("Unknown joint metric")
    units["baseline_usable"],units["perturbed_usable"]=bok,pok
    units["baseline_value"],units["perturbed_value"],units["signed_change"]=np.nan,np.nan,np.nan
    units["absolute_change"]=difference
    return units.reset_index()


def strata(frame):
    return [("all_cached",frame),("recorded_development_exposed",frame[frame.development_exposed.eq(1)]),
            ("no_recorded_development_exposure",frame[frame.development_exposed.eq(0)])]


def run(perturbation,measurement,dependence,out):
    perturbation,measurement,dependence,out=[p.resolve() for p in (perturbation,measurement,dependence,out)]
    if out==ROOT or (ROOT in out.parents and not any(out.is_relative_to(ROOT/p) for p in ("local_data","outputs"))):
        raise ValueError("Use an external or ignored local output directory")
    if out.exists():
        raise ValueError("Preserve earlier summary versions")
    report=json.loads((perturbation/"perturbation_report.json").read_text(encoding="utf-8"))
    if not report["status"].startswith("TECHNICAL_PERTURBATION_COMPLETED") or report["job_states"].get("pending"):
        raise ValueError("Incomplete or stopped perturbation grid cannot become a completed summary")
    paths=[perturbation/"perturbations.sqlite",measurement/"measurements.sqlite",dependence/"dependence_groups.sqlite"]
    expected=[report["output_database_sha256"],report["execution_contract"]["input_database_sha256"]["measurement"],report["execution_contract"]["input_database_sha256"]["dependence_groups"]]
    for path,sha in zip(paths,expected):
        if digest(path)!=sha:
            raise ValueError("Summary input identity changed: "+path.name)
    conns=[sqlite3.connect(p.as_uri()+"?mode=ro",uri=True) for p in paths]
    try:
        probes,base,groups=conns
        baseline=pd.read_sql_query("SELECT head_id,endpoint_id,value AS baseline_value,status AS baseline_status FROM endpoints",base)
        heads=pd.read_sql_query("SELECT head_id,sha256,component_id FROM jobs",probes)
        exposure=pd.read_sql_query("SELECT component_id,development_pool_exposed AS development_exposed FROM components WHERE n_cached_objects>0",groups)
        heads=heads.merge(exposure,on="component_id",how="left",validate="many_to_one")
        if heads.development_exposed.isna().any():
            raise ValueError("Missing component exposure annotation")
        baseline=baseline.merge(heads,on="head_id",how="left",validate="many_to_one")
        if baseline.component_id.isna().any() or baseline.duplicated(["head_id","endpoint_id"]).any() or len(baseline)!=len(heads)*27:
            raise ValueError("Baseline head/endpoint denominator changed")
        records=[]
        conditions=report["execution_contract"]["contract"]["conditions"]
        expected_rows=len(heads)*27*len(conditions)
        if probes.execute("SELECT COUNT(*) FROM endpoints").fetchone()[0]!=expected_rows:
            raise ValueError("Incomplete perturbation endpoint denominator")
        for condition in conditions:
            changed=pd.read_sql_query("SELECT head_id,endpoint_id,value AS perturbed_value,status AS perturbed_status FROM endpoints WHERE condition=?",probes,params=[condition["id"]])
            table=baseline.merge(changed,on=["head_id","endpoint_id"],how="left",validate="one_to_one")
            if len(changed)!=len(baseline) or table.perturbed_status.isna().any():
                raise ValueError("Incomplete paired condition membership")
            for endpoint in registry():
                frame=scalar_frame(table[table.endpoint_id.eq(endpoint["endpoint_id"])])
                for label,subset in strata(frame):
                    records.append({"condition":condition["id"],"metric_id":endpoint["endpoint_id"],"metric_kind":"registered_endpoint","unit":endpoint["unit"],"exposure_stratum":label,**summarize_frame(subset)})
            joint_specs=[("joint_hue_angular_change",["corolla_hue_sin","corolla_hue_cos"],"hue","degrees"),
                         ("joint_colour_composition_total_variation",[f"corolla_{s}_pixel_fraction" for s in ("white","redmagenta","purple","yellow")],"composition","fraction")]
            for name,endpoints,kind,unit in joint_specs:
                for label,subset in strata(joint_frame(table,endpoints,kind)):
                    records.append({"condition":condition["id"],"metric_id":name,"metric_kind":"joint_diagnostic_not_additional_endpoint","unit":unit,"exposure_stratum":label,**summarize_frame(subset,joint=True)})
    finally:
        for conn in conns:
            conn.close()
    if len(records)!=14*29*3:
        raise ValueError("Incomplete condition/endpoint/stratum summary grid")
    out.mkdir(parents=True)
    pd.DataFrame(records).to_csv(out/"technical_sensitivity_summary.csv",index=False,lineterminator="\n")
    result={"status":"COMPONENT_WEIGHTED_TECHNICAL_SENSITIVITY_SUMMARIZED","specification":SPECIFICATION,
            "input_database_sha256":dict(zip(["perturbation","measurement","dependence"],expected)),
            "implementation_sha256_text_lf":text_digest(Path(__file__)),"summary_rows":len(records),
            "summary_sha256":digest(out/"technical_sensitivity_summary.csv"),
            "claim_status":"specified_pixel_and_crop_sensitivity_only","independent_accuracy_estimated":False,"ecological_models_executed":False,
            "limits":["Do not treat hue components, composition parts, exposure strata or perturbation conditions as independent biological tests.",
                       "Pairwise change summaries describe surviving usable pairs; report eligibility loss alongside them.",
                       "Equal dependence-component weighting is not a probability-sampling correction or proof of biological independence.",
                       "No recorded development exposure is not a guarantee of fresh independent evaluation.",
                       "No automated pass/fail accuracy threshold or ecological claim is inferred from these descriptive summaries."]}
    (out/"technical_sensitivity_summary_report.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8",newline="\n")
    print(json.dumps(result,indent=2))
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ("perturbation","measurement","dependence","out-dir"):
        parser.add_argument("--"+name,type=Path,required=True)
    args=parser.parse_args()
    run(args.perturbation,args.measurement,args.dependence,args.out_dir)


if __name__=="__main__":
    main()
