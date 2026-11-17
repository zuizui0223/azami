"""V4 S2: native/introduced colour predictors chosen separately, spatially held out.

Source and 6 frozen coefficients rechecked BEFORE exploratory subset comparison.
One WCVP species x native/introduced x photographic colour endpoint is one
selection problem. Photo/observer/location identifiers never exported.
"""
from __future__ import annotations
import argparse
import itertools
import json
from pathlib import Path
import numpy as np
from sklearn.model_selection import GroupKFold
from analysis.v4 import s2_original_colour_spatial_block_sensitivity as s2
from analysis.v4.s1_s2_extreme_range_models import _formula
from analysis.v4.s1_s2_vif5_si10_alltraits import FOUR_WARM

EXPOSURES=tuple(FOUR_WARM)
SUBSETS=tuple(tuple(a) for n in range(1,5)
               for a in itertools.combinations(range(4),n))
TRAITS=("corolla_lab_lightness","corolla_lab_chroma")
SPECIES=("Cirsium arvense","Cirsium vulgare")
STATUSES=("native","introduced")
assert len(EXPOSURES)==4 and len(SUBSETS)==15

def _design(full,subset):
    if full.ndim!=2 or full.shape[1]<5 or len(subset)==0:
        raise ValueError("S2_ORIGINAL_DESIGN_MISMATCH")
    return full[:,[0]+[1+j for j in subset]+list(range(5,full.shape[1]))]

def _weighted_predict(X,y,w,train,test):
    a=X[train];b=X[test]
    ww=np.sqrt(w[train])
    coef=np.linalg.lstsq(a*ww[:,None],y[train]*ww,rcond=1e-10)[0]
    pred=b@coef
    if not np.isfinite(pred).all():
        raise ValueError("S2_MODEL_NONFINITE")
    return pred

def _weighted_loss(y,pred,w):
    return float(np.dot(w,(y-pred)**2))

def _select_within_train(full,y,w,groups,train):
    x=full[train];yy=y[train];wt=w[train];g=groups[train]
    if len(np.unique(g))<3:raise ValueError("S2_SOURCE_CELL_SUPPORT_TOO_SMALL")
    fold=GroupKFold(n_splits=3)
    totals={}
    for part in SUBSETS:
        A=_design(x,part);err=0.;base=0.
        for j,(t,v) in enumerate(fold.split(A,yy,groups=g)):
            pred=_weighted_predict(A,yy,wt,t,v)
            mu=float(np.average(yy[t],weights=wt[t]))
            err+=_weighted_loss(yy[v],pred,wt[v])
            base+=_weighted_loss(yy[v],np.full(len(v),mu),wt[v])
        totals[part]=float(err/base) if base>0 else float("inf")
    best=min(totals,key=lambda t:(totals[t],len(t),t))
    return best

def one_panel(frame,species,state,trait):
    balanced,coverage=s2.original_sample(frame,species,state,trait)
    full,_,_,_,_=_formula(balanced,list(EXPOSURES))
    y=balanced[trait].to_numpy(float)
    group=balanced._cell.to_numpy()
    if not np.isfinite(y).all():raise ValueError("NONFINITE_PHOTO_COLOUR")
    counts=balanced.groupby("_cell").size()
    w=balanced._cell.map(1/counts).to_numpy(float)
    outer=GroupKFold(n_splits=4)
    errors={"traitwise_selected":0.,"full4":0.,"climate3_no_elev":0.,
            "radiation_only":0.}
    baseline=0.;selection_counts={k:0 for k in EXPOSURES};choices={}
    fixed=(0,1,2);sun=(2,)
    for train,test in outer.split(full,y,groups=group):
        part=_select_within_train(full,y,w,group,train)
        key=";".join(EXPOSURES[j] for j in part)
        choices[key]=choices.get(key,0)+1
        for j in part:selection_counts[EXPOSURES[j]]+=1
        mu=float(np.average(y[train],weights=w[train]))
        baseline+=_weighted_loss(y[test],np.full(len(test),mu),w[test])
        for model,subset in (("traitwise_selected",part),
                             ("full4",(0,1,2,3)),
                             ("climate3_no_elev",fixed),
                             ("radiation_only",sun)):
            A=_design(full,subset)
            pred=_weighted_predict(A,y,w,train,test)
            errors[model]+=_weighted_loss(y[test],pred,w[test])
    return {"species":species,"range_state":state,"photographic_trait":trait,
            "n_original_balanced_photo_rows":int(len(balanced)),
            "n_unique_original_geographic_cells":int(len(counts)),
            "source_support":coverage,
            "nested_geo_cell_holdout_skill":{
                k:float(1-v/baseline) for k,v in errors.items()},
            "outer_selected_feature_frequency":{
                k:float(v/4) for k,v in selection_counts.items()},
            "outer_selected_subsets":choices}

def run(focal,extra,measured,output):
    frame,provenance=s2.source_cohort(focal,extra,measured)
    if len(frame)!=297622:raise ValueError("SOURCE_S2_ROSTER_DRIFT")
    reference=[]
    for species,state,trait,exposure,n,expected in s2.EXPECTED:
        part,_=s2.original_sample(frame,species,state,trait)
        if len(part)!=n:raise ValueError("FROZEN_S2_COHORT_DRIFT")
        _,_,_,_,fit=s2._original_fit(part,trait,exposure)
        if not np.isclose(fit["beta"],expected,rtol=0,atol=2e-7):
            raise ValueError("FROZEN_S2_COEFFICIENT_CHANGED")
        reference.append({"species":species,"range_state":state,"trait":trait,
                          "exposure":exposure,"verified_n":n,
                          "reference_beta":float(expected)})
    panels=[one_panel(frame,s,state,t) for s in SPECIES
            for state in STATUSES for t in TRAITS]
    out={"schema":"azami_v4_s2_native_intro_traitwise_subset_geo_cv_v1",
         "status":"original_s2_source_verified_exploratory",
         "original_wild_source_n":297622,
         "source_integrity":provenance,
         "preserves_frozen_S2_128_Holm_and_S1_H3c":True,
         "frozen_reference_rows_checked":reference,
         "candidate_exposure_columns":list(EXPOSURES),
         "fifteen_subsets_tested_with_inner_geo_cell_threefold":True,
         "outer_geographic_cell_group_kfold":4,
         "selection_repeated_only_in_outer_training":True,
         "inference":"predictive_comparison_not_range_interaction_or_evolution",
         "panels":panels}
    prohibited={"obs_id","taxon_name","accepted_key","latitude","longitude",
                "user_id","photographer_id","photo_id","individual_prediction"}
    def audit(obj):
        if isinstance(obj,dict):
            for k,v in obj.items():
                if k in prohibited:raise ValueError("PRIVATE_KEY_IN_PUBLIC_S2")
                audit(v)
        elif isinstance(obj,list):
            for v in obj:audit(v)
    audit(out)
    Path(output).parent.mkdir(parents=True,exist_ok=True)
    Path(output).write_text(json.dumps(out,ensure_ascii=False,indent=2,
                                      allow_nan=False)+"\n",encoding="utf-8")
    return out

def main():
    ap=argparse.ArgumentParser()
    for flag in ("focal","extra","measured","out"):
        ap.add_argument("--"+flag,required=True)
    args=ap.parse_args()
    out=run(args.focal,args.extra,args.measured,args.out)
    print(json.dumps({"status":out["status"],"n_panels":len(out["panels"]),
                      "skills":[[x["species"],x["range_state"],
                                 x["photographic_trait"],
                                 x["nested_geo_cell_holdout_skill"]]
                                for x in out["panels"]]}),flush=True)

if __name__=="__main__":main()
