"""V4 source-verified S1: trait-wise climate predictor COUNT, not one shared set.

Selection is fully nested inside taxon-held-out 5-fold CV. All eight photographic
constructs remain distinct; circular hue is one 2-dimensional response. Original
S1 observations and identifiers never leave the trusted workflow.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold
from analysis.v4 import species_environment_rf_screen as rf
from analysis.v4.species_environment_existing_inputs import load_existing_inputs

CLIMATE=tuple(x for x in rf.PREDICTORS if x!="topo_elevation")
SIZES=(3,6,10,22)
SEED=20261107
TREES=55
assert len(CLIMATE)==22

def _impute(a,b):
    med=np.nanmedian(a,axis=0)
    if not np.isfinite(med).all():
        raise ValueError("CLIMATE_SOURCE_COLUMN_NOT_PRESENT")
    return np.where(np.isfinite(a),a,med),np.where(np.isfinite(b),b,med)

def _forest(x,y,seed):
    f=RandomForestRegressor(n_estimators=TREES,min_samples_leaf=4,
                            max_features=.75,random_state=int(seed),n_jobs=1)
    f.fit(x,y[:,0] if y.shape[1]==1 else y)
    return f

def _predict(model,x,dimensions):
    return np.asarray(model.predict(x),float).reshape(-1,dimensions)

def _loss(y,yp,hue):
    if hue:
        y=rf._unit_hue(y);yp=rf._unit_hue(yp)
    return float(np.square(y-yp).sum())

def _rank(x,y,seed):
    a,_=_impute(x,x)
    f=_forest(a,y,seed)
    return np.argsort(-f.feature_importances_,kind="stable")

def _choose(x,y,hue,seed):
    fold=KFold(n_splits=3,shuffle=True,random_state=seed)
    errors={k:[] for k in SIZES}
    for j,(tr,te) in enumerate(fold.split(x)):
        xt,xv=_impute(x[tr],x[te])
        ranked=np.argsort(-_forest(xt,y[tr],seed+j).feature_importances_,
                          kind="stable")
        for k in SIZES:
            est=_forest(xt[:,ranked[:k]],y[tr],seed+100+j+k)
            pred=_predict(est,xv[:,ranked[:k]],y.shape[1])
            errors[k].append(_loss(y[te],pred,hue)/len(te))
    means={k:float(np.mean(v)) for k,v in errors.items()}
    winner=min(SIZES,key=lambda k:(means[k],k))
    # Conservative one-SE choice avoids adding predictors without clear gains.
    tolerance=float(np.std(errors[winner],ddof=1)/np.sqrt(3))
    pick=min(k for k in SIZES if means[k]<=means[winner]+tolerance)
    return pick,means

def _score_one(x,y,train,test,cols,hue,seed):
    xt,xv=_impute(x[train][:,cols],x[test][:,cols])
    model=_forest(xt,y[train],seed)
    pred=_predict(model,xv,y.shape[1])
    return _loss(y[test],pred,hue)

def run_cohort(cells,trait):
    hue=trait=="hue"
    x=cells[list(CLIMATE)].to_numpy(float)
    y=cells[[f"y{i}" for i in range(len(rf.FEATURES[trait]))]].to_numpy(float)
    if len(cells)<40 or not np.isfinite(y).all():
        raise ValueError("SOURCE_SPECIES_ELIGIBILITY_FAILURE")
    if not np.isfinite(x).any(axis=0).all():
        raise ValueError("SOURCE_MISSING_ENVIRONMENT_AXES")
    outer=KFold(n_splits=5,shuffle=True,random_state=SEED)
    all_cols=np.arange(len(CLIMATE))
    benchmark=np.array([CLIMATE.index(s) for s in
                       ("chelsa_bio08","chelsa_bio18","chelsa_rsds_mean")])
    errors={name:0. for name in ("selected","full22","fixed3")}
    baseline=0.;selection=np.zeros(len(CLIMATE),int);choices=[]
    sector=cells.get("_sector",pd.Series([-1]*len(cells))).to_numpy(int)
    for j,(tr,te) in enumerate(outer.split(x)):
        k,_=_choose(x[tr],y[tr],hue,SEED+1000+j)
        xt,_=_impute(x[tr],x[tr])
        ranked=np.argsort(-_forest(xt,y[tr],SEED+2000+j).feature_importances_,
                          kind="stable")
        cols=ranked[:k]
        selection[cols]+=1;choices.append(k)
        base=np.repeat(y[tr].mean(axis=0)[None,:],len(te),axis=0)
        baseline+=_loss(y[te],base,hue)
        for name,part in (("selected",cols),("full22",all_cols),("fixed3",benchmark)):
            errors[name]+=_score_one(x,y,tr,te,part,hue,SEED+3000+j)
    # Descriptive full-fit RF impurity order only; outer frequencies are separate.
    xt,_=_impute(x,x)
    imp=_forest(xt,y,SEED+9000).feature_importances_
    rows=[{"predictor":CLIMATE[i],
           "full_fit_impurity_importance":float(imp[i]),
           "outer_selected_fraction":float(selection[i]/5),
           "outer_selected_n":int(selection[i])} for i in range(len(CLIMATE))]
    rows.sort(key=lambda v:(-v["outer_selected_fraction"],
                             -v["full_fit_impurity_importance"],v["predictor"]))
    return {"trait":trait,"n_taxa":int(len(cells)),
            "source_sampled_longitude_sector_counts":{str(i):int((sector==i).sum())
                for i in (0,1,2)},
            "selection_size_outer_counts":{str(k):int(choices.count(k)) for k in SIZES},
            "taxon_holdout_skill":{name:float(1-v/baseline) for name,v in errors.items()},
            "variable_ranking":rows}

def run_existing(native,decoded_receipt,seven,nine,nine_manifest,source_receipt,
                 source_sha,output):
    frame,rec=load_existing_inputs(native,decoded_receipt,seven,nine,
                                   source_receipt,nine_manifest,
                                   source_receipt_sha256=source_sha,
                                   expected_n=313350)
    loc=pd.read_csv(native,usecols=["obs_id","longitude"],dtype={"obs_id":str})
    if len(loc)!=313350 or not loc.obs_id.is_unique:
        raise ValueError("SOURCE_GEOGRAPHY_DRIFT")
    ids=frame.obs_id.copy()
    frame=frame.merge(loc,on="obs_id",validate="one_to_one",how="left",sort=False)
    if not frame.obs_id.equals(ids):raise ValueError("SOURCE_JOIN_CHANGED")
    result=[]
    for trait in rf.TRAITS:
        sample=rf.aggregate_species(frame,trait)
        result.append(run_cohort(sample,trait))
    document={"schema":"azami_v4_s1_traitwise_nested_rf_variable_count_v1",
              "status":"source_authenticated_posthoc_exploratory",
              "original_source_n":len(frame),
              "source_native_sha256":rec["input_digests"]["native"],
              "source_receipt_sha256":source_sha,
              "climate_predictors":list(CLIMATE),
              "elevation_excluded_from_selection":True,
              "same_cohort_nested_taxon_fivefold":True,
              "selection_inner_threefold_sizes":list(SIZES),
              "selection_tuned_on_training_only":True,
              "models":["per_trait_selected","all22","fixed_BIO08_BIO18_RSDS"],
              "no_new_confirmatory_p_values":True,
              "original_frozen_GEB_H3c_untouched":True,
              "results":result}
    banned={"obs_id","accepted_key","taxon","taxon_name",
            "latitude","longitude","user_id","photographer_id","photo_id"}
    def audit(z):
        if isinstance(z,dict):
            for k,v in z.items():
                if k in banned:raise ValueError("PROTECTED_KEY_EXPORTED")
                audit(v)
        elif isinstance(z,list):
            for v in z:audit(v)
    audit(document)
    Path(output).parent.mkdir(parents=True,exist_ok=True)
    Path(output).write_text(json.dumps(document,indent=2,allow_nan=False)+"\n")
    return document

def main():
    p=argparse.ArgumentParser()
    for name in ("native","decoded-receipt","seven","nine","nine-manifest",
                 "source-receipt","source-sha","output"):
        p.add_argument("--"+name,required=True)
    a=p.parse_args()
    result=run_existing(a.native,a.decoded_receipt,a.seven,a.nine,
                        a.nine_manifest,a.source_receipt,a.source_sha,a.output)
    print(json.dumps({"status":result["status"],"n_traits":len(result["results"]),
                      "skills":[(x["trait"],x["taxon_holdout_skill"])
                                for x in result["results"]]}),flush=True)

if __name__=="__main__":main()
