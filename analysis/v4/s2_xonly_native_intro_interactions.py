"""V4 S2 native-versus-introduced PHOTO-climate interaction reanalysis.

For each species use the SAME climate adjustment set in both ranges, selected
BEFORE any Y is read by the original-source outcome-blind X-only analysis.
Standardise climate on the joint matched original photo cohort. Fit full
environment × introduction terms with original observer-cell dedup, equal
cell/range weights, dates and spatial-cell clustered robust covariance.

This is a posthoc observational extension, not adaptation or invasion causality.
Frozen S1 80/S2 128/H3c are not changed.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import chi2
from statsmodels.stats.multitest import multipletests
from analysis.v4.s1_s2_extreme_lego import ALL_ENV,S2_RENAMING
from analysis.v4.native_integration import FEATURES
from analysis.v4.s1_s2_extreme_range_models import _get_dated_balanced

SPECIES=("Cirsium arvense","Cirsium vulgare")
TR=tuple(FEATURES)
SUPPORT_Q=(.01,.99)
MIN_PHOTOS_EACH=100
MIN_CELLS_EACH=12

def _selected_sets(blind):
    if (blind.get("schema")!="azami_s2_outcome_blind_common_climate_selection_v1"
        or blind.get("original_source_observation_n")!=297622
        or not blind.get("trait_values_not_loaded")
        or not blind.get("past_results_P_values_and_RF_rankings_not_loaded")):
        raise ValueError("S2_NOT_OUTCOME_BLIND_SELECTION")
    ret={}
    for s in blind["results"]:
        names=s["selected_shared_covariate_set_X_only"]["selected_variables"]
        if not names or "topo_elevation" in names or len(set(names))!=len(names):
            raise ValueError("S2_INVALID_BLIND_CLIMATE_COLUMNS")
        ret[s["species"]]=list(names)
    if set(ret)!=set(SPECIES):raise ValueError("S2_SELECTION_SPECIES_MISMATCH")
    return ret

def _reader(focal,extra,measurements):
    a=pd.read_csv(focal,dtype={"obs_id":str},low_memory=False)
    b=pd.read_csv(extra,dtype={"obs_id":str},low_memory=False)
    m=pd.read_csv(measurements,dtype={"obs_id":str},low_memory=False)
    if any(len(x)!=297622 for x in (a,b)) or any(not x.obs_id.is_unique for x in (a,b,m)):
        raise ValueError("ORIGINAL_S2_SOURCE_ID_NOT_AUTHENTICATED")
    if set(a.obs_id)!=set(b.obs_id):raise ValueError("EXTERNAL_CLIMATE_ID_MISMATCH")
    if set(a.taxon_name)!=set(SPECIES) or set(a.range_state)!={"native","introduced"}:
        raise ValueError("ORIGINAL_S2_TAXON_OR_RANGE_DRIFT")
    append=[col for col in ALL_ENV if col not in a]
    if not set(append)<=set(b):raise ValueError("ORIGINAL_S2_CLIMATE_COLUMNS_MISSING")
    a=a.merge(b[["obs_id",*append]],on="obs_id",how="left",validate="one_to_one",sort=False)
    aliases={}
    for name in [v for vals in FEATURES.values() for v in vals]:
        available=[name,*S2_RENAMING.get(name,[])]
        aliases[name]=next((x for x in available if x in m.columns),None)
    if any(v is None for v in aliases.values()):
        raise ValueError("ORIGINAL_S2_ALL_8_PHOTOGRAPHIC_CONSTRUCTS_UNAVAILABLE")
    fields=list(dict.fromkeys(aliases.values()))
    measured=m[["obs_id",*fields]].rename(columns={
        alias:canonical for canonical,alias in aliases.items()})
    result=a.merge(measured,on="obs_id",how="left",validate="one_to_one",sort=False)
    if len(result)!=297622:raise ValueError("S2_IDENTITY_JOIN_CHANGED")
    return result,aliases

def _common_rectangle(part,names):
    pairs=[]
    for state in ("native","introduced"):
        arr=part.loc[part.range_state.eq(state),names].to_numpy(float)
        if len(arr)<10:raise ValueError("S2_RANGE_NO_SOURCE_ROWS")
        pairs.append(np.quantile(arr,SUPPORT_Q,axis=0))
    lo=np.maximum(pairs[0][0],pairs[1][0])
    hi=np.minimum(pairs[0][1],pairs[1][1])
    return lo,hi

def _prepare(original,species,traits,names):
    # Outcome-blind joint climate support bounds derived from ORIGINAL source
    # before looking at eligible images or their phenotype values.
    subset=original.loc[original.taxon_name.eq(species)].copy()
    climate=subset[names].apply(pd.to_numeric,errors="coerce").to_numpy(float)
    finite=np.isfinite(climate).all(axis=1)
    available=subset.loc[finite].copy()
    if len(available)<1000:raise ValueError("S2_COMMON_CLIMATE_TOO_FEW")
    low,high=_common_rectangle(available,names)
    # A non-overlapping rectangle makes inference unavailable, not a null.
    if (low>=high).any():
        return None,{"status":"no_joint_common_environmental_support"}
    candidate=available[names].to_numpy(float)
    within=np.all((candidate>=low)&(candidate<=high),axis=1)
    rows=available.loc[within].copy()
    for t in traits:rows[t]=pd.to_numeric(rows[t],errors="coerce")
    observed=np.isfinite(rows[list(traits)].to_numpy(float)).all(axis=1)
    rows=rows.loc[observed].copy()
    if len(traits)==2:
        hue=rows[list(traits)].to_numpy(float)
        norm=np.linalg.norm(hue,axis=1)
        rows=rows.loc[norm>=.20].copy()
        vec=rows[list(traits)].to_numpy(float)
        vec=vec/np.linalg.norm(vec,axis=1,keepdims=True)
        rows.loc[:,list(traits)]=vec
    complete=[]
    for state in ("native","introduced"):
        block=rows.loc[rows.range_state.eq(state)].copy()
        reduced=_get_dated_balanced(block)
        complete.append(reduced)
    df=pd.concat(complete,ignore_index=True)
    n_native=int(df.range_state.eq("native").sum())
    n_intro=int(df.range_state.eq("introduced").sum())
    cells={s:int(df.loc[df.range_state.eq(s),"_cell"].nunique()) for s in ("native","introduced")}
    meta={"status":"eligible","n_source_wild_species":int(len(subset)),
          "n_source_common_X_rectangle":int(len(available.loc[within])),
          "n_native_balanced_QC_photos":n_native,
          "n_introduced_balanced_QC_photos":n_intro,
          "cells":cells,
          "common_climate_support_from_source_only":True}
    if (min(n_native,n_intro)<MIN_PHOTOS_EACH or
        min(cells.values())<MIN_CELLS_EACH):
        meta["status"]="insufficient_joint_range_phenotype_or_climate_support"
        return None,meta
    return df,meta

def _cluster_multiy(X,Y,weight,labels):
    n,p=X.shape;q=Y.shape[1]
    groups,codes=np.unique(labels,return_inverse=True)
    if len(groups)<p+4:raise ValueError("S2_TOO_FEW_INDEPENDENT_CELLS")
    gram=X.T@(weight[:,None]*X)
    cond=float(np.linalg.cond(gram))
    if not np.isfinite(cond) or cond>1e12:raise ValueError("S2_COLLINEAR_COMMON_CLIMATE_MODEL")
    inv=np.linalg.inv(gram)
    beta=inv@(X.T@(weight[:,None]*Y))
    resid=Y-X@beta
    scores=np.empty((len(groups),p,q))
    for g in range(len(groups)):
        idx=np.flatnonzero(codes==g)
        scores[g]=inv@(X[idx].T@(weight[idx,None]*resid[idx]))
    correction=len(groups)/(len(groups)-1)*(n-1)/(n-p)
    return beta,scores,correction,cond,len(groups)

def _coefficient(beta,scores,corr,index):
    b=beta[index]
    mat=corr*np.einsum("gq,gr->qr",scores[:,index,:],scores[:,index,:])
    se=np.sqrt(np.maximum(np.diag(mat),0))
    p=float(chi2.sf(float(b@np.linalg.pinv(mat,rcond=1e-9)@b),len(b)))
    return {"beta_components":[float(v) for v in b],
            "SE_components":[float(v) for v in se],
            "pointwise_95CI_components":[[float(v-1.96*s),float(v+1.96*s)]
                                         for v,s in zip(b,se)],
            "joint_cluster_robust_p":p,
            "n_response_dimensions":int(len(b))}

def _fit(data,names,traits):
    n=len(data);q=len(traits)
    intro=data.range_state.eq("introduced").to_numpy(float)
    counts=data.groupby(["range_state","_cell"]).size()
    key=list(zip(data.range_state,data._cell))
    w=np.asarray([.5/counts[item]/
       sum(1 for k in counts.index if k[0]==item[0]) for item in key],float)
    x=data[names].to_numpy(float)
    mu=np.average(x,axis=0,weights=w)
    sd=np.sqrt(np.average((x-mu)**2,weights=w,axis=0))
    if (sd<1e-10).any():raise ValueError("S2_ZERO_VARIANCE_CLIMATE")
    z=(x-mu)/sd
    raw=data[list(traits)].to_numpy(float)
    my=np.average(raw,axis=0,weights=w)
    sy=np.sqrt(np.sum(np.average((raw-my)**2,weights=w,axis=0)))
    if sy<=1e-12:raise ValueError("S2_ZERO_VARIANCE_RESPONSE")
    Y=(raw-my)/sy
    month=data._month.to_numpy(float)
    year=data._year.to_numpy(float)
    yearz=(year-year.mean())/year.std() if year.std()>0 else np.zeros(n)
    nuisance=np.column_stack([np.sin(2*np.pi*month/12),np.cos(2*np.pi*month/12),yearz])
    X=np.column_stack([np.ones(n),intro,z,z*intro[:,None],nuisance])
    names_full=["intercept","intro",*["environment:"+v for v in names],
                *["intro_interaction:"+v for v in names],
                "month_sin","month_cos","year"]
    label=np.asarray([str(st)+":"+str(c) for st,c in key])
    beta,scores,corr,cond,G=_cluster_multiy(X,Y,w,label)
    effects=[]
    for k,var in enumerate(names):
        native=_coefficient(beta,scores,corr,2+k)
        interaction=_coefficient(beta,scores,corr,2+len(names)+k)
        # Introduced gradient = native + interaction on common pooled SD.
        idx=2+k;cross=2+len(names)+k
        introduced_beta=beta[idx]+beta[cross]
        effects.append({"exposure":var,"native_gradient":native,
                        "intro_minus_native_interaction":interaction,
                        "introduced_gradient_beta_components":
                            [float(v) for v in introduced_beta]})
    return {"status":"estimated","common_pooled_standardisation":True,
            "same_adjustment_covariates_both_ranges":True,
            "cell_cluster_n":int(G),
            "condition_number_XtWX":float(cond),
            "dated_month_sin_cos_year_controls":True,
            "effects":effects}

def analyse(original,blind):
    choose=_selected_sets(blind)
    output=[]
    for species in SPECIES:
        for trait,components in FEATURES.items():
            table,meta=_prepare(original,species,tuple(components),choose[species])
            row={"species":species,"trait":trait,"climate_explanatory_columns":choose[species],
                 "source_support":meta}
            if table is None:
                row["status"]=meta["status"]
            else:
                try:row["fit"]=_fit(table,choose[species],tuple(components));row["status"]="estimated"
                except (ValueError,np.linalg.LinAlgError) as exc:
                    row["status"]="unidentifiable"
                    row["reason"]=str(exc)
            output.append(row)
    valid=[e for row in output if row.get("status")=="estimated"
           for e in row["fit"]["effects"]]
    if valid:
        corrected=multipletests([e["intro_minus_native_interaction"]["joint_cluster_robust_p"]
                                 for e in valid],method="holm")[1]
        for e,p in zip(valid,corrected):
            e["intro_minus_native_interaction"]["posthoc_holm_2species_8traits_K"]=float(p)
    return {"schema":"azami_v4_s2_Xonly_blind_climate_common_range_interactions_v1",
            "status":"original_source_only_if_authenticated",
            "original_source_wild_n":297622,
            "n_registered_species":2,"n_registered_traits":8,
            "n_exploratory_interactions":len(valid),
            "selected_X_fully_outcome_blind":True,
            "original_S2_128_Holm_frozen_untouched":True,
            "not_phylogenetic_or_invasion_causation":True,
            "common_support_quantiles":list(SUPPORT_Q),
            "results":output}

def run_original(focal,extra,measured,selected,output):
    blind=json.loads(Path(selected).read_text())
    original,aliases=_reader(focal,extra,measured)
    r=analyse(original,blind)
    r["measurement_source_columns_available"]={k:bool(v) for k,v in aliases.items()}
    banned={"obs_id","accepted_key","taxon_name","latitude","longitude",
            "user_id","photographer_id","photo_id","individual_prediction"}
    def audit(x):
        if isinstance(x,dict):
            for k,v in x.items():
                if k in banned:raise ValueError("SOURCE_IDENTIFIERS_EXPORTED")
                audit(v)
        elif isinstance(x,list):
            for v in x:audit(v)
    audit(r)
    target=Path(output);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(r,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    return r

def main():
    p=argparse.ArgumentParser()
    for arg in ("focal","extra","measured","selected","output"):
        p.add_argument("--"+arg,required=True)
    a=p.parse_args()
    r=run_original(a.focal,a.extra,a.measured,a.selected,a.output)
    print(json.dumps({"status":r["status"],"n_interactions":r["n_exploratory_interactions"],
                      "qualified":[[x["species"],x["trait"],x["status"]] for x in r["results"]]}),flush=True)
if __name__=="__main__":main()
