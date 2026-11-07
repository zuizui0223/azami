"""V4 original S1 between-/within-taxon climate-gradient reanalysis.

Predictors are retrieved ONLY from the source-authenticated X-only climate
selection JSON created before fitting. Never choose predictors using photo
phenotypes, former P values, RF scores, or sPLS loadings in this script.
Inference uses equal total taxon weight and taxon-cluster sandwich covariance.
All results exploratory; original V4 frozen 80/128 Holm and H3c untouched.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import chi2
from statsmodels.stats.multitest import multipletests
from analysis.v4.species_environment_existing_inputs import load_existing_inputs
from analysis.v4.species_environment_scan import FEATURES

MIN_TAXA=35
MIN_PHOTOS=5
YENDPOINTS=tuple(FEATURES)

def selected_axes(choice:dict)->tuple[list[str],list[str]]:
    if (choice.get("status")!="original_S1_authenticated_X_only_blind_selection"
        or choice.get("uses_original_or_previous_flower_trait_results") is not False
        or choice.get("uses_photographic_trait_columns") is not False
        or choice.get("original_native_observations")!=313350):
        raise ValueError("S1_SELECTION_NOT_X_ONLY_AUTHENTICATED")
    section=choice["section31"]
    between=section["between_taxa_X_only"]["selected_variables"]
    within=section["within_taxa_X_only"]["selected_variables"]
    if not between or not within or len(set(between))!=len(between) or len(set(within))!=len(within):
        raise ValueError("S1_X_ONLY_SELECTION_EMPTY_OR_DUPLICATED")
    if "topo_elevation" in between+within:raise ValueError("S1_SELECTION_CONTAINS_ELEVATION")
    return list(between),list(within)

def _prepare(frame:pd.DataFrame,trait:str,between:list[str],within:list[str]):
    ycols=tuple(FEATURES[trait])
    p=tuple(dict.fromkeys([*between,*within]))
    required={"accepted_key",*p,*ycols,*(q+"__eligible" for q in ycols)}
    if not required<=set(frame):raise ValueError("S1_PRIVATE_TRAIT_OR_CLIMATE_MISSING")
    # Photo QC mirrors the original exact endpoint quality gates; all X finite.
    y=frame[list(ycols)].apply(pd.to_numeric,errors="coerce").to_numpy(float)
    x=frame[list(p)].apply(pd.to_numeric,errors="coerce").to_numpy(float)
    valid=(np.isfinite(y).all(axis=1)&np.isfinite(x).all(axis=1)
           &frame.accepted_key.notna().to_numpy())
    for q in ycols:valid &= frame[q+"__eligible"].eq(True).fillna(False).to_numpy(bool)
    if trait=="hue":
        norm=np.linalg.norm(y,axis=1)
        valid &= norm>=.2
        y=np.divide(y,norm[:,None],out=np.zeros_like(y),where=norm[:,None]>1e-12)
    keep=frame.loc[valid,["accepted_key",*p]].copy()
    keep["accepted_key"]=keep.accepted_key.astype(str)
    group=keep.groupby("accepted_key",sort=True,observed=True)
    counts=group.size()
    qualified=counts[counts>=MIN_PHOTOS].index
    mask=keep.accepted_key.isin(qualified).to_numpy(bool)
    keep=keep.loc[mask].copy()
    y=y[valid][mask]
    if keep.accepted_key.nunique()<MIN_TAXA:raise ValueError("S1_TOO_FEW_TAXA")
    group=keep.groupby("accepted_key",sort=True,observed=True)
    taxon_sizes=group.size()
    w=1/keep.accepted_key.map(taxon_sizes).to_numpy(float)
    # Each taxon has weight one; no source photo has independent-species weight.
    if trait=="hue":
        v=y.mean(axis=0)
        if np.linalg.norm(v)<1e-12:raise ValueError("HUE_MEAN_UNIDENTIFIABLE")
    mean=np.average(y,axis=0,weights=w)
    scale=np.sqrt(np.sum(np.average((y-mean)**2,axis=0,weights=w),axis=0))
    if not np.isfinite(scale) or scale<1e-12:raise ValueError("S1_TRAIT_ZERO_VARIANCE")
    y=(y-mean)/scale
    return keep,y,w,{"n_source_trait_photos_before_finite_X":int(valid.sum()),
                     "n_same_photo_eligible_and_complete_X":int(len(keep)),
                     "n_taxa":int(len(taxon_sizes))}

def _columns(frame,between,within,w,geo=None):
    grp=frame.groupby("accepted_key",observed=True,sort=True)
    blocks=[]
    names=["intercept"]
    for role,variables in (("between",between),("within",within)):
        total=[]
        for v in variables:
            m=grp[v].transform("mean").to_numpy(float)
            raw=m if role=="between" else frame[v].to_numpy(float)-m
            center=float(np.average(raw,weights=w))
            sd=float(np.sqrt(np.average((raw-center)**2,weights=w)))
            if not np.isfinite(sd) or sd<1e-10:
                raise ValueError("S1_SELECTED_CLIMATE_ZERO_VARIANCE_"+v)
            total.append((raw-center)/sd)
            names.append(role+":"+v)
        blocks.extend(total)
    x=np.column_stack([np.ones(len(frame)),*blocks])
    if geo is not None:
        for j in range(geo.shape[1]):
            arr=np.asarray(geo[:,j],float)
            mu=float(np.average(arr,weights=w))
            sd=float(np.sqrt(np.average((arr-mu)**2,weights=w)))
            if sd>1e-10:
                x=np.column_stack([x,(arr-mu)/sd])
                names.append("geo_"+str(j))
    return x,names

def _sandwich(X,Y,w,groups):
    n,p=X.shape;q=Y.shape[1]
    uniq,codes=np.unique(groups,return_inverse=True)
    G=len(uniq)
    if n<=p+2 or G<=p+2:raise ValueError("S1_DEGREES_OF_FREEDOM_INSUFFICIENT")
    mat=X.T@(w[:,None]*X)
    cond=float(np.linalg.cond(mat))
    if not np.isfinite(cond) or cond>1e12:
        raise ValueError("S1_DESIGN_ILL_CONDITIONED")
    inverse=np.linalg.inv(mat)
    beta=inverse@(X.T@(w[:,None]*Y))
    resid=Y-X@beta
    correction=(G/(G-1))*((n-1)/(n-p))
    influences=np.empty((G,p,q))
    for g in range(G):
        ix=np.flatnonzero(codes==g)
        influences[g]=inverse@(X[ix].T@(w[ix,None]*resid[ix]))
    return beta,influences,correction,cond

def _inference(beta,influences,correction,names,focus):
    output=[]
    for i,name in enumerate(names):
        if name not in focus:continue
        b=np.asarray(beta[i],float)
        V=correction*np.einsum("gq,gr->qr",influences[:,i,:],influences[:,i,:])
        se=np.sqrt(np.maximum(np.diag(V),0))
        try:
            stat=float(b@np.linalg.pinv(V,rcond=1e-9)@b)
            pv=float(chi2.sf(stat,df=len(b)))
        except (ValueError,np.linalg.LinAlgError):pv=None
        if pv is not None and (not np.isfinite(pv)):pv=None
        output.append({"term":name,"beta_std_components":[float(z) for z in b],
                       "cluster_SE_components":[float(z) for z in se],
                       "pointwise_normal_95CI_components":[[float(t-1.96*s),float(t+1.96*s)]
                                                            for t,s in zip(b,se)],
                       "cluster_robust_joint_p":pv,
                       "joint_response_dimensions":int(len(b))})
    return output

def one_trait(frame,trait,between,within,geography=None):
    selected,y,w,support=_prepare(frame,trait,between,within)
    x,names=_columns(selected,between,within,w)
    groups=selected.accepted_key.to_numpy(str)
    focus=set("between:"+v for v in between)|set("within:"+v for v in within)
    out={"trait":trait,"support":support,
         "selected_between_X":between,"selected_within_X":within,
         "has_source_photo_equal_taxon_weight":True}
    try:
        beta,scores,corr,cond=_sandwich(x,y,w,groups)
        out["status"]="estimated"
        out["condition_number_XtWX"]=cond
        out["raw_inference"]=_inference(beta,scores,corr,names,focus)
    except ValueError as err:
        out["status"]="unidentifiable"
        out["reason"]=str(err)
        return out
    # No geographic covariates in the primary model; broad 8 spherical terms
    # are a distinct sensitivity, not evidence of climate-specific causal effects.
    if geography is not None:
        coords=geography.reindex(selected.index)
        if coords.isna().any().any():raise ValueError("S1_MISSING_GEO_BY_SELECTED_SOURCE")
        # Average photo-level spherical basis per taxon; coordinate effects only
        # on between-taxon component, not within-species pseudo-identities.
        tax_means=coords.copy()
        tax_means["taxon"]=groups
        for v in geography:
            tax_means[v]=tax_means.groupby("taxon",observed=True)[v].transform("mean")
        g=tax_means[list(geography)].to_numpy(float)
        x2,n2=_columns(selected,between,within,w,g)
        try:
            bb,sc,cc,cond2=_sandwich(x2,y,w,groups)
            out["geographic_sensitivity"]={
                "status":"estimated","condition_number_XtWX":cond2,
                "terms":_inference(bb,sc,cc,n2,focus),
                "eight_spherical_photo_mean_basis_not_phylogenetic_correction":True}
        except ValueError as e:
            out["geographic_sensitivity"]={"status":"unidentifiable","reason":str(e)}
    return out

def _spherical(lat,lon):
    a=np.deg2rad(np.asarray(lat,float));b=np.deg2rad(np.asarray(lon,float))
    x=np.cos(a)*np.cos(b);y=np.cos(a)*np.sin(b);z=np.sin(a)
    return np.column_stack([x,y,z,x*y,x*z,y*z,x*x-y*y,3*z*z-1])

def analyse(frame,choice,geo=None):
    between,within=selected_axes(choice)
    if geo is not None and len(geo)!=len(frame):raise ValueError("S1_GEO_LENGTH_CHANGED")
    output=[one_trait(frame,t,between,within,geography=geo) for t in YENDPOINTS]
    # One 8 traits x (6 between + 9 within) posthoc family. Null results remain.
    flat=[row for panel in output if panel["status"]=="estimated"
          for row in panel["raw_inference"] if row["cluster_robust_joint_p"] is not None]
    if flat:
        holm=multipletests([r["cluster_robust_joint_p"] for r in flat],
                          method="holm")[1]
        for r,v in zip(flat,holm):
            r["exploratory_holm_p_in_full_fitted_family"]=float(v)
    return {"schema":"azami_v4_s1_xonly_selected_between_within_multitrait_v1",
            "status":"real_only_if_original_S1_receipt_verified",
            "inference":"exploratory_photographic_traits_not_causal_or_confirmatory",
            "n_original_native_S1":313350,
            "selected_without_traits":True,
            "selected_between_X":between,"selected_within_X":within,
            "n_tests_familywise_original_pvalues":len(flat),
            "not_original_S1_80_tests_or_H3c":True,
            "no_original_photo_taxon_observer_coordinates_exported":True,
            "results":output}

def run_original(native,decoded_receipt,seven,nine,nine_manifest,source_receipt,
                 source_sha,selected,output):
    d,receipt=load_existing_inputs(native,decoded_receipt,seven,nine,
                                   source_receipt,nine_manifest,
                                   source_receipt_sha256=source_sha)
    # Source lat/lon are used only for geographic sensitivity and never exported.
    loc=pd.read_csv(native,usecols=["obs_id","latitude","longitude"],
                    dtype={"obs_id":str},low_memory=False)
    if not loc.obs_id.is_unique or len(loc)!=313350:raise ValueError("SOURCE_COORDINATE_CHECK_FAILED")
    ids=d.obs_id.copy()
    d=d.merge(loc,on="obs_id",how="left",sort=False,validate="one_to_one")
    if not d.obs_id.equals(ids):raise ValueError("SOURCE_GEO_JOIN_DRIFT")
    source=pd.DataFrame(_spherical(d.latitude,d.longitude),
                        columns=["sphere_"+str(k) for k in range(8)],index=d.index)
    choice=json.loads(Path(selected).read_text())
    if choice.get("climate_source_receipt_sha256")!=source_sha:
        raise ValueError("BLIND_CHOICE_SOURCE_RECEIPT_MISMATCH")
    data=analyse(d,choice,geo=source)
    data["source_receipt_sha256"]=source_sha
    data["original_native_file_sha256"]=receipt["input_digests"]["native"]
    banned={"obs_id","accepted_key","taxon","taxon_name","user_id",
            "photographer_id","photo_id","latitude","longitude",
            "individual_prediction"}
    def audit(obj):
        if isinstance(obj,dict):
            for k,v in obj.items():
                if k in banned:raise ValueError("PRIVATE_S1_DATA_EXPORTED")
                audit(v)
        elif isinstance(obj,list):
            for v in obj:audit(v)
    audit(data)
    out=Path(output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(data,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    return data

def main():
    p=argparse.ArgumentParser()
    for n in ("native","decoded-receipt","seven","nine","nine-manifest",
              "source-receipt","source-sha","selected","output"):
        p.add_argument("--"+n,required=True,type=str)
    args=p.parse_args()
    r=run_original(args.native,args.decoded_receipt,args.seven,args.nine,
                   args.nine_manifest,args.source_receipt,args.source_sha,
                   args.selected,args.output)
    print(json.dumps({"status":r["status"],"n_traits":len(r["results"]),
                      "n_tests":r["n_tests_familywise_original_pvalues"],
                      "estimated":[x["trait"] for x in r["results"] if x["status"]=="estimated"]}),flush=True)
if __name__=="__main__":main()
