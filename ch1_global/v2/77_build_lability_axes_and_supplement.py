#!/usr/bin/env python3
"""Final two-axis lability analysis for Chapter 1.

Separates species-level within-species phenotypic variation from species-specific
climate responsiveness. Circular hue is handled as a joint sine/cosine endpoint.
"""
from __future__ import annotations
import argparse, json, math
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PREDICTORS = ["chelsa_bio01", "chelsa_bio04", "chelsa_bio12", "chelsa_bio15"]
LINEAR_TRAITS = [
    "orientation_angle_degrees_median", "corolla_lab_lightness_median",
    "corolla_lab_chroma_median", "shape_aspect_ratio_median",
    "shape_circularity_median", "shape_solidity_median", "shape_width_cv_median",
]
HUE_SIN, HUE_COS, HUE = "corolla_hue_sin_median", "corolla_hue_cos_median", "corolla_hue_angle"
MODULE = {
    "orientation_angle_degrees_median": "orientation",
    "corolla_lab_lightness_median": "colour", "corolla_lab_chroma_median": "colour",
    HUE: "colour", "shape_aspect_ratio_median": "shape",
    "shape_circularity_median": "shape", "shape_solidity_median": "shape",
    "shape_width_cv_median": "shape",
}

def parse_args():
    p=argparse.ArgumentParser()
    p.add_argument("--observations", required=True)
    p.add_argument("--coefficients", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--min-species-n", type=int, default=10)
    p.add_argument("--min-traits", type=int, default=6)
    p.add_argument("--min-predictors", type=int, default=3)
    return p.parse_args()

def numeric(s): return pd.to_numeric(s, errors="coerce")
def mad_scale(s):
    x=numeric(s); med=x.median(); sc=1.4826*(x-med).abs().median()
    if not np.isfinite(sc) or sc<=0: sc=x.std(ddof=1)
    return float(sc) if np.isfinite(sc) and sc>0 else np.nan

def hue_angle(df): return np.degrees(np.arctan2(numeric(df[HUE_SIN]), numeric(df[HUE_COS])))%360
def circ_dist(x,c): return np.abs((x-c+180)%360-180)
def circ_center(x):
    r=np.radians(np.asarray(x,float)); return float(np.degrees(np.arctan2(np.nanmean(np.sin(r)),np.nanmean(np.cos(r))))%360)

def variation_detail(obs,min_n):
    obs=obs.copy(); obs[HUE]=hue_angle(obs); rows=[]
    for trait in LINEAR_TRAITS:
        sc=mad_scale(obs[trait])
        if not np.isfinite(sc): continue
        for taxon,g in obs.groupby("taxon_name",sort=True):
            x=numeric(g[trait]).dropna()
            if len(x)<min_n: continue
            raw=float(np.median(np.abs(x-np.median(x))))
            rows.append(dict(taxon_name=taxon,trait=trait,module=MODULE[trait],n_variation=len(x),within_variation_raw=raw,within_variation_standardized=raw/sc))
    allh=numeric(obs[HUE]).dropna(); denom=np.nan
    if len(allh):
        denom=float(np.median(circ_dist(allh,circ_center(allh))))
        if not np.isfinite(denom) or denom<=0: denom=float(np.std(circ_dist(allh,circ_center(allh)),ddof=1))
    if np.isfinite(denom) and denom>0:
        for taxon,g in obs.groupby("taxon_name",sort=True):
            x=numeric(g[HUE]).dropna()
            if len(x)<min_n: continue
            raw=float(np.median(circ_dist(x,circ_center(x))))
            rows.append(dict(taxon_name=taxon,trait=HUE,module="colour",n_variation=len(x),within_variation_raw=raw,within_variation_standardized=raw/denom))
    return pd.DataFrame(rows)

def linear_slope(y,x):
    w=pd.DataFrame({"y":numeric(y),"x":numeric(x)}).dropna()
    if len(w)<3 or w.x.nunique()<2 or w.y.nunique()<2: return None
    xs=float(w.x.std(ddof=1)); ys=float(w.y.std(ddof=1))
    if xs<=0 or ys<=0: return None
    xz=(w.x-w.x.mean())/xs; yz=(w.y-w.y.mean())/ys
    beta=float(np.dot(xz,yz)/np.dot(xz,xz)); resid=yz-beta*xz
    se=float(np.sqrt(np.dot(resid,resid)/(len(w)-1)/np.dot(xz,xz))) if len(w)>2 else np.nan
    z=beta/se if np.isfinite(se) and se>0 else np.nan
    p=math.erfc(abs(z)/math.sqrt(2)) if np.isfinite(z) else np.nan
    return beta,se,p,len(w),float(w.x.max()-w.x.min())

def hue_slope(g,pred,min_n):
    w=g[[HUE_SIN,HUE_COS,pred]].apply(numeric).dropna()
    if len(w)<min_n or w[pred].nunique()<2: return None
    xs=float(w[pred].std(ddof=1)); joint=math.sqrt(float(w[HUE_SIN].var(ddof=1)+w[HUE_COS].var(ddof=1)))
    if xs<=0 or joint<=0: return None
    xz=(w[pred]-w[pred].mean())/xs
    ys=(w[HUE_SIN]-w[HUE_SIN].mean())/joint; yc=(w[HUE_COS]-w[HUE_COS].mean())/joint
    den=float(np.dot(xz,xz)); bs=float(np.dot(xz,ys)/den); bc=float(np.dot(xz,yc)/den)
    return dict(beta_sin=bs,beta_cos=bc,effect_magnitude=math.hypot(bs,bc),effect_direction_degrees=math.degrees(math.atan2(bs,bc))%360,n_observations=len(w),predictor_range=float(w[pred].max()-w[pred].min()))

def species_slopes(obs,min_n):
    rows=[]
    for taxon,g in obs.groupby("taxon_name",sort=True):
        for trait in LINEAR_TRAITS:
            for pred in PREDICTORS:
                res=linear_slope(g[trait],g[pred])
                if res is None or res[3]<min_n: continue
                b,se,p,n,rng=res
                rows.append(dict(taxon_name=taxon,trait=trait,module=MODULE[trait],predictor=pred,endpoint_type="linear",beta_std=b,se=se,p_value=p,effect_magnitude=abs(b),effect_direction_degrees=0.0 if b>=0 else 180.0,n_observations=n,predictor_range=rng))
        for pred in PREDICTORS:
            res=hue_slope(g,pred,min_n)
            if res:
                rows.append(dict(taxon_name=taxon,trait=HUE,module="colour",predictor=pred,endpoint_type="circular_joint",beta_std=np.nan,se=np.nan,p_value=np.nan,**res))
    return pd.DataFrame(rows)

def bh(p):
    p=np.asarray(p,float); out=np.full(len(p),np.nan); ok=np.isfinite(p)
    vals=p[ok]; order=np.argsort(vals); ranked=vals[order]; q=ranked*len(vals)/np.arange(1,len(vals)+1); q=np.minimum.accumulate(q[::-1])[::-1]; q=np.clip(q,0,1)
    tmp=np.empty(len(vals)); tmp[order]=q; out[np.where(ok)[0]]=tmp; return out

def summarize_axes(var,slopes,min_traits,min_preds):
    vr=var.groupby("taxon_name",as_index=False).agg(species_within_variation_index=("within_variation_standardized","mean"),n_traits_variation=("trait","nunique"),min_trait_n=("n_variation","min"))
    tr=slopes.groupby(["taxon_name","trait","module"],as_index=False).agg(trait_environmental_responsiveness=("effect_magnitude",lambda x:float(np.sqrt(np.mean(np.square(x))))),n_predictors_response=("predictor","nunique"),min_response_n=("n_observations","min"))
    tr=tr[tr.n_predictors_response>=min_preds].copy()
    sr=tr.groupby("taxon_name",as_index=False).agg(species_environmental_responsiveness_index=("trait_environmental_responsiveness",lambda x:float(np.sqrt(np.mean(np.square(x))))),n_traits_response=("trait","nunique"),min_predictors_per_trait=("n_predictors_response","min"))
    axes=vr.merge(sr,on="taxon_name",how="inner")
    axes["primary_complete"]=(axes.n_traits_variation>=min_traits)&(axes.n_traits_response>=min_traits)
    axes["quadrant_eligible"]=axes.primary_complete
    modv=var.groupby(["taxon_name","module"],as_index=False).agg(module_within_variation=("within_variation_standardized","mean"),n_traits_variation=("trait","nunique"))
    modr=tr.groupby(["taxon_name","module"],as_index=False).agg(module_environmental_responsiveness=("trait_environmental_responsiveness",lambda x:float(np.sqrt(np.mean(np.square(x))))),n_traits_response=("trait","nunique"))
    modules=modv.merge(modr,on=["taxon_name","module"],how="inner")
    return axes,tr,modules

def combine_global(coef):
    c=coef.copy(); ok=c[c.status.eq("ok")].copy() if "status" in c else c.copy(); non=ok[~ok.trait.isin([HUE_SIN,HUE_COS])].copy()
    non["trait_endpoint"]=non.trait; non["module"]=non.trait_endpoint.map(MODULE); non["effect_magnitude"]=non.beta_std_within.abs(); non["effect_direction_degrees"]=np.where(non.beta_std_within>=0,0.,180.); non["endpoint_type"]="linear"
    hrs=[]
    for pred,g in ok[ok.trait.isin([HUE_SIN,HUE_COS])].groupby("predictor"):
        v=g.set_index("trait").beta_std_within
        if HUE_SIN in v and HUE_COS in v:
            bs,bc=float(v[HUE_SIN]),float(v[HUE_COS]); hrs.append(dict(trait=HUE,trait_endpoint=HUE,module="colour",predictor=pred,status="ok",endpoint_type="circular_joint",beta_std_within=np.nan,effect_magnitude=math.hypot(bs,bc),effect_direction_degrees=math.degrees(math.atan2(bs,bc))%360,hue_beta_sin=bs,hue_beta_cos=bc,n_observations=int(g.n_observations.min()),n_species=int(g.n_species.min()),p_value=np.nan,q_fdr_bh=np.nan,fdr_significant_0_05=False))
    return pd.concat([non,pd.DataFrame(hrs)],ignore_index=True,sort=False)

def plot_heatmap(globalc,out):
    order=["orientation_angle_degrees_median","corolla_lab_lightness_median","corolla_lab_chroma_median",HUE,"shape_aspect_ratio_median","shape_circularity_median","shape_solidity_median","shape_width_cv_median"]
    mat=globalc.pivot_table(index="trait_endpoint",columns="predictor",values="beta_std_within",aggfunc="first").reindex(order).reindex(columns=PREDICTORS)
    mag=globalc.pivot_table(index="trait_endpoint",columns="predictor",values="effect_magnitude",aggfunc="first").reindex(order).reindex(columns=PREDICTORS)
    arr=mat.to_numpy(float); hi=order.index(HUE); arr[hi,:]=mag.loc[HUE].to_numpy(float)
    vmax=np.nanmax(np.abs(arr)); fig,ax=plt.subplots(figsize=(8.8,6.6)); im=ax.imshow(arr,aspect="auto",vmin=-vmax,vmax=vmax,cmap="coolwarm")
    ax.set_xticks(range(4),PREDICTORS,rotation=30,ha="right"); ax.set_yticks(range(8),[x.replace("_median","") for x in order]); ax.set_title("Figure S1. Standardized within-species climate effects")
    for i,t in enumerate(order):
        for j,p in enumerate(PREDICTORS):
            v=arr[i,j]
            if np.isfinite(v):
                row=globalc[(globalc.trait_endpoint==t)&(globalc.predictor==p)].iloc[0]; star="*" if bool(row.get("fdr_significant_0_05",False)) else ""
                txt=f"{v:.3f}{star}" if t!=HUE else f"|β|={v:.3f}"
                ax.text(j,i,txt,ha="center",va="center",fontsize=8)
    cb=fig.colorbar(im,ax=ax); cb.set_label("signed standardized β; hue row shows joint magnitude")
    fig.tight_layout(); fig.savefig(out.with_suffix('.png'),dpi=300); fig.savefig(out.with_suffix('.pdf')); plt.close(fig)

def plot_quadrants(axes,out):
    p=axes[axes.quadrant_eligible].copy(); xm=p.species_within_variation_index.median(); ym=p.species_environmental_responsiveness_index.median()
    p["quadrant"]=np.select([(p.species_within_variation_index>=xm)&(p.species_environmental_responsiveness_index>=ym),(p.species_within_variation_index<xm)&(p.species_environmental_responsiveness_index>=ym),(p.species_within_variation_index>=xm)&(p.species_environmental_responsiveness_index<ym)], ["high variation / high response","low variation / high response","high variation / low response"],default="low variation / low response")
    fig,ax=plt.subplots(figsize=(9,7)); ax.scatter(p.species_within_variation_index,p.species_environmental_responsiveness_index,s=24,alpha=.65)
    score=((p.species_within_variation_index-p.species_within_variation_index.median())/p.species_within_variation_index.std()).abs()+((p.species_environmental_responsiveness_index-p.species_environmental_responsiveness_index.median())/p.species_environmental_responsiveness_index.std()).abs()
    for idx in score.nlargest(min(12,len(score))).index:
        r=p.loc[idx]; ax.annotate(r.taxon_name.replace("Cirsium ","C. "),(r.species_within_variation_index,r.species_environmental_responsiveness_index),xytext=(3,3),textcoords="offset points",fontsize=7)
    ax.axvline(xm,ls="--",lw=1); ax.axhline(ym,ls="--",lw=1); ax.set_xlabel("Species within-variation index"); ax.set_ylabel("Species environmental responsiveness index"); ax.set_title("Species-level decomposition of capitulum trait lability")
    fig.tight_layout(); fig.savefig(out.with_suffix('.png'),dpi=300); fig.savefig(out.with_suffix('.pdf')); plt.close(fig)
    return p

def run_threshold(obs,min_n,min_traits,min_preds):
    v=variation_detail(obs,min_n); s=species_slopes(obs,min_n); a,t,m=summarize_axes(v,s,min_traits,min_preds); return v,s,a,t,m

def main():
    a=parse_args(); out=Path(a.out_dir); out.mkdir(parents=True,exist_ok=True)
    obs=pd.read_csv(a.observations,low_memory=False); coef=pd.read_csv(a.coefficients,low_memory=False)
    req={"taxon_name",*LINEAR_TRAITS,HUE_SIN,HUE_COS,*PREDICTORS}; miss=sorted(req-set(obs.columns))
    if miss: raise ValueError(f"Missing observation columns: {miss}")
    v,s,axes,tr,mods=run_threshold(obs,a.min_species_n,a.min_traits,a.min_predictors)
    s["q_within_species_bh"]=bh(s.p_value); s["fdr_significant_0_05"]=s.q_within_species_bh<.05
    axes=axes.merge(s.groupby("taxon_name",as_index=False).agg(n_species_specific_effects=("effect_magnitude","count"),n_fdr_significant_linear_effects=("fdr_significant_0_05","sum")),on="taxon_name",how="left")
    quad=plot_quadrants(axes,out/"figure_species_lability_quadrants"); axes=axes.merge(quad[["taxon_name","quadrant"]],on="taxon_name",how="left")
    globalc=combine_global(coef); plot_heatmap(globalc,out/"figure_s1_effect_size_heatmap")
    v.to_csv(out/"species_trait_within_variation.csv",index=False); s.to_csv(out/"table_s2_species_specific_coefficients.csv",index=False); axes.to_csv(out/"species_lability_axes.csv",index=False); tr.to_csv(out/"species_trait_environmental_responsiveness.csv",index=False); mods.to_csv(out/"species_module_lability.csv",index=False); globalc.to_csv(out/"table_s1_all_global_coefficients.csv",index=False)
    sens=[]; base=axes.set_index("taxon_name")
    for n in [5,10,20]:
        _,_,ax,_,_=run_threshold(obs,n,a.min_traits,a.min_predictors); ax.to_csv(out/f"species_lability_axes_min_n_{n}.csv",index=False)
        j=base[["species_within_variation_index","species_environmental_responsiveness_index"]].join(ax.set_index("taxon_name")[["species_within_variation_index","species_environmental_responsiveness_index"]],lsuffix="_primary",rsuffix="_sensitivity",how="inner")
        sens.append(dict(min_species_n=n,n_species_axes=len(ax),n_primary_complete=int(ax.primary_complete.sum()),variation_spearman=j.iloc[:,[0,2]].corr(method="spearman").iloc[0,1] if len(j)>2 else np.nan,responsiveness_spearman=j.iloc[:,[1,3]].corr(method="spearman").iloc[0,1] if len(j)>2 else np.nan))
    pd.DataFrame(sens).to_csv(out/"sensitivity_minimum_sample_size.csv",index=False)
    modsum=mods.groupby("module",as_index=False).agg(median_species_within_variation=("module_within_variation","median"),median_species_environmental_responsiveness=("module_environmental_responsiveness","median"),n_species=("taxon_name","nunique")); modsum.to_csv(out/"module_lability_summary.csv",index=False)
    meta=dict(primary_definition={"min_species_n":a.min_species_n,"min_traits":a.min_traits,"min_predictors_per_trait":a.min_predictors,"variation":"mean trait-standardized MAD; circular median absolute angular deviation for hue","responsiveness":"RMS of species-specific standardized slopes across predictors and eligible traits","hue":"joint sine/cosine vector; no scalar-angle OLS"},counts={"input_observations":len(obs),"input_species":int(obs.taxon_name.nunique()),"species_with_variation":int(v.taxon_name.nunique()),"species_with_response":int(s.taxon_name.nunique()),"species_primary_complete":int(axes.primary_complete.sum())},interpretation_boundary="Indices describe image-derived within-species variation and climate association, not demonstrated plasticity or local adaptation.")
    (out/"analysis_metadata.json").write_text(json.dumps(meta,indent=2),encoding="utf-8")
    print(json.dumps(meta,indent=2))
if __name__=="__main__": main()
