#!/usr/bin/env python3
"""Reapply the frozen-v2 sensitivity sequence to v3 biological constructs.

The endpoint-level v2 analysis remains canonical. This script starts only from
FDR-supported rows of the biological-construct reanalysis, then repeats the v2
sequence: sampling-composition stability -> broad/residual spatial sensitivity
-> 52-tree historical-placement sensitivity.

For scalar constructs the v2 sign-preservation rule is unchanged. For constructs
represented by >1 response component, the circular-hue vector logic is generalized
to Euclidean coefficient vectors: sampling stability is positive cosine alignment.
Spatial and historical pass rules otherwise remain the v2 rules; vector alignment
is additionally recorded but is not used to make a joint construct pass.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from analysis.v3 import run_biological_axis_reanalysis as bio
import legacy.v2.analysis.run_geb_v2_full27_spatial_sensitivity as spatial_v2
import legacy.v2.analysis.run_geb_v2_full27_historical_sensitivity as hist_v2


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--traits", type=Path, required=True)
    p.add_argument("--environment", type=Path, required=True)
    p.add_argument("--axis-among", type=Path, required=True)
    p.add_argument("--axis-within", type=Path, required=True)
    p.add_argument("--regions", type=Path, required=True)
    p.add_argument("--native-status", type=Path, required=True)
    p.add_argument("--tree-dir", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--spatial-permutations", type=int, default=999)
    p.add_argument("--moran-permutations", type=int, default=999)
    p.add_argument("--minimum-taxa-historical", type=int, default=30)
    p.add_argument("--seed", type=int, default=20260910)
    return p.parse_args()


def truth(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def standardize(a: np.ndarray, weights: np.ndarray | None = None) -> np.ndarray:
    a = np.asarray(a, float)
    if weights is None:
        mean = float(np.mean(a)); var = float(np.mean((a - mean) ** 2))
    else:
        weights = np.asarray(weights, float)
        mean = float(np.average(a, weights=weights)); var = float(np.average((a - mean) ** 2, weights=weights))
    if not np.isfinite(a).all() or not np.isfinite(var) or var <= 0:
        raise ValueError("no finite variation")
    return (a - mean) / math.sqrt(var)


def weighted_slope(y: np.ndarray, x: np.ndarray, weights: np.ndarray | None = None) -> float:
    if weights is None:
        weights = np.ones(len(x), float)
    den = float(np.dot(weights * x, x))
    if not np.isfinite(den) or den <= 0:
        raise ValueError("predictor has no finite variation")
    return float(np.dot(weights * x, y) / den)


def demean_standardize(frame: pd.DataFrame, col: str, weights: np.ndarray | None = None) -> np.ndarray:
    s = pd.to_numeric(frame[col], errors="coerce")
    c = s - s.groupby(frame["taxon_name"]).transform("mean")
    return standardize(c.to_numpy(float), weights)


def equal_taxon_weights(taxa: pd.Series) -> np.ndarray:
    counts = taxa.groupby(taxa).transform("size").to_numpy(float)
    return 1.0 / counts


def vector_alignment(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, float); b = np.asarray(b, float)
    den = float(np.linalg.norm(a) * np.linalg.norm(b))
    return float(np.dot(a, b) / den) if den > 0 and np.isfinite(den) else float("nan")


def component_beta(row: pd.Series) -> np.ndarray:
    if str(row["kind"]) == "scalar":
        return np.asarray([float(row["beta_std"])])
    obj = json.loads(str(row["component_betas"]))
    members = bio.CONSTRUCTS[str(row["construct_id"])]["members"]
    return np.asarray([float(obj[m]) for m in members], float)


def load_inputs(args: argparse.Namespace):
    t, e = bio.load_data(args.traits, args.environment)
    r = pd.read_csv(args.regions, usecols=["obs_id", "broad_region"], low_memory=False)
    n = pd.read_csv(args.native_status, usecols=["obs_id", "taxon_name", "native_range_status"], low_memory=False)
    for x in (e, r, n): x["obs_id"] = x["obs_id"].astype(str)
    if r.obs_id.duplicated().any() or n.obs_id.duplicated().any():
        raise ValueError("sampling lookup not unique by obs_id")
    if set(r.obs_id) != set(e.obs_id) or set(n.obs_id) != set(e.obs_id):
        raise ValueError("sampling lookup does not cover exact environment cohort")
    meta = e.merge(r, on="obs_id", validate="one_to_one")
    meta = meta.merge(n.rename(columns={"taxon_name":"native_taxon_name"}), on="obs_id", validate="one_to_one")
    if not meta.taxon_name.astype(str).eq(meta.native_taxon_name.astype(str)).all():
        raise ValueError("native taxon names disagree")
    meta = meta.drop(columns="native_taxon_name")
    return t, e, meta


def construct_obs(t: pd.DataFrame, e: pd.DataFrame, name: str) -> pd.DataFrame:
    d = bio.CONSTRUCTS[name]
    if d["kind"] == "scalar":
        return bio.scalar_obs_table(t, e, d)
    return bio.joint_obs_table(t, e, d)


def stable_top_taxa(e: pd.DataFrame, n: int = 10) -> list[str]:
    x = e.groupby("taxon_name", as_index=False).size().rename(columns={"size":"n"})
    x = x.sort_values(["n","taxon_name"], ascending=[False, True])
    return x.head(n).taxon_name.astype(str).tolist()


def scenarios(top: list[str], regions: list[str], scale: str) -> list[dict[str, Any]]:
    out=[]
    if scale == "within_taxon":
        out.append(dict(family="equal_taxon_weight", scenario="equal_total_weight_per_taxon", excluded_taxa=[], excluded_region="", native_only=False, equal_weight=True))
    for taxon in top:
        out.append(dict(family="dominant_taxon_omission", scenario=f"omit:{taxon}", excluded_taxa=[taxon], excluded_region="", native_only=False, equal_weight=False))
    out.append(dict(family="dominant_taxon_omission", scenario="omit_top2_joint", excluded_taxa=top[:2], excluded_region="", native_only=False, equal_weight=False))
    for region in regions:
        out.append(dict(family="leave_one_broad_region_out", scenario=f"omit_region:{region}", excluded_taxa=[], excluded_region=region, native_only=False, equal_weight=False))
    out.append(dict(family="native_only", scenario="native_only", excluded_taxa=[], excluded_region="", native_only=True, equal_weight=False))
    return out


def apply_scenario(obs: pd.DataFrame, meta: pd.DataFrame, sc: dict[str, Any]) -> tuple[pd.DataFrame,pd.DataFrame]:
    o = obs.merge(meta[["obs_id","broad_region","native_range_status"]], on="obs_id", how="left", validate="one_to_one")
    m = meta.copy()
    if sc["excluded_taxa"]:
        o=o[~o.taxon_name.isin(sc["excluded_taxa"])]; m=m[~m.taxon_name.isin(sc["excluded_taxa"])]
    if sc["excluded_region"]:
        o=o[~o.broad_region.eq(sc["excluded_region"])]; m=m[~m.broad_region.eq(sc["excluded_region"])]
    if sc["native_only"]:
        o=o[o.native_range_status.eq("native")]; m=m[m.native_range_status.eq("native")]
    return o.copy(), m.copy()


def fit_sampling_within(obs: pd.DataFrame, name: str, predictor: str, equal_weight: bool) -> dict[str, Any]:
    d=bio.CONSTRUCTS[name]; members=["trait"] if d["kind"]=="scalar" else d["members"]
    x=obs.dropna(subset=[*members,predictor]).copy()
    counts=x.groupby("taxon_name").size(); x=x[x.taxon_name.isin(counts[counts>=2].index)].copy()
    base={"n_observations":len(x),"n_taxa":x.taxon_name.nunique()}
    if len(x)<100 or x.taxon_name.nunique()<10: return {**base,"status":"insufficient_support"}
    try:
        w=equal_taxon_weights(x.taxon_name) if equal_weight else None
        xp=demean_standardize(x,predictor,w)
        betas=np.asarray([weighted_slope(demean_standardize(x,m,w),xp,w) for m in members])
        return {**base,"status":"ok","beta_vector":json.dumps(betas.tolist()),"effect_magnitude":float(np.linalg.norm(betas))}
    except ValueError as exc:
        return {**base,"status":f"failed:{exc}"}


def fit_sampling_among(obs: pd.DataFrame, meta: pd.DataFrame, name: str, predictor: str) -> dict[str, Any]:
    d=bio.CONSTRUCTS[name]; members=["trait"] if d["kind"]=="scalar" else d["members"]
    x=obs.dropna(subset=members).copy(); counts=x.groupby("taxon_name").size().rename("n_trait_observations")
    med=x.groupby("taxon_name")[members].median().join(counts); med=med[med.n_trait_observations>=5]
    env=meta.groupby("taxon_name")[[predictor]].median(numeric_only=True)
    dat=med.join(env,how="inner").dropna(subset=[*members,predictor]); base={"n_taxa":len(dat)}
    if len(dat)<20: return {**base,"status":"insufficient_support"}
    try:
        xp=standardize(dat[predictor].to_numpy(float))
        betas=np.asarray([weighted_slope(standardize(dat[m].to_numpy(float)),xp) for m in members])
        return {**base,"status":"ok","beta_vector":json.dumps(betas.tolist()),"effect_magnitude":float(np.linalg.norm(betas))}
    except ValueError as exc:
        return {**base,"status":f"failed:{exc}"}


def run_sampling(t,e,meta,among_atlas,within_atlas,out_dir:Path):
    sa=among_atlas[truth(among_atlas.fdr_0_05)].copy(); sw=within_atlas[truth(within_atlas.fdr_0_05)].copy()
    top=stable_top_taxa(e,10); regions=sorted(x for x in meta.broad_region.dropna().astype(str).unique() if x!="UNMAPPED")
    rows=[]
    for scale,selected in (("within_taxon",sw),("among_taxon",sa)):
        for _,r in selected.iterrows():
            name=str(r.construct_id); pred=str(r.predictor); baseline=component_beta(r); obs=construct_obs(t,e,name)
            for sc in scenarios(top,regions,scale):
                oscope,mscope=apply_scenario(obs,meta,sc)
                fit=fit_sampling_within(oscope,name,pred,sc["equal_weight"]) if scale=="within_taxon" else fit_sampling_among(oscope,mscope,name,pred)
                row={"scale":scale,"construct_id":name,"kind":r.kind,"predictor":pred,"baseline_beta_vector":json.dumps(baseline.tolist()),"baseline_effect_magnitude":float(np.linalg.norm(baseline)),"sensitivity_family":sc["family"],"scenario":sc["scenario"],**fit}
                if fit.get("status")=="ok":
                    b=np.asarray(json.loads(fit["beta_vector"]),float); align=vector_alignment(baseline,b)
                    row["vector_alignment_cosine"]=align; row["direction_stable"]=bool(np.isfinite(align) and align>0)
                    row["effect_magnitude_ratio_to_baseline"]=float(np.linalg.norm(b)/np.linalg.norm(baseline))
                else: row["direction_stable"]=False
                rows.append(row)
    detail=pd.DataFrame(rows); sums=[]
    for keys,p in detail.groupby(["scale","construct_id","predictor"],sort=True):
        ok=p.status.eq("ok"); stable=p.loc[ok,"direction_stable"].astype(bool)
        sums.append({"scale":keys[0],"construct_id":keys[1],"predictor":keys[2],"n_scenarios":len(p),"n_evaluable":int(ok.sum()),"all_declared_scenarios_evaluable":bool(ok.all()),"all_directions_stable_where_evaluable":bool(len(stable) and stable.all()),"minimum_effect_magnitude_ratio":float(pd.to_numeric(p.loc[ok,"effect_magnitude_ratio_to_baseline"],errors="coerce").min()) if ok.any() else np.nan,"sampling_composition_stability_class":"stable_all_declared_scenarios" if ok.all() and len(stable) and stable.all() else ("direction_unstable" if len(stable) and (~stable).any() else "stable_where_evaluable_incomplete")})
    summary=pd.DataFrame(sums); out_dir.mkdir(parents=True,exist_ok=True)
    detail.to_csv(out_dir/"biological_axis_sampling_scenarios.csv",index=False); summary.to_csv(out_dir/"biological_axis_sampling_summary.csv",index=False)
    return summary


def taxon_construct_data(t,e,name:str,pred:str) -> tuple[pd.DataFrame,list[str]]:
    d=bio.CONSTRUCTS[name]; members=["trait"] if d["kind"]=="scalar" else d["members"]
    obs=construct_obs(t,e,name).dropna(subset=members)
    counts=obs.groupby("taxon_name").size().rename("n_trait_observations")
    med=obs.groupby("taxon_name")[members].median().join(counts); med=med[med.n_trait_observations>=5]
    env=e.copy(); basis=spatial_v2.spherical_basis(env.latitude,env.longitude)
    for j in range(8): env[f"spatial_basis_{j}"]=basis[:,j]
    cols=[pred,"latitude","longitude",*[f"spatial_basis_{j}" for j in range(8)]]
    envmed=env.groupby("taxon_name")[cols].median(numeric_only=True)
    return med.join(envmed,how="inner").reset_index(),members


def spatial_one(t,e,r:pd.Series,seed:int,perms:int,moran_perms:int,within:bool) -> dict[str,Any]:
    name=str(r.construct_id); pred=str(r.predictor); d=bio.CONSTRUCTS[name]; baseline=component_beta(r)
    if within:
        dat=construct_obs(t,e,name); members=["trait"] if d["kind"]=="scalar" else d["members"]
        dat=dat.merge(e[["obs_id","latitude","longitude"]],on="obs_id",how="left",validate="one_to_one").dropna(subset=[*members,pred,"latitude","longitude"])
        counts=dat.groupby("taxon_name").size(); dat=dat[dat.taxon_name.isin(counts[counts>=2].index)].copy(); taxa=dat.taxon_name.astype(str).to_numpy()
        basis=spatial_v2.spherical_basis(dat.latitude,dat.longitude); design,reduced=spatial_v2.prepare_design(dat[pred].to_numpy(float),basis,taxa)
        if d["kind"]=="scalar": responses=spatial_v2.standardize(spatial_v2.demean_by_taxon(dat[["trait"]].to_numpy(float),taxa)[:,0])
        else:
            raw=spatial_v2.demean_by_taxon(dat[members].to_numpy(float),taxa); responses=np.column_stack([spatial_v2.standardize(raw[:,j]) for j in range(raw.shape[1])])
        groups=taxa
    else:
        dat,members=taxon_construct_data(t,e,name,pred); bcols=[f"spatial_basis_{j}" for j in range(8)]; dat=dat.dropna(subset=[*members,pred,*bcols]); design,reduced=spatial_v2.prepare_design(dat[pred].to_numpy(float),dat[bcols].to_numpy(float),None)
        responses=spatial_v2.standardize(dat[members[0]].to_numpy(float)) if d["kind"]=="scalar" else np.column_stack([spatial_v2.standardize(dat[m].to_numpy(float)) for m in members]); groups=None
    rng=spatial_v2.stable_rng(seed,"bio_axis_spatial","within" if within else "among",name,pred)
    if d["kind"]=="scalar":
        beta,pv,resid=spatial_v2.linear_freedman_lane(responses,design,reduced,groups,perms,rng); betavec=np.asarray([beta]); same=bool(np.sign(beta)==np.sign(baseline[0]))
    else:
        betavec,pv=spatial_v2.multivariate_freedman_lane(responses,design,reduced,groups,perms,rng); resid=responses-design@np.linalg.lstsq(design,responses,rcond=None)[0]; resid=np.linalg.norm(resid,axis=1); same=True
    mrng=spatial_v2.stable_rng(seed,"bio_axis_moran","within" if within else "among",name,pred); mi,mp,mn=spatial_v2.moran_test(resid,dat.latitude.to_numpy(float),dat.longitude.to_numpy(float),moran_perms,5000,mrng)
    align=vector_alignment(baseline,betavec)
    return {"scale":"within_taxon" if within else "among_taxon","construct_id":name,"kind":d["kind"],"predictor":pred,"n":len(dat),"baseline_beta_vector":json.dumps(baseline.tolist()),"spatial_beta_vector":json.dumps(np.asarray(betavec,float).tolist()),"spatial_effect_magnitude":float(np.linalg.norm(betavec)),"spatial_permutation_p_value":pv,"vector_alignment_cosine":align,"residual_morans_i":mi,"residual_morans_p_value":mp,"residual_morans_n":mn,"broad_spatial_sensitivity_pass":bool(pv<.05 and same and np.isfinite(mp) and mp>=.05)}


def run_spatial(t,e,among_atlas,within_atlas,out_dir:Path,seed:int,perms:int,moran_perms:int):
    sa=among_atlas[truth(among_atlas.fdr_0_05)].copy(); sw=within_atlas[truth(within_atlas.fdr_0_05)].copy(); rows=[]
    for _,r in sw.iterrows(): rows.append(spatial_one(t,e,r,seed,perms,moran_perms,True))
    for _,r in sa.iterrows(): rows.append(spatial_one(t,e,r,seed,perms,moran_perms,False))
    df=pd.DataFrame(rows); out_dir.mkdir(parents=True,exist_ok=True); df.to_csv(out_dir/"biological_axis_spatial.csv",index=False)
    return df


def run_historical(t,e,spatial_df,tree_dir:Path,out_dir:Path,minimum_taxa:int):
    selected=spatial_df[(spatial_df.scale=="among_taxon") & truth(spatial_df.broad_spatial_sensitivity_pass)].copy(); trees=[]
    for scenario,rep,tree in hist_v2.load_trees(tree_dir):
        names,cov=hist_v2.tree_covariance(tree); trees.append((scenario,rep,names,cov))
    rows=[]
    for _,r in selected.iterrows():
        name=str(r.construct_id); pred=str(r.predictor); d=bio.CONSTRUCTS[name]; dat,members=taxon_construct_data(t,e,name,pred); dat["taxon_key"]=dat.taxon_name.str.strip().str.replace(" ","_",regex=False); dat=dat.set_index("taxon_key",drop=False); spatial_beta=np.asarray(json.loads(r.spatial_beta_vector),float)
        for scenario,rep,names,cov in trees:
            lookup={x:i for i,x in enumerate(names)}; common=[x for x in names if x in dat.index]; base={"construct_id":name,"kind":d["kind"],"predictor":pred,"scenario":scenario,"replicate":rep,"n_taxa":len(common)}
            if len(common)<minimum_taxa: rows.append({**base,"status":"insufficient_taxa"}); continue
            ordered=dat.loc[common]; pos=[lookup[x] for x in common]; pcov=cov[np.ix_(pos,pos)]
            try:
                fit=hist_v2.fit_pagel(ordered[members].to_numpy(float),ordered[pred].to_numpy(float),pcov); beta=np.asarray(fit.pop("beta"),float); se=fit.pop("standard_error"); rows.append({**base,"status":"ok","pgls_beta_vector":json.dumps(beta.tolist()),"pgls_effect_magnitude":float(np.linalg.norm(beta)),"pgls_standard_error":se if np.isfinite(se) else np.nan,"vector_alignment_with_spatial":vector_alignment(spatial_beta,beta),**fit})
            except Exception as exc: rows.append({**base,"status":f"failed:{type(exc).__name__}:{exc}"})
    models=pd.DataFrame(rows); sums=[]
    if len(models):
        for keys,p in models.groupby(["construct_id","kind","predictor"],sort=True):
            ok=p[p.status.eq("ok")]; linear=keys[1]=="scalar"; dir_ok=bool((pd.to_numeric(ok.vector_alignment_with_spatial,errors="coerce")>0).all()) if linear and len(ok) else (not linear); sig=bool(len(ok) and pd.to_numeric(ok.p_value,errors="coerce").lt(.05).all()); sums.append({"construct_id":keys[0],"kind":keys[1],"predictor":keys[2],"n_successful_placement_trees":len(ok),"n_placement_trees_p_lt_0_05":int(pd.to_numeric(ok.p_value,errors="coerce").lt(.05).sum()) if len(ok) else 0,"minimum_p_value":float(pd.to_numeric(ok.p_value,errors="coerce").min()) if len(ok) else np.nan,"maximum_p_value":float(pd.to_numeric(ok.p_value,errors="coerce").max()) if len(ok) else np.nan,"lambda_min":float(pd.to_numeric(ok['lambda'],errors="coerce").min()) if len(ok) else np.nan,"lambda_max":float(pd.to_numeric(ok['lambda'],errors="coerce").max()) if len(ok) else np.nan,"direction_stable_across_trees":dir_ok if linear else np.nan,"historical_placement_sensitivity_pass":bool(len(ok)==52 and sig and dir_ok)})
    summary=pd.DataFrame(sums); out_dir.mkdir(parents=True,exist_ok=True); models.to_csv(out_dir/"biological_axis_historical_models.csv",index=False); summary.to_csv(out_dir/"biological_axis_historical_summary.csv",index=False)
    return summary


def main() -> int:
    args=parse_args(); args.out_dir.mkdir(parents=True,exist_ok=True); t,e,meta=load_inputs(args); among=pd.read_csv(args.axis_among,low_memory=False); within=pd.read_csv(args.axis_within,low_memory=False)
    sampling=run_sampling(t,e,meta,among,within,args.out_dir/"sampling")
    spatial=run_spatial(t,e,among,within,args.out_dir/"spatial",args.seed,args.spatial_permutations,args.moran_permutations)
    historical=run_historical(t,e,spatial,args.tree_dir,args.out_dir/"historical",args.minimum_taxa_historical)
    report={"analysis_id":"ch1_v3_biological_construct_v2_sensitivity_chain_20260910","starting_constructs_inferential":10,"among_fdr_rows_entered_sampling":int(truth(among.fdr_0_05).sum()),"within_fdr_rows_entered_sampling":int(truth(within.fdr_0_05).sum()),"sampling_pairs_stable_all":int((sampling.sampling_composition_stability_class=="stable_all_declared_scenarios").sum()),"spatial_among_entered":int(((spatial.scale=="among_taxon")).sum()),"spatial_among_passed":int(((spatial.scale=="among_taxon") & truth(spatial.broad_spatial_sensitivity_pass)).sum()),"spatial_within_entered":int(((spatial.scale=="within_taxon")).sum()),"spatial_within_passed":int(((spatial.scale=="within_taxon") & truth(spatial.broad_spatial_sensitivity_pass)).sum()),"historical_pairs_entered":int(len(historical)),"historical_pairs_passed":int(truth(historical.historical_placement_sensitivity_pass).sum()) if len(historical) else 0,"claim_boundary":"construct-level robustness layer only; frozen v2 endpoint-level conclusions remain canonical"}
    (args.out_dir/"biological_axis_sensitivity_chain_report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n",encoding="utf-8"); print(json.dumps(report,indent=2)); return 0

if __name__=="__main__": raise SystemExit(main())
