#!/usr/bin/env python3
from __future__ import annotations
import argparse, hashlib, json, math
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm

PREDICTORS=[
    'chelsa_bio01','chelsa_bio04','chelsa_bio12','chelsa_bio15',
    'chelsa_rsds_mean','chelsa_vpd_mean','chelsa_sfcwind_mean','chelsa_gsp','chelsa_npp'
]
RESTORED5={
    'visible_floret_fraction','corolla_white_pixel_fraction','corolla_redmagenta_pixel_fraction',
    'corolla_purple_pixel_fraction','corolla_yellow_pixel_fraction'
}
CONSTRUCTS={
    'presentation_angle': dict(kind='scalar', members=['orientation_image_vertical_angle'], signs=[1], tier='primary'),
    'floral_lightness': dict(kind='scalar', members=['corolla_lab_lightness'], signs=[1], tier='primary'),
    'floral_chroma': dict(kind='scalar', members=['corolla_lab_chroma'], signs=[1], tier='primary'),
    'floral_hue': dict(kind='joint', members=['corolla_hue_sin','corolla_hue_cos'], tier='primary'),
    'head_elongation': dict(kind='scalar', members=['capitulum_outline_aspect_ratio'], signs=[1], tier='primary'),
    'head_compactness': dict(kind='scalar', members=['capitulum_outline_circularity','capitulum_outline_solidity','capitulum_width_profile_cv'], signs=[1,1,-1], tier='primary'),
    'involucre_form': dict(kind='joint', members=['involucre_length_width_ratio','involucre_apical_taper_ratio','involucre_basal_taper_ratio'], tier='candidate'),
    'projection_prominence': dict(kind='scalar', members=['bract_projection_roughness','bract_projection_p95','bract_projection_maximum','bract_spread_fraction'], signs=[1,1,1,1], tier='candidate'),
    'projection_pattern': dict(kind='joint', members=['bract_projection_peak_density','bract_projection_asymmetry'], tier='candidate'),
    'surface_texture': dict(kind='joint', members=['involucre_surface_edge_density','involucre_surface_lbp_entropy','involucre_surface_high_frequency_energy'], tier='validation_only'),
    'surface_specularity': dict(kind='scalar', members=['involucre_surface_specular_fraction'], signs=[1], tier='validation_only'),
}

def stable_rng(seed,*parts):
    d=hashlib.sha256('|'.join([str(seed),*map(str,parts)]).encode()).digest()
    return np.random.default_rng(int.from_bytes(d[:8],'little'))

def z(a):
    a=np.asarray(a,float); sd=a.std(ddof=0)
    if not np.isfinite(a).all() or sd<=0: raise ValueError('no finite variation')
    return (a-a.mean())/sd

def bh(s):
    p=np.asarray(s,float); n=len(p)
    order=np.argsort(p); ranked=p[order]
    adj=ranked*n/np.arange(1,n+1); adj=np.minimum.accumulate(adj[::-1])[::-1]; adj=np.minimum(adj,1)
    out=np.empty(n); out[order]=adj
    return out

def load_data(traits_path, env_path):
    t=pd.read_csv(traits_path,low_memory=False)
    t=t[~t.endpoint_id.isin(RESTORED5)].copy()
    t['obs_id']=t.obs_id.astype(str); t['taxon_name']=t.taxon_name.astype(str)
    t['value']=pd.to_numeric(t.value,errors='coerce'); t=t[t.value.notna()].copy()
    e=pd.read_csv(env_path,low_memory=False)
    e['obs_id']=e.obs_id.astype(str); e['taxon_name']=e.taxon_name.astype(str)
    for p in PREDICTORS: e[p]=pd.to_numeric(e[p],errors='coerce')
    t=t[t.obs_id.isin(set(e.obs_id))].copy()
    return t,e

def taxon_endpoint_table(t,members,min_n=5):
    x=t[t.endpoint_id.isin(members)].groupby(['taxon_name','endpoint_id']).value.agg(['median','count']).reset_index()
    x=x[x['count']>=min_n]
    med=x.pivot(index='taxon_name',columns='endpoint_id',values='median')
    cnt=x.pivot(index='taxon_name',columns='endpoint_id',values='count')
    med=med.dropna(subset=members)
    return med[members],cnt.reindex(med.index)

def scaling_from_taxon_medians(t,members,min_n=5):
    med,_=taxon_endpoint_table(t,members,min_n)
    means=med.mean(0); sds=med.std(0,ddof=0)
    if (sds<=0).any(): raise ValueError('constant taxon-median member')
    return means,sds

def scalar_among_table(t,e,defn,min_n=5):
    members=defn['members']; signs=np.asarray(defn['signs'],float)
    med,_=taxon_endpoint_table(t,members,min_n)
    means=med.mean(0); sds=med.std(0,ddof=0)
    if (sds<=0).any(): return pd.DataFrame(), {'status':'constant_response'}
    score=((med-means)/sds).to_numpy()@signs/len(signs)
    out=pd.DataFrame({'trait':score},index=med.index)
    env=e.groupby('taxon_name')[PREDICTORS].median()
    return out.join(env,how='inner').dropna(), {'means':means.to_dict(),'sds':sds.to_dict()}

def joint_among_table(t,e,defn,min_n=5):
    members=defn['members']; med,_=taxon_endpoint_table(t,members,min_n)
    env=e.groupby('taxon_name')[PREDICTORS].median()
    return med.join(env,how='inner').dropna(),{}

def scalar_obs_table(t,e,defn):
    members=defn['members']; signs=np.asarray(defn['signs'],float)
    means,sds=scaling_from_taxon_medians(t,members,5)
    part=t[t.endpoint_id.isin(members)][['obs_id','taxon_name','endpoint_id','value']]
    wide=part.pivot(index=['obs_id','taxon_name'],columns='endpoint_id',values='value').reset_index().dropna(subset=members)
    wide['trait']=((wide[members]-means)/sds).to_numpy()@signs/len(signs)
    return wide[['obs_id','taxon_name','trait']].merge(e[['obs_id',*PREDICTORS]],on='obs_id',how='inner')

def joint_obs_table(t,e,defn):
    members=defn['members']; part=t[t.endpoint_id.isin(members)][['obs_id','taxon_name','endpoint_id','value']]
    wide=part.pivot(index=['obs_id','taxon_name'],columns='endpoint_id',values='value').reset_index().dropna(subset=members)
    return wide.merge(e[['obs_id',*PREDICTORS]],on='obs_id',how='inner')

def perm_scalar_among(frame,predictor,perms,rng):
    yz=z(frame.trait.to_numpy(float)); xz=z(frame[predictor].to_numpy(float)); den=float(np.dot(xz,xz))
    obs=float(np.dot(xz,yz)/den); exceed=0; rem=perms
    while rem:
        size=min(256,rem); mat=np.broadcast_to(xz,(size,len(xz))).copy(); xp=rng.permuted(mat,axis=1)
        sims=xp@yz/den; exceed+=int(np.sum(np.abs(sims)>=abs(obs)-1e-15)); rem-=size
    return obs,(exceed+1)/(perms+1)

def perm_joint_among(frame,members,predictor,perms,rng):
    Y=np.column_stack([z(frame[m].to_numpy(float)) for m in members]); xz=z(frame[predictor].to_numpy(float)); den=float(np.dot(xz,xz))
    beta=xz@Y/den; mag=float(np.linalg.norm(beta)); exceed=0; rem=perms
    while rem:
        size=min(256,rem); mat=np.broadcast_to(xz,(size,len(xz))).copy(); xp=rng.permuted(mat,axis=1)
        sims=xp@Y/den; exceed+=int(np.sum(np.linalg.norm(sims,axis=1)>=mag-1e-15)); rem-=size
    return beta,mag,(exceed+1)/(perms+1)

def demean_std(frame,col):
    s=pd.to_numeric(frame[col],errors='coerce'); c=s-s.groupby(frame.taxon_name).transform('mean'); sd=float(c.std(ddof=0))
    if not np.isfinite(sd) or sd<=0: raise ValueError('no within variation')
    return (c/sd).to_numpy(float)

def fit_scalar_within(frame,predictor):
    y=demean_std(frame,'trait'); x=demean_std(frame,predictor)
    fit=sm.OLS(y,x[:,None]).fit(cov_type='cluster',cov_kwds={'groups':frame.taxon_name.to_numpy()})
    return float(fit.params[0]),float(fit.bse[0]),float(fit.conf_int()[0,0]),float(fit.conf_int()[0,1]),float(fit.pvalues[0])

def perm_joint_within(frame,members,predictor,perms,rng):
    Y=np.column_stack([demean_std(frame,m) for m in members]); x=demean_std(frame,predictor); den=float(np.dot(x,x)); beta=x@Y/den; mag=float(np.linalg.norm(beta))
    taxa=frame.taxon_name.astype(str).to_numpy(); groups=[np.flatnonzero(taxa==u) for u in np.unique(taxa)]; exceed=0; rem=perms
    while rem:
        size=min(128,rem); nums=np.zeros((size,Y.shape[1]))
        for idx in groups:
            mat=np.broadcast_to(x[idx],(size,len(idx))).copy(); xp=rng.permuted(mat,axis=1); nums += xp@Y[idx]
        sims=nums/den; exceed += int(np.sum(np.linalg.norm(sims,axis=1)>=mag-1e-15)); rem-=size
    return beta,mag,(exceed+1)/(perms+1)

def eligible_within(frame,min_obs=100,min_taxa=10,min_per_taxon=2):
    counts=frame.groupby('taxon_name').size(); taxa=counts[counts>=min_per_taxon].index; x=frame[frame.taxon_name.isin(taxa)].copy()
    return x if len(x)>=min_obs and x.taxon_name.nunique()>=min_taxa else pd.DataFrame()

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--traits',type=Path,required=True); ap.add_argument('--environment',type=Path,required=True); ap.add_argument('--out-dir',type=Path,required=True); ap.add_argument('--permutations',type=int,default=9999); ap.add_argument('--seed',type=int,default=20260910)
    args=ap.parse_args(); args.out_dir.mkdir(parents=True,exist_ok=True); t,e=load_data(args.traits,args.environment)
    measured=sorted(t.endpoint_id.unique()); expected=sorted({m for d in CONSTRUCTS.values() for m in d['members']})
    if measured!=expected: raise ValueError(f'22 endpoint universe mismatch: measured={len(measured)} expected={len(expected)} extra={set(measured)-set(expected)} missing={set(expected)-set(measured)}')
    among=[]; within=[]; status=[]
    for name,d in CONSTRUCTS.items():
        aframe,_=(scalar_among_table(t,e,d,5) if d['kind']=='scalar' else joint_among_table(t,e,d,5))
        if aframe.empty:
            status.append({'construct_id':name,'scale':'among_taxon_min5','status':'no_finite_variation_or_insufficient','n':0})
        else:
            status.append({'construct_id':name,'scale':'among_taxon_min5','status':'ok','n':len(aframe)})
            for p in PREDICTORS:
                rng=stable_rng(args.seed,'among',name,p)
                if d['kind']=='scalar':
                    beta,pv=perm_scalar_among(aframe,p,args.permutations,rng); among.append({'construct_id':name,'tier':d['tier'],'kind':'scalar','predictor':p,'n_taxa':len(aframe),'beta_std':beta,'effect_magnitude':abs(beta),'component_betas':'','p_value':pv})
                else:
                    beta,mag,pv=perm_joint_among(aframe,d['members'],p,args.permutations,rng); among.append({'construct_id':name,'tier':d['tier'],'kind':'joint','predictor':p,'n_taxa':len(aframe),'beta_std':np.nan,'effect_magnitude':mag,'component_betas':json.dumps({m:float(b) for m,b in zip(d['members'],beta)},sort_keys=True),'p_value':pv})
        if name == 'surface_specularity':
            status.append({'construct_id':name,'scale':'within_taxon','status':'descriptive_only_not_inferential','n':0}); continue
        wframe=(scalar_obs_table(t,e,d) if d['kind']=='scalar' else joint_obs_table(t,e,d)); wframe=eligible_within(wframe)
        if wframe.empty:
            status.append({'construct_id':name,'scale':'within_taxon','status':'insufficient_or_no_variation','n':0}); continue
        status.append({'construct_id':name,'scale':'within_taxon','status':'ok','n':len(wframe),'n_taxa':wframe.taxon_name.nunique()})
        for p in PREDICTORS:
            wf=wframe.dropna(subset=[p]).copy()
            if len(wf)<100 or wf.taxon_name.nunique()<10: continue
            rng=stable_rng(args.seed,'within',name,p)
            if d['kind']=='scalar':
                try:
                    b,se,lo,hi,pv=fit_scalar_within(wf,p); within.append({'construct_id':name,'tier':d['tier'],'kind':'scalar','predictor':p,'n_observations':len(wf),'n_taxa':wf.taxon_name.nunique(),'beta_std':b,'standard_error':se,'ci_low':lo,'ci_high':hi,'effect_magnitude':abs(b),'component_betas':'','p_value':pv,'p_method':'taxon_clustered_ols'})
                except ValueError: pass
            else:
                try:
                    beta,mag,pv=perm_joint_within(wf,d['members'],p,args.permutations,rng); within.append({'construct_id':name,'tier':d['tier'],'kind':'joint','predictor':p,'n_observations':len(wf),'n_taxa':wf.taxon_name.nunique(),'beta_std':np.nan,'standard_error':np.nan,'ci_low':np.nan,'ci_high':np.nan,'effect_magnitude':mag,'component_betas':json.dumps({m:float(b) for m,b in zip(d['members'],beta)},sort_keys=True),'p_value':pv,'p_method':'within_taxon_predictor_permutation'})
                except ValueError: pass
    a=pd.DataFrame(among); w=pd.DataFrame(within); st=pd.DataFrame(status)
    if len(a): a['q_bh']=bh(a.p_value); a['fdr_0_05']=a.q_bh<.05
    if len(w): w['q_bh']=bh(w.p_value); w['fdr_0_05']=w.q_bh<.05
    a.to_csv(args.out_dir/'biological_axes_among_min5.csv',index=False); w.to_csv(args.out_dir/'biological_axes_within.csv',index=False); st.to_csv(args.out_dir/'biological_axes_status.csv',index=False)
    report={'analysis_id':'ch1_v3_biological_axis_reanalysis_20260910','starting_measured_endpoints':22,'constructs':len(CONSTRUCTS),'testable_among_constructs':int(st.query("scale=='among_taxon_min5' and status=='ok'").shape[0]),'testable_within_constructs':int(st.query("scale=='within_taxon' and status=='ok'").shape[0]),'among_tests':len(a),'within_tests':len(w),'among_fdr_hits':int(a.fdr_0_05.sum()) if len(a) else 0,'within_fdr_hits':int(w.fdr_0_05.sum()) if len(w) else 0,'among_hits':a.loc[a.fdr_0_05,['construct_id','predictor','beta_std','effect_magnitude','p_value','q_bh']].replace({np.nan:None}).to_dict('records') if len(a) else [],'within_hits':w.loc[w.fdr_0_05,['construct_id','predictor','beta_std','effect_magnitude','p_value','q_bh']].replace({np.nan:None}).to_dict('records') if len(w) else [],'claim_boundary':'exploratory biological aggregation of frozen v2 measured endpoints; v2 endpoint-level family remains canonical'}
    (args.out_dir/'biological_axes_report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n'); print(json.dumps(report,indent=2,allow_nan=False))

if __name__=='__main__': main()
