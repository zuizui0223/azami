#!/usr/bin/env python3
from __future__ import annotations
import argparse, hashlib, json, math
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
DEFAULT_REPO=Path(__file__).resolve().parents[2]
P=['chelsa_bio01','chelsa_bio04','chelsa_bio12','chelsa_bio15','chelsa_rsds_mean','chelsa_vpd_mean','chelsa_sfcwind_mean','chelsa_gsp','chelsa_npp']
SEED=20260827
B=9999
BOOT=1000

def parse_args():
 p=argparse.ArgumentParser(description='Post-hoc VIFstep + simultaneous multivariable sensitivity for frozen GEB-v2 among-taxon atlas.')
 p.add_argument('--repo-root',type=Path,default=DEFAULT_REPO)
 p.add_argument('--traits-long',type=Path,required=True)
 p.add_argument('--environment',type=Path,required=True)
 p.add_argument('--out-dir',type=Path,required=True)
 p.add_argument('--frozen-among',type=Path,default=None)
 p.add_argument('--permutations',type=int,default=9999)
 p.add_argument('--circular-bootstraps',type=int,default=1000)
 return p.parse_args()

ARGS=parse_args()
REPO=ARGS.repo_root.resolve()
if ARGS.frozen_among is None:
 ARGS.frozen_among=REPO/'analysis_outputs/v2_full27_environment_atlas_2026-08-27/v2_full27_environment_among.csv'
TRAITS=ARGS.traits_long
ENV=ARGS.environment
OUT=ARGS.out_dir
OUT.mkdir(parents=True,exist_ok=True)
B=ARGS.permutations
BOOT=ARGS.circular_bootstraps

def standardize(v):
 v=np.asarray(v,float); sd=np.std(v,ddof=0)
 if not np.isfinite(sd) or sd<=0: raise ValueError('No finite variation')
 return (v-np.mean(v))/sd

def bh(v):
 v=np.asarray(v,float); out=np.full(len(v),np.nan); ok=np.isfinite(v)
 if not ok.any(): return out
 p=v[ok]; order=np.argsort(p); r=p[order]*len(p)/np.arange(1,len(p)+1); r=np.minimum.accumulate(r[::-1])[::-1]; r=np.clip(r,0,1); z=np.empty(len(p)); z[order]=r; out[np.flatnonzero(ok)]=z; return out

def vifs(X):
 X=X.dropna().astype(float); Z=(X-X.mean())/X.std(ddof=0); out={}
 for c in Z.columns:
  y=Z[c].to_numpy(); others=[d for d in Z.columns if d!=c]; D=np.column_stack([np.ones(len(Z)),Z[others].to_numpy()]) if others else np.ones((len(Z),1)); b=np.linalg.lstsq(D,y,rcond=None)[0]; e=y-D@b; den=float(y@y)
  if not np.isfinite(den) or den<=0: out[c]=float('inf'); continue
  R2=1-float(e@e)/den; out[c]=float('inf') if R2>=1 else 1/(1-R2)
 return out

def vifstep(X,thr):
 cur=list(X.columns); hist=[]
 while len(cur)>1:
  vv=vifs(X[cur]); mx=max(vv,key=vv.get)
  if vv[mx] <= thr: return cur,hist,vv
  hist.append((mx,vv[mx],cur.copy())); cur.remove(mx)
 return cur,hist,vifs(X[cur])

def rng(*parts):
 d=hashlib.sha256('|'.join(map(str,(SEED,*parts))).encode()).digest(); return np.random.default_rng(int.from_bytes(d[:8],'little'))

def rx_target(X,j):
 oth=[k for k in range(X.shape[1]) if k!=j]; x=X[:,j]
 if not oth:return x.copy()
 D=np.column_stack([np.ones(len(X)),X[:,oth]]); return x-D@np.linalg.lstsq(D,x,rcond=None)[0]

def fl_linear(y,X,j,R):
 oth=[k for k in range(X.shape[1]) if k!=j]; D=np.column_stack([np.ones(len(X)),X[:,oth]]) if oth else np.ones((len(X),1)); fit=D@np.linalg.lstsq(D,y,rcond=None)[0]; e=y-fit; rx=rx_target(X,j); den=float(rx@X[:,j]); obs=float(rx@y/den); exc=0
 for st in range(0,B,512):
  n=min(512,B-st); mat=np.vstack([R.permutation(e) for _ in range(n)]); sim=mat@rx/den; exc+=int(np.sum(np.abs(sim)>=abs(obs)-1e-15))
 return obs,(exc+1)/(B+1)

def fl_circ(ys,yc,X,j,R):
 oth=[k for k in range(X.shape[1]) if k!=j]; D=np.column_stack([np.ones(len(X)),X[:,oth]]) if oth else np.ones((len(X),1)); es=ys-D@np.linalg.lstsq(D,ys,rcond=None)[0]; ec=yc-D@np.linalg.lstsq(D,yc,rcond=None)[0]; rx=rx_target(X,j); den=float(rx@X[:,j]); bs=float(rx@ys/den); bc=float(rx@yc/den); mag=math.hypot(bs,bc); exc=0; nobs=len(ys)
 for st in range(0,B,512):
  n=min(512,B-st); a=np.empty((n,nobs)); b=np.empty((n,nobs))
  for i in range(n): idx=R.permutation(nobs); a[i]=es[idx]; b[i]=ec[idx]
  exc+=int(np.sum(np.hypot(a@rx/den,b@rx/den)>=mag-1e-15))
 return bs,bc,mag,(exc+1)/(B+1)

def boot_circ(ys,yc,X,j,R,N=1000):
 n=len(ys); vals=[]
 for _ in range(N):
  idx=R.integers(0,n,n); D=sm.add_constant(X[idx])
  try:
   bs=np.linalg.lstsq(D,ys[idx],rcond=None)[0][1+j]; bc=np.linalg.lstsq(D,yc[idx],rcond=None)[0][1+j]; vals.append((bs,bc,math.hypot(bs,bc)))
  except Exception: pass
 if len(vals)<100:return (np.nan,)*6
 q=np.quantile(np.array(vals),[.025,.975],axis=0); return q[0,0],q[1,0],q[0,1],q[1,1],q[0,2],q[1,2]

contract=pd.read_csv(REPO/'ch1_global/v2/ontology/ch1_continuous_trait_contract.csv',dtype=str,keep_default_na=False); contract['circular_group']=contract.circular_group.fillna('').str.strip(); units=[]
for _,r in contract[contract.circular_group.eq('')].iterrows(): units.append(dict(unit_id=r.endpoint_id,members=[r.endpoint_id],module=r.module,tier=r.analysis_tier,validation=r.validation_status,kind='linear'))
for g,p in contract[contract.circular_group.ne('')].groupby('circular_group'): units.append(dict(unit_id=g,members=p.endpoint_id.tolist(),module=p.iloc[0].module,tier=p.iloc[0].analysis_tier,validation=p.iloc[0].validation_status,kind='circular'))
units=sorted(units,key=lambda x:(x['module'],x['unit_id']))
env=pd.read_csv(ENV,low_memory=False); env['obs_id']=env.obs_id.astype(str)
for c in P: env[c]=pd.to_numeric(env[c],errors='coerce')
envt=env.groupby('taxon_name')[P].median()
tr=pd.read_csv(TRAITS,usecols=['obs_id','taxon_name','endpoint_id','measurement_available','value'],low_memory=False); tr['obs_id']=tr.obs_id.astype(str); tr['value']=pd.to_numeric(tr.value,errors='coerce'); tr=tr[tr.measurement_available.astype(str).str.lower().isin(['true','1','yes']) & tr.value.notna()]

def unit_tax(u,minobs):
 p=tr[tr.endpoint_id.isin(u['members'])][['obs_id','taxon_name','endpoint_id','value']]
 if p.empty or set(p.endpoint_id)!=set(u['members']):return pd.DataFrame()
 w=p.pivot(index=['obs_id','taxon_name'],columns='endpoint_id',values='value').reset_index(); w.columns.name=None; w=w.dropna(subset=u['members']); c=w.groupby('taxon_name').size().rename('n_trait_observations'); m=w.groupby('taxon_name')[u['members']].median().join(c).reset_index(); return m[m.n_trait_observations>=minobs]

rows=[]; vrows=[]
for minobs,scope in [(5,'among_taxon_min5'),(2,'among_taxon_min2')]:
 for u in units:
  tt=unit_tax(u,minobs)
  if tt.empty:
   for thr in [10,5]: rows.append(dict(scope=scope,vif_threshold=thr,unit_id=u['unit_id'],member_endpoint_ids='|'.join(u['members']),module=u['module'],analysis_tier=u['tier'],validation_status=u['validation'],inferential_unit=u['kind'],predictor='',status='unexecuted_no_measurement',n_taxa=0))
   continue
  dat=tt.merge(envt.reset_index(),on='taxon_name',how='inner',validate='one_to_one').dropna(subset=[*u['members'],*P]).copy(); n=len(dat)
  if n<20:
   for thr in [10,5]: rows.append(dict(scope=scope,vif_threshold=thr,unit_id=u['unit_id'],member_endpoint_ids='|'.join(u['members']),module=u['module'],analysis_tier=u['tier'],validation_status=u['validation'],inferential_unit=u['kind'],predictor='',status='insufficient_support',n_taxa=n))
   continue
  for thr in [10,5]:
   keep,hist,final=vifstep(dat[P],thr)
   for step,(rem,vv,before) in enumerate(hist,1): vrows.append(dict(scope=scope,vif_threshold=thr,unit_id=u['unit_id'],n_taxa=n,step=step,action='remove',predictor=rem,vif=vv,predictors_before='|'.join(before)))
   for p,v in final.items(): vrows.append(dict(scope=scope,vif_threshold=thr,unit_id=u['unit_id'],n_taxa=n,step=len(hist)+1,action='retain_final',predictor=p,vif=v,predictors_before='|'.join(keep)))
   try:
    X=np.column_stack([standardize(dat[p].to_numpy(float)) for p in keep])
    if u['kind']=='linear': y=standardize(dat[u['members'][0]].to_numpy(float))
   except ValueError as e:
    for pred in keep: rows.append(dict(scope=scope,vif_threshold=thr,unit_id=u['unit_id'],member_endpoint_ids='|'.join(u['members']),module=u['module'],analysis_tier=u['tier'],validation_status=u['validation'],inferential_unit=u['kind'],predictor=pred,status=f'failed:ValueError:{e}',n_taxa=n,retained_predictors='|'.join(keep),final_max_vif=max(final.values())))
    continue
   if u['kind']=='linear':
    fit=sm.OLS(y,sm.add_constant(X)).fit(cov_type='HC3')
    for j,pred in enumerate(keep):
     beta,pv=fl_linear(y,X,j,rng('endpoint_vif',scope,u['unit_id'],pred,'|'.join(keep))); ci=fit.conf_int()[1+j]
     rows.append(dict(scope=scope,vif_threshold=thr,unit_id=u['unit_id'],member_endpoint_ids='|'.join(u['members']),module=u['module'],analysis_tier=u['tier'],validation_status=u['validation'],inferential_unit='linear_endpoint',predictor=pred,status='ok',n_taxa=n,retained_predictors='|'.join(keep),n_retained_predictors=len(keep),final_max_vif=max(final.values()),beta_std_adjusted=beta,standard_error_hc3=fit.bse[1+j],confidence_low_95=ci[0],confidence_high_95=ci[1],p_value=pv,p_value_method=f'Freedman-Lane residual permutation {B} + HC3 CI'))
   else:
    sine=next(m for m in u['members'] if 'sin' in m); cosine=next(m for m in u['members'] if 'cos' in m); ys=standardize(dat[sine]); yc=standardize(dat[cosine])
    for j,pred in enumerate(keep):
     bs,bc,mag,pv=fl_circ(ys,yc,X,j,rng('endpoint_vif_circ',scope,u['unit_id'],pred,'|'.join(keep))); q=boot_circ(ys,yc,X,j,rng('endpoint_vif_boot',scope,u['unit_id'],pred,'|'.join(keep)),BOOT)
     rows.append(dict(scope=scope,vif_threshold=thr,unit_id=u['unit_id'],member_endpoint_ids='|'.join(u['members']),module=u['module'],analysis_tier=u['tier'],validation_status=u['validation'],inferential_unit='circular_joint',predictor=pred,status='ok',n_taxa=n,retained_predictors='|'.join(keep),n_retained_predictors=len(keep),final_max_vif=max(final.values()),beta_sine_std_adjusted=bs,beta_cosine_std_adjusted=bc,effect_magnitude_adjusted=mag,effect_direction_degrees_adjusted=math.degrees(math.atan2(bs,bc))%360,beta_sine_ci_low_95=q[0],beta_sine_ci_high_95=q[1],beta_cosine_ci_low_95=q[2],beta_cosine_ci_high_95=q[3],effect_magnitude_ci_low_95=q[4],effect_magnitude_ci_high_95=q[5],p_value=pv,p_value_method=f'paired Freedman-Lane residual permutation {B}; taxon bootstrap {BOOT} CI'))

res=pd.DataFrame(rows); pd.DataFrame(vrows).to_csv(OUT/'endpoint_cohort_vif_step_history.csv',index=False,float_format='%.12g'); res['q_fdr_bh_global_family']=np.nan
for (_,thr),idx in res.groupby(['scope','vif_threshold']).groups.items(): res.loc[list(idx),'q_fdr_bh_global_family']=bh(pd.to_numeric(res.loc[list(idx),'p_value'],errors='coerce'))
res['fdr_significant_0_05']=pd.to_numeric(res.q_fdr_bh_global_family,errors='coerce')<.05; res.to_csv(OUT/'endpoint_cohort_multivariable_vif_adjusted.csv',index=False,float_format='%.12g'); res[(res.status=='ok')&res.fdr_significant_0_05].to_csv(OUT/'endpoint_cohort_multivariable_fdr_signals.csv',index=False,float_format='%.12g')
fro=pd.read_csv(ARGS.frozen_among); orig=fro[(fro.scope=='among_taxon_min5')&(pd.to_numeric(fro.q_fdr_bh_global_family,errors='coerce')<.05)]
cmp=[]
for _,o in orig.iterrows():
 for thr in [10,5]:
  q=res[(res.scope=='among_taxon_min5')&(res.vif_threshold==thr)&(res.unit_id==o.unit_id)&(res.predictor==o.predictor)]
  if len(q):
   r=q.iloc[0]; z={'unit_id':o.unit_id,'predictor':o.predictor,'vif_threshold':thr,'original_beta_std':o.get('beta_std',np.nan),'original_effect_magnitude':o.get('effect_magnitude',np.nan),'original_effect_direction_degrees':o.get('effect_direction_degrees',np.nan),'original_p_value':o.p_value,'original_q_fdr':o.q_fdr_bh_global_family,'adjusted_status':r.status,'adjusted_n_taxa':r.n_taxa,'retained_predictors':r.get('retained_predictors',''),'final_max_vif':r.get('final_max_vif',np.nan),'adjusted_beta_std':r.get('beta_std_adjusted',np.nan),'adjusted_ci_low_95':r.get('confidence_low_95',np.nan),'adjusted_ci_high_95':r.get('confidence_high_95',np.nan),'adjusted_effect_magnitude':r.get('effect_magnitude_adjusted',np.nan),'adjusted_effect_direction_degrees':r.get('effect_direction_degrees_adjusted',np.nan),'adjusted_p_value':r.get('p_value',np.nan),'adjusted_q_fdr':r.get('q_fdr_bh_global_family',np.nan),'adjusted_fdr_significant_0_05':bool(r.get('fdr_significant_0_05',False))}
   if o.inferential_unit=='linear_endpoint' and np.isfinite(z['adjusted_beta_std']): z['sign_concordant']=np.sign(float(o.beta_std))==np.sign(float(z['adjusted_beta_std']))
   elif o.inferential_unit=='circular_joint' and np.isfinite(z['adjusted_effect_direction_degrees']):
    d=abs(float(z['adjusted_effect_direction_degrees'])-float(o.effect_direction_degrees))%360; z['direction_difference_degrees']=min(d,360-d)
   cmp.append(z)
  else:
   uq=res[(res.scope=='among_taxon_min5')&(res.vif_threshold==thr)&(res.unit_id==o.unit_id)]
   st='predictor_removed_by_endpoint_vifstep' if len(uq) and (uq.status=='ok').any() else (uq.status.iloc[0] if len(uq) else 'not_modelled')
   cmp.append({'unit_id':o.unit_id,'predictor':o.predictor,'vif_threshold':thr,'original_beta_std':o.get('beta_std',np.nan),'original_effect_magnitude':o.get('effect_magnitude',np.nan),'original_effect_direction_degrees':o.get('effect_direction_degrees',np.nan),'original_p_value':o.p_value,'original_q_fdr':o.q_fdr_bh_global_family,'adjusted_status':st,'adjusted_fdr_significant_0_05':False})
cmp=pd.DataFrame(cmp); cmp.to_csv(OUT/'original_10_endpoint_cohort_vif_comparison.csv',index=False,float_format='%.12g')
report={'analysis_id':'geb_v2_full27_among_endpoint_cohort_vifstep_multivariable_posthoc_v1','status':'strict_posthoc_sensitivity_not_primary_replacement','inputs':{'traits_long':str(TRAITS),'environment':str(ENV),'environment_rows':int(len(env)),'environment_taxa':int(env.taxon_name.nunique()),'frozen_among':str(ARGS.frozen_among)},'vif_design':'VIFstep is recomputed inside each response-specific taxon cohort before fitting the simultaneous model. This guarantees final VIF <= threshold for every successful model cohort.','thresholds':[10,5],'inference':{'linear':f'simultaneous standardized OLS; Freedman-Lane residual permutation {B}; HC3 CI','circular':f'simultaneous sine/cosine; paired Freedman-Lane residual permutation {B}; {BOOT} taxon bootstrap CI','multiplicity':'BH across all successful tested unit-predictor rows separately by among scope and VIF threshold'},'counts':{},'original_min5_fdr_signals':len(orig),'claim_boundary':'Retrospective response-cohort-specific collinearity sensitivity; not a replacement for frozen marginal atlas and not causal.'}
for scope in ['among_taxon_min5','among_taxon_min2']:
 for thr in [10,5]:
  q=res[(res.scope==scope)&(res.vif_threshold==thr)]; ok=q[q.status=='ok']; report['counts'][f'{scope}_vif{thr}']={'successful_test_rows':int(len(ok)),'successful_units':int(ok.unit_id.nunique()),'fdr_signals':int(ok.fdr_significant_0_05.sum()),'max_final_vif':float(ok.final_max_vif.max()) if len(ok) else None,'median_predictors_retained':float(ok.groupby('unit_id').n_retained_predictors.first().median()) if len(ok) else None}
for thr in [10,5]:
 q=cmp[cmp.vif_threshold==thr]; report['counts'][f'original_10_vif{thr}']={'removed_by_endpoint_vifstep':int((q.adjusted_status=='predictor_removed_by_endpoint_vifstep').sum()),'tested':int((q.adjusted_status=='ok').sum()),'fdr_survivors':int(q.adjusted_fdr_significant_0_05.fillna(False).sum())}
(OUT/'endpoint_cohort_report.json').write_text(json.dumps(report,indent=2)+'\n'); print(json.dumps(report,indent=2))
