"""Observation-level within/between models with joint species-cluster inference.

Each response has a Gaussian random species intercept and the same eight climate
terms. Response-specific missingness is retained. Stacked cluster scores provide
joint colour/projection tests; this is marginal working-covariance inference,
not an unrestricted multivariate random-effects likelihood.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import f, t
from statsmodels.stats.multitest import multipletests

ENV = ['chelsa_bio01','chelsa_rsds_mean','chelsa_vpd_mean','chelsa_bio12']
TERMS = ['intercept']+[v+s for v in ENV for s in ('__between','__within')]
MODULES = {
    'colour': ['corolla_lab_lightness','corolla_lab_chroma','corolla_hue_sin','corolla_hue_cos'],
    'orientation': ['orientation_image_vertical_angle'],
    'bract_projection': ['bract_projection_p95','bract_spread_fraction'],
}


def fit_response(x, y, groups):
    """Exact random-intercept REML using sums of crossproducts, not mean responses."""
    x,y,groups = np.asarray(x,float),np.asarray(y,float),np.asarray(groups)
    n,p=x.shape
    labels,inverse=np.unique(groups,return_inverse=True)
    g=len(labels)
    if g<=p+2 or np.linalg.matrix_rank(x)<p or np.std(y)<=1e-12:
        return {'status':'not_identifiable','n':n,'n_species':g}
    sizes=np.bincount(inverse).astype(float)
    sx=np.zeros((g,p));sy=np.zeros(g)
    np.add.at(sx,inverse,x);np.add.at(sy,inverse,y)
    xx=np.zeros((g,p,p));xy=np.zeros((g,p))
    for j in range(p):
        np.add.at(xy[:,j],inverse,x[:,j]*y)
        for k in range(p):
            np.add.at(xx[:,j,k],inverse,x[:,j]*x[:,k])
    yy=float(y@y)
    def evaluate(ratio, details=False):
        a=ratio/(1+ratio*sizes)
        info=xx.sum(axis=0)-np.einsum('g,gi,gj->ij',a,sx,sx)
        rhs=xy.sum(axis=0)-sx.T@(a*sy)
        beta=np.linalg.solve(info,rhs)
        sse=yy-np.sum(a*sy*sy)-beta@rhs
        sign,logdet=np.linalg.slogdet(info)
        if sign<=0 or sse<=0:
            return np.inf
        objective=(n-p)*np.log(sse/(n-p))+np.log1p(ratio*sizes).sum()+logdet
        if not details:return objective
        inv=np.linalg.inv(info)
        score=xy-np.einsum('gij,j->gi',xx,beta)-sx*(a*(sy-sx@beta))[:,None]
        influence=score@inv
        influence*=np.sqrt(g/(g-1)*(n-1)/(n-p))
        return beta,influence,sse/(n-p)
    opt=minimize_scalar(lambda log_ratio:evaluate(np.exp(log_ratio)),
                        bounds=(-16,16),method='bounded',options={'xatol':1e-8})
    if not opt.success or not np.isfinite(opt.fun):
        return {'status':'not_converged','n':n,'n_species':g}
    ratio=float(np.exp(opt.x)) if opt.fun<evaluate(0.) else 0.
    beta,influence,sigma2=evaluate(ratio,True)
    return {'status':'estimated','n':n,'n_species':g,'labels':labels,
            'beta':beta,'influence':influence,'variance_ratio':ratio,
            'residual_variance':float(sigma2),
            'variance_boundary':bool(ratio==0 or abs(opt.x)>15.99),
            'objective':float(evaluate(ratio))}


def fit_module(frame, outcomes, scales):
    finite_environment=np.isfinite(frame[ENV]).all(axis=1)
    labels=np.sort(frame.accepted_key.unique())
    label_index={label:i for i,label in enumerate(labels)}
    fits=[];influences=[];components={}
    for outcome in outcomes:
        eligible=finite_environment & frame[outcome+'__eligible'].eq(True) & np.isfinite(frame[outcome])
        d=frame.loc[eligible]
        x=np.column_stack([np.ones(len(d))]+[
            ((d[v+s]-(scales[v]['mean'] if s=='__between' else 0))/scales[v]['sd']).to_numpy()
            for v in ENV for s in ('__between','__within')])
        y=d[outcome].to_numpy(float)
        center=float(y.mean()) if len(y) else 0.
        scale=float(y.std()) if len(y) else 0.
        fit=fit_response(x,(y-center)/scale if scale>0 else y,d.accepted_key.to_numpy())
        components[outcome]={k:v for k,v in fit.items() if k not in ('labels','beta','influence')}
        if fit['status']!='estimated':
            return {'status':'not_identifiable','components':components}
        influence=np.zeros((len(labels),len(TERMS)))
        influence[[label_index[a] for a in fit['labels']]]=fit['influence']
        influences.append(influence);fits.append(fit)
        se=np.sqrt(np.sum(fit['influence']**2,axis=0))
        critical=t.ppf(.975,fit['n_species']-1)
        components[outcome]['effects']={
            term:{'standardized_beta':float(fit['beta'][j]),
                  'raw_units_per_environment_sd':float(fit['beta'][j]*scale),
                  'raw_ci95':[float((fit['beta'][j]-critical*se[j])*scale),
                              float((fit['beta'][j]+critical*se[j])*scale)]}
            for j,term in enumerate(TERMS) if j}
        components[outcome]['response_sd']=scale
        components[outcome]['response_mean']=center
    scores=np.concatenate(influences,axis=1)
    covariance=scores.T@scores
    coefficients=np.concatenate([a['beta'] for a in fits])
    tests={}
    denominator=min(a['n_species'] for a in fits)-1
    for j,term in enumerate(TERMS[1:],1):
        ix=[i*len(TERMS)+j for i in range(len(outcomes))]
        b=coefficients[ix];c=covariance[np.ix_(ix,ix)]
        sd=np.sqrt(np.maximum(np.diag(c),0))
        corr=c/np.outer(sd,sd) if np.all(sd>0) else np.zeros_like(c)
        if np.linalg.matrix_rank(corr)<len(outcomes):
            tests[term]={'status':'not_identifiable'}
        else:
            statistic=float((b/sd)@np.linalg.solve(corr,b/sd)/len(outcomes))
            tests[term]={'status':'estimated','F':statistic,'df1':len(outcomes),
                         'df2':denominator,'p':float(f.sf(statistic,len(outcomes),denominator))}
    return {'status':'estimated','outcomes':outcomes,'components':components,'tests':tests,
            'inference':'approximate stacked marginal REML species-cluster Wald F',
            'random_effects':'response-specific species intercept; no phylogenetic correction',
            'pointwise_component_intervals_not_joint_significance_tests':True}


def execute(table, receipt, output):
    report=json.loads(Path(receipt).read_text())
    payload=Path(table).read_bytes()
    if hashlib.sha256(payload).hexdigest()!=report['prepared_table_sha256']:
        raise ValueError('Prepared input SHA mismatch')
    frame=pd.read_csv(table,dtype={'obs_id':str},low_memory=False)
    if len(frame)!=report['native_observations'] or not frame.obs_id.is_unique:
        raise ValueError('Observation population changed')
    scales={v:{'mean':float(frame[v].mean()),'sd':float(frame[v].std(ddof=0))} for v in ENV}
    if any(a['sd']<=0 for a in scales.values()):raise ValueError('Constant climate predictor')
    ranking=frame.accepted_key.value_counts().index
    result={'input':report,'scales':scales,'models':{},
            'model_code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'primary_family':'24 joint tests: 3 modules x 4 predictors x within/between; Holm',
            'scope':'retrospective observational associations; not adaptation or measured defence',
            'limitations':['working random-intercept covariance and conditional cluster inference',
                          'no broad spatial or phylogenetic adjustment in this primary stage',
                          'legacy image colour is not calibrated reflectance or anthocyanin',
                          'bract projection and spread fraction are not anatomical opening angles',
                          'missing measurements are not imputed or inferred to be absent traits']}
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    for drop in (0,2,5):
        d=frame.loc[~frame.accepted_key.isin(ranking[:drop])]
        for module,outcomes in MODULES.items():
            print(f'Fitting {module}, exclude {drop} most represented native taxa',flush=True)
            fit=fit_module(d,outcomes,scales)
            fit['excluded_taxa_count']=drop
            result['models'][f'{module}/drop{drop}']=fit
    family=[test for name,model in result['models'].items() if name.endswith('/drop0')
            for test in model.get('tests',{}).values()]
    if len(family)!=24 or any(test['status']!='estimated' for test in family):
        result['primary_status']='incomplete_inference'
    else:
        for test,q in zip(family,multipletests([test['p'] for test in family],method='holm')[1]):
            test['holm_p']=float(q)
        result['primary_status']='estimated'
    (output/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    rows=[]
    for name,model in result['models'].items():
        for term,test in model.get('tests',{}).items():
            rows.append({'model':name,'term':term,**test})
    pd.DataFrame(rows).to_csv(output/'joint_tests.csv',index=False)
    rows=[]
    for name,model in result['models'].items():
        for outcome,component in model.get('components',{}).items():
            for term,effect in component.get('effects',{}).items():
                rows.append({'model':name,'outcome':outcome,'term':term,
                             'n':component['n'],'n_species':component['n_species'],**effect})
    pd.DataFrame(rows).to_csv(output/'component_effects.csv',index=False)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('table','receipt','output'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();execute(a.table,a.receipt,a.output)
