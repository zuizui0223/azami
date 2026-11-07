"""Original V4 S1 joint eight-trait sparse PLS2 climate exploration.

One jointly photo-QC-eligible taxon-mean row per taxon. X: 22 existing climate
axes (BIO01-19, RSDS, VPD, wind), NO elevation or terrain/soil in this fit.
Y: eight constructs in nine response columns; circular hue is a unit sin/cos
pair and shares one unit of standardized response variance. Taxon-, sector-
and inner-CV holdouts are entirely fitted within their training partitions.

Uses explicit *X-sparse PLS2*, an alternating covariance singular-vector
algorithm with a soft-thresholded X weight vector and a shared unpenalized
multivariate Y score (not the R mixOmics implementation). Deflation is PLS2
regression deflation. This is exploratory observational prediction, not causal
inference, hypothesis confirmation, or evidence of biological adaptation.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import linalg
from sklearn.model_selection import KFold, RepeatedKFold

from analysis.v4.species_environment_scan import FEATURES, PREDICTORS
from analysis.v4.species_environment_existing_inputs import load_existing_inputs

CLIMATE = tuple(p for p in PREDICTORS if p != 'topo_elevation')
TRAITS = tuple(FEATURES)
YCOLS = tuple(c for t in TRAITS for c in FEATURES[t])
OUTCOMES = {t: tuple(YCOLS.index(c) for c in FEATURES[t]) for t in TRAITS}
assert len(CLIMATE) == 22 and len(YCOLS) == 9 and len(OUTCOMES) == 8
MIN_PHOTOS = 5
MIN_TAXA = 40
SEED = 20261107
K_GRID = (2, 4, 7, 11, 16, 22)
H_GRID = (1, 2)


def aggregate_species_joint(d:pd.DataFrame)->tuple[pd.DataFrame,dict]:
    """Private taxon means, preserving original source-photo joint 8-trait QC."""
    required={'obs_id','accepted_key','longitude',*YCOLS,
              *(x+'__eligible' for x in YCOLS)}
    if not required <= set(d):raise ValueError('JOINT_ORIGINAL_S1_COLUMNS_MISSING')
    if d.obs_id.isna().any() or not d.obs_id.is_unique or d.accepted_key.isna().any():
        raise ValueError('ORIGINAL_S1_IDENTITIES_INVALID')
    d=d.sort_values('obs_id',key=lambda s:s.astype(str),kind='mergesort').reset_index(drop=True)
    y=d[list(YCOLS)].apply(pd.to_numeric,errors='coerce').to_numpy(float)
    eligible=np.isfinite(y).all(axis=1)
    for col in YCOLS:eligible &= d[col+'__eligible'].eq(True).fillna(False).to_numpy(bool)
    hue=OUTCOMES['hue']
    norm=np.linalg.norm(y[:,hue],axis=1)
    eligible &= norm>=.20
    good=d.loc[eligible,['accepted_key','longitude',*CLIMATE] if set(CLIMATE)<=set(d) else ['accepted_key','longitude',*[p for p in CLIMATE if p in d]]].copy()
    # Do not infer missing external climate columns as zero.
    for p in CLIMATE:
        if p not in good:good[p]=np.nan
    for i,col in enumerate(YCOLS):good[col]=y[eligible,i]
    z=good[list(YCOLS)].to_numpy(float)
    z[:,hue]=z[:,hue]/np.linalg.norm(z[:,hue],axis=1,keepdims=True)
    good[list(YCOLS)]=z
    good['taxon']=good['accepted_key'].astype(str)
    for p in CLIMATE:good[p]=pd.to_numeric(good[p],errors='coerce').replace([np.inf,-np.inf],np.nan)
    lon=pd.to_numeric(good.longitude,errors='coerce').to_numpy(float)
    good['_sector']=np.where(~np.isfinite(lon)|(np.abs(lon)>180),-1,np.where(lon< -30,0,np.where(lon<60,1,2))).astype(int)
    grouped=good.groupby('taxon',sort=True,observed=True)
    means=grouped[list(YCOLS)+list(CLIMATE)].mean()
    counts=grouped.size()
    means=means.loc[counts>=MIN_PHOTOS].copy()
    # Each exposure needs >=5 actual photo locations to supply its taxon mean.
    exposure_counts=grouped[list(CLIMATE)].count().reindex(means.index)
    for p in CLIMATE:means.loc[exposure_counts[p]<MIN_PHOTOS,p]=np.nan
    means['n_photos']=counts.reindex(means.index)
    means['sector']=grouped['_sector'].agg(lambda x:int(x[x>=0].mode().iloc[0]) if (x>=0).any() else -1).reindex(means.index)
    huevec=means[list(FEATURES['hue'])].to_numpy(float)
    hue_norm=np.linalg.norm(huevec,axis=1)
    means=means.loc[hue_norm>=.2].copy()
    if len(means):
        v=means[list(FEATURES['hue'])].to_numpy(float)
        means.loc[:,list(FEATURES['hue'])]=v/np.linalg.norm(v,axis=1,keepdims=True)
    means=means.reset_index().rename(columns={'taxon':'_taxon_private'})
    if len(means)<MIN_TAXA:raise ValueError('JOINT_SPARSE_PLS2_INSUFFICIENT_TAXA')
    cov=[{'predictor':p,'n_taxa_with_finite_species_mean':int(np.isfinite(means[p].to_numpy(float)).sum()),
          'fraction_supported':float(np.isfinite(means[p].to_numpy(float)).mean())} for p in CLIMATE]
    return means,{'original_source_observations':int(len(d)),'all_eight_joint_photo_eligible':int(eligible.sum()),
                  'n_joint_supported_taxa':int(len(means)),
                  'source_photo_contributions_to_taxon_means':int(means.n_photos.sum()),
                  'sector_taxa_counts':{str(int(k)):int(v) for k,v in means.sector.value_counts().sort_index().items()},
                  'predictor_support':cov,'min_joint_QC_photos_per_taxon':MIN_PHOTOS,
                  'not_equal_to_individual_plants':True}


def _sparse_weight(cov:np.ndarray, keep:int, init:np.ndarray|None=None)->np.ndarray:
    """Rank-one X-L1 soft-threshold, <=keep nonzero features, stable tie breaking."""
    p,q=cov.shape
    if not 1<=keep<=p:raise ValueError('INVALID_KEEP_X')
    if init is None:
        _,_,vt=np.linalg.svd(cov,full_matrices=False)
        u=vt[0]
    else:u=np.asarray(init,float).copy()
    if np.linalg.norm(u)<1e-12:u=np.eye(q)[0]
    u=u/np.linalg.norm(u)
    w=np.zeros(p)
    for _ in range(100):
        z=cov@u
        if keep==p:v=z
        else:
            order=np.lexsort((np.arange(p),-np.abs(z)))
            a=np.abs(z[order[keep-1]])
            next_a=np.abs(z[order[keep]])
            threshold=(a+next_a)/2 if a>next_a+1e-12 else a*0.99
            v=np.zeros_like(z)
            support=order[:keep]
            v[support]=np.sign(z[support])*np.maximum(np.abs(z[support])-threshold,0)
        norm=np.linalg.norm(v)
        if norm<1e-12:raise ValueError('SPARSE_COMPONENT_NO_COVARIANCE')
        w=v/norm
        u_new=cov.T@w
        u_norm=np.linalg.norm(u_new)
        if u_norm<1e-12:raise ValueError('PLS2_RESPONSE_SCORE_DEGENERATE')
        u_new=u_new/u_norm
        if np.dot(u_new,u)<0:u_new=-u_new
        if np.linalg.norm(u_new-u)<1e-9:
            break
        u=u_new
    return w


def _components(x:np.ndarray,y:np.ndarray,keep:int,components:int):
    X=np.asarray(x,float).copy();Y=np.asarray(y,float).copy()
    p=X.shape[1]
    W=[];P=[];Q=[]
    for _ in range(components):
        cov=X.T@Y
        if np.linalg.norm(cov)<1e-11:break
        try:w=_sparse_weight(cov,keep)
        except ValueError:break
        t=X@w
        ss=float(t@t)
        if ss<1e-12:break
        load_x=X.T@t/ss;load_y=Y.T@t/ss
        W.append(w);P.append(load_x);Q.append(load_y)
        X-=np.outer(t,load_x);Y-=np.outer(t,load_y)
    if not W:raise ValueError('SPARSE_PLS2_NO_IDENTIFIABLE_COMPONENT')
    W=np.column_stack(W);P=np.column_stack(P);Q=np.column_stack(Q)
    coeff=W@np.linalg.pinv(P.T@W,rcond=1e-10)@Q.T
    return coeff,W,P,Q


class SparsePLS2:
    def __init__(self,keep:int,components:int):
        self.keep=int(keep);self.components=int(components)
    def fit(self,X:np.ndarray,Y:np.ndarray):
        X=np.asarray(X,float);Y=np.asarray(Y,float)
        if X.ndim!=2 or Y.ndim!=2 or X.shape[0]!=Y.shape[0] or X.shape[0]<10:
            raise ValueError('SPLS2_FIT_GRAIN_INVALID')
        if not np.isfinite(Y).all() or np.isinf(X).any():raise ValueError('SPLS2_NONFINITE')
        if not 1<=self.keep<=X.shape[1] or not 1<=self.components<=3:
            raise ValueError('SPLS2_HYPERPARAM_INVALID')
        finite=np.isfinite(X)
        if not finite.any(axis=0).all():raise ValueError('TRAIN_FOLD_ALL_MISSING_PREDICTOR')
        self.med=np.nanmedian(np.where(finite,X,np.nan),axis=0)
        xi=np.where(finite,X,self.med[None,:]);self.xmean=xi.mean(0)
        self.xstd=xi.std(0,ddof=1);self.xstd[self.xstd<1e-10]=1.
        xs=(xi-self.xmean)/self.xstd
        self.ymean=Y.mean(0);self.yscale=Y.std(0,ddof=1)
        hue=OUTCOMES['hue'];self.yscale[list(hue)]=np.sqrt(np.sum(self.yscale[list(hue)]**2))
        self.yscale[self.yscale<1e-10]=1.
        ys=(Y-self.ymean)/self.yscale
        self.beta,self.weights,self.xloadings,self.yloadings=_components(xs,ys,self.keep,self.components)
        self.selected=np.flatnonzero((np.abs(self.weights)>1e-12).any(axis=1))
        return self
    def predict(self,X:np.ndarray)->np.ndarray:
        X=np.asarray(X,float)
        if X.ndim!=2 or X.shape[1]!=len(self.med):raise ValueError('SPLS2_PREDICT_DIMENSION_INVALID')
        xi=np.where(np.isfinite(X),X,self.med[None,:])
        return ((xi-self.xmean)/self.xstd)@self.beta*self.yscale+self.ymean
    def response_loss(self,X:np.ndarray,Y:np.ndarray)->tuple[float,float,np.ndarray,np.ndarray]:
        pred=self.predict(X)
        sq=np.square((pred-Y)/self.yscale)
        baseline=np.square((Y-self.ymean)/self.yscale)
        return float(sq.sum()),float(baseline.sum()),sq.sum(axis=0),baseline.sum(axis=0)


def tune(X:np.ndarray,Y:np.ndarray,seed:int,*,k_grid=K_GRID,h_grid=H_GRID,inner_folds=4):
    inner=KFold(n_splits=min(inner_folds,len(X)//8),shuffle=True,random_state=seed)
    records=[]
    for h in h_grid:
        for k in k_grid:
            l=0.;base=0.
            for train,test in inner.split(X):
                try:
                    model=SparsePLS2(k,h).fit(X[train],Y[train])
                    a,b,_,_=model.response_loss(X[test],Y[test]);l+=a;base+=b
                except (ValueError,np.linalg.LinAlgError):
                    l=np.inf;break
            skill=float(1-l/base) if np.isfinite(l) and base>1e-12 else -np.inf
            records.append({'keepX':k,'n_components':h,'inner_skill':skill})
    allowed=[r for r in records if np.isfinite(r['inner_skill'])]
    if not allowed:raise ValueError('NO_INNER_CV_MODEL_ESTIMABLE')
    selected=sorted(allowed,key=lambda x:(-x['inner_skill'],x['keepX'],x['n_components']))[0]
    return selected,records


def outer_cv(X:np.ndarray,Y:np.ndarray,*,seed=SEED,folds=None,k_grid=K_GRID,h_grid=H_GRID):
    if folds is None:
        iterator=RepeatedKFold(n_splits=5,n_repeats=2,random_state=seed)
        folds=[('taxon_holdout',tr,te) for tr,te in iterator.split(X)]
    score=0.;baseline=0.;subscore=np.zeros(Y.shape[1]);subbaseline=np.zeros(Y.shape[1]);
    weights=Counter();choices=Counter();details=[]
    for fold,(label,tr,te) in enumerate(folds):
        chosen,_=tune(X[tr],Y[tr],seed+fold*101,k_grid=k_grid,h_grid=h_grid)
        model=SparsePLS2(chosen['keepX'],chosen['n_components']).fit(X[tr],Y[tr])
        a,b,aa,bb=model.response_loss(X[te],Y[te]);score+=a;baseline+=b;subscore+=aa;subbaseline+=bb
        for f in model.selected:weights[int(f)]+=1
        choices[(chosen['keepX'],chosen['n_components'])]+=1
        details.append({'validation_partition':str(label),'n_heldout_taxa':int(len(te)),
                        'keepX':int(chosen['keepX']),'n_components':int(chosen['n_components']),
                        'skill':float(1-a/b) if b>1e-12 else None})
    skills={t:float(1-subscore[list(idx)].sum()/subbaseline[list(idx)].sum())
            if subbaseline[list(idx)].sum()>1e-12 else None for t,idx in OUTCOMES.items()}
    return {'skill_all_eight':float(1-score/baseline),'skill_by_construct':skills,
            'baseline_scaled_loss':float(baseline),'model_scaled_loss':float(score),
            'n_outer_partitions':len(folds),'selected_frequency':[weights[i]/len(folds) for i in range(X.shape[1])],
            'selection_n_folds':{f'keep{k}_components{h}':int(n) for (k,h),n in sorted(choices.items())},
            'partition_details':details,'nonindependent_repeated_CV_partitions':True}


def run_spls2(frame:pd.DataFrame,*,seed=SEED,bootstrap=100,k_grid=K_GRID,h_grid=H_GRID)->dict:
    taxa,coverage=aggregate_species_joint(frame)
    X=taxa[list(CLIMATE)].to_numpy(float)
    Y=taxa[list(YCOLS)].to_numpy(float)
    if (np.isinf(X).any() or not np.isfinite(Y).all()):raise ValueError('NONFINITE_JOINT_SPECIES_MEANS')
    # Retain all 22 registered predictors, training-only imputation in each fold.
    if not np.isfinite(X).any(axis=0).all():raise ValueError('ORIGINAL_SOURCE_CLIMATE_AXIS_MISSING')
    cv=outer_cv(X,Y,seed=seed,k_grid=k_grid,h_grid=h_grid)
    final_choice,grid=tune(X,Y,seed+999,k_grid=k_grid,h_grid=h_grid)
    fitted=SparsePLS2(final_choice['keepX'],final_choice['n_components']).fit(X,Y)
    rng=np.random.default_rng(seed+841);bootstrap_freq=np.zeros(X.shape[1]);usable=0
    for _ in range(bootstrap):
        ix=rng.integers(0,len(X),size=len(X))
        try:
            m=SparsePLS2(final_choice['keepX'],final_choice['n_components']).fit(X[ix],Y[ix])
            bootstrap_freq[m.selected]+=1;usable+=1
        except ValueError:continue
    if usable==0:raise ValueError('NO_BOOTSTRAP_FIT')
    # Correlation clusters remain DESCRIPTIVE (no test leakage); selection done inside CV.
    reg=taxa['sector'].to_numpy(int)
    geo={'status':'not_estimable_geography'}
    if set(reg)=={0,1,2} and min(np.bincount(reg,minlength=3))>=10:
        folds=[(f'sector_{k}',np.flatnonzero(reg!=k),np.flatnonzero(reg==k)) for k in (0,1,2)]
        geo={'status':'estimated_three_dominant_longitude_sectors',
             **outer_cv(X,Y,seed=seed+222,folds=folds,k_grid=k_grid,h_grid=h_grid)}
    axes=[{'component':int(j+1),'selected_X':int(np.sum(np.abs(fitted.weights[:,j])>1e-12)),
           'x_weights':[{'predictor':CLIMATE[i],'weight':float(fitted.weights[i,j])}
                        for i in range(len(CLIMATE)) if abs(fitted.weights[i,j])>1e-12],
           'joint_Y_loadings':{col:float(fitted.yloadings[i,j]) for i,col in enumerate(YCOLS)}}
          for j in range(fitted.weights.shape[1])]
    importance=[{'predictor':name,
                 'outer_selection_fraction':float(cv['selected_frequency'][i]),
                 'bootstrap_selection_fraction_conditional_on_full_data_tuning':float(bootstrap_freq[i]/usable),
                 'full_fit_component_count':int(np.sum(np.abs(fitted.weights[i,:])>1e-12))}
                for i,name in enumerate(CLIMATE)]
    importance.sort(key=lambda r:(-r['outer_selection_fraction'],
                                  -r['bootstrap_selection_fraction_conditional_on_full_data_tuning'],r['predictor']))
    out={'schema':'azami_v4_original_s1_species_joint_spls2_climate22_v1',
         'status':'genuine_source_only_if_authenticated_runner','n_registered_constructs':8,
         'response_columns':list(YCOLS),'predictors':list(CLIMATE),
         'uses_soil':False,'includes_elevation_or_terrain':False,
         'estimator':'X_sparse_PLS2_soft_threshold_iterative_covariance_multivariate_Y_deflation',
         'not_miXOmics_software':True,
         'tuning':'inner_four_fold_keepX_2_4_7_11_16_22_components_1_2_outer_fivefold_twice',
         'climate_is_observation_site_mean_not_germination_or_flowering_microclimate':True,
         'photo_trait_orientation_is_camera_frame_not_gravitational_nodding':True,
         'phenotype_color_is_photo_proxy_not_pigment_chemistry':True,
         'no_source_taxon_or_observer_location_ids_in_output':True,
         'no_new_confirmatory_p_values':True,'causality_or_adaptation_established':False,
         'species_support':coverage,
         'taxon_holdout_nested_cv':cv,'three_sector_holdout_nested_cv':geo,
         'final_fit':{'hyperparameters':final_choice,
                      'in_sample_loss':fitted.response_loss(X,Y)[0],
                      'axes':axes},
         'environment_selection_stability':importance,
         'bootstrap_success_count':usable,'bootstrap_draws_attempted':bootstrap,
         'sparsity_tuning_search':grid,
         'prior_GEB_S1_S2_H3c_and_RF_frozen':True}
    return out


def run_existing_spls2(native:Path,decoded_receipt:Path,seven:Path,nine:Path,
                         source_receipt:Path,nine_manifest:Path,output:Path,*,
                         source_receipt_sha256:str,expected_n=313350,
                         bootstrap=100)->dict:
    inputs=[native,decoded_receipt,seven,nine,source_receipt,nine_manifest]
    if Path(output).resolve() in {Path(p).resolve() for p in inputs}:
        raise ValueError('OUTPUT_OVERWRITES_SOURCE')
    frame,receipt=load_existing_inputs(native,decoded_receipt,seven,nine,
                                       source_receipt,nine_manifest,
                                       source_receipt_sha256=source_receipt_sha256,expected_n=expected_n)
    loc=pd.read_csv(native,usecols=['obs_id','longitude'],dtype={'obs_id':str},low_memory=False)
    if len(loc)!=len(frame) or not loc.obs_id.is_unique:raise ValueError('ORIGINAL_S1_LONGITUDE_ID_MISMATCH')
    ids=frame.obs_id.copy()
    frame=frame.merge(loc,on='obs_id',validate='one_to_one',how='left',sort=False)
    if not frame.obs_id.equals(ids):raise ValueError('ORIGINAL_SOURCE_JOIN_ORDER_CHANGED')
    result=run_spls2(frame,bootstrap=bootstrap)
    result['status']='V4_SOURCE_VERIFIED_EXPLORATORY_SPLS2'
    result['authenticated_source_receipt_sha256']=source_receipt_sha256
    result['original_source_S1_sha256']=receipt['input_digests']['native']
    result['original_source_observation_n']=int(len(frame))
    forbidden={'obs_id','accepted_key','taxon','taxon_name','user_id',
               'photographer_id','photo_id','latitude','longitude','individual_prediction'}
    def privacy(x):
        if isinstance(x,dict):
            for key,val in x.items():
                if key in forbidden:raise ValueError('PUBLIC_SPLS2_EXPOSES_SOURCE_KEY_'+key)
                privacy(val)
        elif isinstance(x,list):
            for v in x:privacy(v)
    privacy(result)
    out=Path(output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for p in ('native','decoded-receipt','seven','nine','nine-manifest','source-receipt','output'):
        ap.add_argument('--'+p,type=Path,required=True)
    ap.add_argument('--source-receipt-sha256',required=True)
    ap.add_argument('--bootstrap',type=int,default=100)
    a=ap.parse_args()
    result=run_existing_spls2(native=a.native,decoded_receipt=a.decoded_receipt,
         seven=a.seven,nine=a.nine,nine_manifest=a.nine_manifest,
         source_receipt=a.source_receipt,source_receipt_sha256=a.source_receipt_sha256,
         output=a.output,bootstrap=a.bootstrap)
    print(json.dumps({'status':result['status'],'taxa':result['species_support']['n_joint_supported_taxa'],
          'joint_skill':result['taxon_holdout_nested_cv']['skill_all_eight'],
          'sector_skill':result['three_sector_holdout_nested_cv'].get('skill_all_eight')}))

if __name__=='__main__':main()
