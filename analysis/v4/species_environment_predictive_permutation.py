"""Exploratory V4 S1 original-source RF and sparse-PLS2 predictive permutation tests.

Data privacy boundary: source table and taxon identities stay in the protected
workspace. Only anonymous null distributions, scores, P values and calibration
metadata can be exported. RF group tests use a GROUP permutation-importance
statistic, tested against shuffled *outcome labels with RF retraining*; they
are NOT conditional independence tests for that group given correlated climate.
Both unrestricted and within-longitude-sector label shuffles are reported.
Posthoc group selection means P values are not prospective confirmation.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
from collections import OrderedDict

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold

from analysis.v4 import species_environment_rf_screen as rf
from analysis.v4 import species_environment_spls2 as spls
from analysis.v4.species_environment_existing_inputs import load_existing_inputs

SEED=20261107
CLIMATE=tuple(spls.CLIMATE)
CANDIDATES=OrderedDict([
    ('precipitation__angle', ('angle', tuple(p for p in CLIMATE if rf.FAMILIES[p]=='precipitation'))),
    ('radiation__chroma', ('chroma', ('chelsa_rsds_mean',))),
    ('BIO08__hue', ('hue', ('chelsa_bio08',))),
])
assert len(CLIMATE)==22 and len(CANDIDATES)==3
assert all(len(v[1])>=1 for v in CANDIDATES.values())


def holm(vals:dict[str,float])->dict[str,float]:
    """Monotone Holm over the explicitly selected three posthoc RF candidates."""
    pairs=sorted(vals.items(),key=lambda t:(t[1],t[0]))
    out={};current=0.;m=len(pairs)
    for j,(name,p) in enumerate(pairs):
        current=max(current,min(1.,float(p)*(m-j)))
        out[name]=float(current)
    return out


def shuffled_labels(Y:np.ndarray, sectors:np.ndarray, seed:int, *,within_sector:bool)->np.ndarray:
    """Permute WHOLE multivariate response vectors to preserve hue and trait covariance."""
    Y=np.asarray(Y,float);sect=np.asarray(sectors,int)
    if len(Y)!=len(sect) or not np.isfinite(Y).all():raise ValueError('INVALID_SHUFFLE_GRAIN')
    rng=np.random.default_rng(seed);ix=np.arange(len(Y))
    if within_sector:
        if set(np.unique(sect))!={0,1,2}:raise ValueError('REGIONAL_SUPPORT_NOT_THREE_SECTORS')
        for g in (0,1,2):
            pos=np.flatnonzero(sect==g)
            if len(pos)<10:raise ValueError('INSUFFICIENT_SECTOR_TAXA')
            ix[pos]=rng.permutation(pos)
    else:ix=rng.permutation(ix)
    return Y[ix]


def prep_rf_species(frame:pd.DataFrame)->dict[str,dict]:
    out={}
    for name,(trait,pair) in CANDIDATES.items():
        c=rf.aggregate_species(frame,trait).sort_values('taxon',kind='mergesort').reset_index(drop=True)
        if len(c)<40:raise ValueError('INSUFFICIENT_RF_TAXA')
        available=rf._coverage(c)
        predictors=[r['name'] for r in available if r['name'] in CLIMATE and r['status']=='used']
        if len(predictors)!=22 or not set(pair)<=set(predictors):
            raise ValueError('REQUIRED_V4_CLIMATE_SOURCE_NOT_COVERED')
        if '_sector' not in c or set(c._sector.to_numpy(int))!={0,1,2}:
            raise ValueError('SOURCE_GEOGRAPHY_UNAVAILABLE')
        Y=c[[f'y{i}' for i in range(len(rf.FEATURES[trait]))]].to_numpy(float)
        X=c[predictors].to_numpy(float)
        out[name]={'x':X,'y':Y,'sector':c._sector.to_numpy(int),
                  'features':predictors,'target':trait,
                  'group':pair,'taxa':len(c),'n_photo_contributions':int(c.n_source_photos.sum())}
    return out


def rf_group_stat(cohort:dict, Y:np.ndarray, *, seed:int=SEED, trees:int=80,
                  n_folds:int=5, permutations_per_fold:int=3)->dict:
    """Heldout RF group-permutation delta-loss / training mean reference loss.

    Global training/holdout outcome permutation is reapplied at every null draw.
    All imputation is train-only, fixed-seed KFold is reused. The group itself
    is shuffled JOINTLY on test rows, preserving its internal correlation but
    not conditioning on the other 22 predictors.
    """
    X=cohort['x'];Y=np.asarray(Y,float)
    if Y.shape!=cohort['y'].shape:raise ValueError('RF_SHUFFLED_Y_SHAPE_MISMATCH')
    names=cohort['features'];idx=[names.index(p) for p in cohort['group']]
    kfold=KFold(n_splits=n_folds,shuffle=True,random_state=SEED)
    total_ref=0.;total_model=0.;total_permuted=0.;details=[]
    for f,(train,test) in enumerate(kfold.split(X)):
        trainx,testx,_=rf.impute_train_only(X[train],X[test])
        forest=RandomForestRegressor(n_estimators=trees,min_samples_leaf=4,
                   max_features=.75,n_jobs=1,random_state=SEED+f)
        target=Y[train,0] if Y.shape[1]==1 else Y[train]
        forest.fit(trainx,target)
        base=np.repeat(Y[train].mean(axis=0)[None,:],len(test),axis=0)
        pred=rf._prediction(forest,testx,Y.shape[1])
        circular=cohort['target']=='hue'
        ref=rf._loss(Y[test],base,circular=circular)
        loss=rf._loss(Y[test],pred,circular=circular)
        permuted=0.
        for k in range(permutations_per_fold):
            shuffled=testx.copy()
            rng=np.random.default_rng(SEED+f*213+k*103+int(hashlib.sha256('+'.join(cohort['group']).encode()).hexdigest()[:8],16))
            order=rng.permutation(len(test))
            shuffled[:,idx]=testx[order][:,idx]
            xp=rf._prediction(forest,shuffled,Y.shape[1])
            permuted+=rf._loss(Y[test],xp,circular=circular)/permutations_per_fold
        total_ref+=ref;total_model+=loss;total_permuted+=permuted
        details.append({'heldout_taxa':int(len(test)),'group_delta_loss':float(permuted-loss),
                        'skill':float(1-loss/ref) if ref>1e-12 else None})
    if total_ref<=1e-10:raise ValueError('DEGENERATE_RF_BASELINE_LOSS')
    return {'group_importance_normalized':float((total_permuted-total_model)/total_ref),
            'pooled_skill':float(1-total_model/total_ref),
            'baseline_loss':float(total_ref),'model_loss':float(total_model),
            'group_permuted_loss':float(total_permuted),'folds':details}


def pls_stats(model:spls.SparsePLS2)->dict[str,float|list[float]]:
    """Coefficients in standardized X/Y units, invariant to PLS component sign."""
    beta=model.beta
    col=list(CLIMATE)
    get=lambda k:col.index(k)
    angle=spls.OUTCOMES['angle'][0];chroma=spls.OUTCOMES['chroma'][0]
    hue=list(spls.OUTCOMES['hue'])
    rain=np.array([beta[get(k),angle] for k in CANDIDATES['precipitation__angle'][1]],float)
    return {'precipitation__angle':float(np.linalg.norm(rain)),
            'precipitation__angle_sum_signed':float(np.sum(rain)),
            'radiation__chroma':float(abs(beta[get('chelsa_rsds_mean'),chroma])),
            'radiation__chroma_signed':float(beta[get('chelsa_rsds_mean'),chroma]),
            'BIO08__hue':float(np.linalg.norm(beta[get('chelsa_bio08'),hue])),
            'BIO08__hue_vector':beta[get('chelsa_bio08'),hue].astype(float).tolist()}


def pls_cv_and_coeffs(X:np.ndarray,Y:np.ndarray, *, seed:int=SEED,
                      nested:bool=True, k_grid=None,h_grid=None)->dict:
    if k_grid is None:k_grid=spls.K_GRID
    if h_grid is None:h_grid=spls.H_GRID
    opt,_=spls.tune(X,Y,seed+999,k_grid=k_grid,h_grid=h_grid)
    fit=spls.SparsePLS2(opt['keepX'],opt['n_components']).fit(X,Y)
    result={'n_components':opt['n_components'],'keepX':opt['keepX'],
            'coefficient_stats':pls_stats(fit)}
    if nested:
        cv=spls.outer_cv(X,Y,seed=seed,k_grid=k_grid,h_grid=h_grid)
        result['nested_cv_skill']=cv['skill_all_eight']
        result['trait_skills']=cv['skill_by_construct']
    return result


def observed_and_null(frame:pd.DataFrame,*,unrestricted:int=99,regional:int=99,
                      trees:int=80,pls_cv:bool=True,
                      seed:int=SEED)->dict:
    if unrestricted<19 or regional<19 or trees<10:raise ValueError('INVALID_EXPLORATORY_NULL_BUDGET')
    prepared=prep_rf_species(frame)
    m,cov=spls.aggregate_species_joint(frame)
    X=m[list(CLIMATE)].to_numpy(float)
    Y=m[list(spls.YCOLS)].to_numpy(float)
    sectors=m.sector.to_numpy(int)
    if len(m)<40 or not np.isfinite(Y).all():raise ValueError('SPLS2_COHORT_UNAVAILABLE')
    if not np.isfinite(X).any(axis=0).all():raise ValueError('SPLS2_CLIMATE_COVERAGE_INVALID')
    original_rf={name:rf_group_stat(cohort,cohort['y'],trees=trees) for name,cohort in prepared.items()}
    original_pls=pls_cv_and_coeffs(X,Y,nested=pls_cv)
    # Unrestricted null addresses ANY climate–trait relationship; within-sector
    # shuffles retain longitude composition and address association BEYOND it.
    null={}
    for mode,n in [('unrestricted_taxon_labels',unrestricted),('labels_within_3_longitude_sectors',regional)]:
        result={'RF':{n:[] for n in CANDIDATES},
                'sPLS2_nested_cv_skill':[],
                'sPLS2_coefficient_magnitudes':{n:[] for n in CANDIDATES},
                'null_count':n}
        for k in range(n):
            nonce=int(hashlib.sha256(f'{seed}-{mode}-{k}'.encode()).hexdigest()[:12],16)
            local_seed=nonce % (2**32-1)
            isregional=(mode=='labels_within_3_longitude_sectors')
            for name,cohort in prepared.items():
                Yp=shuffled_labels(cohort['y'],cohort['sector'],local_seed,within_sector=isregional)
                score=rf_group_stat(cohort,Yp,trees=trees)
                result['RF'][name].append(score['group_importance_normalized'])
            yp=shuffled_labels(Y,sectors,local_seed,within_sector=isregional)
            stat=pls_cv_and_coeffs(X,yp,nested=pls_cv)
            if pls_cv:result['sPLS2_nested_cv_skill'].append(stat['nested_cv_skill'])
            for name in CANDIDATES:
                result['sPLS2_coefficient_magnitudes'][name].append(stat['coefficient_stats'][name])
            if (k+1)%25==0:print(json.dumps({'permutation':mode,'done':k+1,'total':n}),flush=True)
        null[mode]=result
    # All empirical P use +1 convention. Tie >=, one-sided importance / skill.
    tests={}
    for mode,res in null.items():
        raw_rf={}
        for name in CANDIDATES:
            obs=original_rf[name]['group_importance_normalized'];vals=res['RF'][name]
            raw_rf[name]=(1+sum(v>=obs-1e-12 for v in vals))/(1+len(vals))
        raw_pls={}
        for name in CANDIDATES:
            obs=original_pls['coefficient_stats'][name];vals=res['sPLS2_coefficient_magnitudes'][name]
            raw_pls[name]=(1+sum(v>=obs-1e-12 for v in vals))/(1+len(vals))
        target=original_pls.get('nested_cv_skill')
        arr=res['sPLS2_nested_cv_skill']
        tests[mode]={'RF_group_importance':{
                name:{'observed_importance':original_rf[name]['group_importance_normalized'],
                      'null_p_one_sided':float(raw_rf[name]),
                      'posthoc_family3_holm':holm(raw_rf)[name],
                      'null_q95':float(np.quantile(res['RF'][name],.95))} for name in CANDIDATES},
               'sPLS2_selected_coefficient':{
                   name:{'observed_magnitude':float(original_pls['coefficient_stats'][name]),
                         'null_p_one_sided':float(raw_pls[name]),
                         'posthoc_family3_holm':holm(raw_pls)[name],
                         'null_q95':float(np.quantile(res['sPLS2_coefficient_magnitudes'][name],.95))} for name in CANDIDATES},
               'sPLS2_global_nested_cv':{
                   'observed_skill':float(target) if target is not None else None,
                   'null_p_one_sided':float((1+sum(x>=target-1e-12 for x in arr))/(1+len(arr))) if arr else None,
                   'null_q95':float(np.quantile(arr,.95)) if arr else None}}
    return {'schema':'V4_S1_original_photo_RF_PLS2_posthoc_permutation_v1',
            'status':'source_real_exploratory_permutation_calibration',
            'data':{'original_source_observation_n':int(len(frame)),
                'sPLS2_joint_taxa':int(len(m)),'RF_trait_taxa':{key:int(d['taxa']) for key,d in prepared.items()},
                'X_variables':list(CLIMATE),'response_constructs':list(spls.TRAITS)},
            'hypotheses':{name:{'construct':trait,'variables':list(group)} for name,(trait,group) in CANDIDATES.items()},
            'methods':{'seed':seed,'RF_trees':trees,'RF_folds':5,'RF_fold_repeats':1,
                       'RF_heldout_group_permutations_per_fold':3,'RF_leaf_min':4,
                       'sPLS2_hyperparam_selection':'nested 5fold*2 outer / fourfold inner and all full-fit coefficient refits',
                       'permutation_unrestricted_n':unrestricted,'permutation_within_sector_n':regional,
                       'P_values_ties_included_and_plus_one':True,
                       'correlated_group_importance_is_NOT_unique_conditional_effect':True,
                       'three_candidates_selected_after_inspecting_original_data':True,
                       'Holm3_corrects_only_three_tested_candidates_NOT_prior_search':True,
                       'sector_null_preserves_region_level_composition_not_full_geographic_or_phylogenetic_structure':True,
                       'legacy_frozen_Holm_results_unchanged':True,
                       'no_individual_taxon_photo_coordinate_observer_prediction_export':True},
            'observed_RF':{key:{'group_importance_normalized':score['group_importance_normalized'],
                                'pooled_skill':score['pooled_skill']} for key,score in original_rf.items()},
            'observed_sPLS2':original_pls,
            'permutation_inference':tests,
            'null_distributions':null}


def run_existing(native:Path,decoded_receipt:Path,seven:Path,nine:Path,
                 source_receipt:Path,nine_manifest:Path,output:Path,*,
                 source_receipt_sha256:str,unrestricted:int=99,regional:int=99,trees:int=80)->dict:
    inputs=[native,decoded_receipt,seven,nine,source_receipt,nine_manifest]
    if Path(output).resolve() in {Path(v).resolve() for v in inputs}:
        raise ValueError('OUTPUT_OVERWRITES_SOURCE')
    d,receipt=load_existing_inputs(native,decoded_receipt,seven,nine,source_receipt,nine_manifest,
        source_receipt_sha256=source_receipt_sha256,expected_n=313350)
    lon=pd.read_csv(native,usecols=['obs_id','longitude'],dtype={'obs_id':str},low_memory=False)
    if len(lon)!=len(d) or not lon.obs_id.is_unique or set(lon.obs_id)!=set(d.obs_id):
        raise ValueError('SOURCE_LONGITUDE_IDENTITY_INVALID')
    ids=d.obs_id.copy();d=d.merge(lon,on='obs_id',validate='one_to_one',how='left',sort=False)
    if not d.obs_id.equals(ids):raise ValueError('SOURCE_LONGITUDE_ORDER_DRIFT')
    result=observed_and_null(d,unrestricted=unrestricted,regional=regional,trees=trees)
    result.update(original_source_receipt_sha256=source_receipt_sha256,
                  source_receipt_validation_status='already_authenticated_SHA_original_S1_313350_rows',
                  original_S1_source_row_count=313350,
                  secure_local_only_no_data_or_id_export=True)
    forbidden={'obs_id','accepted_key','taxon','taxon_name','latitude','longitude','photographer_id',
               'user_id','photo_id','individual_prediction','_taxon_private'}
    def audit(node):
        if isinstance(node,dict):
            for k,v in node.items():
                if k in forbidden:raise ValueError('PRIVATE_KEY_EXPORT_REFUSED_'+k)
                audit(v)
        elif isinstance(node,list):
            for v in node:audit(v)
    audit(result)
    dest=Path(output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('native','decoded-receipt','seven','nine','nine-manifest','source-receipt','output'):
        parser.add_argument('--'+key,type=Path,required=True)
    parser.add_argument('--source-receipt-sha256',required=True)
    parser.add_argument('--unrestricted',type=int,default=99)
    parser.add_argument('--regional',type=int,default=99)
    parser.add_argument('--trees',type=int,default=80)
    a=parser.parse_args()
    r=run_existing(native=a.native,decoded_receipt=a.decoded_receipt,
        seven=a.seven,nine=a.nine,source_receipt=a.source_receipt,
        nine_manifest=a.nine_manifest,output=a.output,
        source_receipt_sha256=a.source_receipt_sha256,
        unrestricted=a.unrestricted,regional=a.regional,trees=a.trees)
    print(json.dumps({'status':r['status'],
                      'RF_p_unrestricted':{k:v['null_p_one_sided'] for k,v in r['permutation_inference']['unrestricted_taxon_labels']['RF_group_importance'].items()},
                      'sPLS2_p_unrestricted':r['permutation_inference']['unrestricted_taxon_labels']['sPLS2_global_nested_cv']['null_p_one_sided']}),flush=True)

if __name__=='__main__':main()
