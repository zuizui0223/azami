"""Cross-fitted, source-verified V4 S1 climate sPLS2 axes versus geography and elevation.

One original native-S1 taxon per row, all eight photographic constructs in a
jointly QC-eligible original-photo panel. This is NOT a confirmation of the
supervised sPLS2 axes, a causal effect estimate, a phylogenetic correction, or
an alteration of the frozen GEB S1/S2/H3c analyses.

All original identifiers and locations remain confined to the authorised
source environment; only nonidentifying model/coverage summaries are written.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.model_selection import RepeatedKFold

from analysis.v4.species_environment_spls2 import (
    CLIMATE, YCOLS, OUTCOMES, MIN_PHOTOS, MIN_TAXA, SEED, SparsePLS2,
    aggregate_species_joint, tune,
)
from analysis.v4.species_environment_partial_rda import (
    GEO_COLS, aggregate_joint_geo,
)
from analysis.v4.species_environment_existing_inputs import load_existing_inputs

ELEVATION = 'topo_elevation'
KEEP_GRID = (2,4,7,11,16,22)
RIDGE_ALPHA = 5.0
MODELS = {
    'climate_axes': ('climate',),
    'geography_only': ('geography',),
    'climate_axes_plus_geography': ('climate','geography'),
    'geography_plus_elevation': ('geography','elevation'),
    'climate_axes_plus_geography_plus_elevation': ('climate','geography','elevation'),
}
FORBIDDEN = {'obs_id','accepted_key','taxon','taxon_name','user_id',
             'photographer_id','photo_id','latitude','longitude',
             'individual_prediction','_taxon_private'}


def aggregate_joint_elevation(d:pd.DataFrame)->tuple[pd.DataFrame,dict]:
    """Use SAME all-eight-QC photographs as existing joint sPLS2 and geo panel.

    Elevation mean is measured over finite pixels at those photos, but all
    outcome/climate values retain their original same-trait-QC sample.
    Require >=5 authentic elevation pixels per supported taxon, then use the
    identical retained taxa in every predictive comparison.
    """
    if not {'obs_id','accepted_key',ELEVATION,*YCOLS,
            *(x+'__eligible' for x in YCOLS)}.issubset(d.columns):
        raise ValueError('S1_ELEVATION_OR_JOINT_PHOTO_FIELDS_MISSING')
    geo,receipt=aggregate_joint_geo(d)
    y=d[list(YCOLS)].apply(pd.to_numeric,errors='coerce').to_numpy(float)
    mask=np.isfinite(y).all(axis=1)
    for col in YCOLS:
        mask &= d[col+'__eligible'].eq(True).fillna(False).to_numpy(bool)
    hue=list(OUTCOMES['hue'])
    mask &= np.linalg.norm(y[:,hue],axis=1)>=.20
    if int(mask.sum())!=receipt['all_eight_joint_photo_eligible']:
        raise ValueError('PHOTO_QC_COHORT_DIVERGED_FROM_JOINT_SPLS2')
    used=d.loc[mask,['accepted_key',ELEVATION]].copy()
    used[ELEVATION]=pd.to_numeric(used[ELEVATION],errors='coerce').replace([np.inf,-np.inf],np.nan)
    t=used.groupby('accepted_key',sort=True,observed=True)[ELEVATION].agg(['mean','count'])
    t.index=t.index.astype(str)
    merged=geo.merge(t,left_on='_taxon_private',right_index=True,how='left',validate='one_to_one',sort=False)
    if len(merged)!=receipt['n_joint_supported_taxa']:
        raise ValueError('ELEVATION_JOIN_MAY_HAVE_LOST_JOINT_TAXA')
    matched=merged.loc[(merged['count']>=MIN_PHOTOS)&np.isfinite(merged['mean'])].copy()
    matched[ELEVATION]=matched['mean'].astype(float)
    matched=matched.drop(columns=['mean','count'])
    if len(matched)<MIN_TAXA:
        raise ValueError('S1_TOO_FEW_EXACT_CLIMATE_GEO_ELEV_MATCHED_TAXA')
    if not np.isfinite(matched[[*GEO_COLS,ELEVATION]].to_numpy(float)).all():
        raise ValueError('S1_GEO_ELEV_NONFINITE')
    if matched['_taxon_private'].duplicated().any():
        raise ValueError('DUPLICATE_SOURCE_TAXA')
    receipt['n_source_species_before_elevation_gate']=int(len(merged))
    receipt['n_source_species_matching_elevation_gate']=int(len(matched))
    receipt['matched_joint_source_photo_contributions']=int(matched.n_photos.sum())
    receipt['elevation_source_photo_QC_is_same_as_all_eight_joint']=True
    receipt['min_elevation_finite_photo_pixels_per_taxon']=MIN_PHOTOS
    return matched.sort_values('_taxon_private',kind='stable').reset_index(drop=True),receipt


def project_scores(model:SparsePLS2,X:np.ndarray)->np.ndarray:
    """Projection W(P'W)^(-1) accounts for PLS2 X deflation correctly."""
    a=np.asarray(X,float)
    if a.ndim!=2 or a.shape[1]!=len(model.med) or np.isinf(a).any():
        raise ValueError('PLS2_PROJECTION_COLUMNS_INVALID')
    imp=np.where(np.isfinite(a),a,model.med[None,:])
    z=(imp-model.xmean)/model.xstd
    proj=model.weights@np.linalg.pinv(model.xloadings.T@model.weights,rcond=1e-10)
    projected=z@proj
    if projected.ndim!=2 or projected.shape[1]!=model.weights.shape[1] or not np.isfinite(projected).all():
        raise ValueError('PLS2_DEFLATION_PROJECTION_FAILED')
    return projected


def _standardized_y(Ytr,Yte):
    avg=np.mean(Ytr,axis=0);sd=np.std(Ytr,axis=0,ddof=1)
    idx=list(OUTCOMES['hue']);sd[idx]=float(np.sqrt(np.sum(sd[idx]**2)))
    if np.any(sd<1e-8) or not np.isfinite(sd).all():
        raise ValueError('CONSTANT_TRAIN_RESPONSE')
    return (Ytr-avg)/sd,(Yte-avg)/sd


def _scale(train,test):
    train=np.asarray(train,float);test=np.asarray(test,float)
    if train.ndim!=2 or test.ndim!=2 or train.shape[1]!=test.shape[1] or not np.isfinite(train).all() or not np.isfinite(test).all():
        raise ValueError('GEO_ELEV_OR_SCORE_NONFINITE')
    avg=train.mean(0);std=train.std(0,ddof=1)
    # Near-invariant spherical polynomial columns should not generate spikes.
    std=np.where(std>1e-8,std,1.)
    return (train-avg)/std,(test-avg)/std


def _make_design(tr,te,which):
    parts=[];test_parts=[]
    for key in which:
        a,b=tr[key],te[key]
        x,y=_scale(a,b)
        parts.append(x);test_parts.append(y)
    return np.column_stack(parts),np.column_stack(test_parts)


def _evaluate_fold(X,Y,G,E,tr,te,*,seed:int)->dict:
    # Training-only 4-fold tuning of sparsity; 2 sPLS2 latent axes fixed by the
    # previously published exploratory 2-component structure. Test taxa never
    # influence component weights, missing-value imputation or ridge fitting.
    hyper,_=tune(X[tr],Y[tr],seed=seed,k_grid=KEEP_GRID,h_grid=(2,))
    pls=SparsePLS2(hyper['keepX'],2).fit(X[tr],Y[tr])
    cltr=project_scores(pls,X[tr]);clte=project_scores(pls,X[te])
    train={'climate':cltr,'geography':G[tr],'elevation':E[tr,None]}
    test={'climate':clte,'geography':G[te],'elevation':E[te,None]}
    ytr,yte=_standardized_y(Y[tr],Y[te])
    baseline=np.square(yte - ytr.mean(0)).sum(axis=0)
    result={'n_train_taxa':int(len(tr)),'n_test_taxa':int(len(te)),
            'sparse_pls_selected_predictors':int(hyper['keepX']),
            'sparse_pls_components':2,'baseline_loss_by_response':baseline.tolist(),
            'models':{}}
    for name,terms in MODELS.items():
        sx,tx=_make_design(train,test,terms)
        reg=Ridge(alpha=RIDGE_ALPHA,fit_intercept=True,solver='svd')
        reg.fit(sx,ytr)
        pred=reg.predict(tx)
        d=np.sum(np.square(yte-pred),axis=0)
        result['models'][name]=d.tolist()
    return result


def compare_cv(frame:pd.DataFrame,*,repeats:int=2,region_holdout:bool=True,seed:int=SEED)->dict:
    if repeats<1 or repeats>5:raise ValueError('REPEATS_OUTSIDE_FIXED_EXPLORATORY_BUDGET')
    n=len(frame)
    if n<MIN_TAXA:raise ValueError('S1_MATCHED_TAXA_INSUFFICIENT')
    X=frame[list(CLIMATE)].to_numpy(float)
    Y=frame[list(YCOLS)].to_numpy(float)
    G=frame[list(GEO_COLS)].to_numpy(float)
    E=frame[ELEVATION].to_numpy(float)
    sectors=frame['sector'].to_numpy(int)
    if np.isinf(X).any() or not np.isfinite(Y).all() or not np.isfinite(G).all() or not np.isfinite(E).all():
        raise ValueError('S1_BAD_JOINT_PANEL')
    def evaluate(partitions,source_seed):
        sums={k:np.zeros(Y.shape[1]) for k in MODELS};base=np.zeros(Y.shape[1])
        choices=Counter();holdout=[]
        for i,(name,tr,te) in enumerate(partitions):
            result=_evaluate_fold(X,Y,G,E,tr,te,seed=source_seed+101*i)
            base+=np.array(result['baseline_loss_by_response'])
            for model,losses in result['models'].items():sums[model]+=np.array(losses)
            choices[result['sparse_pls_selected_predictors']]+=1
            foldscore={k:float(1-sum(result['models'][k])/sum(result['baseline_loss_by_response']))
                       for k in MODELS}
            holdout.append({'split':name,'heldout_taxa':int(len(te)),
                            'skill_against_training_mean':foldscore})
        den=float(base.sum())
        if den<=1e-10:raise ValueError('REFERENCE_LOSS_DEGENERATE')
        output={'n_partitions':len(partitions),
                'model_skill_against_training_mean':{k:float(1-sums[k].sum()/den) for k in MODELS},
                'model_skill_by_construct':{},
                'incremental_skill_relative_to_train_mean':{
                    'climate_after_geography':float((sums['geography_only'].sum()-sums['climate_axes_plus_geography'].sum())/den),
                    'elevation_after_climate_geography':float((sums['climate_axes_plus_geography'].sum()-sums['climate_axes_plus_geography_plus_elevation'].sum())/den),
                    'elevation_after_geography':float((sums['geography_only'].sum()-sums['geography_plus_elevation'].sum())/den),
                    'climate_after_geography_elevation':float((sums['geography_plus_elevation'].sum()-sums['climate_axes_plus_geography_plus_elevation'].sum())/den)},
                'selected_sparsity_frequency':{str(k):int(v) for k,v in sorted(choices.items())},
                'holdout_partition_scores':holdout,
                'repeat_partitions_are_not_independent_replications':True}
        for trait,idx in OUTCOMES.items():
            selection=np.array(idx,int);td=float(base[selection].sum())
            output['model_skill_by_construct'][trait]={k:float(1-sums[k][selection].sum()/td)
                        for k in MODELS} if td>1e-10 else {k:None for k in MODELS}
        return output
    split=[(f'taxon_fold_{i}',tr,te) for i,(tr,te) in enumerate(
        RepeatedKFold(n_splits=5,n_repeats=repeats,random_state=seed).split(X))]
    cv=evaluate(split,seed+41)
    if region_holdout:
        if sorted(set(sectors))==[0,1,2] and min(np.bincount(sectors,minlength=3))>=10:
            folds=[(f'sector_{k}',np.flatnonzero(sectors!=k),np.flatnonzero(sectors==k)) for k in (0,1,2)]
            cv['geographic_holdout']={'status':'estimated',**evaluate(folds,seed+399)}
        else:cv['geographic_holdout']={'status':'insufficient_region_support'}
    else:cv['geographic_holdout']={'status':'not_requested'}
    return cv


def run(frame:pd.DataFrame,*,repeats:int=2)->dict:
    units,receipt=aggregate_joint_elevation(frame)
    comparison=compare_cv(units,repeats=repeats)
    if (comparison['n_partitions']!=5*repeats or
            len(comparison['model_skill_by_construct'])!=8 or
            len(comparison['model_skill_against_training_mean'])!=5):
        raise ValueError('S1_AXIS_COMPARISON_OUTPUT_CONTRACT_INVALID')
    out={'schema':'azami_v4_s1_joint_spls2_axes_geo_elevation_crossfit_v1',
         'status':'source_authentication_required_outside_this_function',
         'original_input_observations':int(len(frame)),
         'species_support':receipt,
         'environment_axes':{'number_of_original_bioclim_climate_predictors':len(CLIMATE),
                             'number_of_spls2_axes':2,
                             'sparsity_tuned_only_within_training_taxa':True,
                             'climate_projection':'Xstd @ W @ pinv(P.T @ W)',
                             'response_y_influences_training_axis':True},
         'geography':{'n_spherical_basis_terms':len(GEO_COLS),
                      'basis':'8 spherical degree-1/degree-2 terms, mean of exact QC photo coordinates',
                      'not_a_phylogenetic_or_causal_correction':True},
         'elevation':{'model_feature':'WorldClim original source pixel topo_elevation',
                      'original_QC_pixel_support_minimum':MIN_PHOTOS,
                      'extra_information_not_causal_altitude_effect':True},
         'regression':{'estimator':'fixed ridge multivariate prediction',
                       'ridge_alpha':RIDGE_ALPHA,'baseline':'training taxon response mean',
                       'circular_hue':'two response coordinates, joint standardized variance one'},
         'heldout':comparison,
         'no_confirmatory_p_values':True,
         'all_models_use_exact_same_matched_source_taxa_and_fold_indices':True,
         'no_S2_or_H3c_reanalysis_claimed':True,
         'original_frozen_S1_S2_H3c_unchanged':True,
         'no_individual_taxa_observations_locations_or_predictions_exported':True}
    return out


def run_existing(native:Path,decoded_receipt:Path,seven:Path,nine:Path,
                 source_receipt:Path,nine_manifest:Path,output:Path,*,
                 source_receipt_sha256:str,repeats:int=2)->dict:
    sources=[native,decoded_receipt,seven,nine,source_receipt,nine_manifest]
    if Path(output).resolve() in set(Path(x).resolve() for x in sources):
        raise ValueError('OUTPUT_OVERWRITES_PROTECTED_INPUT')
    frame,receipt=load_existing_inputs(native,decoded_receipt,seven,nine,
                       source_receipt,nine_manifest,
                       source_receipt_sha256=source_receipt_sha256,expected_n=313350)
    # Original coordinates are read from the same decoded source and validated
    # against exact obs IDs, never fetched from a different copy or provider.
    loc=pd.read_csv(native,usecols=['obs_id','latitude','longitude'],dtype={'obs_id':str},low_memory=False)
    if len(loc)!=len(frame) or not loc.obs_id.is_unique:
        raise ValueError('ORIGINAL_S1_GEOGRAPHY_GRAIN_INVALID')
    ids=frame.obs_id.copy()
    frame=frame.merge(loc,on='obs_id',validate='one_to_one',how='left',sort=False)
    if not frame.obs_id.equals(ids):
        raise ValueError('SOURCE_GEO_JOIN_ORDER_CHANGED')
    result=run(frame,repeats=repeats)
    result['status']='V4_ORIGINAL_S1_SOURCE_VERIFIED_CROSSFITTED_SPLS2_GEO_ELEV'
    result['original_source_s1_sha256']=receipt['input_digests']['native']
    result['source_receipt_sha256']=source_receipt_sha256
    def audit(obj):
        if isinstance(obj,dict):
            for k,v in obj.items():
                if k in FORBIDDEN:raise ValueError('PUBLIC_OUTPUT_PRIVATE_KEY_'+str(k))
                audit(v)
        elif isinstance(obj,list):
            for v in obj:audit(v)
    audit(result)
    output=Path(output)
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('native','decoded-receipt','seven','nine','nine-manifest','source-receipt','output'):
        p.add_argument('--'+k,required=True,type=Path)
    p.add_argument('--source-receipt-sha256',required=True)
    p.add_argument('--repeats',type=int,default=2)
    a=p.parse_args()
    data=run_existing(a.native,a.decoded_receipt,a.seven,a.nine,
                      a.source_receipt,a.nine_manifest,a.output,
                      source_receipt_sha256=a.source_receipt_sha256,repeats=a.repeats)
    print(json.dumps({'status':data['status'],
                      'n_taxa':data['species_support']['n_source_species_matching_elevation_gate'],
                      'n_climate_axes':data['environment_axes']['number_of_spls2_axes'],
                      'n_model_comparisons':len(data['heldout']['model_skill_against_training_mean'])}))

if __name__=='__main__':main()
