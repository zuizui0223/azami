"""Exploratory RF environmental-variable ranking for Cirsium species means.

Different estimand from the 26-taxon within-species geographic RF: a single
mean phenotype per eligible taxon, without an eight-cell/40-photo support gate.
All computations are exploratory. Held-out permutation importance means reliance
of this RF on a source environment axis; it is NOT an independent causal effect.
Only anonymous scores/rank frequencies and source provenance may be exported.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import RepeatedKFold, LeaveOneGroupOut

from analysis.v4.species_environment_scan import FEATURES, FAMILIES, PREDICTORS
from analysis.v4.species_environment_existing_inputs import load_existing_inputs

TRAITS = tuple(FEATURES)
MIN_OBS = 5
MIN_TAXA = 30
MIN_PREDICTOR_COVERAGE = .70
CLUSTER_ABS_SPEARMAN = .80
SEED = 20261107


def _as_float(s: pd.Series) -> np.ndarray:
    return pd.to_numeric(s, errors='coerce').to_numpy(float)


def aggregate_species(frame: pd.DataFrame, construct: str, *, min_obs: int = MIN_OBS) -> pd.DataFrame:
    """Private one-row/taxon summary; response cohort identical for all X.

    Environment averages use only source observations eligible for the response,
    separately for each X's finite rows; missing X does not discard a taxon or
    change its response mean. No geographic-cell support or 40-photo threshold.
    """
    if construct not in FEATURES or min_obs < 2:
        raise ValueError('UNREGISTERED_RF_CONSTRUCT_OR_SUPPORT')
    if 'obs_id' in frame:
        # Stable floating-point group aggregation and deterministic taxon panels.
        if frame.obs_id.isna().any() or not frame.obs_id.is_unique:
            raise ValueError('RF_SOURCE_OBSERVATION_IDS_INVALID')
        frame=frame.sort_values('obs_id',key=lambda s:s.astype(str),kind='mergesort')
    ycols = list(FEATURES[construct])
    missing = {'accepted_key', *ycols, *(s+'__eligible' for s in ycols)} - set(frame)
    cols=['taxon','n_source_photos',*[f'y{i}' for i in range(len(ycols))],*PREDICTORS]
    if missing:
        d=pd.DataFrame(columns=cols)
        d.attrs['source_status']='source_trait_columns_unavailable'
        return d
    y = frame[ycols].apply(pd.to_numeric,errors='coerce').to_numpy(float)
    valid = (np.isfinite(y).all(axis=1) & frame.accepted_key.notna().to_numpy(bool))
    for c in ycols:
        valid &= frame[c+'__eligible'].eq(True).fillna(False).to_numpy(bool)
    if construct=='hue':
        norms=np.linalg.norm(y,axis=1)
        valid &= norms >= .20
        y = np.divide(y,norms[:,None],out=np.full_like(y,np.nan),where=norms[:,None]>1e-12)
    good=pd.DataFrame({'taxon':frame.loc[valid,'accepted_key'].astype(str).to_numpy()})
    for i in range(len(ycols)):
        good[f'y{i}'] = y[valid,i]
    for p in PREDICTORS:
        good[p] = _as_float(frame.loc[valid,p]) if p in frame else np.nan
    if 'longitude' in frame:
        lon = _as_float(frame.loc[valid,'longitude'])
        good['_sector'] = np.where(~np.isfinite(lon)| (abs(lon)>180),-1,
                                   np.where(lon < -30,0,np.where(lon < 60,1,2)))
    grouped=good.groupby('taxon',sort=True,observed=True)
    out = grouped[[*[f'y{i}' for i in range(len(ycols))],*PREDICTORS]].mean()
    out['n_source_photos']=grouped.size()
    out=out.loc[out.n_source_photos>=min_obs].copy()
    if '_sector' in good:
        # Most-observed sector, a sampling-composition proxy, not ancestry.
        mode=grouped['_sector'].agg(lambda s:int(s[s>=0].mode().iloc[0]) if (s>=0).any() else -1)
        out['_sector']=mode.reindex(out.index).astype(int)
    if construct=='hue' and len(out):
        z=out[['y0','y1']].to_numpy(float)
        norms=np.linalg.norm(z,axis=1)
        out=out.loc[norms>=.20].copy()
        z=out[['y0','y1']].to_numpy(float)
        norm=np.linalg.norm(z,axis=1)
        out[['y0','y1']]=z/norm[:,None]
    out=out.reset_index()
    out.attrs['source_status']='qualified_species_means'
    return out


def impute_train_only(x_train:np.ndarray,x_test:np.ndarray)->tuple[np.ndarray,np.ndarray,np.ndarray]:
    """Impute environmental features with TRAINING taxon medians only."""
    train=np.asarray(x_train,float);test=np.asarray(x_test,float)
    if train.ndim!=2 or test.ndim!=2 or train.shape[1]!=test.shape[1]:
        raise ValueError('INCOMPATIBLE_RF_FEATURE_DIMENSIONS')
    if not np.isfinite(train).any(axis=0).all():
        raise ValueError('ALL_MISSING_TRAINING_FEATURE')
    median=np.nanmedian(train,axis=0)
    if not np.isfinite(median).all():raise ValueError('NONFINITE_TRAINING_IMPUTER')
    return (np.where(np.isfinite(train),train,median[None,:]),
            np.where(np.isfinite(test),test,median[None,:]),median)


def _unit_hue(y:np.ndarray)->np.ndarray:
    y=np.asarray(y,float)
    norm=np.linalg.norm(y,axis=1,keepdims=True)
    return np.divide(y,norm,out=np.zeros_like(y),where=norm>1e-10)


def _loss(y:np.ndarray,pred:np.ndarray,*,circular:bool)->float:
    a=np.asarray(y,float);b=np.asarray(pred,float)
    if circular:
        a=_unit_hue(a);b=_unit_hue(b)
    return float(np.sum((a-b)**2))


def _corr_groups(x:np.ndarray,names:list[str],threshold:float=CLUSTER_ABS_SPEARMAN)->list[list[int]]:
    """Outcome-blind correlated-variable complete-linkage groups in train only."""
    if x.ndim!=2 or x.shape[1]!=len(names):raise ValueError('CLUSTER_DIMENSION_MISMATCH')
    if len(names)==1:return [[0]]
    c=np.eye(len(names));ranks=np.apply_along_axis(stats.rankdata,0,x)
    for i in range(len(names)):
        for j in range(i+1,len(names)):
            if np.std(ranks[:,i])<1e-12 or np.std(ranks[:,j])<1e-12:
                r=0.
            else:
                r=float(np.corrcoef(ranks[:,i],ranks[:,j])[0,1])
            c[i,j]=c[j,i]=abs(r) if np.isfinite(r) else 0.
    d=np.clip(1-c,0,1)
    labels=fcluster(linkage(squareform(d,checks=True),method='complete'),
                    t=1-threshold+1e-12,criterion='distance')
    return [list(np.flatnonzero(labels==k)) for k in np.unique(labels)]


def _prediction(est:RandomForestRegressor,x:np.ndarray,n_dim:int)->np.ndarray:
    v=np.asarray(est.predict(x),float)
    return v.reshape(-1,n_dim)


def _percentile(x:list[float])->dict:
    if not x:return {'mean_skill_loss':None,'q10':None,'q50':None,'q90':None,
                     'fraction_positive':None,'n_folds':0}
    v=np.asarray(x,float)
    return {'mean_skill_loss':float(v.mean()),'q10':float(np.quantile(v,.1)),
            'q50':float(np.median(v)),'q90':float(np.quantile(v,.9)),
            'fraction_positive':float((v>0).mean()),'n_folds':int(len(v))}


def _coverage(cells:pd.DataFrame)->list[dict]:
    n=len(cells);out=[]
    for p in PREDICTORS:
        z=_as_float(cells[p]) if p in cells else np.full(n,np.nan)
        count=int(np.isfinite(z).sum())
        valid=z[np.isfinite(z)]
        if p not in cells or count==0:status='unavailable_source_column'
        elif count<max(20,int(np.ceil(n*MIN_PREDICTOR_COVERAGE))):status='insufficient_species_coverage'
        elif float(valid.std())<1e-12:status='constant_predictor'
        else:status='used'
        out.append({'name':p,'environment_family':FAMILIES[p],
                    'n_taxa_finite':count,'coverage_fraction':float(count/n) if n else 0.,
                    'status':status})
    return out


def _folds(cells:pd.DataFrame,*,seed:int,repeats:int,n_splits:int,geographic:bool):
    if not geographic:
        k=min(n_splits,len(cells))
        if k<3:return []
        return [('taxon_holdout',train,test) for train,test in RepeatedKFold(
            n_splits=k,n_repeats=repeats,random_state=seed).split(np.arange(len(cells)))]
    if '_sector' not in cells or (cells['_sector']<0).any():return []
    groups=cells['_sector'].to_numpy(int)
    if set(np.unique(groups))!={0,1,2} or (pd.Series(groups).value_counts()<10).any():return []
    return [(int(groups[test[0]]),train,test) for train,test in LeaveOneGroupOut().split(
        np.arange(len(cells)),groups=groups)]


def _fold_evaluate(cells:pd.DataFrame,trait:str,names:list[str],folds:list,
                   *,seed:int,trees:int,permutation_repeats:int,
                   do_importances:bool=True)->tuple[dict,dict,dict,dict]:
    circular=trait=='hue'; ycols=[f'y{i}' for i in range(len(FEATURES[trait]))]
    x=cells[names].to_numpy(float);y=cells[ycols].to_numpy(float)
    var_sums=defaultdict(list);fam_sums=defaultdict(list);cluster_sums=defaultdict(list)
    top5=defaultdict(int);details=[];total_base=0.;total_rf=0.
    for fold,(label,train,test) in enumerate(folds):
        xt,xv,_=impute_train_only(x[train],x[test])
        if not np.isfinite(xt).all() or not np.isfinite(xv).all():raise ValueError('IMPUTATION_NONFINITE')
        forest=RandomForestRegressor(n_estimators=trees,min_samples_leaf=4,
            max_features=.75,n_jobs=1,random_state=seed+fold)
        forest.fit(xt,y[train,0] if y.shape[1]==1 else y[train])
        pred=_prediction(forest,xv,y.shape[1]);base=np.repeat(y[train].mean(axis=0)[None,:],len(test),axis=0)
        ref_loss=_loss(y[test],base,circular=circular);model_loss=_loss(y[test],pred,circular=circular)
        total_base+=ref_loss;total_rf+=model_loss
        details.append({'heldout_region_or_fold':label,'heldout_taxa':len(test),
                        'skill_vs_training_mean':float(1-model_loss/ref_loss) if ref_loss>1e-10 else None})
        if not do_importances or ref_loss<1e-10:continue
        group_types={'family':{},'cluster':{}}
        for family in sorted(set(FAMILIES[k] for k in names)):
            group_types['family'][family]=[i for i,p in enumerate(names) if FAMILIES[p]==family]
        for indices in _corr_groups(xt,names):
            signature='+'.join(sorted(names[i] for i in indices))
            group_types['cluster'][signature]=indices
        group_types['variable']={name:[i] for i,name in enumerate(names)}
        imp_this={}
        for kind,groups in group_types.items():
            for name,indices in groups.items():
                rng=np.random.default_rng(seed+fold*1000+int.from_bytes(hashlib.sha256((kind+name).encode()).digest()[:4],'big'))
                perm_l=[]
                for rep in range(permutation_repeats):
                    xp=xv.copy();ix=rng.permutation(len(test))
                    xp[:,indices]=xv[ix][:,indices]
                    score=(_loss(y[test],_prediction(forest,xp,y.shape[1]),circular=circular)-model_loss)/ref_loss
                    perm_l.append(float(score))
                d=float(np.mean(perm_l))
                if kind=='variable':var_sums[name].append(d);imp_this[name]=d
                elif kind=='family':fam_sums[name].append(d)
                else:cluster_sums[name].append(d)
        for name in sorted(imp_this,key=lambda n:(-imp_this[n],n))[:5]:top5[name]+=1
    n_important=sum(len(v) for v in var_sums.values())
    scored_folds=max([len(v) for v in var_sums.values()],default=0)
    variables=[{'name':name,'environment_family':FAMILIES[name],
                **_percentile(var_sums.get(name,[])),
                'fraction_top_five':float(top5[name]/scored_folds) if scored_folds else None}
               for name in names]
    variables.sort(key=lambda r:(-(r['mean_skill_loss'] if r['mean_skill_loss'] is not None else -1e9),r['name']))
    families=[{'family':name,**_percentile(v)} for name,v in fam_sums.items()]
    families.sort(key=lambda r:-r['mean_skill_loss'])
    clusters=[{'members':name.split('+'),**_percentile(v)} for name,v in cluster_sums.items()]
    clusters.sort(key=lambda r:-(r['mean_skill_loss'] if r['mean_skill_loss'] is not None else -1e9))
    status='estimated_exploratory' if total_base>1e-10 else 'reference_loss_degenerate'
    summ={'status':status,'n_folds':len(folds),
          'pooled_skill_vs_training_mean':float(1-total_rf/total_base) if total_base>1e-10 else None,
          'fraction_folds_skill_positive':float(np.mean([r['skill_vs_training_mean']>0 for r in details if r['skill_vs_training_mean'] is not None]))
          if any(r['skill_vs_training_mean'] is not None for r in details) else None,
          'folds':details}
    return summ,variables,families,clusters


def evaluate_trait(cells:pd.DataFrame,trait:str,*,trees:int=120,repeats:int=3,
                   permutation_repeats:int=2,n_splits:int=5,seed:int=SEED,
                   geography:bool=True)->dict:
    if trait not in TRAITS or trees<5 or repeats<1 or permutation_repeats<1 or n_splits<3:
        raise ValueError('INVALID_EXPLORATORY_RF_SETUP')
    coverage=_coverage(cells)
    base={'construct':trait,'n_taxa':int(len(cells)),
          'n_source_photo_contributions':int(cells['n_source_photos'].sum()) if 'n_source_photos' in cells else 0,
          'response_dimensions':len(FEATURES[trait]),
          'loss':'joint_unit_hue_chord_squared' if trait=='hue' else 'squared_trait_difference',
          'predictor_coverage':coverage,'registered_predictor_count':23,
          'variable_importance':[],'family_importance':[],
          'correlation_cluster_importance':[],
          'geographic_holdout':{'status':'not_requested'} if not geography else {'status':'not_evaluated_no_rf_fit'},
          'exploratory_not_causal':True,'no_confirmatory_p_values':True}
    if cells.attrs.get('source_status')=='source_trait_columns_unavailable':
        return {**base,'status':'source_trait_columns_unavailable'}
    if len(cells)<MIN_TAXA:return {**base,'status':'insufficient_species_support'}
    yy=cells[[f'y{i}' for i in range(len(FEATURES[trait]))]].to_numpy(float)
    if not np.isfinite(yy).all():return {**base,'status':'nonfinite_response'}
    if np.linalg.norm(np.std(yy,axis=0))<1e-10:return {**base,'status':'insufficient_response_variation'}
    names=[r['name'] for r in coverage if r['status']=='used']
    if not names:return {**base,'status':'no_available_environment_predictors'}
    ordered=cells.sort_values('taxon',kind='stable').reset_index(drop=True)
    folds=_folds(ordered,seed=seed,repeats=repeats,n_splits=n_splits,geographic=False)
    if not folds:return {**base,'status':'insufficient_folds'}
    scores,variable,family,clusters=_fold_evaluate(ordered,trait,names,folds,
              trees=trees,seed=seed,permutation_repeats=permutation_repeats)
    base.update(status='exploratory_fitted_synthetic_or_authorized_inputs',
      n_predictors_used=len(names),predictors_used=names,
      pooled_skill_vs_training_mean=scores['pooled_skill_vs_training_mean'],
      taxon_holdout=scores,variable_importance=variable,
      family_importance=family,correlation_cluster_importance=clusters,
      n_repeated_taxon_folds=scores['n_folds'],
      variable_rank_is_conditional_on_model_and_other_predictors=True,
      negative_overall_skill_does_not_establish_absence_of_association=True)
    if geography:
        sectors=_folds(ordered,seed=seed,repeats=repeats,n_splits=n_splits,geographic=True)
        if sectors:
            geo,_,geo_fam,_ = _fold_evaluate(ordered,trait,names,sectors,
                trees=trees,seed=seed+351,permutation_repeats=permutation_repeats,do_importances=False)
            base['geographic_holdout']={'status':'estimated_exploratory_three_sectors',
                'pooled_skill_vs_training_mean':geo['pooled_skill_vs_training_mean'],
                'sectors':[{'sector_label':row['heldout_region_or_fold'],
                            'heldout_taxa':row['heldout_taxa'],'skill':row['skill_vs_training_mean']}
                           for row in geo['folds']],
                'method':'leave_one_dominant_original_photo_longitude_sector_out',
                'geographic_label_is_sampling_composition_not_lineage':True}
        else:
            base['geographic_holdout']={'status':'unavailable_no_source_locations' if '_sector' not in ordered else 'insufficient_3sector_support'}
    return base


def run_screen(frame:pd.DataFrame,*,trees:int=120,repeats:int=3,permutation_repeats:int=2,
               seed:int=SEED,geography:bool=True)->dict:
    r=[evaluate_trait(aggregate_species(frame,trait),trait,trees=trees,repeats=repeats,
          permutation_repeats=permutation_repeats,seed=seed,geography=geography) for trait in TRAITS]
    return {'schema':'azami_species_mean_environment_rf_screen_v1',
        'status':'exploratory_random_forest_variable_screen',
        'n_original_input_observations':int(len(frame)),
        'registered_constructs':list(TRAITS),'registered_predictors':list(PREDICTORS),
        'number_of_registered_predictors':23,
        'sampling_unit':'one_QC_eligible_mean_per_taxon',
        'minimum_source_observations_per_taxon':MIN_OBS,
        'minimum_taxa_for_rf':MIN_TAXA,
        'minimum_species_level_predictor_coverage':MIN_PREDICTOR_COVERAGE,
        'minimum_geographic_cells_per_taxon':None,
        'rf_configuration':{'trees':trees,'leaf_min':4,'max_features':.75,'seed':seed,
                            'taxon_folds':5,'repeats':repeats,'permutation_repeats':permutation_repeats,
                            'corr_cluster_abs_spearman':CLUSTER_ABS_SPEARMAN},
        'importance':'heldout permutation delta loss / heldout train-mean-reference loss; correlated group permuted jointly',
        'new_species_mean_measurement_reliability_estimated':False,
        'causal_or_confirmatory_discovery':False,
        'frozen_GEB_S1_S2_H3c_untouched':True,
        'no_photo_taxon_id_location_predictions_in_output':True,
        'results':r}


def run_existing_rf(native:Path,decoded_receipt:Path,seven:Path,nine:Path,
                    source_receipt:Path,nine_manifest:Path,output:Path,*,
                    source_receipt_sha256:str,trees:int=120,repeats:int=3,
                    permutation_repeats:int=2,geography:bool=True,
                    expected_n:int=313350)->dict:
    inputs=[native,decoded_receipt,seven,nine,source_receipt,nine_manifest]
    if Path(output).resolve() in {Path(v).resolve() for v in inputs}:
        raise ValueError('OUTPUT_OVERWRITES_SOURCE')
    d,receipt=load_existing_inputs(native,decoded_receipt,seven,nine,source_receipt,nine_manifest,
        source_receipt_sha256=source_receipt_sha256,expected_n=expected_n)
    # Re-read source longitude only for a broad, optional holdout, never export rows.
    if geography:
        try:
            loc=pd.read_csv(native,usecols=['obs_id','longitude'],dtype={'obs_id':str},low_memory=False)
        except ValueError:
            loc=None
        if loc is not None:
            if len(loc)!=len(d) or not loc.obs_id.is_unique:
                raise ValueError('SOURCE_GEOGRAPHY_IDENTITY_INVALID')
            ids=d.obs_id.copy()
            d=d.merge(loc,on='obs_id',validate='one_to_one',how='left',sort=False)
            if not d.obs_id.equals(ids) or d.longitude.isna().all():
                raise ValueError('SOURCE_GEOGRAPHY_MISMATCH')
    r=run_screen(d,trees=trees,repeats=repeats,
                 permutation_repeats=permutation_repeats,geography=geography)
    r.update(source_provenance=receipt,
             scope='already-authorized local source files only; no source recovery or raster fetch',
             genuine_original_source_checked_by_existing_receipt=True)
    out=Path(output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(r,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return r


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for key in ('native','decoded-receipt','seven','nine','nine-manifest','source-receipt','output'):
        ap.add_argument('--'+key,type=Path,required=True)
    ap.add_argument('--source-receipt-sha256',required=True)
    ap.add_argument('--trees',type=int,default=120)
    ap.add_argument('--repeats',type=int,default=3)
    ap.add_argument('--permutation-repeats',type=int,default=2)
    ap.add_argument('--no-geography',action='store_true')
    a=ap.parse_args()
    r=run_existing_rf(native=a.native,decoded_receipt=a.decoded_receipt,
        seven=a.seven,nine=a.nine,nine_manifest=a.nine_manifest,
        source_receipt=a.source_receipt,source_receipt_sha256=a.source_receipt_sha256,
        output=a.output,trees=a.trees,repeats=a.repeats,
        permutation_repeats=a.permutation_repeats,geography=not a.no_geography)
    print(json.dumps({'status':r['status'],'n_constructs':len(r['results']),
        'n_fitted':sum(z['status'].startswith('exploratory_fitted') for z in r['results'])}))

if __name__=='__main__':main()
