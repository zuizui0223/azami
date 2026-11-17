"""All-eight photographic traits: observer-disjoint, spatially held-out climate prediction.

Protected original S1 only. This is a new seven-predictor exploratory comparison,
not a repeat of the 12/23/28/34-variable RF or a change to frozen GEB findings.
The estimand is a known taxon's phenotype in a held-out geographic block as
measured by observers never used in training. All IDs and predictions stay private.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import sklearn
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler

FEATURES = {
    'angle': ['orientation_image_vertical_angle'],
    'lightness': ['corolla_lab_lightness'],
    'chroma': ['corolla_lab_chroma'],
    'hue': ['corolla_hue_sin', 'corolla_hue_cos'],
    'aspect': ['capitulum_outline_aspect_ratio'],
    'circularity': ['capitulum_outline_circularity'],
    'solidity': ['capitulum_outline_solidity'],
    'width_profile': ['capitulum_width_profile_cv'],
}
PREDICTORS = ('chelsa_bio01', 'chelsa_bio04', 'chelsa_bio12', 'chelsa_bio15',
              'chelsa_rsds_mean', 'chelsa_vpd_mean', 'chelsa_sfcwind_mean')
SEEDS = (20261101, 20261103, 20261105)
KEYS = ['taxon', 'lat_cell', 'lon_cell']
EXPECTED_N, EXPECTED_USERS = 313350, 92839


def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8 << 20), b''): h.update(block)
    return h.hexdigest()


def unit(x: np.ndarray) -> np.ndarray:
    x=np.asarray(x,float);norm=np.linalg.norm(x,axis=1,keepdims=True)
    return np.divide(x,norm,out=np.zeros_like(x),where=norm>1e-12)


def side(user: str, seed: int) -> int:
    return hashlib.sha256(f'{seed}\0{user}'.encode()).digest()[-1] & 1


def distribution(x) -> dict:
    x=np.asarray(x,float);x=x[np.isfinite(x)]
    if not len(x): return {'n':0,'q05':None,'q50':None,'q95':None}
    q=np.quantile(x,[.05,.5,.95])
    return dict(n=int(len(x)),q05=float(q[0]),q50=float(q[1]),q95=float(q[2]))


def make_pairs(frame: pd.DataFrame, construct: str, degrees: int, seed: int):
    """Return PRIVATE paired cells and public counts; no model fitting here."""
    if construct not in FEATURES or degrees not in (1,5):
        raise ValueError('INVALID_ENDPOINT_OR_GRAIN')
    cols=FEATURES[construct];hue=construct=='hue'
    needed={'obs_id','user_id','accepted_key','latitude','longitude',*PREDICTORS,
            *cols,*(c+'__eligible' for c in cols)}
    if not needed <= set(frame):raise ValueError('SOURCE_COLUMNS_MISSING')
    if frame.obs_id.isna().any() or not frame.obs_id.is_unique:
        raise ValueError('DUPLICATE_OR_MISSING_OBSERVATION')
    y=frame[cols].apply(pd.to_numeric,errors='coerce').to_numpy(float)
    x=frame[list(PREDICTORS)].apply(pd.to_numeric,errors='coerce').to_numpy(float)
    lat=pd.to_numeric(frame.latitude,errors='coerce').to_numpy(float)
    lon=pd.to_numeric(frame.longitude,errors='coerce').to_numpy(float)
    valid=(np.isfinite(x).all(1)&np.isfinite(y).all(1)&np.isfinite(lat)&np.isfinite(lon)&
           (np.abs(lat)<=90)&(np.abs(lon)<=180)&frame.accepted_key.notna().to_numpy()&
           frame.user_id.notna().to_numpy())
    for col in cols:valid &= frame[col+'__eligible'].eq(True).to_numpy(bool)
    if hue:valid &= np.linalg.norm(y,axis=1)>=.2
    audit={'candidate_records':int(valid.sum()),'selected_records':0,'observer_overlap':0}
    if not valid.any():return pd.DataFrame(),audit
    d=frame.loc[valid,['obs_id','user_id','accepted_key']].copy()
    d.rename(columns={'accepted_key':'taxon'},inplace=True)
    d[['obs_id','user_id','taxon']]=d[['obs_id','user_id','taxon']].astype(str)
    d['lat_cell']=np.minimum(np.floor((lat[valid]+90)/degrees).astype(int),180//degrees-1)
    d['lon_cell']=np.minimum(np.floor((lon[valid]+180)/degrees).astype(int),360//degrees-1)
    for j,c in enumerate(cols):d[c]=(unit(y[valid]) if hue else y[valid])[:,j]
    for j,c in enumerate(PREDICTORS):d[c]=x[valid,j]
    d['_priority']=[hashlib.sha256(v.encode()).hexdigest() for v in d.obs_id]
    d=d.sort_values(['_priority','obs_id']).drop_duplicates(KEYS+['user_id'])
    audit['selected_records']=len(d)
    allocation={u:side(u,seed) for u in sorted(d.user_id.unique())}
    d['_half']=d.user_id.map(allocation).astype(int)
    audit['observer_overlap']=len(set(d.loc[d._half==0,'user_id']) & set(d.loc[d._half==1,'user_id']))
    if audit['observer_overlap']:raise ValueError('OBSERVER_LEAKAGE')
    if set(d._half.unique())!={0,1}:return pd.DataFrame(),audit
    g=d.groupby(KEYS+['_half'],observed=True,sort=True)
    med=g[cols+list(PREDICTORS)].median();med['n']=g.size()
    p=med.xs(0,level='_half').add_suffix('_a').join(med.xs(1,level='_half').add_suffix('_b'),how='inner').reset_index()
    p=p.loc[(p.n_a>=3)&(p.n_b>=3)].copy()
    if hue and len(p):
        keep=np.ones(len(p),bool)
        for s in ('a','b'):keep &= np.linalg.norm(p[[c+'_'+s for c in cols]],axis=1)>=.2
        p=p.loc[keep].copy()
        for s in ('a','b'):p.loc[:,[c+'_'+s for c in cols]]=unit(p[[c+'_'+s for c in cols]].to_numpy())
    if len(p):
        n=p.groupby('taxon').agg(cells=('n_a','size'),a=('n_a','sum'),b=('n_b','sum'))
        keep=n.index[(n.cells>=8)&((n.a+n.b)>=40)]
        p=p.loc[p.taxon.isin(keep)].copy()
    return p.sort_values(KEYS).reset_index(drop=True),audit


def fit_fold(train: pd.DataFrame, test: pd.DataFrame, cols: list[str],
             source: str, target: str, hue: bool, *, trees: int, seed: int) -> dict:
    """Private predictions; all learned quantities depend on training half only."""
    yc=[c+'_'+source for c in cols];xc=[c+'_'+source for c in PREDICTORS]
    means=train.groupby('taxon')[yc+xc].mean()
    support=train.groupby('taxon').size()
    keep=support.index[support>=3]
    if hue:keep=keep[np.linalg.norm(means.loc[keep,yc],axis=1)>1e-12]
    te=test.loc[test.taxon.isin(keep)].copy()
    if te.empty:return {}
    tr=train.loc[train.taxon.isin(keep)].copy()
    y0=means.loc[tr.taxon,yc].to_numpy();x0=means.loc[tr.taxon,xc].to_numpy()
    yt=tr[yc].to_numpy()-y0;xt=tr[xc].to_numpy()-x0
    base=means.loc[te.taxon,yc].to_numpy()
    xe=te[[c+'_'+target for c in PREDICTORS]].to_numpy()-means.loc[te.taxon,xc].to_numpy()
    obs=te[[c+'_'+target for c in cols]].to_numpy()
    counts=tr.groupby('taxon').size()
    w=1/counts.loc[tr.taxon].to_numpy(float);w=w*len(w)/w.sum()
    scaler=StandardScaler().fit(xt,sample_weight=w)
    linear=Ridge(alpha=1.0).fit(scaler.transform(xt),yt,sample_weight=w)
    rf=RandomForestRegressor(n_estimators=trees,min_samples_leaf=5,max_features=.75,
                             random_state=int(seed),n_jobs=1)
    rf.fit(xt,yt.ravel() if len(cols)==1 else yt,sample_weight=w)
    pred={'baseline':base,'ridge':base+np.asarray(linear.predict(scaler.transform(xe))).reshape(-1,len(cols)),
          'rf':base+np.asarray(rf.predict(xe)).reshape(-1,len(cols))}
    if hue:
        fallback=unit(base)
        for name,v in pred.items():
            tiny=np.linalg.norm(v,axis=1)<=1e-12
            pred[name]=unit(v);pred[name][tiny]=fallback[tiny]
    pred.update(observed=obs,groups=te.taxon.to_numpy())
    return pred


def _score(losses: pd.DataFrame, seed: int) -> dict:
    # Equal taxon weighting after collecting all held-out folds/directions.
    a=losses.groupby('_taxon')[['baseline','ridge','rf']].mean().to_numpy()
    rng=np.random.default_rng(seed);ix=rng.integers(0,len(a),size=(500,len(a)))
    boot=a[ix].mean(axis=1);mean=a.mean(axis=0);out={}
    for j,name in enumerate(('ridge','rf'),start=1):
        finite=boot[:,0]>1e-15
        bs=1-boot[finite,j]/boot[finite,0]
        out[name]={'skill':float(1-mean[j]/mean[0]) if mean[0]>1e-15 else None,
                   'taxon_bootstrap_95_conditional_on_fits':[float(v) for v in np.quantile(bs,[.025,.975])] if len(bs) else None,
                   'taxon_equal_MSE':float(mean[j]),'baseline_taxon_equal_MSE':float(mean[0]),
                   'n_taxa':len(a)}
    # Pair models within the same held-out rows, taxon weights and bootstrap draws.
    # A positive difference means RF has lower loss than Ridge, not necessarily
    # that either model improves on the training-taxon-mean baseline.
    finite = boot[:, 0] > 1e-15
    difference = (boot[finite, 1] - boot[finite, 2]) / boot[finite, 0]
    out['rf_minus_ridge'] = {
        'skill_difference': float((mean[1] - mean[2]) / mean[0]) if mean[0] > 1e-15 else None,
        'taxon_equal_MSE_reduction': float(mean[1] - mean[2]),
        'taxon_bootstrap_95_conditional_on_fits':
            [float(v) for v in np.quantile(difference, [.025, .975])] if len(difference) else None,
        'same_evaluation_rows_taxon_weights_and_bootstrap_draws': True,
    }
    out['evaluation_support'] = {
        'n_evaluated_taxa': int(len(a)),
        'at_least_five_evaluated_taxa': bool(len(a) >= 5),
        'note': 'Scored taxa can be fewer than paired-input taxa after training-support exclusions; '
                'a narrow conditional interval is not evidence of broad generality.',
    }
    return out


def evaluate(pairs: pd.DataFrame, construct: str, degrees: int, *, seed: int, trees: int=80) -> dict:
    if pairs.empty or pairs.taxon.nunique()<5:
        return {'status':'insufficient_paired_taxa','n_pairs':len(pairs)}
    p=pairs.reset_index(drop=True);cols=FEATURES[construct];hue=construct=='hue'
    blocks=(p.lat_cell*degrees//5).astype(str)+'_'+(p.lon_cell*degrees//5).astype(str)
    if blocks.nunique()<5:return {'status':'insufficient_spatial_blocks','n_pairs':len(p)}
    folds=list(GroupKFold(n_splits=5).split(p,groups=blocks))
    records=[];directions={};overlap=0;excluded=0
    for src,dst in (('a','b'),('b','a')):
        direction=[]
        for fold,(tr,te) in enumerate(folds):
            overlap+=len(set(blocks.iloc[tr])&set(blocks.iloc[te]))
            fitted=fit_fold(p.iloc[tr],p.iloc[te],cols,src,dst,hue,trees=trees,seed=seed+fold)
            if not fitted:excluded+=len(te);continue
            excluded+=len(te)-len(fitted['observed'])
            d=pd.DataFrame({'_taxon':fitted['groups']})
            for name in ('baseline','ridge','rf'):d[name]=np.sum((fitted[name]-fitted['observed'])**2,axis=1)
            direction.append(d)
        if not direction:return {'status':'insufficient_training_taxon_support','n_pairs':len(p)}
        d=pd.concat(direction,ignore_index=True);directions[src+'_to_'+dst]=_score(d,seed)
        records.append(d)
    if overlap:raise ValueError('GEOGRAPHIC_BLOCK_LEAKAGE')
    losses=pd.concat(records,ignore_index=True)
    alignment=[]
    for _,d in p.groupby('taxon'):
        a=d[[c+'_a' for c in cols]].to_numpy();b=d[[c+'_b' for c in cols]].to_numpy()
        a=a-a.mean(0);b=b-b.mean(0);den=np.sqrt(np.sum(a*a)*np.sum(b*b))
        if den>1e-12:alignment.append(np.sum(a*b)/den)
    return {'status':'estimated','n_pairs':len(p),'n_paired_taxa':int(p.taxon.nunique()),
            'n_contributing_observer_selected_records':int((p.n_a+p.n_b).sum()),
            'n_spatial_blocks':int(blocks.nunique()),'n_folds_per_direction':5,
            'n_prediction_cells_both_directions':len(losses),
            'n_test_cells_excluded_training_support':excluded,'spatial_overlap_count':overlap,
            'matched_cohort_photo_geographic_alignment':distribution(alignment),
            'pooled':_score(losses,seed),'directions':directions}


def load_source(native: Path, receipt_path: Path, roster: Path, *,
                expected_n: int=EXPECTED_N,expected_users: int=EXPECTED_USERS):
    receipt=json.loads(Path(receipt_path).read_text());digest=sha256(native)
    if (digest!=receipt.get('prepared_table_sha256') or receipt.get('native_observations')!=expected_n or
        receipt.get('environment_decoding_status')!='physical_units_gsp_invalid_encoding_masked_v1'):
        raise ValueError('SOURCE_AUTHENTICATION_FAILED')
    decode=receipt.get('environment_decoding',{})
    if any(not isinstance(decode.get(e,{}).get('source_sha256'),str) or
           len(decode[e]['source_sha256'])!=64 for e in PREDICTORS):
        raise ValueError('CLIMATE_SOURCE_RECEIPTS_MISSING')
    columns=['obs_id','accepted_key','latitude','longitude',*PREDICTORS]
    for cols in FEATURES.values():
        for c in cols:columns.extend([c,c+'__eligible'])
    d=pd.read_csv(native,usecols=columns,dtype={'obs_id':str},low_memory=False)
    u=pd.read_csv(roster,usecols=['obs_id','user_id'],dtype=str,low_memory=False)
    if (len(d)!=expected_n or len(u)!=expected_n or not d.obs_id.is_unique or not u.obs_id.is_unique or
        d.obs_id.isna().any() or u.obs_id.isna().any() or u.user_id.isna().any() or
        u.user_id.nunique()!=expected_users or set(d.obs_id)!=set(u.obs_id)):
        raise ValueError('SOURCE_OR_ROSTER_IDENTITY_DRIFT')
    d=d.merge(u,on='obs_id',how='left',sort=False,validate='one_to_one')
    return d,{'source_observations':expected_n,'source_observers':expected_users,
              'decoded_sha256':digest,'roster_sha256':sha256(roster),
              'environment_sources':{e:decode[e]['source_sha256'] for e in PREDICTORS}}


def privacy_check(obj):
    forbidden={'obs_id','user_id','accepted_key','taxon','taxon_name','latitude','longitude',
               'lat_cell','lon_cell','photo_id','predictions','model_weights','photographer_id'}
    if isinstance(obj,dict):
        if forbidden&set(obj):raise ValueError('PRIVATE_OUTPUT_FIELD')
        for value in obj.values():privacy_check(value)
    elif isinstance(obj,list):
        for value in obj:privacy_check(value)


def run(native: Path, receipt: Path, roster: Path, output: Path) -> dict:
    d,source=load_source(native,receipt,roster)
    rows=[]
    for c in FEATURES:
        for degree in (1,5):
            for seed in SEEDS:
                p,audit=make_pairs(d,c,degree,seed)
                rows.append({'construct':c,'degrees':degree,'seed':seed,'selection':audit,
                             **evaluate(p,c,degree,seed=seed)})
            print('Completed aggregate condition:',c,degree,flush=True)
    out={'schema':'azami_observer_disjoint_climate_all8_v1_1','source':source,
         'status':'source_authenticated_exploratory_prediction',
         'design':{'constructs':list(FEATURES),'climate_predictors':list(PREDICTORS),
                   'seeds':list(SEEDS),'degrees':[1,5],'holdout_block_degrees':5,
                   'min_observers_per_half_cell':3,'min_cells_per_taxon':8,'min_records_per_taxon':40,
                   'trees':80,'min_leaf':5,'max_features':.75,'ridge_alpha':1.0,
                   'global_observer_disjoint':True,'training_only_taxon_centering':True,
                   'taxon_weighting':'equal after aggregating all held-out losses',
                   'estimand':'known taxa, new 5degree blocks, independent observer group',
                   'loss':'scalar squared error; normalized joint hue chord loss'},
         'software':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'sklearn':sklearn.__version__},
         'limits':['No spatial buffer; neighboring blocks can remain close.',
                   'Bootstrap resamples taxa conditional on fitted models, not phylogenetic or refit uncertainty.',
                   'Three overlapping seeded splits are not independent experiments.',
                   'Seven native-source predictors; not a direct 23/28/34-variable comparison.',
                   'Photographic measures, not physical gravity angle, calibrated pigment or anatomy.',
                   'Unknown taxa, causal adaptation and future temporal prediction are not tested.'],
         'original_GEB_S1_S2_H3c_unchanged':True,'private_rows_and_predictions_not_exported':True,
         'results':rows}
    privacy_check(out)
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(out,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    return out


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('native','receipt','roster','output'):p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args();run(a.native,a.receipt,a.roster,a.output)
