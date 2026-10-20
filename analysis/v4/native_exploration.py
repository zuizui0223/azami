"""Native-only diagnostic of inherited upstream measurements, NOT v4 remeasurement.

The all-source native roster is determined before joining phenotype rows. These
exploratory models do not validate image phenotypes or identify causal effects.
"""
from __future__ import annotations
import argparse, hashlib, itertools, json, math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.cluster.vq import kmeans2
from analysis.v4.native_full import native_measurement_join, digest

VARIABLES=['chelsa_rsds_mean','chelsa_bio01','chelsa_bio04','chelsa_bio12','chelsa_bio15','chelsa_vpd_mean','chelsa_gsp','chelsa_sfcwind_mean','chelsa_npp']
BLOCKS={'radiation':[0],'thermal':[1,2],'water_related':[3,4,5,6],'wind':[7],'productivity':[8]}
ANGLE=['orientation_angle_degrees_median']
COLOUR=['corolla_lab_lightness_median','corolla_lab_chroma_median','corolla_hue_sin_median','corolla_hue_cos_median']
LEGACY_SHA='b49a093c849928c456d6747ec51772b5a6b78debb2ae5eb5dbe2a91b08f66dd7'


def standardize(X,Y,w=None):
    X=np.asarray(X,float);Y=np.asarray(Y,float)
    w=np.ones(len(Y)) if w is None else np.asarray(w,float);w=w/w.sum()
    def z(a):
        b=a-np.average(a,axis=0,weights=w);s=np.sqrt(np.average(b*b,axis=0,weights=w))
        return b/np.where(s>1e-12,s,1.)
    return z(X),z(Y),w


def moments(X,Y,w=None,response_weights=None):
    X,Y,w=standardize(X,Y,w)
    if response_weights is not None: Y=Y*np.sqrt(response_weights)
    return X.T@(w[:,None]*X),X.T@(w[:,None]*Y),Y.T@(w[:,None]*Y)


def r2(m,cols):
    if not cols:return 0.
    xx,xy,yy=m;cross=xy[cols]
    return float(np.trace(cross.T@np.linalg.pinv(xx[np.ix_(cols,cols)])@cross)/np.trace(yy))


def environment_summary(X,Y,w=None,response_weights=None):
    m=moments(X,Y,w,response_weights);space=list(range(9,X.shape[1]));env=list(range(9));names=list(BLOCKS)
    def alloc(base):
        scores={}
        for k in range(32):
            ids=base+[j for i,n in enumerate(names) if k&(1<<i) for j in BLOCKS[n]]
            scores[k]=r2(m,ids)
        ans={}
        for i,n in enumerate(names):
            ans[n]=sum((scores[k|(1<<i)]-scores[k])/(5*math.comb(4,k.bit_count())) for k in range(32) if not k&(1<<i))
        return ans,{n:scores[31]-scores[31^(1<<i)] for i,n in enumerate(names)}
    re,rs,rf=r2(m,env),r2(m,space),r2(m,env+space)
    sh,un=alloc([]);shs,uns=alloc(space)
    return {'r2_environment':re,'r2_space':rs,'r2_joint':rf,'partition':{'unique_environment':rf-rs,'unique_space':rf-re,'shared':re+rs-rf,'unexplained':1-rf},'shapley_environment':sh,'shapley_given_space':shs,'unique_blocks_given_space':uns,'scope':'in_sample_descriptive; correlated unequal-size blocks; no causal ranking'}


def spatial(frame):
    lat=np.deg2rad(frame.latitude.to_numpy());lon=np.deg2rad(frame.longitude.to_numpy())
    x=np.cos(lat)*np.cos(lon);y=np.cos(lat)*np.sin(lon);z=np.sin(lat)
    return np.column_stack([x,y,z,x*y,x*z,y*z,x*x-y*y,3*z*z-1])


def crossvalidate(X,Y,xyz,response_weights=None):
    _,labels=kmeans2(np.asarray(xyz),6,minit='++',seed=20261020)
    rows=[]
    for held in sorted(set(labels)):
        train=labels!=held;test=~train
        if train.sum()<25 or test.sum()<3:continue
        ym=Y[train].mean(0);ys=Y[train].std(0);ys=np.where(ys>1e-12,ys,1)
        a=(Y-ym)/ys
        if response_weights is not None:a=a*np.sqrt(response_weights)
        fold={'fold':int(held),'train_taxa':int(train.sum()),'test_taxa':int(test.sum())}
        for name,cols in [('environment',list(range(9))),('space',list(range(9,17))),('joint',list(range(17)))]:
            xx=X[:,cols];xm=xx[train].mean(0);xs=xx[train].std(0);xs=np.where(xs>1e-12,xs,1)
            xx=np.c_[np.ones(len(X)),(xx-xm)/xs]
            beta=np.linalg.lstsq(xx[train],a[train],rcond=None)[0]
            fold[name+'_sse']=float(np.sum((a[test]-xx[test]@beta)**2))
        fold['train_mean_baseline_sse']=float(np.sum(a[test]**2));rows.append(fold)
    denom=sum(r['train_mean_baseline_sse'] for r in rows)
    return {'folds':rows,'pooled_skill':{n:1-sum(r[n+'_sse'] for r in rows)/denom for n in ['environment','space','joint']},'definition':'six spherical-geography clusters of taxon native centroids; all standardization fitted on training taxa; score relative to training taxon mean'}


def headline(d,env,column,predictor,rng):
    a=d.dropna(subset=[column]);counts=a.groupby('analysis_taxon').size();a=a[a.analysis_taxon.isin(counts[counts>=5].index)]
    med=a.groupby('analysis_taxon')[column].median().to_frame().join(env.groupby('analysis_taxon')[predictor].median()).dropna()
    x=med[predictor].to_numpy();y=med[column].to_numpy();b=float(np.corrcoef(x,y)[0,1])
    draws=[float(np.corrcoef(x[ix],y[ix])[0,1]) for ix in rng.integers(0,len(x),(2000,len(x)))]
    exclusions={}
    for name,bad in [('exclude_vulgare',['Cirsium vulgare']),('exclude_arvense',['Cirsium arvense']),('exclude_both',['Cirsium vulgare','Cirsium arvense'])]:
        ids=set(d.loc[d.taxon_name.isin(bad),'analysis_taxon']);m=med.drop(index=list(ids),errors='ignore')
        exclusions[name]={'n_taxa':len(m),'beta':float(m[column].corr(m[predictor]))}
    return {'n_observations_with_measurement':len(a),'n_accepted_taxa':len(med),'beta_standardized':b,'taxon_bootstrap_95':np.quantile(draws,[.025,.975]).tolist(),'exclusions':exclusions}


def rv(X,Y):
    X,Y,_=standardize(X,Y)
    if Y.shape[1]==4:Y=Y*np.sqrt([1,1,.5,.5])
    return float(np.sum((X.T@Y)**2)/np.sqrt(np.sum((X.T@X)**2)*np.sum((Y.T@Y)**2)))


def run(directory,legacy_path,outdir):
    directory=Path(directory);outdir=Path(outdir);outdir.mkdir(exist_ok=True,parents=True)
    if digest(legacy_path)!=LEGACY_SHA:raise ValueError('Wrong legacy source')
    native=pd.read_pickle(directory/'native_observations.pkl');env=pd.read_pickle(directory/'native_environment.pkl')
    t=pd.read_csv(legacy_path,dtype={'obs_id':str},low_memory=False)
    d=native_measurement_join(native,t)
    if (d.taxon_name!=d.taxon_name_legacy).any():raise ValueError('Source taxon labels disagree between observation versions')
    d=d.merge(env[['obs_id',*VARIABLES]],on='obs_id',validate='one_to_one')
    d=d.dropna(subset=VARIABLES).copy()
    print('Native legacy measurements with complete environment',len(d),flush=True)
    sp=[f'space{i}' for i in range(8)];d[sp]=spatial(d);expo=env.dropna(subset=VARIABLES).copy();expo[sp]=spatial(expo)
    expo=expo.groupby('analysis_taxon')[[*VARIABLES,*sp]].median()
    result={'analysis_scope':'retrospective native-only diagnostic using inherited 406582 upstream measurements, NOT new pixel measurements','legacy_source_sha256':LEGACY_SHA,'native_input_sha256':digest(directory/'native_observations.csv.gz'),'environment_input_sha256':digest(directory/'native_environment.csv.gz'),'matched_native_observations':len(d),'modules':{}}
    rng=np.random.default_rng(20261020)
    result['headline_chroma_radiation']=headline(d,env,COLOUR[1],VARIABLES[0],rng)
    result['headline_angle_precipitation']=headline(d,env,ANGLE[0],'chelsa_bio12',rng)
    for name,cols,weights in [('presentation',ANGLE,None),('colour',COLOUR,[1,1,.5,.5])]:
        a=d.dropna(subset=cols).copy();cnt=a.groupby('analysis_taxon').size();a=a[a.analysis_taxon.isin(cnt[cnt>=5].index)]
        med=a.groupby('analysis_taxon')[cols].median().join(expo).dropna()
        X=med[[*VARIABLES,*sp]].to_numpy();Y=med[cols].to_numpy()
        values=a[[*VARIABLES,*sp,*cols]]-a.groupby('analysis_taxon')[[*VARIABLES,*sp,*cols]].transform('mean')
        w=1/a.groupby('analysis_taxon').obs_id.transform('size').to_numpy(float)
        result['modules'][name]={'n_observations':len(a),'n_taxa':len(med),'among':environment_summary(X,Y,response_weights=weights),'within_equal_taxon':environment_summary(values[[*VARIABLES,*sp]].to_numpy(),values[cols].to_numpy(),w,weights),'spatial_holdout_among':crossvalidate(X,Y,med[sp[:3]].to_numpy(),weights)}
        print(name,'fitted',len(a),len(med),flush=True)
    common=d.dropna(subset=ANGLE+COLOUR).copy();ct=common.groupby('analysis_taxon').size();common=common[common.analysis_taxon.isin(ct[ct>=5].index)].reset_index(drop=True)
    med=common.groupby('analysis_taxon')[ANGLE+COLOUR].median();center=common[ANGLE+COLOUR]-common.groupby('analysis_taxon')[ANGLE+COLOUR].transform('mean')
    groups=common.groupby('analysis_taxon').indices;ar=rv(med[ANGLE].to_numpy(),med[COLOUR].to_numpy());diff=[]
    for _ in range(1000):
        ix=[rng.choice(v) for v in groups.values()];s=center.iloc[ix];diff.append(ar-rv(s[ANGLE].to_numpy(),s[COLOUR].to_numpy()))
    result['presentation_colour_covariance']={'common_observations':len(common),'common_taxa':len(med),'among_RV':ar,'matched_row_count_diff_median':float(np.median(diff)),'conditional_resampling_95':np.quantile(diff,[.025,.975]).tolist(),'positive_fraction':float(np.mean(np.array(diff)>0)),'scope':'taxa and among-taxon medians fixed; not a taxon-bootstrap CI or P value; hue coordinates share one outcome weight'}
    result['new_anatomical_measurement_inference']=False
    (outdir/'native_legacy_exploration.json').write_text(json.dumps(result,indent=2)+'\n')
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--native-dir',type=Path,required=True);p.add_argument('--legacy',type=Path,required=True);p.add_argument('--out-dir',type=Path,required=True);a=p.parse_args()
    print(json.dumps(run(a.native_dir,a.legacy,a.out_dir),indent=2))
