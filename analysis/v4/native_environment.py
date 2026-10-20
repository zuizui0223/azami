"""Joint environment and space associations in the full native upstream cohort.

Exploratory multivariate RDA-style fits, grouped attribution and honest geography
holdouts. No phylogenetic correction or causal importance is asserted.
"""
import itertools
import math
import numpy as np
import pandas as pd

GROUPS={'radiation':['chelsa_rsds_mean'],'thermal':['chelsa_bio01','chelsa_bio04'],
        'water_demand':['chelsa_bio12','chelsa_bio15','chelsa_vpd_mean','chelsa_gsp'],
        'wind':['chelsa_sfcwind_mean'],'productivity':['chelsa_npp']}
VARS=sum(GROUPS.values(),[])
MODULES={'presentation':['orientation_image_vertical_angle'],
         'colour':['corolla_lab_lightness','corolla_lab_chroma','corolla_hue_sin','corolla_hue_cos']}


def weights(f):
    a=f.groupby(['taxon_name','cell']).taxon_name.transform('size').to_numpy(float)
    k=f.groupby('taxon_name').cell.transform('nunique').to_numpy(float);w=1/(a*k)
    return w/w.sum()


def normalize(a,w=None,hue=False):
    a=np.asarray(a,float)
    if a.ndim==1:a=a[:,None]
    w=np.full(len(a),1/len(a)) if w is None else np.asarray(w)/np.sum(w)
    mu=np.average(a,axis=0,weights=w);var=np.average((a-mu)**2,axis=0,weights=w)
    scale=np.sqrt(var)
    if hue:scale[-2:]=np.sqrt(var[-2:].sum())
    scale=np.where(scale>1e-12,scale,1.)
    return (a-mu)/scale,mu,scale


def r2(x,y,w=None):
    y=np.asarray(y,float)
    if y.ndim==1:y=y[:,None]
    w=np.full(len(y),1/len(y)) if w is None else np.asarray(w)/np.sum(w)
    x=normalize(x,w)[0];X=np.column_stack([np.ones(len(x)),x]);sw=np.sqrt(w)[:,None]
    b=np.linalg.lstsq(X*sw,y*sw,rcond=None)[0];mu=np.average(y,axis=0,weights=w)
    den=np.sum(w[:,None]*(y-mu)**2)
    if den<=1e-15:raise ValueError('Constant response')
    return float(1-np.sum(w[:,None]*(y-X@b)**2)/den)


def partition(e,s,y,w=None):
    re=r2(e,y,w);rs=r2(s,y,w);both=r2(np.column_stack([e,s]),y,w)
    return {'environment_total':re,'space_total':rs,'combined':both,'pure_environment':both-rs,
            'pure_space':both-re,'shared':re+rs-both,'unexplained':1-both}


def spherical(lat,lon):
    a=np.radians(lat);b=np.radians(lon);x=np.cos(a)*np.cos(b);y=np.cos(a)*np.sin(b);z=np.sin(a)
    return np.column_stack([x,y,z,x*y,x*z,y*z,x*x-y*y,3*z*z-1])


def taxon_summary(f,cols,w):
    x=f[cols].mul(w,axis=0);x['taxon_name']=f.taxon_name.to_numpy()
    den=pd.Series(w).groupby(f.taxon_name).sum()
    return x.groupby('taxon_name').sum().div(den,axis=0)


def demean(f,cols,w):
    a=taxon_summary(f,cols,w)
    return f[cols].to_numpy()-a.loc[f.taxon_name,cols].to_numpy()


def attribution(e,s,y):
    group_columns={k:[VARS.index(x) for x in v] for k,v in GROUPS.items()}
    names=list(GROUPS);cache={}
    for mask in range(1<<len(names)):
        cols=[c for i,g in enumerate(names) if mask&(1<<i) for c in group_columns[g]]
        cache[mask]=r2(e[:,cols],y) if cols else 0.
    full=(1<<len(names))-1;shapley={};unique={};joint=r2(np.column_stack([e,s]),y)
    for i,g in enumerate(names):
        shapley[g]=sum((cache[m|(1<<i)]-cache[m])/(len(names)*math.comb(len(names)-1,m.bit_count())) for m in cache if not m&(1<<i))
        keep=[c for c in range(e.shape[1]) if c not in group_columns[g]]
        unique[g]=joint-r2(np.column_stack([e[:,keep],s]),y)
    return {'shapley_in_sample':shapley,'unique_delta_given_other_environment_and_space':unique,
            'full_environment_r2':cache[full]}


def fold_ids(groups,n_folds=5):
    g=np.asarray(groups);unique,counts=np.unique(g,return_counts=True);load=np.zeros(n_folds,int);mapping={}
    for k,n in sorted(zip(unique,counts),key=lambda x:(-x[1],str(x[0]))):
        fold=int(np.argmin(load));mapping[k]=fold;load[fold]+=n
    return np.array([mapping[x] for x in g])


def ridge_predict(x,y,test,w=None,alpha=.1):
    z,mu,sd=normalize(x,w);test=(np.asarray(test)-mu)/sd;y=np.asarray(y)
    if y.ndim==1:y=y[:,None]
    w=np.full(len(y),1/len(y)) if w is None else np.asarray(w)/np.sum(w)
    ym=np.average(y,axis=0,weights=w);yc=y-ym
    coef=np.linalg.lstsq(z.T@(z*w[:,None])+alpha*np.eye(z.shape[1]),z.T@(yc*w[:,None]),rcond=None)[0]
    return ym+test@coef


def among_cv(e,y,geo,alpha=.1):
    ids=fold_ids(geo);err=base=0.;folds=[]
    for fold in np.unique(ids):
        tr=ids!=fold;te=ids==fold
        yt,ym,ys=normalize(y[tr],hue=y.shape[1]==4);yv=(y[te]-ym)/ys
        pred=ridge_predict(e[tr],yt,e[te],alpha=alpha)
        err+=float(np.sum((yv-pred)**2));base+=float(np.sum(yv*yv));folds.append(int(te.sum()))
    return {'r2':float(1-err/base),'fold_sizes':folds,'alpha':alpha,
            'scope':'whole taxon-centroid 10-degree blocks; species geographic ranges may overlap folds'}


def within_cv(f,cols,alpha=.1):
    # Hold out actual observation locations, not just taxon centroids. Only taxa
    # represented in training are scored. Training means are never recomputed
    # using the held-out response. This evaluates new locations of seen taxa.
    f=f.copy();f['_cv_weight']=weights(f)
    geo=np.floor((f.latitude+90)/10).astype(int)*36+np.floor((f.longitude+180)/10).astype(int)
    ids=fold_ids(geo);variants={'full':list(range(len(VARS)))}
    variants.update({'without_'+g:[i for i,v in enumerate(VARS) if v not in members] for g,members in GROUPS.items()})
    errs={k:0. for k in variants};base=0.;tested=0;excluded=0;folds=[]
    for fold in np.unique(ids):
        tr=f[ids!=fold].reset_index(drop=True);te=f[ids==fold].reset_index(drop=True)
        ns=tr.groupby('taxon_name').size();keep=ns[ns>=3].index;tr=tr[tr.taxon_name.isin(keep)].reset_index(drop=True)
        excluded+=int((~te.taxon_name.isin(keep)).sum());te=te[te.taxon_name.isin(keep)].reset_index(drop=True)
        if len(te)==0:continue
        w=weights(tr);wv=te._cv_weight.to_numpy();means=taxon_summary(tr,cols+VARS,w)
        X=tr[VARS].to_numpy()-means.loc[tr.taxon_name,VARS].to_numpy();Y=tr[cols].to_numpy()-means.loc[tr.taxon_name,cols].to_numpy()
        XV=te[VARS].to_numpy()-means.loc[te.taxon_name,VARS].to_numpy();YV=te[cols].to_numpy()-means.loc[te.taxon_name,cols].to_numpy()
        Y,ym,ys=normalize(Y,w,hue=len(cols)==4);YV=(YV-ym)/ys
        base+=float(np.sum(wv[:,None]*YV*YV))
        for name,index in variants.items():
            pred=ridge_predict(X[:,index],Y,XV[:,index],w=w,alpha=alpha)
            errs[name]+=float(np.sum(wv[:,None]*(YV-pred)**2))
        tested+=len(te);folds.append({'test_observations':len(te),'train_observations':len(tr),'tested_taxa':te.taxon_name.nunique()})
    scores={k:float(1-v/base) for k,v in errs.items()}
    return {'full_r2':scores['full'],'drop_block_delta_r2':{g:scores['full']-scores['without_'+g] for g in GROUPS},
            'tested_observations':tested,'excluded_unseen_or_sparse_taxon_rows':excluded,'folds':folds,'alpha':alpha,
            'scope':'5 folds of actual 10-degree observation blocks; training-only species means; known-taxon transfer; no buffer'}


def analyse(frame,run_within_cv=True):
    result={}
    for name,cols in MODULES.items():
        f=frame.dropna(subset=cols+VARS).copy();n=f.groupby('taxon_name').size()
        f=f[f.taxon_name.isin(n[n>=5].index)].reset_index(drop=True);w=weights(f)
        S=spherical(f.latitude,f.longitude);sc=['s'+str(i) for i in range(S.shape[1])];f[sc]=S
        a=taxon_summary(f,cols+VARS+sc,w);E=a[VARS].to_numpy();A=a[sc].to_numpy();Y=normalize(a[cols].to_numpy(),hue=name=='colour')[0]
        WYE=normalize(demean(f,cols,w),w,hue=name=='colour')[0]
        pa=partition(E,A,Y);pw=partition(demean(f,VARS,w),demean(f,sc,w),WYE,w)
        lat=np.degrees(np.arctan2(A[:,2],np.hypot(A[:,0],A[:,1])));lon=np.degrees(np.arctan2(A[:,1],A[:,0]));geo=np.floor((lat+90)/10).astype(int)*36+np.floor((lon+180)/10).astype(int)
        cv=among_cv(E,a[cols].to_numpy(),geo);drop={}
        for g,members in GROUPS.items():
            keep=[i for i,v in enumerate(VARS) if v not in members]
            drop[g]=cv['r2']-among_cv(E[:,keep],a[cols].to_numpy(),geo)['r2']
        result[name]={'n_observations':len(f),'n_taxa':len(a),'among_partition':pa,'within_partition':pw,
                      'attribution':attribution(E,A,Y),'among_cv':cv,'among_cv_drop_block_delta_r2':drop,
                      'within_cv':within_cv(f,cols) if run_within_cv else None}
        print(name,result[name]['n_observations'],result[name]['n_taxa'],cv['r2'],flush=True)
    return result
