"""Native-only within/among covariance with taxon and matched-n uncertainty.

Functions receive already native-filtered observations. Neither the source
eligibility rule nor the phenotype definition is estimated from model outcomes.
"""
import numpy as np
import pandas as pd


def moments(frame,columns,min_n=5,weighting='taxon_cell'):
    f=frame.dropna(subset=columns).copy()
    n=f.groupby('taxon_name').size();f=f[f.taxon_name.isin(n[n>=min_n].index)].reset_index(drop=True)
    if f.taxon_name.nunique()<3:raise ValueError('Fewer than three eligible taxa')
    Y=f[columns].to_numpy(float)
    if not np.isfinite(Y).all():raise ValueError('Non-finite trait')
    if weighting=='taxon_cell':
        ncell=f.groupby(['taxon_name','cell']).taxon_name.transform('size').to_numpy(float)
        k=f.groupby('taxon_name').cell.transform('nunique').to_numpy(float);raw=1/(ncell*k)
    elif weighting=='taxon':raw=1/f.groupby('taxon_name').taxon_name.transform('size').to_numpy(float)
    else:raise ValueError('Use taxon_cell or taxon')
    w=raw/raw.sum();mu=[];cov=[];groups=[];R=np.empty_like(Y)
    taxa=[]
    for name,idx in f.groupby('taxon_name',sort=True).indices.items():
        p=w[idx]/w[idx].sum();mean=np.average(Y[idx],axis=0,weights=p);r=Y[idx]-mean
        mu.append(mean);cov.append((r*p[:,None]).T@r);R[idx]=r
        groups.append((idx,np.cumsum(p)));taxa.append(str(name))
    M=np.asarray(mu);V=np.asarray(cov);W=V.mean(0);c=M-M.mean(0);B=c.T@c/len(M)
    return {'n_observations':len(f),'n_taxa':len(M),'taxa':taxa,'Y':Y,'weights':w,
            'means':M,'taxon_covariances':V,'within':W,'among':B,'residuals':R,'groups':groups}


def rv(C,left,right):
    A=C[np.ix_(left,left)];B=C[np.ix_(right,right)];AB=C[np.ix_(left,right)]
    d=np.sqrt(np.sum(A*A)*np.sum(B*B))
    return float(np.sum(AB*AB)/d) if d>1e-20 else float('nan')


def scaling(W,B,colour_metric):
    v=np.diag(W+B);s=np.sqrt(v)
    if colour_metric=='joint_hue':s[-2:]=np.sqrt(v[-2:].sum())
    elif colour_metric!='components':raise ValueError('Unknown colour metric')
    return np.where(s>1e-12,s,1.)


def integration(frame,columns,replicates=1000,seed=20261020,min_n=5,
                weighting='taxon_cell',colour_metric='components'):
    if len(columns)!=5:raise ValueError('Expected angle, L*, C*, hue sin, hue cos in order')
    m=moments(frame,columns,min_n,weighting);M=m['means'];V=m['taxon_covariances'];W=m['within'];B=m['among']
    s=scaling(W,B,colour_metric);SS=np.outer(s,s);left=[0];right=[1,2,3,4]
    wa=rv(W/SS,left,right);ba=rv(B/SS,left,right);rng=np.random.default_rng(seed)
    nested=[];raw=[];n=m['n_taxa'];R=m['residuals'];groups=m['groups']
    for _ in range(replicates):
        selected=rng.integers(n,size=n);mb=M[selected];mc=mb-mb.mean(0)
        b=mc.T@mc/n;w=V[selected].mean(0);s=scaling(w,b,colour_metric);ss=np.outer(s,s)
        ar=rv(b/ss,left,right);raw.append(ar-rv(w/ss,left,right))
        picks=[]
        for j in selected:
            idx,cdf=groups[j];picks.append(idx[min(len(idx)-1,np.searchsorted(cdf,rng.random(),side='right'))])
        r=R[picks];r=r-r.mean(0);c=r.T@r/n
        nested.append(ar-rv(c/ss,left,right))
    nested=np.asarray(nested);raw=np.asarray(raw)
    if not np.isfinite(nested).all():raise ValueError('Degenerate bootstrap replicate; no silent omission')
    return {'n_observations':m['n_observations'],'n_taxa':n,'columns':columns,'weighting':weighting,
            'colour_metric':colour_metric,'within_rv':wa,'among_rv':ba,
            'raw_taxon_bootstrap_delta_interval':np.quantile(raw,[.025,.975]).tolist(),
            'nested_delta_median':float(np.median(nested)),
            'nested_delta_interval':np.quantile(nested,[.025,.975]).tolist(),
            'nested_positive_fraction':float(np.mean(nested>0)),
            'replicates':replicates,'seed':seed,
            'interval_definition':'taxon bootstrap plus one weighted centred observation per sampled taxon; not a permutation P value',
            'limitations':'photographic measurement error remains; resampling cannot identify genetic covariance'}
