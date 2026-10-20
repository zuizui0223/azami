"""Observer-balanced within-versus-among presentation/colour covariance (v4 exploration).
Only aggregate RV summaries are written; observational frequency is not biological
abundance. Taxon bootstrap plus one weighted observation per sampled taxon is a
descriptive finite-sample contrast, not a permutation p value.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from analysis.v4.observer_bias import (
    ALL_TRAITS, ENV, build_inputs, weights
)

def rv(matrix,left,right):
    a=matrix[np.ix_(left,left)]
    b=matrix[np.ix_(right,right)]
    ab=matrix[np.ix_(left,right)]
    den=np.sqrt(np.sum(a*a)*np.sum(b*b))
    return float(np.sum(ab*ab)/den) if den>1e-20 else float('nan')

def integration(frame, weight="equal_taxon_cell_observer", omit=(), seed=20261020, reps=400):
    f=frame.loc[~frame.taxon.isin(omit)].dropna(subset=ALL_TRAITS+ENV).copy()
    f=f.loc[f.groupby('taxon').obs_id.transform('size')>=5].reset_index(drop=True)
    w=weights(f,weight)
    Y=f[ALL_TRAITS].to_numpy(float)
    codes,taxa=pd.factorize(f.taxon,sort=True)
    n=len(taxa)
    if n<20:raise ValueError('Too few full-colour native taxa')
    means=[];covs=[];choices=[];cdf=[]
    for j in range(n):
        indices=np.flatnonzero(codes==j)
        p=w[indices]/w[indices].sum()
        center=np.average(Y[indices],axis=0,weights=p)
        means.append(center)
        r=Y[indices]-center
        covs.append((r*p[:,None]).T@r)
        choices.append(indices)
        cdf.append(np.cumsum(p))
    means=np.asarray(means);covs=np.asarray(covs)
    W=covs.mean(axis=0)
    M=means-means.mean(axis=0)
    B=M.T@M/n
    sd=np.sqrt(np.diag(W+B));sd=np.where(sd>1e-12,sd,1.)
    denom=np.outer(sd,sd)
    within=rv(W/denom,[0],[1,2,3,4])
    among=rv(B/denom,[0],[1,2,3,4])
    rng=np.random.default_rng(seed)
    differences=[]
    for _ in range(reps):
        taxon_samples=rng.integers(n,size=n)
        current=means[taxon_samples]
        mc=current-current.mean(axis=0)
        between=mc.T@mc/n
        average_within=covs[taxon_samples].mean(axis=0)
        scale=np.sqrt(np.maximum(np.diag(between+average_within),1e-20))
        mat=np.outer(scale,scale)
        sampled=np.stack([
            Y[choices[k][min(len(choices[k])-1,int(np.searchsorted(cdf[k],rng.random(),side='right')))]]-means[k]
            for k in taxon_samples])
        sampled-=sampled.mean(axis=0)
        matched=sampled.T@sampled/n
        differences.append(rv(between/mat,[0],[1,2,3,4])-rv(matched/mat,[0],[1,2,3,4]))
    differences=np.asarray(differences)
    if not np.isfinite(differences).all():raise ValueError('Invalid bootstrap RV')
    return {'method':weight,'n_observations':len(f),'n_taxa':n,
            'within_rv':within,'among_rv':among,
            'raw_difference':among-within,
            'matched_taxon_bootstrap_median':float(np.median(differences)),
            'matched_taxon_bootstrap_95pct':[float(v) for v in np.quantile(differences,[.025,.975])],
            'fraction_positive':float((differences>0).mean()),
            'replicates':reps,
            'limits':'Taxon-bootstrap + matched-n within reference. In-sample photo phenotypes; not genetic or adaptive covariance.'}

def main():
    p=argparse.ArgumentParser()
    for key in ['native-roster','observer-private-csv','archived-traits','native-climate','out-dir']:
        p.add_argument('--'+key,required=True)
    p.add_argument('--replicates',type=int,default=400)
    a=p.parse_args()
    native=pd.read_csv(a.native_roster,dtype={'obs_id':str})
    observer=pd.read_csv(a.observer_private_csv,dtype={'obs_id':str,'user_id':str})
    measured=pd.read_csv(a.archived_traits,dtype={'obs_id':str},
                         usecols=['obs_id','taxon_name',*ALL_TRAITS])
    climate=pd.read_csv(a.native_climate,dtype={'obs_id':str},
                        usecols=['obs_id',*ENV])
    source,frame=build_inputs(native,observer,measured,climate)
    top=source.taxon.value_counts()
    rows=[]
    for method in ['equal_taxon_cell','equal_taxon_cell_observer']:
        for omit, names in [('none',[]),('top2',list(top.index[:2])),
                            ('top5',list(top.index[:5]))]:
            row=integration(frame,method,names,reps=a.replicates)
            row['omitted_source_taxa']=omit
            rows.append(row)
    out=Path(a.out_dir);out.mkdir(exist_ok=True,parents=True)
    (out/'observer_covariance_summary.json').write_text(
         json.dumps({'analysis':'observer-balanced exploratory within/among RV',
                     'results':rows,
                     'source_scope':{'native_observations':len(source),
                                     'native_taxa':int(source.taxon.nunique())}},
                    ensure_ascii=False,indent=2)+'\n')
    print('Six observer/cell RV conditions completed',len(rows))
if __name__=='__main__':
    main()
