"""V5 source-matched comparative analysis. All exported tables are aggregates."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.optimize import minimize_scalar
from statsmodels.stats.multitest import multipletests
from sklearn.model_selection import GroupKFold
from analysis.v4.s1_s2_extreme_lego import read_tables, ALL_ENV
from analysis.v4.s1_s2_vif5_si10_alltraits import selected_source
from analysis.v4.s1_s2_extreme_range_models import _get_dated_balanced, _formula
from analysis.v4.s1_climate_trait_hypervolume import load_source, geo_aggregate, standardise_cells
from analysis.v4 import s1_h3_climate_analogue_alignment as analog

ROOT=Path('work/submission')
PRIVATE=ROOT/'private'
OUT=ROOT/'v5'
TRAITS={'angle':'orientation_image_vertical_angle','chroma':'corolla_lab_chroma','aspect':'capitulum_outline_aspect_ratio'}
SPECS={'warm':['bio18','bio05','rsds_mean'],'warm_elev':['bio18','bio05','rsds_mean','elevation'],
       'annual':['bio12','bio01','rsds_mean'],'dry':['bio18','bio05','rsds_mean','bio14'],
       'seasonality':['bio18','bio01','rsds_mean','bio04','bio15'],
       'extremes':['bio13','bio17','bio06','rsds_mean'],
       'water':['bio12','bio14','bio15'], 'thermal':['bio01','bio04','bio05','bio06']}
FOCUS={'angle':'bio18','chroma':'rsds_mean','aspect':'bio18'}
SEED=20261109

def rename_env(frame):
    return frame.rename(columns={k:k.removeprefix('chelsa_').removeprefix('topo_') for k in ALL_ENV})

def adjust(df, groups):
    df['holm']=df.groupby(groups).p.transform(lambda a:multipletests(a,method='holm')[1])
    return df

def blocks(frame,deg):
    lat=np.clip(np.floor((frame.latitude+90)/deg),0,180/deg-1).astype(int)
    lon=np.floor(((frame.longitude+180)%360)/deg).astype(int)
    return lat*int(360/deg)+lon

def cluster_fit(y,x,w,groups):
    assert np.linalg.matrix_rank(x)==x.shape[1]
    return sm.WLS(y,x,weights=w).fit(cov_type='cluster',cov_kwds={'groups':np.asarray(groups)},use_t=True)

def stats_row(fit,j):
    lo,hi=np.asarray(fit.conf_int())[j]
    return dict(beta=float(np.asarray(fit.params)[j]),lower=lo,upper=hi,p=float(np.asarray(fit.pvalues)[j]))

def heterogeneity(cells,axes,yname,focal=None):
    """REML synthesis of partial slopes of well-sampled taxa, common units."""
    estimates=[]; variances=[]
    for _,g in cells.groupby('accepted_key'):
        if len(g)<12: continue
        x=sm.add_constant(g[axes].to_numpy(),has_constant='add')
        if np.linalg.matrix_rank(x)<x.shape[1] or np.linalg.cond(x)>100: continue
        fit=sm.OLS(g[yname],x).fit(cov_type='HC3')
        j=axes.index(focal or FOCUS[yname])+1
        v=fit.bse.iloc[j]**2 if hasattr(fit.bse,'iloc') else fit.bse[j]**2
        if np.isfinite(v) and v>0:
            estimates.append(np.asarray(fit.params)[j]);variances.append(v)
    b=np.array(estimates);v=np.array(variances)
    if len(b)<3:return dict(n_taxa=len(b),status='insufficient')
    def objective(t):
        w=1/(v+t)
        return np.log(v+t).sum()+np.log(w.sum())+np.sum(w*(b-np.average(b,weights=w))**2)
    opt=minimize_scalar(objective,bounds=(0,max(1,5*np.var(b))),method='bounded')
    tau2=opt.x if objective(opt.x)<objective(0) else 0
    loo=[np.average(np.delete(b,k),weights=1/np.delete(v+tau2,k)) for k in range(len(b))]
    return dict(loo_mean_min=min(loo),loo_mean_max=max(loo),n_taxa=len(b),mean_partial_slope=np.average(b,weights=1/(v+tau2)),tau=np.sqrt(tau2),
                observed_positive=int((b>0).sum()),observed_negative=int((b<0).sum()))

def s1(native):
    slopes=[];scores=[];variance=[];random=[]
    for label,trait in TRAITS.items():
        frame,_=selected_source(native,[trait])
        frame=rename_env(frame); frame['_cell']=blocks(frame,5)
        axes=sorted(set(v for spec in SPECS.values() for v in spec))
        cells=frame.groupby(['accepted_key','_cell'])[[trait,'latitude','longitude',*axes]].mean().reset_index().rename(columns={trait:label})
        species=cells.accepted_key; region=blocks(cells,20)
        w=(1/species.map(species.value_counts())).to_numpy()
        for k in [label,*axes]:
            val=cells[k];mean=np.average(val,weights=w)
            cells[k]=(val-mean)/np.sqrt(np.average((val-mean)**2,weights=w))
        y=cells[label].to_numpy()
        means=cells.groupby('accepted_key')[axes+[label]].transform('mean')
        variance.append(dict(trait=label,n_photo=len(frame),n_taxa=species.nunique(),taxon_cells=len(cells),
                             blocks=region.nunique(),observed_between_share=np.average(means[label]**2,weights=w)))
        geo=pd.get_dummies(region,drop_first=True,dtype=float).to_numpy()
        for spec,axis in SPECS.items():
            between=means[axis].to_numpy(); within=cells[axis].to_numpy()-between
            design=np.column_stack([np.ones(len(y)),between,within]); terms=[(s,k) for s in ['between','within'] for k in axis]
            for geo_adjust in [False,True]:
                x=np.column_stack([design,geo]) if geo_adjust else design
                # Cluster by geographic block for all fits, keep point estimates identical across SE choices.
                fit=cluster_fit(y,x,w,region)
                for j,(scale,key) in enumerate(terms,1):
                    slopes.append(dict(trait=label,spec=spec,geo_adjust=geo_adjust,scale=scale,axis=key,**stats_row(fit,j)))
            err=0.;null=0.
            for tr,te in GroupKFold(n_splits=5).split(design,groups=region):
                coef=np.linalg.lstsq(design[tr]*np.sqrt(w[tr,None]),y[tr]*np.sqrt(w[tr]),rcond=None)[0]
                err+=np.sum(w[te]*(y[te]-design[te]@coef)**2)
                null+=np.sum(w[te]*(y[te]-np.average(y[tr],weights=w[tr]))**2)
            covariance=np.cov(between.T,aweights=w)
            sd=np.sqrt(np.diag(covariance))
            corr=covariance/np.outer(sd,sd)
            scores.append(dict(trait=label,spec=spec,spatial_Q2=1-err/null,
                               max_VIF_between=np.diag(np.linalg.inv(corr)).max()))
            if spec in ['warm','warm_elev']:
                random.append(dict(trait=label,spec=spec,axis=FOCUS[label],**heterogeneity(cells,axis,label)))
                if spec=='warm' and label=='chroma':
                    random.append(dict(trait=label,spec=spec,axis='bio05',**heterogeneity(cells,axis,label,'bio05')))
    adjust(pd.DataFrame(slopes),['spec','geo_adjust']).to_csv(OUT/'s1_slopes.csv',index=False)
    for file,rows in [('s1_spatial_cv',scores),('s1_variance',variance),('s1_slope_heterogeneity',random)]:
        pd.DataFrame(rows).to_csv(OUT/(file+'.csv'),index=False)

def s2(focal):
    results=[];support=[]
    for (species,g) in focal.groupby('taxon_name'):
        basis=rename_env(g.loc[np.isfinite(g[ALL_ENV]).all(axis=1)])
        for label,trait in TRAITS.items():
            usable=basis[np.isfinite(basis[trait])]
            if label=='angle': usable=usable[usable[trait].between(0,180)]
            balanced=pd.concat([_get_dated_balanced(h) for _,h in usable.groupby('range_state')],ignore_index=True)
            for spec,axes in SPECS.items():
                # Use climate-only overlap even for elevation sensitivity.
                overlap_axes=SPECS['warm'] if spec=='warm_elev' else axes
                keep=np.ones(len(balanced),bool)
                for a in overlap_axes:
                    bounds=balanced.groupby('range_state')[a].quantile([.05,.95]).unstack()
                    keep &= balanced[a].between(bounds[.05].max(),bounds[.95].min()).to_numpy()
                sub=balanced.loc[keep].copy()
                counts=sub.groupby('range_state')._cell.nunique()
                support.append(dict(species=species,trait=label,spec=spec,n=len(sub),original=len(balanced),
                                    cells_native=counts.get('native',0),cells_introduced=counts.get('introduced',0)))
                if len(counts)!=2 or counts.min()<8:continue
                design,_,_,_,_=_formula(sub,axes)
                status=(sub.range_state=='introduced').to_numpy(float)
                x=np.column_stack([design,design*status[:,None]])
                # Remove redundant date terms in a stratum if unidentifiable; preserve climatic contrasts.
                if np.linalg.matrix_rank(x)<x.shape[1]:continue
                w=1/sub.groupby(['range_state','_cell'])[trait].transform('size')
                w=w/w.groupby(sub.range_state).transform('sum')
                y=(sub[trait]-sub[trait].mean())/sub[trait].std(ddof=0)
                fit=cluster_fit(y,x,w,sub._cell)
                for j,a in enumerate(axes,1):
                    difference=stats_row(fit,design.shape[1]+j)
                    native=float(np.asarray(fit.params)[j])
                    results.append(dict(species=species,trait=label,spec=spec,axis=a,n=len(sub),native=native,
                                        introduced=native+difference['beta'],**difference))
    adjust(pd.DataFrame(results),['spec']).to_csv(OUT/'s2_interactions.csv',index=False)
    pd.DataFrame(support).to_csv(OUT/'s2_support.csv',index=False)

def s3():
    source=PRIVATE/'decoded'
    frame,_=load_source(source/'observations.csv.gz',source/'measurement_receipt.json')
    rows=[]
    for deg in [1,5]:
        cells,coverage=geo_aggregate(frame,deg)
        z,_=standardise_cells(cells);taxa=analog.eligible_taxa(z)
        for reverse in [False,True]:
            for threshold in [.7,.9,1.1]:
                rng=np.random.default_rng(SEED)
                nullrng=np.random.default_rng(SEED+1)
                obs=[];null=[];pairs=[]
                original=analog.CLIMATE_ANALOG_Z
                try:
                    analog.CLIMATE_ANALOG_Z=threshold
                    for _ in range(24):
                        x,y,sector=analog._sample_taxa(taxa,rng)
                        if reverse:x,y=y,x
                        pair,meta=analog.matched_pairs(x,y,sector)
                        if meta['status']!='qualified':continue
                        value,permuted=analog.observed_and_permutations(pair,len(taxa),nullrng,199)
                        obs.append(value);null.append(permuted);pairs.append(meta['n_pairs'])
                finally:analog.CLIMATE_ANALOG_Z=original
                if len(obs)<12:
                    rows.append(dict(deg=deg,reverse=reverse,threshold=threshold,n_taxa=len(taxa),draws=len(obs)))
                    continue
                value=np.mean(obs);perm=np.mean(null,axis=0);avg=perm.mean()
                p=(1+sum(abs(perm-avg)>=abs(value-avg)))/200
                rows.append(dict(deg=deg,reverse=reverse,threshold=threshold,n_taxa=len(taxa),draws=len(obs),
                                 mean_pairs=np.mean(pairs),reduction=value,null_mean=avg,excess=value-avg,p=p))
    result=pd.DataFrame(rows);chosen=result.threshold.eq(.9)&result.p.notna()
    result.loc[chosen,'holm_primary']=multipletests(result.loc[chosen,'p'],method='holm')[1]
    result.to_csv(OUT/'s3_analogues.csv',index=False)

def main():
    OUT.mkdir(exist_ok=True)
    native,focal=read_tables(PRIVATE/'decoded/observations.csv.gz',PRIVATE/'focal_wild_five.csv.gz',
                            PRIVATE/'extended/environment_ext.csv.gz',PRIVATE/'original-measured.csv.gz')
    s1(native);print('S1 complete',flush=True)
    s2(focal);print('S2 complete',flush=True)
    s3();print('S3 complete',flush=True)

if __name__=='__main__':
    main()
