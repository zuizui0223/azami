"""Held-out geography: does known species identity add predictive information?"""
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.model_selection import GroupKFold
from analysis.v5.run import TRAITS, FOCUS, blocks, cluster_fit, stats_row, adjust

def run():
    root=Path("work/bio18_geo")
    f=pd.read_csv(root/"private/decoded/observations.csv.gz",low_memory=False,dtype={"obs_id":str})
    e=pd.read_csv(root/"private/extended-source.csv.gz",dtype={"obs_id":str})
    f=f.merge(e,on="obs_id",validate="one_to_one")
    assert len(f)==313350
    base=['chelsa_bio18','chelsa_bio05','chelsa_rsds_mean']
    results=[]; slopes=[]
    for trait,column in TRAITS.items():
        ok=f[column+'__eligible'].fillna(False).eq(True)&np.isfinite(f[base+[column,'latitude','longitude']]).all(axis=1)
        if trait=='angle':ok &= f[column].between(0,180)
        d=f.loc[ok].copy()
        counts=d.accepted_key.value_counts()
        d=d[d.accepted_key.isin(counts[counts>=5].index)]
        d['_cell']=blocks(d,5)
        c=d.groupby(['accepted_key','_cell'])[base+[column,'latitude','longitude']].mean().reset_index()
        c=c.rename(columns={a:a.removeprefix('chelsa_') for a in base})
        t=c.accepted_key; blk=blocks(c,20).to_numpy(); w=(1/t.map(t.value_counts())).to_numpy()
        y=c[column].to_numpy(); ysd=np.sqrt(np.average((y-np.average(y,weights=w))**2,weights=w))
        y=(y-np.average(y,weights=w))/ysd
        lat,lon=np.radians(c.latitude),np.radians(c.longitude)
        a,b,h=np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)
        geography=np.column_stack([a,b,h,a*b,a*h,b*h,a*a-b*b,3*h*h-1])
        exposure=FOCUS[trait]
        z=c[['bio18','bio05','rsds_mean']].copy()
        for v in z:
            m=np.average(z[v],weights=w)
            z[v]=(z[v]-m)/np.sqrt(np.average((z[v]-m)**2,weights=w))
        means=z.groupby(t).transform('mean');within=z-means
        for group,axes in [('focal_only',[exposure]),('conditional',list(z))]:
            X=np.column_stack([np.ones(len(c)),means[axes],within[axes]])
            for geo in [False,True]:
                design=np.column_stack([X,geography]) if geo else X
                fit=cluster_fit(y,design,w,blk)
                focal_axes=([exposure] if group=='focal_only' or trait=='angle' else ['bio05','rsds_mean'])
                if trait=='aspect' and group=='focal_only':continue
                for axis in focal_axes:
                    for scale,j in [('between',axes.index(axis)+1),('within',axes.index(axis)+1+len(axes))]:
                        slopes.append(dict(trait=trait,axis=axis,scale=scale,spec=group,geo=geo,n_photo=len(d),
                                           n_taxa=t.nunique(),**stats_row(fit,j)))
        # Identical held-out cells and weights for all competing predictors.
        err={k:0. for k in ['mean','species_only','climate_geo','climate_geo_species']}
        ntest=0;seen=set(); fold_rows=[]
        X=np.column_stack([np.ones(len(c)),means,within,geography])
        for tr,te in GroupKFold(n_splits=5).split(c,groups=blk):
            te=te[t.iloc[te].isin(t.iloc[tr]).to_numpy()]
            if len(te)==0:continue
            ntest+=len(te);seen.update(t.iloc[te])
            raw=y[tr]; ww=(1/t.iloc[tr].map(t.iloc[tr].value_counts())).to_numpy()
            avg=np.average(raw,weights=ww)
            tab=pd.DataFrame({'species':t.iloc[tr].to_numpy(),'y':raw})
            ybar=tab.groupby('species').y.mean()
            dtr=pd.DataFrame(X[tr]);dtr['species']=t.iloc[tr].to_numpy()
            xbar=dtr.groupby('species').mean()
            centered=X[tr]-xbar.loc[t.iloc[tr]].to_numpy()
            def solve(matrix,rhs):
                return np.linalg.lstsq(matrix*np.sqrt(ww[:,None]),rhs*np.sqrt(ww),rcond=1e-9)[0]
            pooled=solve(X[tr],raw)
            # Identical between-species columns cancel after species demeaning.
            with_id=solve(centered,raw-ybar.loc[t.iloc[tr]].to_numpy())
            pred={'mean':avg,'species_only':ybar.loc[t.iloc[te]].to_numpy(),'climate_geo':X[te]@pooled,
                  'climate_geo_species':ybar.loc[t.iloc[te]].to_numpy()+(X[te]-xbar.loc[t.iloc[te]].to_numpy())@with_id}
            for name,value in pred.items():
                err[name]+=np.sum(w[te]*(y[te]-value)**2)
        for name in err:
            results.append(dict(trait=trait,model=name,Q2=1-err[name]/err['mean'],
                                tested_cells=ntest,all_cells=len(c),tested_taxa=len(seen),n_photo=len(d)))
    pd.DataFrame(results).to_csv(root/'public/v5_species_identity.csv',index=False)
    adjust(pd.DataFrame(slopes),['spec','geo']).to_csv(root/'public/v5_focal_no_elevation_gate.csv',index=False)
if __name__=='__main__':
    run()
