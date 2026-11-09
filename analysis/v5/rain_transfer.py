"""Native climate models of photographed capitulum angle transferred across ranges."""
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from sklearn.model_selection import GroupKFold

TARGET='orientation_angle_degrees_median'
AXES=['bio18','bio05','rsds_mean']
MATCH=['bio01','bio04','bio12','bio15','rsds_mean']
SEED=20261109

def cell_means(d):
    d=d.loc[d[TARGET].between(0,180)].dropna(subset=AXES+MATCH).copy()
    # Keep novelty groups separate in cells with both types of climate.
    d['block']=(np.minimum(8,((d.latitude+90)//20).astype(int))*18
                +(((d.longitude+180)%360)//20).astype(int))
    keys=['_cell','block','novel']
    cols=list(dict.fromkeys(AXES+MATCH+[TARGET]))
    g=d.groupby(keys,as_index=False)[cols].mean()
    g['photos']=d.groupby(keys).size().values
    return g

def design(d):
    return np.column_stack([np.ones(len(d)),d[AXES]])

def fit(x,y):
    if np.linalg.matrix_rank(x)!=x.shape[1]:raise ValueError('Climate design lacks rank')
    return np.linalg.lstsq(x,y,rcond=None)[0]

def metrics(y,pred,baseline):
    return dict(observed=y.mean(),predicted=pred.mean(),bias=(y-pred).mean(),
                rmse=np.sqrt(np.mean((y-pred)**2)),q2=1-np.sum((y-pred)**2)/np.sum((y-baseline)**2))

def resample_blocks(labels,rng):
    groups=np.unique(labels)
    return np.concatenate([np.flatnonzero(labels==b) for b in rng.choice(groups,len(groups),replace=True)])

def run_transfer(sp,sn,si,out):
    native=sn.rename(columns=lambda c:c.removeprefix('chelsa_')).copy()
    native['novel']=-1
    n=cell_means(native);i=cell_means(si)
    x=design(n); y=n[TARGET].to_numpy(); coef=fit(x,y)
    pred=np.zeros(len(n));base=np.zeros(len(n))
    for train,test in GroupKFold(5).split(x,groups=n.block):
        pred[test]=x[test]@fit(x[train],y[train]);base[test]=y[train].mean()
    rows=[];direction=[];rain_model=[]
    rows.append(dict(species=sp,group='native_spatial_cv',n_cells=len(n),n_blocks=n.block.nunique(),n_photos=int(n.photos.sum()),**metrics(y,pred,base)))
    native_blocks=n.block.to_numpy()
    rng=np.random.default_rng(SEED)
    bs=[resample_blocks(native_blocks,rng) for _ in range(999)]
    coefs=np.array([fit(x[t],y[t]) for t in bs])
    scale=n[AXES].std(ddof=0).to_numpy()
    for j,a in enumerate(AXES,1):
        l,h=np.percentile(coefs[:,j]*scale[j-1],[2.5,97.5])
        rain_model.append(dict(species=sp,axis=a,native_sd=scale[j-1],deg_per_sd=coef[j]*scale[j-1],lower=l,upper=h))
    # Nearest native cells in the original five climatic coordinates (no head traits or BIO18).
    ref=n[MATCH].to_numpy(); sd=ref.std(axis=0); ref=ref/sd
    for status,name in [(0,'introduced_analogous'),(1,'introduced_novel')]:
        g=i[i.novel==status].reset_index(drop=True); xx=design(g); yy=g[TARGET].to_numpy();bb=g.block.to_numpy()
        pp=xx@coef
        row=dict(species=sp,group=name,n_cells=len(g),n_blocks=g.block.nunique(),n_photos=int(g.photos.sum()),**metrics(yy,pp,np.full(len(g),y.mean())))
        dist=cdist(g[MATCH].to_numpy()/sd,ref)
        rain=g.bio18.to_numpy()
        lo,hi=n.bio18.quantile([.05,.95])
        row['fraction_bio18_outside_native_central90']=np.mean((rain<lo)|(rain>hi))
        row['median_5th_neighbor_distance']=np.median(np.sort(dist,axis=1)[:,4])
        sims=[];diffs=[]
        for rep,draw in enumerate(bs):
            t=resample_blocks(bb,rng)
            prediction=xx[t]@coefs[rep]; residual=yy[t]-prediction
            sims.append([residual.mean(),1-np.sum(residual**2)/np.sum((yy[t]-y[draw].mean())**2)])
            # Replace native reference and recompute matches for every bootstrap draw.
            nn=draw[np.argsort(dist[t][:,draw],axis=1)[:,:5]]
            dr=rain[t]-n.bio18.to_numpy()[nn].mean(1)
            da=yy[t]-y[nn].mean(1)
            diffs.append([da[dr>0].mean() if (dr>0).any() else np.nan,
                          da[dr<0].mean() if (dr<0).any() else np.nan,
                          np.cov(dr,da,ddof=0)[0,1]/np.var(dr) if np.var(dr)>0 else np.nan])
        bounds=np.percentile(sims,[2.5,97.5],axis=0)
        row.update(bias_lo=bounds[0,0],bias_hi=bounds[1,0],q2_lo=bounds[0,1],q2_hi=bounds[1,1])
        rows.append(row)
        nn=np.argsort(dist,axis=1)[:,:5]
        dr=rain-n.bio18.to_numpy()[nn].mean(1);da=yy-y[nn].mean(1)
        bounds=np.nanpercentile(diffs,[2.5,97.5],axis=0)
        for j,tag in enumerate(['wetter','drier','slope_deg_per_mm_month']):
            ix=(dr>0 if j==0 else dr<0) if j<2 else np.ones(len(g),bool)
            val=da[ix].mean() if j<2 else np.cov(dr,da,ddof=0)[0,1]/np.var(dr)
            direction.append(dict(species=sp,group=name,comparison=tag,n=int(ix.sum()),rain_difference=dr[ix].mean(),angle_estimate=val,
                                  lower=bounds[0,j],upper=bounds[1,j]))
    for filename,table in [('transfer',rows),('rain_displacement',direction),('native_coefficients',rain_model)]:
        path=out/(filename+'.csv')
        pd.DataFrame(table).to_csv(path,mode='a',header=not path.exists(),index=False,columns=list(dict.fromkeys(k for r in table for k in r)))

