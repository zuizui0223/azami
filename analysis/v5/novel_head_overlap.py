"""Does introduced novelty reduce head-trait overlap with the native distribution?"""
from pathlib import Path
from itertools import product
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist

T=['orientation_angle_degrees_median','corolla_lab_chroma_median','shape_aspect_ratio_median']
SEED=20261109

def make_cells(frame, native=False):
    # Use a single unit per 1-degree cell, preventing paired pseudo-replication.
    d=frame.loc[frame[T[0]].between(0,180)].dropna(subset=T).copy()
    if native:d['novel']=0
    d['i']=np.minimum(179,np.floor(d.latitude+90)).astype(int)
    d['j']=np.floor((d.longitude+180)%360).astype(int)
    g=d.groupby(['i','j'])
    c=g[T+['latitude','longitude','novel']].median().reset_index()
    c['photos']=g.size().values
    c['block']=c.i.floordiv(20)*18+c.j.floordiv(20)
    return c

def xyz(c):
    lat,lon=np.radians(c.latitude),np.radians(c.longitude)
    return np.column_stack([np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)])

def match_cells(c):
    result=[]
    for block,g in c.groupby('block'):
        base=g[g.novel==0]; novel=g[g.novel==1]
        if len(base) and len(novel):
            a,b=linear_sum_assignment(cdist(xyz(base),xyz(novel)))
            result.append((base.index[a].to_numpy(),novel.index[b].to_numpy()))
    return result

def analysis(sp,native,intro,level):
    n=make_cells(native,True); c=make_cells(intro)
    ties=int((c.novel==.5).sum());c=c[c.novel!=.5].reset_index(drop=True)
    pairs=match_cells(c)
    result=[]
    subsets=[T,T[:2]]
    if level==.95:
        subsets.extend([[T[0],T[2]],T[1:],*[T[j:j+1] for j in range(3)]])
    for cols in subsets:
        x=n[cols].to_numpy();v=c[cols].to_numpy(); mean=x.mean(0);sd=x.std(0)
        x=(x-mean)/sd;v=(v-mean)/sd
        # Robust bounds, with extreme observations retained in tail bins.
        lim=np.quantile(np.vstack([x,v]),[.005,.995],axis=0)
        a=np.clip(x,lim[0],lim[1]);b=np.clip(v,lim[0],lim[1])
        for h in ([.35,.5,.75] if level==.95 else [.5]):
            edges=[np.linspace(lo-3*h,hi+3*h,41) for lo,hi in lim.T]
            def surface(pts):
                hist=np.histogramdd(pts,bins=edges)[0]
                s=gaussian_filter(hist,[h/(e[1]-e[0]) for e in edges],mode='constant').ravel()
                return s/s.sum()
            ref=surface(a)
            def D(u,v):return np.minimum(u,v).sum()
            sim=c.novel==0;new=c.novel==1
            u,v=surface(b[sim]),surface(b[new])
            base=dict(traits="+".join(["angle","chroma","outline"][T.index(v)] for v in cols),species=sp,native95_mass=level,n_traits=len(cols),bandwidth=h,
                      n_native=len(n),n_analogous=int(sim.sum()),n_novel=int(new.sum()),mixed_ties=ties,
                      photo_native=n.photos.sum(),photo_analogous=c[sim].photos.sum(),photo_novel=c[new].photos.sum())
            result.append(dict(**base,design='all_cells',D_native_analogous=D(ref,u),D_native_novel=D(ref,v),
                               D_analogous_novel=D(u,v),drop_in_overlap=D(ref,u)-D(ref,v)))
            if len(pairs)<2:continue
            total=sum(len(p) for p,q in pairs)
            left=np.array([surface(b[p])*len(p)/total for p,q in pairs])
            right=np.array([surface(b[q])*len(q)/total for p,q in pairs])
            pool=left.sum(0)+right.sum(0); u=left.sum(0);v=right.sum(0)
            effect=D(ref,u)-D(ref,v)
            row=dict(**base,design='geographic_pairs',pairs=total,blocks=len(pairs),
                     D_native_analogous=D(ref,u),D_native_novel=D(ref,v),D_analogous_novel=D(u,v),drop_in_overlap=effect)
            # Swap all matched cell labels together within each geographic block.
            if len(pairs)<=11:
                flips=np.array(list(product([0.,1.],repeat=len(pairs))));exact=True
            else:
                flips=np.random.default_rng(SEED).integers(0,2,size=(999,len(pairs))).astype(float);exact=False
            extreme=0
            for start in range(0,len(flips),16):
                uu=u+flips[start:start+16]@(right-left);vv=pool-uu
                delta=np.minimum(ref,uu).sum(axis=1)-np.minimum(ref,vv).sum(axis=1)
                extreme+=(np.abs(delta)>=abs(effect)-1e-12).sum()
            row['p_swap']=(extreme+int(not exact))/(len(flips)+int(not exact));row['n_maps']=len(flips)
            result.append(row)
    return result
