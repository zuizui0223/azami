"""Do native/introduced regional pairs with more similar climates have more similar heads?"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
from scipy.spatial.distance import cdist
from scipy.stats import rankdata
from statsmodels.stats.multitest import multipletests
from analysis.v4.native_introduced_opportunity import angular_cells
from analysis.v4.native_cached_climate import VARIABLES

TRAITS=['orientation_angle_degrees_median','corolla_lab_chroma_median','shape_aspect_ratio_median']
ROOT=Path('work/submission')
OUT=ROOT/'paired_overlap'
SEED=20261109

def project(cells,cols):
    a=cells[cols].to_numpy()
    w=(1/cells.groupby('range_state').range_state.transform('size')).to_numpy()
    m=np.average(a,axis=0,weights=w)
    sd=np.sqrt(np.average((a-m)**2,axis=0,weights=w))
    z=(a-m)/sd
    val,vec=np.linalg.eigh((z.T*w)@z/w.sum())
    order=np.argsort(val)[::-1];val=val[order];vec=vec[:,order]
    return z@vec[:,:2]/np.sqrt(val[:2]),val[:2].sum()/val.sum()

def grid(points):
    return [np.linspace(points[:,j].min()-2.5,points[:,j].max()+2.5,65) for j in [0,1]]

def density(points,edges,h):
    hist=np.histogram2d(*points.T,bins=edges)[0]
    d=gaussian_filter(hist,[h/(v[1]-v[0]) for v in edges],mode='constant').ravel()
    return d/d.sum()

def overlap(a,b):
    return 1-cdist(a,b,metric='cityblock')/2

def infer(a,b):
    ra=rankdata(a).reshape(a.shape);rb=rankdata(b).reshape(b.shape)
    ra=(ra-ra.mean())/np.linalg.norm(ra-ra.mean());rb=(rb-rb.mean())/np.linalg.norm(rb-rb.mean())
    rho=(ra*rb).sum()
    rng=np.random.default_rng(SEED+9); null=[]
    for _ in range(4999):
        null.append((ra*rb[rng.permutation(len(rb))][:,rng.permutation(rb.shape[1])]).sum())
    return rho,(1+np.count_nonzero(np.abs(null)>=abs(rho)-1e-12))/5000

def calculate(c,sp,traits,h,summary,pairs,full):
    x,xvar=project(c,list(VARIABLES)); y,yvar=project(c,traits)
    ex,ey=grid(x),grid(y); g=c.range_state.to_numpy()
    # Whole-range overlap uses every occupied one-degree cell once.
    full.append(dict(species=sp,shape=len(traits)==3,h=h,climate_D=overlap([density(x[g=='native'],ex,h)],[density(x[g=='introduced'],ex,h)]).item(),
                     head_D=overlap([density(y[g=='native'],ey,h)],[density(y[g=='introduced'],ey,h)]).item(),n_native=(g=='native').sum(),n_intro=(g=='introduced').sum(),
                     n_photos=c.photos.sum(),climate_PC12=xvar,head_PC12=yvar))
    rs=[[(int(k),v.index.to_numpy()) for k,v in c[c.range_state==a].groupby('block') if len(v)>=8] for a in ['native','introduced']]
    n,m=map(len,rs)
    if min(n,m)<3:raise ValueError('Too few independently permutable regions')
    dx=np.zeros((n,m));dy=dx.copy();rng=np.random.default_rng(SEED)
    for _ in range(32):
        # Identical cells for climate and head within each draw and sensitivity model.
        subset=[[rng.choice(ids,8,replace=False) for _,ids in group] for group in rs]
        for z,edges,acc in [(x,ex,dx),(y,ey,dy)]:
            a,b=[[density(z[idx],edges,h) for idx in grp] for grp in subset]
            acc+=overlap(a,b)/32
    rho,p=infer(dx,dy)
    # Remove one whole native or introduced block, never one pair only.
    leave=[]
    for axis,count in enumerate([n,m]):
        for i in range(count):
            u=np.delete(dx,i,axis).ravel();v=np.delete(dy,i,axis).ravel()
            leave.append(np.corrcoef(rankdata(u),rankdata(v))[0,1])
    summary.append(dict(species=sp,shape=len(traits)==3,h=h,n_native_blocks=n,n_intro_blocks=m,n_pairs=n*m,rho=rho,p=p,
                        leave_min=min(leave),leave_max=max(leave)))
    if h==.5:
        for i in range(n):
            for j in range(m):pairs.append(dict(species=sp,shape=len(traits)==3,native_block=i+1,intro_block=j+1,climate_D=dx[i,j],head_D=dy[i,j]))

def main():
    p=ROOT/'private'
    src=pd.read_csv(p/'focal_wild_five.csv.gz',dtype={'obs_id':str},low_memory=False)
    assert len(src)==297622
    meas=pd.read_csv(p/'original-measured.csv.gz',usecols=['obs_id',*TRAITS],dtype={'obs_id':str})
    src=src.drop(columns=TRAITS,errors='ignore').merge(meas,on='obs_id',how='left',validate='one_to_one')
    src=src[src[TRAITS[0]].between(0,180)].dropna(subset=[*TRAITS,*VARIABLES]).copy()
    src['coarse']=angular_cells(src.latitude,src.longitude)
    src=src.sort_values('obs_id').drop_duplicates(['taxon_name','range_state','coarse','user_id']).copy()
    src['i']=(src.latitude+90).clip(upper=179.99999).astype(int)
    src['j']=((src.longitude+180)%360).astype(int)
    src['block']=(src.i//20)*18+src.j//20
    from analysis.v5 import trait_coupling
    from analysis.v5 import native_region_baseline
    baselines=[];region_values=[]
    components=[];component_pairs=[]
    summary=[];pairs=[];full=[]
    for sp,s in src.groupby('taxon_name'):
        grouped=s.groupby(['range_state','block','i','j'])
        c=grouped[[*VARIABLES,*TRAITS]].median().reset_index()
        c['photos']=grouped.size().to_numpy()
        base,values=native_region_baseline.run(c,sp,TRAITS,project,grid,density,SEED)
        baselines.extend(base);region_values.extend(values)
        rows,points=trait_coupling.analyze(c,sp,TRAITS,project,grid,density,overlap,SEED)
        components.extend(rows);component_pairs.extend(points)
        for t in [TRAITS,TRAITS[:2]]:
            for h in [.35,.5,.75]:
                calculate(c,sp,t,h,summary,pairs,full)
    tab=pd.DataFrame(summary)
    primary=(tab.h==.5)&tab['shape']
    tab.loc[primary,'holm_primary']=multipletests(tab.loc[primary,'p'],method='holm')[1]
    OUT.mkdir(exist_ok=True)
    pd.DataFrame(baselines).to_csv(OUT/"native_baseline.csv",index=False)
    pd.DataFrame(region_values).to_csv(OUT/"native_baseline_regions.csv",index=False)
    trait_coupling.save(components,component_pairs,OUT)
    for file,df in [('regional',tab),('pairs',pd.DataFrame(pairs)),('whole',pd.DataFrame(full))]:
        df.to_csv(OUT/(file+'.csv'),index=False)
if __name__=='__main__': main()
