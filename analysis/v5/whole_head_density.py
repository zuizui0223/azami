"""Full-range head density overlap and novelty strata on one 3D grid."""
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
from analysis.v5.novel_head_overlap import T,make_cells


def run(sp,native,introduced,labeled,out):
    n=make_cells(native,True)
    # All introduced cells including ties or cells lacking a novelty assignment.
    i=make_cells(introduced,True)
    c=make_cells(labeled)
    data=[n,i,c[c.novel==0],c[c.novel==1]]
    rows=[]
    export={}
    for cols in [T,T[:2]]:
        v=[g[cols].to_numpy() for g in data]
        center=v[0].mean(0);sd=v[0].std(0)
        z=[(m-center)/sd for m in v]
        bounds=np.quantile(np.vstack(z[:2]),[.005,.995],axis=0)
        for h in [.35,.5,.75]:
            edges=[np.linspace(lo-3*h,hi+3*h,41) for lo,hi in bounds.T]
            surfaces=[]
            for m in z:
                mass=np.histogramdd(np.clip(m,*bounds),bins=edges)[0]
                s=gaussian_filter(mass,[h/np.diff(e)[0] for e in edges],mode="constant")
                surfaces.append(s/s.sum())
            for j,group in enumerate(["introduced_all","introduced_analogue","introduced_novel"],1):
                D=np.minimum(surfaces[0],surfaces[j]).sum()
                rows.append(dict(species=sp,dimensions=len(cols),bandwidth=h,group=group,n_native=len(n),
                                 n_cells=len(data[j]),records_native=n.photos.sum(),records=data[j].photos.sum(),
                                 ties=sum(c.novel==.5),D=D))
            if len(cols)==3 and h==.5:
                for j,s in enumerate(surfaces):
                    export["density_"+str(j)]=s
                for j,e in enumerate(edges):
                    export["edges_"+str(j)]=e*sd[j]+center[j]
    np.savez_compressed(out/(sp.replace(" ","_")+"_head_densities.npz"),**export)
    path=out/"whole_head_overlap.csv"
    pd.DataFrame(rows).to_csv(path,index=False,mode="a",header=not path.exists())
