"""Relate the non-overlapping parts of two densities through their shared observation cells."""
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
from analysis.v4.native_cached_climate import VARIABLES


class Surface:
    def __init__(self, points, h, bins=64):
        edges=[np.linspace(v.min()-2.5,v.max()+2.5,bins+1) for v in points.T]
        ij=[np.searchsorted(e,v,side="right")-1 for e,v in zip(edges,points.T)]
        self.shape=(bins,)*points.shape[1]
        self.flat=np.ravel_multi_index(ij,self.shape)
        self.sigma=[h/np.diff(e)[0] for e in edges]
        self.kernel_mass=self.smooth(np.ones(self.shape)).ravel()[self.flat]

    def smooth(self,x):
        return gaussian_filter(x,self.sigma,mode="constant")

    def scores(self,native,w):
        """Back-project excess density to its source cells using the same kernel."""
        p=[]
        for subset in [native,~native]:
            x=np.bincount(self.flat[subset],weights=w[subset],minlength=np.prod(self.shape))
            d=self.smooth(x.reshape(self.shape));p.append(d/d.sum())
        nonoverlap=1-np.minimum(*p).sum()
        score=np.zeros(len(w))
        for subset,own,other in [(native,*p),(~native,p[1],p[0])]:
            share=np.divide(np.maximum(own-other,0),own,out=np.zeros_like(own),where=own>0)
            score[subset]=self.smooth(share).ravel()[self.flat[subset]]/self.kernel_mass[subset]
            np.testing.assert_allclose(np.average(score[subset],weights=w[subset]*self.kernel_mass[subset]),
                                       nonoverlap,atol=1e-9)
        return score,1-nonoverlap


def summary(x,y,raw,w):
    a=np.average(y,weights=w*x)
    b=np.average(y,weights=w*(1-x))
    shift=np.average(raw,axis=0,weights=w*x)-np.average(raw,axis=0,weights=w*(1-x))
    return [a,b,a-b,*shift]


def run(c,sp,traits,project):
    native=c.range_state.eq("native").to_numpy()
    w=1/c.groupby("range_state").range_state.transform("size")
    raw=c[traits].to_numpy()
    profiled=c[[*traits,*VARIABLES]].to_numpy()
    mean=np.average(raw,axis=0,weights=w)
    z=(raw-mean)/np.sqrt(np.average((raw-mean)**2,axis=0,weights=w))
    x,_=project(c,list(VARIABLES));y,_=project(c,traits)
    rng=np.random.default_rng(20261123)
    all_weights=[np.ones(len(c))]
    for rep in range(399):
        weights=np.zeros(len(c))
        for subset in [native,~native]:
            groups,idx=np.unique(c.block.to_numpy()[subset],return_inverse=True)
            weights[subset]=rng.multinomial(len(groups),np.full(len(groups),1/len(groups)))[idx]
        all_weights.append(weights)
    result=[]
    for h in [.35,.5,.75]:
        climate=Surface(x,h)
        for name,points in [("joint_PC12",y),("angle",z[:,[0]]),("chroma",z[:,[1]]),("outline",z[:,[2]]),("joint_3D",z)]:
            head=Surface(points,h,bins=40 if name=="joint_3D" else 64)
            rows=[[],[]]
            for k,weight in enumerate(all_weights if h==.5 and name!="joint_3D" else all_weights[:1]):
                cx,dc=climate.scores(native,weight)
                hy,dh=head.scores(native,weight)
                for j,subset in enumerate([native,~native]):
                    rows[j].append(summary(cx[subset],hy[subset],profiled[subset],weight[subset]))
                if k==0: orig=(dc,dh)
            for side,values in zip(["native","introduced"],rows):
                arr=np.array(values)
                row=dict(species=sp,side=side,trait=name,h=h,climate_D=orig[0],trait_D=orig[1],
                         cells=sum(native if side=="native" else ~native))
                for j,key in enumerate(["head_excess_in_climate_excess","head_excess_in_climate_overlap","difference",
                                         "angle_shift","chroma_shift","outline_shift",*[axis+"_shift" for axis in VARIABLES]]):
                    row[key]=arr[0,j]
                    if len(arr)>1:
                        row[key+"_lo"],row[key+"_hi"]=np.quantile(arr[1:,j],[.025,.975])
                result.append(row)
    return result
