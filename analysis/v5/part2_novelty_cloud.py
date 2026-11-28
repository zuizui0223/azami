"""Aggregate-only novelty cloud for the frozen Part 2 support analysis.

No new model, test or threshold search. Same native splits and distances;
pool the 32 reference draws to visualize the exact joint category masses.
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from analysis.v5.part2_hypervolume import (
    ROOT, EXPECTED, VARIABLES, T, authority_source, angular_cells,
    background_mesh, apply_original_five_CHELSA, make_cells,
    native_splits, distances, joint_novelty)
OUT=ROOT/"part2_novelty_cloud"
# Uniform 0.25 log2 intervals; terminal bins include values beyond +/-4.
EDGES=np.linspace(-4,4,33)
EXPECTED_FRACTIONS={
 ("Cirsium arvense",3):[.835938,.028379,.130560,.005124],
 ("Cirsium vulgare",3):[.646236,.041731,.288246,.023787],
}
EXPECTED_RHO={
 ("Cirsium arvense",3):.062913,("Cirsium arvense",2):.077028,
 ("Cirsium vulgare",3):.149545,("Cirsium vulgare",2):.141558,
}

def make_joint_cells(f):
    c=make_cells(f,True)
    d=f.loc[f[T[0]].between(0,180)].dropna(subset=T).copy()
    d["i"]=np.minimum(179,np.floor(d.latitude+90)).astype(int)
    d["j"]=np.floor((d.longitude+180)%360).astype(int)
    return c.merge(d.groupby(["i","j"])[list(VARIABLES)].median().reset_index(),on=["i","j"],validate="one_to_one")

def bins(a,b):
    # Explicitly preserve <=1 versus >1 at the central boundary, even for ties.
    z=np.column_stack((a,b))
    assert np.isfinite(z).all() and (z>=0).all()
    log=np.log2(np.maximum(z,2.**-20))
    log=np.clip(log,EDGES[0],np.nextafter(EDGES[-1],-np.inf))
    idx=np.searchsorted(EDGES,log,side="right")-1
    mid=len(EDGES)//2
    idx=np.clip(idx,0,len(EDGES)-2)
    idx=np.where((z<=1)&(idx>=mid),mid-1,idx)
    idx=np.where((z>1)&(idx<mid),mid,idx)
    np.testing.assert_array_equal(idx>=mid,z>1)
    h=np.zeros((len(EDGES)-1,len(EDGES)-1),float)
    np.add.at(h,(idx[:,0],idx[:,1]),1/len(z))
    np.testing.assert_allclose(h.sum(),1)
    return h

def quadrants(h):
    m=h.shape[0]//2
    return np.array([h[:m,:m].sum(),h[:m,m:].sum(),h[m:,:m].sum(),h[m:,m:].sum()])

def run_species(sp,sn,si,rows,verification):
    n,c=make_joint_cells(sn),make_joint_cells(si)
    assert (len(n),len(c))==EXPECTED[sp][:2]
    splits=native_splits(n.block.to_numpy())
    x=n[list(VARIABLES)].to_numpy();xx=c[list(VARIABLES)].to_numpy()
    check={"joint_support":[],"joint_tests":[]}
    joint_novelty(sp,sn,si,check)
    for dim in (3,2):
        y=n[T[:dim]].to_numpy();yy=c[T[:dim]].to_numpy()
        grids={k:[] for k in ("heldout_native","introduced")}
        for tr,cal in splits:
            xc=distances(x[tr],x[cal]);yc=distances(y[tr],y[cal])
            qx,qy=np.quantile(xc,.95),np.quantile(yc,.95)
            assert qx>0 and qy>0
            xi=distances(x[tr],xx)/qx;yi=distances(y[tr],yy)/qy
            grids["heldout_native"].append(bins(xc/qx,yc/qy))
            grids["introduced"].append(bins(xi,yi))
        for k,values in grids.items():
            h=np.mean(values,axis=0)
            actual=quadrants(h)
            expected=np.array([r["fraction"] for r in check["joint_support"] if r["dimensions"]==dim and r["range"]==k])
            np.testing.assert_allclose(actual,expected,atol=1e-12)
            if k=="introduced" and (sp,dim) in EXPECTED_FRACTIONS:
                np.testing.assert_allclose(actual,EXPECTED_FRACTIONS[(sp,dim)],atol=6e-7,rtol=0)
            for i,j in np.ndindex(h.shape):
                rows.append(dict(species=sp,dimensions=dim,range=k,x_bin=i,y_bin=j,
                    climate_lo=EDGES[i],climate_hi=EDGES[i+1],
                    trait_lo=EDGES[j],trait_hi=EDGES[j+1],mass=h[i,j],
                    terminal_bin_clipped=(i in (0,h.shape[0]-1) or j in (0,h.shape[1]-1))))
            verification.append(dict(species=sp,dimensions=dim,range=k,
                mass_total=h.sum(),max_quadrant_error=np.max(np.abs(actual-expected)),
                n_draws=len(splits),n_cells=len(c) if k=="introduced" else 128))
        a=next(r for r in check["joint_tests"] if r["dimensions"]==dim)
        np.testing.assert_allclose(a["rho_within_region"],EXPECTED_RHO[(sp,dim)],atol=6e-7,rtol=0)

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    auth=authority_source((ROOT/"private/authority350.zip").read_bytes())
    wild=pd.read_csv(ROOT/"private/focal_wild.csv.gz",dtype={"obs_id":str},low_memory=False)
    assert len(wild)==297622
    wild["_cell"]=angular_cells(wild.latitude,wild.longitude)
    wild=wild.sort_values("obs_id").drop_duplicates(["taxon_name","range_state","_cell","user_id"]).copy()
    land=background_mesh(auth)
    wild,_,_=apply_original_five_CHELSA(wild,land,source_meta=Path("analysis/v4/environment_source_metadata.json"),work=ROOT/"private/rasters")
    measured=pd.read_csv(ROOT/"private/original-measured.csv.gz",usecols=["obs_id",*T],dtype={"obs_id":str})
    wild=wild.drop(columns=T,errors="ignore").merge(measured,on="obs_id",how="left",validate="one_to_one")
    rows=[];verification=[]
    for sp,f in wild.groupby("taxon_name"):
        f=f.dropna(subset=list(VARIABLES))
        run_species(sp,f.loc[f.range_state=="native"],f.loc[f.range_state=="introduced"],rows,verification)
    assert len(rows)==8192 and len(verification)==8
    pd.DataFrame(rows).to_csv(OUT/"novelty_density_grid.csv",index=False)
    pd.DataFrame(verification).to_csv(OUT/"verification.csv",index=False)
    print("Verified frozen joint support and rank correlations; aggregate grids only.",flush=True)

def self_test():
    a=np.array([0.,1.,1.000001,1000.])
    b=np.array([1.1,.1,.2,1000.])
    h=bins(a,b)
    np.testing.assert_allclose(quadrants(h),[.25,.25,.25,.25])
    print("Synthetic threshold, clipping and mass checks passed",flush=True)
if __name__=="__main__":
    self_test() if "--self-test" in sys.argv else main()
