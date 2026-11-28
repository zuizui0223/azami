"""Pre-specified geometric occupancy and joint novelty diagnostics.
Only aggregate outputs. No fitness, ancestral or causal inference."""
from pathlib import Path
from itertools import product
import sys,json
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
from scipy.stats import rankdata
from statsmodels.stats.multitest import multipletests
from analysis.v5.phenotype_three_axis_dispersion import (
    ROOT, EXPECTED, VARIABLES, authority_source, angular_cells, background_mesh,
    apply_original_five_CHELSA, load_three_axis, eligible, make_cells, match_cells,
    head_cells, _background_species_audit, _pca_env_2d, _probability_content, T)
from analysis.v5.native_pca_env import _weights_per_5degree_cell
from analysis.v5.phenotype_native_calibration import native_splits, distances

OUT=ROOT/"part2_hypervolume"
SEED=20261128
LEVELS=(.90,.95,.99)
H=(.35,.5,.75)

def hdr(a,q):
    a=np.asarray(a,float)
    assert np.isfinite(a).all() and np.min(a)>=0 and np.isclose(a.sum(),1)
    s=np.sort(a.ravel())[::-1]
    return a>=s[np.searchsorted(np.cumsum(s),q)]

def partition(a,b,voxel,q=.95,ref=None):
    A,B=hdr(a,q),hdr(b,q)
    i=int((A&B).sum());u=int((A&~B).sum());e=int((B&~A).sum())
    va,vb=i+u,i+e
    D=(u+e)/(va+vb)
    turnover=2*min(u,e)/(va+vb);volume=abs(e-u)/(va+vb)
    assert abs(turnover+volume-D)<1e-12
    d=dict(reference_volume=va*voxel,comparison_volume=vb*voxel,
        shared_volume=i*voxel,reference_only_volume=u*voxel,comparison_only_volume=e*voxel,
        log_volume_ratio=np.log(vb/va),volume_ratio=vb/va,
        jaccard=i/(i+u+e),sorensen=2*i/(va+vb),turnover=turnover,
        net_volume_difference=volume,dissimilarity=D,
        comparison_unique_fraction=e/vb,reference_unique_fraction=u/va,
        comparison_mass_outside_reference=float(b[~A].sum()),
        reference_mass_outside_comparison=float(a[~B].sum()),
        density_overlap=float(np.minimum(a,b).sum()),
        reference_hdr_mass=float(a[A].sum()),comparison_hdr_mass=float(b[B].sum()))
    if ref is not None:
        mask=hdr(ref,q)
        d["novelty_mass_difference"]=float(b[~mask].sum()-a[~mask].sum())
        d["analogous_mass_outside_native"]=float(a[~mask].sum())
        d["novel_mass_outside_native"]=float(b[~mask].sum())
    return d

def edges_for(values,bins):
    lo=values.min(0)-3*max(H);hi=values.max(0)+3*max(H)
    return [np.linspace(a,b,bins+1) for a,b in zip(lo,hi)]

def hist_smooth(v,edges,h,weights=None):
    hist=np.histogramdd(v,bins=edges,weights=weights)[0]
    expect=len(v) if weights is None else weights.sum()
    np.testing.assert_allclose(hist.sum(),expect,rtol=1e-8,atol=1e-10)
    sm=gaussian_filter(hist,[h/(e[1]-e[0]) for e in edges],mode="constant")
    assert sm.sum()>0
    return sm/sm.sum()

def climate5(sp,en,ei,sn,si,out):
    # Same source records, occurrence weighting and area-weighted available backgrounds
    # as PCA-env. All five standardized physical axes retained; no PCA truncation.
    ew=[x._area_weight.to_numpy(float) for x in (en,ei)]
    bg=[x[list(VARIABLES)].to_numpy(float) for x in (en,ei)]
    ob=[x[list(VARIABLES)].to_numpy(float) for x in (sn,si)]
    ow=[_weights_per_5degree_cell(x) for x in (sn,si)]
    w=np.r_[ew[0]/ew[0].sum()/2,ew[1]/ew[1].sum()/2]
    x=np.vstack(bg);mu=w@x;sd=np.sqrt(w@((x-mu)**2))
    assert np.all(sd>0)
    bg=[(v-mu)/sd for v in bg];ob=[(v-mu)/sd for v in ob]
    allv=np.vstack(bg+ob)
    for bins,h in [(20,.5),(20,.35),(20,.75),(16,.5),(24,.5)]:
        edges=edges_for(allv,bins)
        dx=np.array([e[1]-e[0] for e in edges]);voxel=float(dx.prod())
        raw=[];corrected=[];supported=[]
        for j in range(2):
            be=hist_smooth(bg[j],edges,h,ew[j]);oc=hist_smooth(ob[j],edges,h,ow[j])
            floor=max(be.max()*1e-4,1e-12)
            supported.append(float(oc[be>=floor].sum()))
            z=np.divide(oc,be,out=np.zeros_like(oc),where=be>=floor)
            assert z.sum()>0
            raw.append(oc);corrected.append(z/z.sum())
            del be,z
        for q in LEVELS:
            for typ,arr in [("background_corrected",corrected),("raw",raw)]:
                d=partition(*arr,voxel,q)
                out["climate5"].append(dict(species=sp,bins=bins,bandwidth=h,coverage=q,
                    density=typ,max_grid_step_over_bandwidth=float(dx.max()/h),
                    resolution_status="coarse" if dx.max()/h>1 else "resolved",
                    native_supported_mass=supported[0],introduced_supported_mass=supported[1],**d))
        del raw,corrected

def flip_patterns(n):
    if n<=10:return np.array(list(product([0.,1.],repeat=n))),True
    return np.vstack([np.zeros(n),np.random.default_rng(SEED).integers(0,2,size=(999,n))]),False

def paired_hv(sp,definition,n,c,pairs,out):
    # Cell-median primary cohort exactly as the existing 3D clouds. Same number
    # of cells per arm and same geographic blocks. Full tails retained.
    a=n[T].to_numpy(float);v=c[T].to_numpy(float)
    mu=a.mean(0);sd=a.std(0);assert np.all(sd>0)
    a=(a-mu)/sd;v=(v-mu)/sd
    counts=np.array([len(l) for l,r in pairs]);total=counts.sum()
    assert len(pairs)>=8 and total>=25
    left=np.concatenate([l for l,r in pairs]);right=np.concatenate([r for l,r in pairs])
    assert len(np.unique(np.r_[left,right]))==2*total
    for dim in (3,2):
        x=a[:,:dim];z=v[:,:dim]
        edges=edges_for(np.vstack([x,z]),40)
        voxel=float(np.prod([e[1]-e[0] for e in edges]))
        for h in H:
            ref=hist_smooth(x,edges,h)
            aa=hist_smooth(z[left],edges,h);bb=hist_smooth(z[right],edges,h)
            for q in LEVELS:
                out["matched_sensitivity"].append(dict(species=sp,definition=definition,
                    dimensions=dim,bandwidth=h,coverage=q,pairs=int(total),blocks=len(counts),
                    **partition(aa,bb,voxel,q,ref)))
            if h!=.5 or dim!=3:continue
            # Fixed primary 3D h=.5 q=.95. Mixture of individually smoothed
            # regional densities; equal cell weights, exact block label reversal.
            # Normalize the pooled mixtures (same small edge-loss rule per region).
            l=np.stack([hist_smooth(z[p],edges,h).ravel()*len(p)/total for p,r in pairs])
            r=np.stack([hist_smooth(z[q],edges,h).ravel()*len(q)/total for p,q in pairs])
            aa=l.sum(0);bb=r.sum(0);pool=aa+bb
            ref=ref.ravel()
            refmask=hdr(ref,.95)
            def metrics(u,v,shift):
                A,B=hdr(u,.95),hdr(v,.95)
                i=(A&B).sum();m=(A&~B).sum();e=(B&~A).sum()
                return np.array([np.log(B.sum()/A.sum()),2*min(m,e)/(A.sum()+B.sum()),
                    np.linalg.norm(shift),v[~refmask].sum()-u[~refmask].sum()])
            g=np.stack([(z[q]-z[p]).sum(0)/total for p,q in pairs])
            observed=metrics(aa,bb,g.sum(0))
            flips,exact=flip_patterns(len(pairs))
            null=[]
            for flip in flips:
                u=aa+flip@(r-l);w=pool-u
                u=np.maximum(u,0);w=np.maximum(w,0);u/=u.sum();w/=w.sum()
                null.append(metrics(u,w,(1-2*flip)@g))
            null=np.stack(null)
            np.testing.assert_allclose(null[0],observed,atol=1e-9)
            for j,name in enumerate(("log_volume_ratio","turnover","centroid_distance","novelty_mass_difference")):
                # For symmetric dissimilarity metrics the upper tail is relevant;
                # signed contrasts use the absolute-value two-sided statistic.
                ob=observed[j];v0=null[:,j]
                p=float(np.mean((v0 if j in (1,2) else abs(v0)) >= (ob if j in (1,2) else abs(ob))-1e-12))
                out["matched_tests"].append(dict(species=sp,definition=definition,metric=name,
                    observed=float(ob),p=p,exact=exact,permutations=len(flips),pairs=int(total),blocks=len(counts),
                    null_median=float(np.median(v0)),null_q025=float(np.quantile(v0,.025)),
                    null_q975=float(np.quantile(v0,.975))))

def joint_novelty(sp,sn,si,out):
    # Independently calibrated multivariate nearest-neighbour support; not
    # background-corrected occupancy, not native phenotype-on-climate prediction.
    def cells(f):
        c=make_cells(f,True)
        d=f.loc[f[T[0]].between(0,180)].dropna(subset=T).copy()
        d["i"]=np.minimum(179,np.floor(d.latitude+90)).astype(int)
        d["j"]=np.floor((d.longitude+180)%360).astype(int)
        c=c.merge(d.groupby(["i","j"])[list(VARIABLES)].median().reset_index(),on=["i","j"],validate="one_to_one")
        return c
    n,c=cells(sn),cells(si)
    assert (len(n),len(c))==EXPECTED[sp][:2]
    split=native_splits(n.block.to_numpy())
    x=n[list(VARIABLES)].to_numpy();xx=c[list(VARIABLES)].to_numpy()
    for dim in (3,2):
        y=n[T[:dim]].to_numpy();yy=c[T[:dim]].to_numpy()
        xs=[];ys=[];joint=[];native_joint=[]
        for tr,cal in split:
            xc=distances(x[tr],x[cal]);yc=distances(y[tr],y[cal])
            qx,qy=np.quantile(xc,.95),np.quantile(yc,.95)
            assert qx>0 and qy>0
            xi=distances(x[tr],xx)/qx;yi=distances(y[tr],yy)/qy
            xs.append(xi);ys.append(yi)
            joint.append(np.bincount((xi>1).astype(int)*2+(yi>1),minlength=4)/len(xi))
            native_joint.append(np.bincount((xc>qx).astype(int)*2+(yc>qy),minlength=4)/len(xc))
        a=np.mean(xs,axis=0);b=np.mean(ys,axis=0);g=c.block.to_numpy()
        for label,z in [("introduced",joint),("heldout_native",native_joint)]:
            for k,value in enumerate(np.mean(z,axis=0)):
                out["joint_support"].append(dict(species=sp,dimensions=dim,range=label,
                    climate_novel=bool(k&2),trait_novel=bool(k&1),fraction=float(value),
                    reference_cells=256,calibration_cells=128,draws=len(split),introduced_cells=len(c)))
        # Region-centred rank covariance; region sign exchange under score symmetry.
        ar=rankdata(a)/len(a);br=rankdata(b)/len(b)
        ss=[];sx=[];sy=[];ids=[]
        for k in np.unique(g):
            mask=g==k
            if mask.sum()<3:continue
            u=ar[mask]-ar[mask].mean();v=br[mask]-br[mask].mean()
            if u@u<=0 or v@v<=0:continue
            ss.append(float(u@v));sx.append(float(u@u));sy.append(float(v@v));ids.append(k)
        den=np.sqrt(sum(sx)*sum(sy))
        if len(ss)<8 or den<=0:
            out["joint_tests"].append(dict(species=sp,dimensions=dim,status="insufficient_support",blocks=len(ss),p=np.nan));continue
        s=np.array(ss);obs=s.sum()/den
        flips,exact=flip_patterns(len(s))
        null=(1-2*flips)@s/den
        p=np.mean(abs(null)>=abs(obs)-1e-12)
        loo=[(s.sum()-s[i])/np.sqrt((sum(sx)-sx[i])*(sum(sy)-sy[i])) for i in range(len(s))]
        out["joint_tests"].append(dict(species=sp,dimensions=dim,status="estimated",blocks=len(ss),
            cells=int(np.isin(g,ids).sum()),rho_within_region=obs,p=p,exact=exact,permutations=len(flips),
            loo_min=min(loo),loo_max=max(loo),draws=len(split)))

def run_species(sp,sn,si,transform,edges,nat,inv,out):
    n=make_cells(sn,True);full=make_cells(si,True);ex=EXPECTED[sp]
    assert (len(n),len(full))==ex[:2]
    coords=transform(si)
    uv=np.column_stack([np.searchsorted(e,coords[:,k],side="right")-1 for k,e in enumerate(edges)])
    for k,e in enumerate(edges):uv[coords[:,k]==e[-1],k]=len(e)-2
    inside=np.all((uv>=0)&(uv<np.array([len(e)-1 for e in edges])),axis=1)
    lab=si.loc[inside].copy();lab["novel"]=(~_probability_content(nat,.95)[uv[inside,0],uv[inside,1]]).astype(int)
    c=make_cells(lab);c=c.loc[c.novel!=.5].reset_index(drop=True)
    pairs=match_cells(c);assert sum(len(l) for l,r in pairs)==ex[2]
    paired_hv(sp,"native_core95",n,c,pairs,out)
    hc=head_cells(si,transform,edges)[["i","j","u","v"]]
    c=full.merge(hc,on=["i","j"],validate="one_to_one")
    u,v=c.u.to_numpy(),c.v.to_numpy()
    c=c.loc[inv[u,v]>0].reset_index(drop=True);u,v=c.u.to_numpy(),c.v.to_numpy()
    c["novel"]=(inv[u,v]>nat[u,v]).astype(int)
    pairs=match_cells(c);assert sum(len(l) for l,r in pairs)==ex[3]
    paired_hv(sp,"density_surplus",n,c,pairs,out)
    joint_novelty(sp,sn,si,out)

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    auth=authority_source((ROOT/"private/authority350.zip").read_bytes())
    wild=pd.read_csv(ROOT/"private/focal_wild.csv.gz",dtype={"obs_id":str},low_memory=False)
    assert len(wild)==297622
    wild["_cell"]=angular_cells(wild.latitude,wild.longitude)
    wild=wild.sort_values("obs_id").drop_duplicates(["taxon_name","range_state","_cell","user_id"]).copy()
    land=background_mesh(auth)
    wild,land,_=apply_original_five_CHELSA(wild,land,source_meta=Path("analysis/v4/environment_source_metadata.json"),work=ROOT/"private/rasters")
    # Match the established cell-median cloud definitions.
    measured=pd.read_csv(ROOT/"private/original-measured.csv.gz",usecols=["obs_id",*T],dtype={"obs_id":str})
    wild=wild.drop(columns=T,errors="ignore").merge(measured,on="obs_id",how="left",validate="one_to_one")
    frozen=json.loads(Path("analysis/v4/results/20261021_unified_abcde/S2_TWO_SPECIES_FULL_WCVP_NICHE_VERIFIED_20261021.json").read_text())
    out={k:[] for k in ("climate5","matched_sensitivity","matched_tests","joint_support","joint_tests")}
    for sp,rows in wild.groupby("taxon_name"):
        env,_=_background_species_audit(auth,land.dropna(subset=list(VARIABLES)),sp,land)
        complete=rows.dropna(subset=list(VARIABLES))
        sn=complete.loc[complete.range_state=="native"];si=complete.loc[complete.range_state=="introduced"]
        en=env[env._range=="native"];ei=env[env._range=="introduced"]
        transform,edges,nat,inv=_pca_env_2d(en,ei,sn,si,return_grid=True)
        assert abs(np.minimum(nat,inv).sum()-frozen["species"][sp]["PCA_env_background_corrected_Schoener_D"])<1e-5
        climate5(sp,en,ei,sn,si,out)
        run_species(sp,sn,si,transform,edges,nat,inv,out)
        print({"species":sp,"aggregate_counts":{k:len(v) for k,v in out.items()}},flush=True)
    for k,rows in out.items():
        d=pd.DataFrame(rows)
        if k=="matched_tests":
            assert len(d)==16
            d["holm16"]=multipletests(d.p,method="holm")[1]
        if k=="joint_tests":
            assert len(d)==4
            d["holm4"]=multipletests(d.p.fillna(1),method="holm")[1]
            d.loc[d.status!="estimated","holm4"]=np.nan
        d.to_csv(OUT/(k+".csv"),index=False)
    print("Complete: aggregate results only",flush=True)

def self_test():
    a=np.zeros((5,5));a[1:3,1:3]=.25
    b=np.zeros((5,5));b[2:4,1:3]=.25
    d=partition(a,b,1)
    assert d["volume_ratio"]==1 and d["turnover"]==.5 and d["net_volume_difference"]==0
    d=partition(a,a,1,ref=a)
    assert d["jaccard"]==1 and d["novelty_mass_difference"]==0
    b=np.zeros((5,5));b[1,1]=1.
    d=partition(a,b,1)
    assert d["volume_ratio"]==.25 and d["turnover"]==0
    rng=np.random.default_rng(42);x=rng.normal(size=(60,3))
    e=edges_for(x,20);z=hist_smooth(x,e,.5)
    assert np.isclose(z.sum(),1)
    assert flip_patterns(9)[0].shape==(512,9)
    assert flip_patterns(14)[0].shape==(1000,14)
    print("Synthetic partition, normalization, permutation checks passed",flush=True)

if __name__=="__main__":
    self_test() if "--self-test" in sys.argv else main()
