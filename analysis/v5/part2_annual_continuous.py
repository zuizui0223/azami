"""Part 2 annual five-axis continuous trait associations.
No quintiles, native phenotype truncation, novelty categories, or seasonal
predictors in the direct climate model. Only aggregate outputs are exported.
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy.stats import t as student_t
from statsmodels.stats.multitest import multipletests
from analysis.v5.phenotype_environment_explanation import (
    ROOT, SPECIES, VARIABLES, NAMES, COLS, authority_source, angular_cells,
    background_mesh, apply_original_five_CHELSA, load_three_axis, eligible,
    CALIPER, MAX_AXIS, MIN_PAIRS, MIN_REGIONS,
)
from analysis.v5.phenotype_redistribution import weighted_quantile

OUT = ROOT / "part2_annual_continuous"
AXES = ("chelsa_bio12", "chelsa_bio01", "chelsa_bio04",
        "chelsa_bio15", "chelsa_rsds_mean")
assert set(AXES) == set(VARIABLES) and len(AXES) == 5
RANGES = ("native", "introduced")
GRIDS = (1, 5)
BW = .25

def cell_data(frame, degrees):
    # Coordinates are used for grouping only; no geographic ID is exported.
    f = eligible(frame)
    q = f[list(AXES)+list(COLS)].to_numpy(float)
    f = f.loc[np.isfinite(q).all(1)].copy()
    f["i"] = np.minimum(180//degrees-1, np.floor((f.latitude+90)/degrees)).astype(int)
    f["j"] = np.floor(((f.longitude+180)%360)/degrees).astype(int)
    g = f.groupby(["i","j"], sort=True)
    c = g[list(AXES)].median()
    c[list(COLS)] = g[list(COLS)].mean()
    c["count"] = g.size()
    c = c.reset_index()
    c["block"] = (c.i*degrees//20)*18+c.j*degrees//20
    index = pd.MultiIndex.from_frame(c[["i","j"]])
    ix = index.get_indexer(pd.MultiIndex.from_frame(f[["i","j"]]))
    assert (ix >= 0).all()
    assert np.array_equal(np.bincount(ix,minlength=len(c)),c["count"].to_numpy())
    return c, f[list(COLS)].to_numpy(float), ix

def support(a,b):
    # Symmetric nearest-neighbour support in all five continuous axes.
    # Not a climate category or a phenotype filter. Fixed, previously used caliper.
    aa=np.zeros(len(a),bool); bb=np.zeros(len(b),bool)
    for i in range(0,len(a),128):
        d=a[i:i+128,None,:]-b[None,:,:]
        valid=(np.sqrt(np.mean(d*d,axis=2))<=CALIPER+1e-12)&(np.max(abs(d),axis=2)<=MAX_AXIS+1e-12)
        aa[i:i+128]=valid.any(1)
        bb|=valid.any(0)
    return aa,bb

def regression(x,y,g):
    # Cell-equal estimates; absorb distinct region intercepts within each range.
    ids=np.unique(g)
    if len(x)<MIN_PAIRS or len(ids)<MIN_REGIONS:
        return None
    x=x.copy(); y=y.copy()
    for block in ids:
        k=g==block
        x[k]-=x[k].mean(0); y[k]-=y[k].mean(0)
    if np.linalg.matrix_rank(x)!=5 or len(x)<=len(ids)+5:
        return None
    beta=np.linalg.lstsq(x,y,rcond=None)[0]
    u=y-x@beta
    bread=np.linalg.inv(x.T@x)
    correction=np.sqrt(len(ids)/(len(ids)-1)*(len(x)-1)/(len(x)-len(ids)-5))
    influence={int(k):correction*bread@(x[g==k].T@u[g==k]) for k in ids}
    return beta,influence,float(np.linalg.cond(x))

def contrast(a,b):
    # Covariance between ranges is retained if they share geographic blocks.
    if a is None or b is None: return None
    x,ia,ca=a; y,ib,cb=b
    keys=set(ia)|set(ib)
    influence={k:ib.get(k,np.zeros_like(y))-ia.get(k,np.zeros_like(x)) for k in keys}
    return y-x,influence,max(ca,cb)

def infer(fit,j,k):
    if fit is None: return dict(status="limited_support_or_rank")
    beta,inf,condition=fit
    se=np.sqrt(sum(v[j,k]**2 for v in inf.values()))
    df=len(inf)-1
    est=float(beta[j,k])
    if se<=0 or df<=0:
        return dict(status="degenerate_variance",beta=est)
    crit=student_t.ppf(.975,df)
    return dict(status="estimated",beta=est,se=se,lower=est-crit*se,
                upper=est+crit*se,p=2*student_t.sf(abs(est/se),df),
                df=df,condition=condition)

def profile(c,v,ix,x,center,scale,meta):
    # Unadjusted local trait distributions along each physical climate axis;
    # every cell contributes equally before the continuous kernel weights.
    out=[]
    counts=c["count"].to_numpy()
    blocks=c.block.to_numpy()
    for j,axis in enumerate(AXES):
        for value in np.linspace(-2,2,81):
            w=np.exp(-.5*((x[:,j]-value)/BW)**2)
            w[abs(x[:,j]-value)>2*BW]=0
            if w.sum()<=0: continue
            ess=w.sum()**2/(w@w)
            gw=np.array([w[blocks==g].sum() for g in np.unique(blocks)])
            gess=gw.sum()**2/(gw@gw)
            rw=w[ix]/counts[ix]
            take=rw>0
            rw=rw[take]/rw[take].sum()
            for k,trait in enumerate(NAMES):
                z=v[take,k]
                m=float(rw@z)
                q=weighted_quantile(z,rw,[.1,.5,.9])
                out.append(dict(**meta,axis=axis,trait=trait,climate_value=center[j]+scale[j]*value,
                    native_climate_z=value,mean=m,sd=np.sqrt(max(float(rw@(z*z))-m*m,0)),
                    q10=q[0],q50=q[1],q90=q[2],effective_cells=ess,effective_regions=gess,
                    status="descriptive" if ess>=20 and gess>=3 else "sparse"))
    return out

def analyze(sp,frames,degrees,out):
    data=[cell_data(frames.loc[frames.range_state==r],degrees) for r in RANGES]
    c0,v0,ix0=data[0]
    mu=c0[list(AXES)].to_numpy().mean(0)
    sd=c0[list(AXES)].to_numpy().std(0)
    assert np.isfinite(sd).all() and (sd>0).all()
    # Common native scale for both ranges. All measured phenotypes retained.
    w=1/c0["count"].to_numpy()[ix0]; w/=w.sum()
    ym=w@v0; ys=np.sqrt(w@((v0-ym)**2))
    assert np.isfinite(ys).all() and (ys>0).all()
    x=[(c[list(AXES)].to_numpy()-mu)/sd for c,_,_ in data]
    shared=support(*x)
    for j,axis in enumerate(AXES):
        out["scales"].append(dict(species=sp,degrees=degrees,axis=axis,center=mu[j],sd=sd[j]))
    for k,trait in enumerate(NAMES):
        out["trait_scales"].append(dict(species=sp,degrees=degrees,trait=trait,center=ym[k],sd=ys[k]))
    for spec,masks in (("all_cells",[np.ones(len(a),bool) for a in x]),("shared_five_axis",shared)):
        fits=[]
        for i,(c,v,ix) in enumerate(data):
            m=masks[i]; g=c.block.to_numpy()
            y=(c[list(COLS)].to_numpy()-ym)/ys
            meta=dict(species=sp,degrees=degrees,spec=spec,range=RANGES[i])
            out["cohort"].append(dict(**meta,cells=int(m.sum()),total_cells=len(c),
                records=int(c.loc[m,"count"].sum()),regions=len(np.unique(g[m]))))
            fits.append(regression(x[i][m],y[m],g[m]))
        fits.append(contrast(*fits))
        for label,fit in zip((*RANGES,"introduced_minus_native"),fits):
            for j,axis in enumerate(AXES):
                for k,trait in enumerate(NAMES):
                    out["coefficients"].append(dict(species=sp,degrees=degrees,spec=spec,
                        contrast=label,axis=axis,trait=trait,**infer(fit,j,k)))
    if degrees==1:
        for i,(c,v,ix) in enumerate(data):
            out["profiles"].extend(profile(c,v,ix,x[i],mu,sd,
                dict(species=sp,degrees=degrees,range=RANGES[i])))

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    auth=authority_source((ROOT/"private/authority350.zip").read_bytes())
    wild=pd.read_csv(ROOT/"private/focal_wild.csv.gz",dtype={"obs_id":str},low_memory=False)
    assert len(wild)==297622
    wild["_cell"]=angular_cells(wild.latitude,wild.longitude)
    wild=wild.sort_values("obs_id").drop_duplicates(["taxon_name","range_state","_cell","user_id"]).copy()
    land=background_mesh(auth)
    wild,_,_=apply_original_five_CHELSA(wild,land,source_meta=Path("analysis/v4/environment_source_metadata.json"),work=ROOT/"private/rasters")
    measured=load_three_axis(ROOT/"private/original-measured.csv.gz")
    wild=wild.drop(columns=[c for c in measured if c!="obs_id"],errors="ignore").merge(measured,on="obs_id",how="left",validate="one_to_one")
    assert set(wild.taxon_name.unique())==set(SPECIES)
    out={k:[] for k in ("scales","trait_scales","cohort","coefficients","profiles")}
    for sp,frame in wild.groupby("taxon_name"):
        frame=frame.dropna(subset=list(AXES))
        for degrees in GRIDS: analyze(sp,frame,degrees,out)
    d=pd.DataFrame(out["coefficients"])
    assert len(d)==360
    # Separate 30 planned comparisons per resolution, support, and contrast.
    # Ineligible estimates keep their status; p=1 only for correction family size.
    d["holm30"]=d.groupby(["degrees","spec","contrast"]).p.transform(
        lambda p:multipletests(p.fillna(1),method="holm")[1])
    d.loc[d.status!="estimated","holm30"]=np.nan
    d.to_csv(OUT/"coefficients.csv",index=False)
    for k in ("scales","trait_scales","cohort","profiles"):
        pd.DataFrame(out[k]).to_csv(OUT/(k+".csv"),index=False)
    print({"cohorts":len(out["cohort"]),"coefficients":len(d),"estimated":int((d.status=="estimated").sum())},flush=True)

def self_test():
    rng=np.random.default_rng(42)
    n=120; g=np.repeat(np.arange(12),10)
    x=rng.normal(size=(n,5))
    beta=rng.normal(size=(5,3))
    y=x@beta+g[:,None]*2
    a=regression(x,y,g); b=regression(x,y+2*x@beta,g)
    np.testing.assert_allclose(a[0],beta,atol=1e-10)
    np.testing.assert_allclose(contrast(a,b)[0],2*beta,atol=1e-10)
    assert support(x,x)[0].all()
    assert not support(x,x+100)[0].any()
    assert regression(x[:5],y[:5],g[:5]) is None
    print("synthetic tests passed",flush=True)

if __name__=="__main__":
    self_test() if "--self-test" in sys.argv else main()
