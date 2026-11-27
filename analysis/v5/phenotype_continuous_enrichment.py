"""Unbinned re-estimation of the existing three-trait enrichment curves.

Same paired cells, same underlying records, native-range inclusion and
pair intercepts. Replace quintile rank by measured value (native-SD scale,
centered on native median). Quadratic form fixed; do not select curvature
or a bandwidth from the results. Descriptive smooth marginal density and
phenotype-conditioned water profiles are not conditional coefficients.
No fitness, selection-gradient or causal interpretation.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.stats import t as student_t
from statsmodels.stats.multitest import multipletests
from analysis.v5.phenotype_curve_mechanism import (
    ROOT,SPECIES,DEFINITIONS,VARIABLES,RAIN,SEASON,NAMES,COLS,EXPECTED,
    EXPECTED_BETA,joint_records,eligible,load_three_axis,make_cells,match_cells,
    head_cells,background_mesh,_background_species_audit,_pca_env_2d,
    _probability_content,apply_original_five_CHELSA,authority_source,
    angular_cells,X,K,native_cuts,rank_codes,cell_hist,fit,climate_cells,
    reference_scale,
)
from analysis.v5.phenotype_redistribution import weighted_quantile
from analysis.v5.phenotype_within_species_enrichment import sign_patterns,group_sums
OUT=ROOT/"phenotype_continuous_enrichment"
BANDWIDTH=.25  # native SD; descriptive, fixed without comparing bandwidths
Q=np.array([.05,.10,.30,.50,.70,.90,.95])

def features(z):
    return np.column_stack([z,z*z])

def raw_pairs(v,ix,valid,left,right,kept,center,scale,blocks):
    n=len(kept);n_pair=len(left)
    assert len(right)==n_pair and len(blocks)==n_pair
    assert len(np.unique(np.r_[left,right]))==2*n_pair
    pidx=np.full(n,-1,dtype=int);group=np.full(n,-1,dtype=int)
    pidx[left]=np.arange(n_pair);pidx[right]=np.arange(n_pair)
    group[left]=0;group[right]=1
    mask=valid&(pidx[ix]>=0)
    ids=ix[mask]
    z=(v[mask]-center)/scale
    p=pidx[ids];y=group[ids].astype(float);w=1./kept[ids]
    for label in (0,1):
        a=np.bincount(p[y==label],weights=w[y==label],minlength=n_pair)
        np.testing.assert_allclose(a,1.,rtol=0,atol=1e-10)
    return dict(z=z,f=features(z),p=p,y=y,w=w,blocks=blocks,
                cells=ids,v=v[mask])

def state(d,beta,hessian=True):
    f,p,y,w=d["f"],d["p"],d["y"],d["w"]
    n=len(d["blocks"])
    eta=f@beta
    if not np.isfinite(eta).all():
        raise ValueError("nonfinite predictor")
    minimum=np.full(n,np.inf);maximum=np.full(n,-np.inf)
    np.minimum.at(minimum,p,eta);np.maximum.at(maximum,p,eta)
    lo=-maximum-40.;hi=-minimum+40.
    for _ in range(52):
        mid=(lo+hi)/2
        pred=expit(eta+mid[p])
        higher=np.bincount(p,weights=w*pred,minlength=n)>1.
        hi=np.where(higher,mid,hi);lo=np.where(higher,lo,mid)
    eta=eta+((lo+hi)/2)[p];pred=expit(eta)
    np.testing.assert_allclose(np.bincount(p,weights=w*pred,minlength=n),
                              1.,rtol=0,atol=1e-10)
    loss=float(w@(np.logaddexp(0.,eta)-y*eta))
    score=np.zeros((n,f.shape[1]))
    np.add.at(score,p,(w*(y-pred))[:,None]*f)
    if not hessian:
        return loss,score,None
    a=w*pred*(1-pred)
    s=np.bincount(p,weights=a,minlength=n)
    if np.any(s<=1e-12):
        raise ValueError("separated pair")
    cross=np.zeros_like(score);np.add.at(cross,p,a[:,None]*f)
    h=f.T@(a[:,None]*f)-cross.T@(cross/s[:,None])
    return loss,score,h

def fit_raw(d,active):
    active=np.asarray(active,int)
    beta=np.zeros(d["f"].shape[1])
    for _ in range(80):
        loss,s,h=state(d,beta)
        g=s.sum(0)[active]
        if not len(active) or np.max(np.abs(g))<1e-7:
            return beta,s,h
        hh=h[np.ix_(active,active)]
        if np.linalg.eigvalsh(hh).min()<=1e-9:
            raise ValueError("singular profiled information")
        step=np.linalg.solve(hh,g)
        for k in range(30):
            b=beta.copy();b[active]+=step/2**k
            if np.max(np.abs(b))>50:
                continue
            try:
                new,_,_=state(d,b,hessian=False)
            except (ValueError,AssertionError):
                continue
            if new<=loss+1e-10:
                beta=b
                break
        else:
            raise ValueError("line search failed")
    raise ValueError("fit did not converge")

def inference(d):
    beta,s,h=fit_raw(d,np.arange(6))
    g=d["blocks"];n=len(np.unique(g));flips=sign_patterns(g)
    if np.linalg.eigvalsh(h).min()<=1e-9:
        raise ValueError("singular full information")
    bs=group_sums(s,g)
    inv=np.linalg.inv(h)
    cov=inv@(bs.T@bs)@inv*n/(n-1)
    crit=student_t.ppf(.975,n-1)
    rows=[]
    for j in range(6):
        nuisance=np.delete(np.arange(6),j)
        _,ss,hh=fit_raw(d,nuisance)
        q=ss[:,j]-ss[:,nuisance]@np.linalg.solve(hh[np.ix_(nuisance,nuisance)],hh[nuisance,j])
        q=group_sums(q[:,None],g).ravel()
        f=flips@q
        p=float(np.mean(np.abs(f)>=abs(f[0])-1e-12))
        se=np.sqrt(max(cov[j,j],0.))
        rows.append(dict(trait=NAMES[j%3],term="directional" if j<3 else "quadratic",
                    status="estimated",coefficient=float(beta[j]),se=float(se),
                    lo=float(beta[j]-crit*se),hi=float(beta[j]+crit*se),p_score=p,
                    positive_score_regions=int((q>0).sum())))
    return beta,cov,crit,rows

def relative_curve(value,j,center,scale,beta,cov,crit):
    z=(value-center[j])/scale[j]
    f=np.zeros(6);f[j]=z;f[j+3]=z*z
    b=float(f@beta);se=float(np.sqrt(max(f@cov@f,0.)))
    return dict(relative_enrichment=float(np.exp(b)),
                lo=float(np.exp(b-crit*se)),hi=float(np.exp(b+crit*se)),
                log_relative_enrichment=b)

def native_reference(nv,ni,nw,n):
    _,scale=reference_scale(nv,ni,n)
    qq=np.stack([weighted_quantile(nv[:,j],nw,Q) for j in range(3)])
    center=qq[:,3]
    assert np.isfinite(scale).all() and (scale>0).all()
    return center,scale,qq

def curve_grid(qq,j):
    return np.linspace(qq[j,0],qq[j,-1],101)

def kde_mass(v,ix,kept,cells,values,h):
    """Equal cell mass, Gaussian kernel in continuous phenotype units."""
    mask=np.isin(ix,cells)
    vv=v[mask];ii=ix[mask]
    if not len(vv):
        raise ValueError("empty group")
    w=1./kept[ii]
    m=(np.exp(-.5*((vv[:,None]-values[None,:])/h)**2)*w[:,None])
    res=np.zeros((len(kept),len(values)))
    np.add.at(res,ii,m)
    return res[cells]

def env_profile(v,ix,kept,cells,climate,blocks,center,scale,qq,meta):
    # P(water | observed continuous trait) with fixed Gaussian smoothing.
    # Cell weights are summed kernel mass; no coordinates or identifiers.
    out=[]
    for j,trait in enumerate(NAMES):
        grid=curve_grid(qq,j)
        mass=kde_mass(v[:,j],ix,kept,cells,grid,BANDWIDTH*scale[j])
        for k,value in enumerate(grid):
            w=mass[:,k];total=w.sum()
            if total<=1e-12:
                continue
            gmass=group_sums(w[:,None],blocks[cells]).ravel()
            ess=total*total/(w@w)
            gess=total*total/(gmass@gmass)
            for axis in (RAIN,SEASON):
                x=climate[cells,axis]
                q=weighted_quantile(x,w,[.1,.5,.9])
                out.append(dict(**meta,trait=trait,value=float(value),
                     native_z=float((value-center[j])/scale[j]),
                     axis=VARIABLES[axis],bandwidth=float(BANDWIDTH*scale[j]),
                     effective_cells=float(ess),effective_regions=float(gess),
                     cells=len(cells),status="descriptive" if ess>=20 and gess>=3 else "sparse",
                     climate_q10=float(q[0]),climate_median=float(q[1]),
                     climate_q90=float(q[2]),climate_mean=float(w@x/total)))
    return out

def analyze(sp,definition,n,c,pairs,fn,fi,out):
    nv,ni,nw=joint_records(fn,n,COLS,False)
    iv,ii,_=joint_records(fi,c,COLS,False)
    assert len(nv)==EXPECTED[sp][4]
    cuts,lo,hi=native_cuts(nv,nw)
    nc,nvalid=rank_codes(nv,cuts,lo,hi)
    assert nvalid.all()
    support=np.bincount(nc,minlength=K)>0
    assert support.sum()==125
    h,tot,kept=cell_hist(iv,ii,len(c),cuts,lo,hi,support)
    code,valid=rank_codes(iv,cuts,lo,hi);valid&=support[code]
    np.testing.assert_array_equal(np.bincount(ii[valid],minlength=len(c)),kept)
    left=np.concatenate([l for l,r in pairs]);right=np.concatenate([r for l,r in pairs])
    g=c.iloc[left].block.to_numpy();np.testing.assert_array_equal(g,c.iloc[right].block.to_numpy())
    original=fit(h[left],h[right],X,np.arange(6))[0]
    np.testing.assert_allclose(original,EXPECTED_BETA[(sp,definition)],atol=1e-6,rtol=0)
    center,scale,qq=native_reference(nv,ni,nw,len(n))
    d=raw_pairs(iv,ii,valid,left,right,kept,center,scale,g)
    meta=dict(species=sp,definition=definition,pairs=len(left),regions=len(np.unique(g)))
    out["cohort"].append(dict(**meta,native_records=len(nv),native_cells=len(n),
              analogous_records=int(kept[left].sum()),novel_records=int(kept[right].sum()),
              excluded_analogous=int((tot[left]-kept[left]).sum()),
              excluded_novel=int((tot[right]-kept[right]).sum())))
    # The aggregate binned objective and the raw-record objective must
    # agree exactly if the raw predictors are replaced by old ranks.
    old=raw_pairs(iv,ii,valid,left,right,kept,center,scale,g)
    keep=valid&np.isin(ii,np.r_[left,right])
    old["f"]=X[code[keep]]
    _,score,_=state(old,original)
    np.testing.assert_allclose(score.sum(0),0.,atol=1e-6,rtol=0)
    try:
        beta,cov,crit,rows=inference(d)
    except (ValueError,np.linalg.LinAlgError,AssertionError) as exc:
        beta=None
        rows=[dict(trait=t,term=k,status="fit_failed",reason=str(exc),p_score=np.nan)
              for t in NAMES for k in ("directional","quadratic")]
    for row in rows:
        j=NAMES.index(row["trait"]);k=j+3*(row["term"]=="quadratic")
        out["coefficients"].append(dict(**meta,old_rank_coefficient=float(original[k]),
             native_center=float(center[j]),native_sd=float(scale[j]),**row))
    climate=climate_cells(fi,c);blocks=c.block.to_numpy()
    for j,trait in enumerate(NAMES):
        grid=curve_grid(qq,j);bandwidth=BANDWIDTH*scale[j]
        # Marginal empirical density ratio, smoothed in the actual trait units.
        # It is intentionally distinct from B's conditional ratio and is not
        # interpreted as a model residual or as an additional selection test.
        ma=kde_mass(iv[valid,j],ii[valid],kept,left,grid,bandwidth)
        mb=kde_mass(iv[valid,j],ii[valid],kept,right,grid,bandwidth)
        central_a=kde_mass(iv[valid,j],ii[valid],kept,left,np.array([center[j]]),bandwidth).mean()
        central_b=kde_mass(iv[valid,j],ii[valid],kept,right,np.array([center[j]]),bandwidth).mean()
        for k,value in enumerate(grid):
            if beta is not None:
                out["curves"].append(dict(**meta,trait=trait,value=float(value),native_z=float((value-center[j])/scale[j]),
                    **relative_curve(value,j,center,scale,beta,cov,crit)))
            wa,wb=ma[:,k],mb[:,k]
            if min(wa.sum(),wb.sum(),central_a,central_b)>1e-12:
                out["marginal_density"].append(dict(**meta,trait=trait,value=float(value),
                    bandwidth=float(bandwidth),relative_marginal_density_ratio=float((wb.mean()/wa.mean())/(central_b/central_a)),
                    analogous_effective_cells=float(wa.sum()**2/(wa@wa)),
                    novel_effective_cells=float(wb.sum()**2/(wb@wb))))
        for rank,pct,value in zip(range(-2,3),[.1,.3,.5,.7,.9],qq[j,1:-1]):
            r=dict(**meta,trait=trait,native_percentile=pct,value=float(value),
                old_rank_relative_enrichment=float(np.exp(original[j]*rank+original[j+3]*rank**2)),
                status="estimated" if beta is not None else "fit_failed")
            if beta is not None:
                r.update(relative_curve(value,j,center,scale,beta,cov,crit))
            out["quantile_comparison"].append(r)
    for label,idx in (("analogous",left),("novel",right)):
        out["environment_profiles"].extend(env_profile(iv[valid],ii[valid],kept,idx,
            climate,blocks,center,scale,qq,dict(**meta,range=label)))

def full_environment(sp,sn,si,out):
    fn,fi=eligible(sn),eligible(si)
    n=make_cells(sn,True);c=make_cells(si,True)
    nv,ni,nw=joint_records(fn,n,COLS,False)
    iv,ii,_=joint_records(fi,c,COLS,False)
    cuts,lo,hi=native_cuts(nv,nw)
    code,valid=rank_codes(nv,cuts,lo,hi)
    assert valid.all()
    support=np.bincount(code,minlength=K)>0
    assert support.sum()==125
    center,scale,qq=native_reference(nv,ni,nw,len(n))
    for label,v,ix,cells,frame in (("native",nv,ni,n,fn),("introduced",iv,ii,c,fi)):
        code,valid=rank_codes(v,cuts,lo,hi);valid&=support[code]
        kept=np.bincount(ix[valid],minlength=len(cells))
        selected=np.flatnonzero(kept>0)
        assert len(selected)>0
        out["environment_profiles"].extend(env_profile(v[valid],ix[valid],kept,selected,
            climate_cells(frame,cells),cells.block.to_numpy(),center,scale,qq,
            dict(species=sp,definition="full_range",range=label,pairs=0,regions=len(np.unique(cells.block.to_numpy()[selected])))))

def run_species(sp,sn,si,transform,edges,nat,inv,out):
    fn,fi=eligible(sn),eligible(si)
    n=make_cells(sn,True);full=make_cells(si,True)
    ex=EXPECTED[sp]
    assert (len(n),len(full))==ex[:2]
    coords=transform(si)
    uv=np.column_stack([np.searchsorted(e,coords[:,k],side="right")-1 for k,e in enumerate(edges)])
    for k,e in enumerate(edges):uv[coords[:,k]==e[-1],k]=len(e)-2
    inside=np.all((uv>=0)&(uv<np.array([len(e)-1 for e in edges])),axis=1)
    labeled=si.loc[inside].copy()
    labeled["novel"]=(~_probability_content(nat,.95)[uv[inside,0],uv[inside,1]]).astype(int)
    c=make_cells(labeled);c=c.loc[c.novel!=.5].reset_index(drop=True)
    pairs=match_cells(c)
    assert sum(len(l) for l,r in pairs)==ex[2]
    analyze(sp,"native_core95",n,c,pairs,fn,fi,out)
    hc=head_cells(si,transform,edges)[["i","j","u","v"]]
    c=full.merge(hc,on=["i","j"],validate="one_to_one")
    u,v=c.u.to_numpy(),c.v.to_numpy()
    c=c.loc[inv[u,v]>0].reset_index(drop=True);u,v=c.u.to_numpy(),c.v.to_numpy()
    c["novel"]=(inv[u,v]>nat[u,v]).astype(int)
    pairs=match_cells(c)
    assert sum(len(l) for l,r in pairs)==ex[3]
    analyze(sp,"density_surplus",n,c,pairs,fn,fi,out)

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    auth=authority_source((ROOT/"private/authority350.zip").read_bytes())
    wild=pd.read_csv(ROOT/"private/focal_wild.csv.gz",dtype={"obs_id":str},low_memory=False)
    assert len(wild)==297622
    wild["_cell"]=angular_cells(wild.latitude,wild.longitude)
    wild=wild.sort_values("obs_id").drop_duplicates(["taxon_name","range_state","_cell","user_id"]).copy()
    land=background_mesh(auth)
    wild,land,_=apply_original_five_CHELSA(wild,land,source_meta=Path("analysis/v4/environment_source_metadata.json"),work=ROOT/"private/rasters")
    measured=load_three_axis(ROOT/"private/original-measured.csv.gz")
    wild=wild.drop(columns=[c for c in measured if c!="obs_id"],errors="ignore").merge(measured,on="obs_id",how="left",validate="one_to_one")
    frozen=json.loads(Path("analysis/v4/results/20261021_unified_abcde/S2_TWO_SPECIES_FULL_WCVP_NICHE_VERIFIED_20261021.json").read_text())
    out={"cohort":[],"coefficients":[],"curves":[],"quantile_comparison":[],"environment_profiles":[],"marginal_density":[]}
    for sp,rows in wild.groupby("taxon_name"):
        env,_=_background_species_audit(auth,land.dropna(subset=list(VARIABLES)),sp,land)
        complete=rows.dropna(subset=list(VARIABLES))
        sn=complete.loc[complete.range_state=="native"];si=complete.loc[complete.range_state=="introduced"]
        transform,edges,nat,inv=_pca_env_2d(env[env._range=="native"],env[env._range=="introduced"],sn,si,return_grid=True)
        assert abs(1-abs(nat-inv).sum()/2-frozen["species"][sp]["PCA_env_background_corrected_Schoener_D"])<1e-5
        run_species(sp,sn,si,transform,edges,nat,inv,out)
        full_environment(sp,sn,si,out)
    d=pd.DataFrame(out["coefficients"])
    assert len(d)==24
    d["holm_all24"]=multipletests(d.p_score.fillna(1.),method="holm")[1]
    d["holm_definition12"]=d.groupby("definition").p_score.transform(
        lambda v:multipletests(v.fillna(1.),method="holm")[1])
    d.to_csv(OUT/"continuous_coefficients.csv",index=False)
    for key in ("cohort","curves","quantile_comparison","environment_profiles","marginal_density"):
        pd.DataFrame(out[key]).to_csv(OUT/(key+".csv"),index=False)
    print("Complete; aggregate outputs only.",flush=True)

if __name__=="__main__":
    main()
