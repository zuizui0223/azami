"""Record-level seasonal exposure support and continuous phenotype distributions.

Freeze continuous B and the annual climate assignments. WorldClim monthly
normals describe typical observation-period environments, not the weather
during the individual photograph or measured developmental exposure.
Only aggregate outputs leave the private runner.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.optimize import linprog
from scipy.sparse import coo_matrix, csr_matrix
from analysis.v5.phenotype_environment_explanation import (
    ROOT,SPECIES,DEFINITIONS,VARIABLES,NAMES,COLS,EXPECTED,EXPECTED_BETA,
    authority_source,angular_cells,background_mesh,apply_original_five_CHELSA,
    load_three_axis,_background_species_audit,_pca_env_2d,_probability_content,
    eligible,make_cells,match_cells,head_cells,joint_records,native_cuts,
    rank_codes,cell_hist,fit,fit_raw,climate_cells,native_reference,raw_pairs,
    X,K,CONT,CALIPER,MAX_AXIS,MIN_PAIRS,MIN_REGIONS,BANDWIDTH,curve_grid,rematch,
)
from analysis.v5.phenotype_redistribution import cellkeys
from analysis.v5.seasonal_exposure import normals
OUT=ROOT/"phenotype_flowering_exposure"
EXPOSURES=("observation_month_prec","observation_month_srad","previous_month_tmax")
SCHEMES={EXPOSURES[j]:(j,) for j in range(3)}
SCHEMES.update(seasonal_joint=(0,1,2),annual_five_and_seasonal=(0,1,2))
MIN_MASS=.5
ENV_BW=.25
GRID=41
MIN_ESS=20.
MIN_GESS=3.

def attach_exposure(wild):
    # Sample only records with all three usable-head-gated measurements and
    # dates. Never delete rows before the frozen annual cohort is reconstructed.
    good=np.isfinite(wild[COLS].to_numpy(float)).all(1)
    month=pd.to_datetime(wild.observed_on,format="mixed",errors="coerce",utc=True).dt.month
    good&=month.notna().to_numpy()
    for field in EXPOSURES:
        wild[field]=np.nan
    if not good.any():
        raise ValueError("No dated joint-trait records")
    d=wild.loc[good]
    folder=ROOT/"temp";folder.mkdir(parents=True,exist_ok=True)
    monthly=normals(d,folder)
    m=month.loc[good].astype(int).to_numpy()-1
    i=np.arange(len(d))
    vals=np.column_stack((monthly["prec"][i,m],monthly["srad"][i,m],
                           monthly["tmax"][i,(m-1)%12]))
    wild.loc[good,list(EXPOSURES)]=vals
    return wild

def align_records(frame,c,values,indices):
    # Exactly the predicate and row order of phenotype_reservoir.records.
    ix=pd.Index(cellkeys(c)).get_indexer(frame.cellkey)
    v=frame[COLS].to_numpy(float)
    ok=(ix>=0)&np.isfinite(v).all(1)
    np.testing.assert_array_equal(ix[ok],indices)
    np.testing.assert_allclose(v[ok],values,rtol=0,atol=0)
    return frame.loc[ok,list(EXPOSURES)].to_numpy(float)

def reference_exposure(x,ix):
    good=np.isfinite(x).all(1)
    count=np.bincount(ix[good])
    if not good.any():
        raise ValueError("No complete native exposures")
    w=1./count[ix[good]]
    mu=np.average(x[good],weights=w,axis=0)
    sd=np.sqrt(np.average((x[good]-mu)**2,weights=w,axis=0))
    if not np.isfinite(sd).all() or np.min(sd)<=0:
        raise ValueError("Degenerate native exposure scales")
    return mu,sd,int(good.sum()),int((count>0).sum())

def flow(a,b,wa,wb,axes):
    """Maximum common continuous exposure mass, then minimum RMS distance.

    Each eligible record carries its original cell mass, including loss
    from missing exposures. Only mass on admissible edges can be retained.
    No phenotype used. Two LPs avoid trading off match quantity vs distance.
    Returns edge endpoints and flows; no observation IDs are exported.
    """
    ga=np.isfinite(a).all(1);gb=np.isfinite(b).all(1)
    ia=np.flatnonzero(ga);ib=np.flatnonzero(gb)
    if not len(ia) or not len(ib):
        return np.empty(0,int),np.empty(0,int),np.empty(0)
    delta=a[ia][:,None,axes]-b[ib][None,:,axes]
    dist=np.sqrt(np.mean(delta*delta,axis=2))
    allowed=(dist<=CALIPER+1e-12)&(np.max(abs(delta),axis=2)<=MAX_AXIS+1e-12)
    l,r=np.nonzero(allowed);m=len(l)
    if not m:
        return np.empty(0,int),np.empty(0,int),np.empty(0)
    l,r=ia[l],ib[r]
    col=np.arange(m)
    mat=coo_matrix((np.ones(2*m),(np.r_[l,len(a)+r],np.r_[col,col])),
                   shape=(len(a)+len(b),m)).tocsr()
    bounds=np.r_[wa,wb]
    fit1=linprog(-np.ones(m),A_ub=mat,b_ub=bounds,bounds=(0,None),method="highs")
    if not fit1.success:
        raise ValueError("max-mass transport failed")
    mass=float(fit1.x.sum())
    if mass<=1e-10:
        return np.empty(0,int),np.empty(0,int),np.empty(0)
    cost=np.sqrt(np.mean((a[l][:,axes]-b[r][:,axes])**2,axis=1))
    fit2=linprog(cost,A_ub=mat,b_ub=bounds,A_eq=csr_matrix(np.ones((1,m))),
                 b_eq=[mass],bounds=(0,None),method="highs")
    if not fit2.success:
        raise ValueError("min-distance transport failed")
    f=np.maximum(fit2.x,0.)
    assert np.max(mat@f-bounds)<1e-7
    np.testing.assert_allclose(f.sum(),mass,atol=1e-7,rtol=0)
    keep=f>1e-12
    return l[keep],r[keep],f[keep]

def curve_rows(out,meta,beta,center,scale,qq,status):
    for j,trait in enumerate(NAMES):
        row=dict(**meta,trait=trait,status=status)
        if beta is not None:
            b,c=beta[j],beta[j+3]
            point=center[j]-b*scale[j]/(2*c) if abs(c)>1e-9 else np.nan
            row.update(linear=float(b),quadratic=float(c),turning_point=float(point),
                turning_point_in_display=bool(qq[j,0]<=point<=qq[j,-1]))
            for value in curve_grid(qq,j):
                z=(value-center[j])/scale[j]
                out["curves"].append(dict(**meta,trait=trait,value=float(value),
                    relative_enrichment=float(np.exp(b*z+c*z*z))))
        out["coefficients"].append(row)

def selected_fit(v,ix,valid,kept,left,right,blocks,center,scale,record_weights):
    d=raw_pairs(v,ix,valid,left,right,kept,center,scale,blocks)
    use=valid&np.isin(ix,np.r_[left,right])
    if record_weights is not None:
        d["w"]=record_weights[use]
    for y in (0,1):
        q=d["y"]==y
        np.testing.assert_allclose(np.bincount(d["p"][q],weights=d["w"][q],minlength=len(left)),1.,
                                  atol=1e-8,rtol=0)
    return fit_raw(d,np.arange(6))[0]

def compare(out,meta,v,ix,valid,kept,left,right,blocks,center,scale,qq,x,z,axes):
    lookup={int(c):np.flatnonzero(valid&(ix==c)) for c in np.r_[left,right]}
    matched=[];mass_list=[];errors=0
    wmatch=np.zeros(len(v));wcomplete=np.zeros(len(v))
    diffs=[];sums=[]
    for i,(l,r) in enumerate(zip(left,right)):
        a,b=lookup[int(l)],lookup[int(r)]
        ga=np.isfinite(z[a]).all(1);gb=np.isfinite(z[b]).all(1)
        if not ga.any() or not gb.any():
            mass_list.append(0.)
            continue
        try:
            ia,ib,f=flow(z[a],z[b],np.full(len(a),1./kept[l]),
                         np.full(len(b),1./kept[r]),axes)
        except ValueError:
            errors+=1;mass_list.append(np.nan);continue
        mass=float(f.sum());mass_list.append(mass)
        if mass<MIN_MASS-1e-9:
            continue
        matched.append(i)
        np.add.at(wmatch,a[ia],f/mass)
        np.add.at(wmatch,b[ib],f/mass)
        wcomplete[a[ga]]=1./ga.sum();wcomplete[b[gb]]=1./gb.sum()
        d=z[b[ib]]-z[a[ia]]
        diffs.append((f/mass)@np.abs(d))
        sums.append((f/mass)@d)
    ind=np.array(matched,int)
    a,b,g=left[ind],right[ind],blocks[ind]
    enough=len(a)>=MIN_PAIRS and len(np.unique(g))>=MIN_REGIONS
    meta=dict(**meta,candidate_pairs=len(left),candidate_regions=len(np.unique(blocks)),
              matched_pairs=len(a),matched_regions=len(np.unique(g)))
    out["coverage"].append(dict(**meta,status=("transport_incomplete" if errors else
          "descriptive" if enough else "limited_support"),
          transport_failures=errors,minimum_mass=MIN_MASS,
          mean_matched_mass_all_pairs=float(np.nanmean(mass_list)) if len(left) and np.isfinite(mass_list).any() else np.nan,
          mean_matched_mass_retained=float(np.asarray(mass_list)[ind].mean()) if len(ind) else np.nan,
          original_records=int(kept[np.r_[left,right]].sum()),
          retained_records=int((wmatch>0).sum())))
    if len(ind):
        for j,axis in enumerate(EXPOSURES):
            # Same retained cell pairs, date-complete records before transport.
            delta=[]
            for l,r in zip(a,b):
                ia=lookup[int(l)];ib=lookup[int(r)]
                ia=ia[wcomplete[ia]>0];ib=ib[wcomplete[ib]>0]
                delta.append(wcomplete[ib]@z[ib,j]-wcomplete[ia]@z[ia,j])
            out["balance"].append(dict(**meta,axis=axis,matched_axis=j in axes,
                complete_signed_difference_sd=float(np.mean(delta)),
                matched_signed_difference_sd=float(np.mean(sums,axis=0)[j]),
                matched_mean_absolute_edge_difference_sd=float(np.mean(diffs,axis=0)[j])))
    if errors or not enough:
        for weighting in ("same_pairs_all_records","same_pairs_dated","exposure_matched"):
            curve_rows(out,dict(**meta,weighting=weighting),None,center,scale,qq,
                       "transport_incomplete" if errors else "limited_support")
        return
    for weighting,w in (("same_pairs_all_records",None),("same_pairs_dated",wcomplete),
                        ("exposure_matched",wmatch)):
        try:
            beta=selected_fit(v,ix,valid,kept,a,b,g,center,scale,w)
            curve_rows(out,dict(**meta,weighting=weighting),beta,center,scale,qq,"descriptive_refit")
        except (ValueError,AssertionError,np.linalg.LinAlgError):
            curve_rows(out,dict(**meta,weighting=weighting),None,center,scale,qq,"fit_failed")

def local_weights(exposure,p,y,w,blocks,at):
    u=(exposure-at)/ENV_BW
    # Missing values must get zero without NaN propagation.
    k=np.zeros_like(w)
    good=np.isfinite(u)&(abs(u)<=2.)
    k[good]=w[good]*np.exp(-.5*u[good]**2)
    n=len(blocks)
    sums=[np.bincount(p[y==j],weights=k[y==j],minlength=n) for j in (0,1)]
    common=np.minimum(*sums)
    total=common.sum()
    if total<=0:
        return np.zeros_like(w),0.,0.,0.
    factor=[np.divide(common,s,out=np.zeros_like(s),where=s>0) for s in sums]
    wnew=k*np.where(y==0,factor[0][p],factor[1][p])/total
    for j in (0,1):
        np.testing.assert_allclose(np.bincount(p[y==j],weights=wnew[y==j],minlength=n),
                                   common/total,atol=1e-10,rtol=0)
    gs=np.array([common[blocks==g].sum()/total for g in np.unique(blocks)])
    return wnew,float(total**2/(common@common)),float(1/(gs@gs)),float(total/n)

def local_profiles(out,meta,d,exposure,mu,sd,qq,scale):
    z=(exposure-mu)/sd
    p,y,w,blocks=d["p"],d["y"],d["w"],d["blocks"]
    complete=np.isfinite(z).all(1)
    z[~complete]=np.nan
    v=d["v"]
    grids=[curve_grid(qq,j)[::2] for j in range(3)]
    density=[np.exp(-.5*((v[:,j,None]-grid[None,:])/(BANDWIDTH*scale[j]))**2)/
             (np.sqrt(2*np.pi)*BANDWIDTH*scale[j]) for j,grid in enumerate(grids)]
    for j,axis in enumerate(EXPOSURES):
        good=complete&np.isfinite(z[:,j])
        if not good.any():
            continue
        bounds=np.quantile(z[good,j],[.05,.95])
        for at in np.linspace(*bounds,GRID):
            ww,ess,gess,mass=local_weights(z[:,j],p,y,w,blocks,at)
            m=dict(**meta,axis=axis,exposure_value=float(mu[j]+sd[j]*at),
                   effective_pairs=ess,effective_regions=gess,common_mass_fraction=mass,
                   status="descriptive" if ess>=MIN_ESS and gess>=MIN_GESS else "sparse")
            if mass<=0:
                out["local_profiles"].append(dict(**m,trait="all",status_empty=True))
                continue
            a,b=y==0,y==1
            for k,trait in enumerate(NAMES):
                ma=float(ww[a]@v[a,k]);mb=float(ww[b]@v[b,k])
                va=float(np.sqrt(max(ww[a]@(v[a,k]**2)-ma*ma,0)))
                vb=float(np.sqrt(max(ww[b]@(v[b,k]**2)-mb*mb,0)))
                out["local_profiles"].append(dict(**m,trait=trait,
                    analogous_mean=ma,novel_mean=mb,mean_difference=mb-ma,
                    analogous_sd=va,novel_sd=vb,sd_difference=vb-va))
                fa=ww[a]@density[k][a];fb=ww[b]@density[k][b]
                for t,aa,bb in zip(grids[k],fa,fb):
                    out["local_density"].append(dict(**m,trait=trait,value=float(t),
                        analogous_density=float(aa),novel_density=float(bb),difference=float(bb-aa)))

def analyze(sp,definition,n,c,pairs,fn,fi,out):
    nv,ni,nw=joint_records(fn,n,COLS,False)
    iv,ii,_=joint_records(fi,c,COLS,False)
    assert len(nv)==EXPECTED[sp][4]
    cuts,lo,hi=native_cuts(nv,nw)
    rc,good=rank_codes(nv,cuts,lo,hi);assert good.all()
    support=np.bincount(rc,minlength=K)>0;assert support.sum()==125
    h,tot,kept=cell_hist(iv,ii,len(c),cuts,lo,hi,support)
    code,valid=rank_codes(iv,cuts,lo,hi);valid&=support[code]
    np.testing.assert_array_equal(np.bincount(ii[valid],minlength=len(c)),kept)
    left=np.concatenate([l for l,r in pairs]);right=np.concatenate([r for l,r in pairs])
    blocks=c.iloc[left].block.to_numpy()
    np.testing.assert_array_equal(blocks,c.iloc[right].block.to_numpy())
    old=fit(h[left],h[right],X,np.arange(6))[0]
    np.testing.assert_allclose(old,EXPECTED_BETA[(sp,definition)],atol=1e-6,rtol=0)
    center,scale,qq=native_reference(nv,ni,nw,len(n))
    d=raw_pairs(iv,ii,valid,left,right,kept,center,scale,blocks)
    beta=fit_raw(d,np.arange(6))[0]
    np.testing.assert_allclose(beta,CONT[(sp,definition)],atol=2e-6,rtol=0)
    xn=align_records(fn,n,nv,ni);xi=align_records(fi,c,iv,ii)
    mu,sd,nrec,ncell=reference_exposure(xn,ni)
    z=(xi-mu)/sd
    meta=dict(species=sp,definition=definition,original_pairs=len(left),
              original_regions=len(np.unique(blocks)))
    for j,axis in enumerate(EXPOSURES):
        out["native_scales"].append(dict(**meta,axis=axis,native_mean=float(mu[j]),
                       native_sd=float(sd[j]),native_dated_records=nrec,native_dated_cells=ncell))
    curve_rows(out,dict(**meta,scheme="original",weighting="all_records",
              candidate_pairs=len(left),candidate_regions=len(np.unique(blocks)),
              matched_pairs=len(left),matched_regions=len(np.unique(blocks))),
              beta,center,scale,qq,"original_reproduced")
    ann=climate_cells(fi,c);ref=climate_cells(fn,n)
    annual_z=(ann-ref.mean(0))/ref.std(0)
    la,ra,ga=rematch(left,right,blocks,annual_z,tuple(range(len(VARIABLES))))
    for name,axes in SCHEMES.items():
        a,b,g=(la,ra,ga) if name=="annual_five_and_seasonal" else (left,right,blocks)
        compare(out,dict(**meta,scheme=name),iv,ii,valid,kept,a,b,g,center,scale,qq,xi,z,axes)
    original_mask=valid&np.isin(ii,np.r_[left,right])
    local_profiles(out,meta,d,xi[original_mask],mu,sd,qq,scale)

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
    wild=attach_exposure(wild)
    frozen=json.loads(Path("analysis/v4/results/20261021_unified_abcde/S2_TWO_SPECIES_FULL_WCVP_NICHE_VERIFIED_20261021.json").read_text())
    out={k:[] for k in ("coverage","balance","coefficients","curves","local_profiles","local_density","native_scales")}
    for sp,rows in wild.groupby("taxon_name"):
        env,_=_background_species_audit(auth,land.dropna(subset=list(VARIABLES)),sp,land)
        complete=rows.dropna(subset=list(VARIABLES))
        sn=complete.loc[complete.range_state=="native"];si=complete.loc[complete.range_state=="introduced"]
        transform,edges,nat,inv=_pca_env_2d(env[env._range=="native"],env[env._range=="introduced"],sn,si,return_grid=True)
        assert abs(1-abs(nat-inv).sum()/2-frozen["species"][sp]["PCA_env_background_corrected_Schoener_D"])<1e-5
        run_species(sp,sn,si,transform,edges,nat,inv,out)
    for key,rows in out.items():
        pd.DataFrame(rows).to_csv(OUT/(key+".csv"),index=False)
    print("Complete; aggregate outputs only.",flush=True)

if __name__=="__main__":
    main()
