"""Native climatic periphery versus introduced novel occupancy.

A comparison of contemporary, continuous phenotype distributions, not a
comparison of ancestral and evolved generations. Freeze the original
annual climate core and B; no new data, thresholds selected by outcomes,
native prediction/transfer, or estimates of fitness or evolutionary rates.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from analysis.v5.phenotype_environment_explanation import (
    ROOT,SPECIES,VARIABLES,NAMES,COLS,EXPECTED,EXPECTED_BETA,
    authority_source,angular_cells,background_mesh,apply_original_five_CHELSA,
    load_three_axis,_background_species_audit,_pca_env_2d,_probability_content,
    eligible,make_cells,match_cells,head_cells,joint_records,native_cuts,
    rank_codes,cell_hist,fit,fit_raw,climate_cells,native_reference,raw_pairs,
    X,K,CONT,CALIPER,MAX_AXIS,MIN_PAIRS,MIN_REGIONS,BANDWIDTH,curve_grid,match_one,kde_mass,
)
OUT=ROOT/"phenotype_native_periphery"

def classify(frame,transform,edges,nat):
    coords=transform(frame)
    uv=np.column_stack([np.searchsorted(e,coords[:,k],side="right")-1 for k,e in enumerate(edges)])
    for k,e in enumerate(edges):
        uv[coords[:,k]==e[-1],k]=len(e)-2
    inside=np.all((uv>=0)&(uv<np.array([len(e)-1 for e in edges])),axis=1)
    labeled=frame.loc[inside].copy()
    labeled["novel"]=(~_probability_content(nat,.95)[uv[inside,0],uv[inside,1]]).astype(int)
    all_cells=make_cells(labeled)
    c=all_cells.loc[all_cells.novel!=.5].reset_index(drop=True)
    return c,len(all_cells),int((all_cells.novel==.5).sum())

def setup(frame,c,cuts,lo,hi,support):
    v,ix,_=joint_records(frame,c,COLS,False)
    h,tot,kept=cell_hist(v,ix,len(c),cuts,lo,hi,support)
    code,valid=rank_codes(v,cuts,lo,hi)
    valid&=support[code]
    np.testing.assert_array_equal(np.bincount(ix[valid],minlength=len(c)),kept)
    x=climate_cells(frame,c)
    return dict(c=c,v=v,ix=ix,h=h,tot=tot,kept=kept,valid=valid,x=x)

def pairs_for(d):
    # Cells without eligible inside-support records cannot enter the model.
    # Match only on geography; retain the original index for record recovery.
    c=d["c"].copy()
    keep=d["kept"]>0
    q=c.loc[keep].copy()
    p=match_cells(q)
    if not p:
        return np.empty(0,int),np.empty(0,int),np.empty(0,int)
    a=np.concatenate([l for l,r in p]);b=np.concatenate([r for l,r in p])
    g=c.iloc[a].block.to_numpy()
    np.testing.assert_array_equal(g,c.iloc[b].block.to_numpy())
    return a,b,g

def subset_pairs(p,ind):
    return tuple(x[ind] for x in p)

def transition_features(x,p,mu,sd):
    a,b,_=p
    return np.column_stack(((x[a]-mu)/sd,(x[b]-mu)/sd))

def comparable_transitions(dn,pn,di,pi,mu,sd):
    a=transition_features(dn["x"],pn,mu,sd)
    b=transition_features(di["x"],pi,mu,sd)
    return match_one(a,b,tuple(range(2*len(mu))))

def summaries(out,meta,d,idx,center,scale,qq):
    v,ix,valid,kept=d["v"],d["ix"],d["valid"],d["kept"]
    idx=np.asarray(idx,int)
    idx=idx[kept[idx]>0]
    c=d["c"];x=d["x"]
    blocks=c.iloc[idx].block.to_numpy()
    info=dict(**meta,cells=len(idx),regions=len(np.unique(blocks)),
              retained_records=int(kept[idx].sum()),
              outside_joint_support=int((d["tot"][idx]-kept[idx]).sum()))
    out["cohort"].append(info)
    if not len(idx):
        return
    for j,axis in enumerate(VARIABLES):
        q=np.quantile(x[idx,j],[.1,.5,.9])
        out["climate_profiles"].append(dict(**info,axis=axis,mean=float(x[idx,j].mean()),
                                           q10=float(q[0]),median=float(q[1]),q90=float(q[2])))
    use=valid&np.isin(ix,idx)
    w=1./kept[ix[use]]/len(idx)
    np.testing.assert_allclose(w.sum(),1.,atol=1e-10)
    for j,trait in enumerate(NAMES):
        y=v[use,j];m=float(w@y)
        out["phenotype_profiles"].append(dict(**info,trait=trait,mean=m,
            sd=float(np.sqrt(max(w@(y*y)-m*m,0))),
            native_center=float(center[j]),native_sd=float(scale[j])))
        grid=curve_grid(qq,j);bw=BANDWIDTH*scale[j]
        a=kde_mass(v[valid,j],ix[valid],kept,idx,grid,bw)/(np.sqrt(2*np.pi)*bw)
        for k,t in enumerate(grid):
            mass=a[:,k]
            ess=mass.sum()**2/(mass@mass) if (mass@mass)>0 else 0.
            out["marginal_density"].append(dict(**info,trait=trait,value=float(t),
                  density=float(mass.mean()),effective_cells=float(ess)))

def fit_contrast(out,meta,d,p,center,scale,qq,require_gate=True):
    a,b,g=p
    status="estimated"
    beta=None
    if not len(a) or (require_gate and (len(a)<MIN_PAIRS or len(np.unique(g))<MIN_REGIONS)):
        status="limited_support"
    else:
        try:
            rd=raw_pairs(d["v"],d["ix"],d["valid"],a,b,d["kept"],center,scale,g)
            beta=fit_raw(rd,np.arange(6))[0]
        except (ValueError,np.linalg.LinAlgError,AssertionError):
            status="fit_failed"
    info=dict(**meta,pairs=len(a),regions=len(np.unique(g)),status=status)
    for j,trait in enumerate(NAMES):
        row=dict(**info,trait=trait)
        if beta is not None:
            row.update(linear=float(beta[j]),quadratic=float(beta[j+3]))
            for t in curve_grid(qq,j):
                z=(t-center[j])/scale[j]
                lg=beta[j]*z+beta[j+3]*z*z
                out["curves"].append(dict(**info,trait=trait,value=float(t),
                    log_relative_enrichment=float(lg),relative_enrichment=float(np.exp(lg))))
        out["coefficients"].append(row)
    return beta

def run_species(sp,sn,si,transform,edges,nat,inv,out):
    fn,fi=eligible(sn),eligible(si)
    n=make_cells(sn,True);full=make_cells(si,True)
    assert (len(n),len(full))==EXPECTED[sp][:2]
    nv,ni,nw=joint_records(fn,n,COLS,False)
    assert len(nv)==EXPECTED[sp][4]
    cuts,lo,hi=native_cuts(nv,nw)
    rc,valid=rank_codes(nv,cuts,lo,hi);assert valid.all()
    support=np.bincount(rc,minlength=K)>0;assert support.sum()==125
    center,scale,qq=native_reference(nv,ni,nw,len(n))
    cn,n_class,n_tie=classify(sn,transform,edges,nat)
    ci,i_class,i_tie=classify(si,transform,edges,nat)
    dn=setup(fn,cn,cuts,lo,hi,support)
    di=setup(fi,ci,cuts,lo,hi,support)
    pn=pairs_for(dn);pi=pairs_for(di)
    a,b,g=pi
    assert len(a)==EXPECTED[sp][2]
    old=fit(di["h"][a],di["h"][b],X,np.arange(6))[0]
    np.testing.assert_allclose(old,EXPECTED_BETA[(sp,"native_core95")],atol=1e-6,rtol=0)
    meta=dict(species=sp,definition="native_core95")
    for side,d,p,source_n,nclass,ntie in (
        ("native",dn,pn,len(n),n_class,n_tie),
        ("introduced",di,pi,len(full),i_class,i_tie)):
        c=d["c"]
        for k,group in ((0,"core"),(1,"outside_core")):
            summaries(out,dict(**meta,range=side,comparison="all_classified",group=group,
                       source_cells=source_n,classified_before_ties=nclass,tied_cells=ntie),
                       d,np.flatnonzero(c.novel.to_numpy()==k),center,scale,qq)
            summaries(out,dict(**meta,range=side,comparison="geographic_pairs",group=group,
                       source_cells=source_n,classified_before_ties=nclass,tied_cells=ntie),
                       d,p[k],center,scale,qq)
    bn=fit_contrast(out,dict(**meta,range="native",comparison="geographic_pairs"),dn,pn,center,scale,qq)
    bi=fit_contrast(out,dict(**meta,range="introduced",comparison="geographic_pairs"),di,pi,center,scale,qq,False)
    np.testing.assert_allclose(bi,CONT[(sp,"native_core95")],atol=2e-6,rtol=0)
    # The two domains may have very different climate transitions. Match
    # transitions only on both endpoints of the same 5 annual axes,
    # retaining each original within-range geographic pair intact.
    x=climate_cells(fn,n);mu=x.mean(0);sd=x.std(0)
    assert np.isfinite(sd).all() and np.min(sd)>0
    l,r=comparable_transitions(dn,pn,di,pi,mu,sd)
    pnn=subset_pairs(pn,l);pii=subset_pairs(pi,r)
    nn,ni2=len(np.unique(pnn[2])),len(np.unique(pii[2]))
    ok=len(l)>=MIN_PAIRS and min(nn,ni2)>=MIN_REGIONS
    out["transition_support"].append(dict(**meta,native_pairs=len(pn[0]),introduced_pairs=len(pi[0]),
         native_regions=len(np.unique(pn[2])),introduced_regions=len(np.unique(pi[2])),
         common_transitions=len(l),native_retained_regions=nn,introduced_retained_regions=ni2,
         status="descriptive" if ok else "limited_support"))
    if len(l):
        za=transition_features(dn["x"],pnn,mu,sd)
        zb=transition_features(di["x"],pii,mu,sd)
        delta=zb-za
        assert np.max(np.sqrt(np.mean(delta**2,axis=1)))<=CALIPER+1e-10
        assert np.max(abs(delta))<=MAX_AXIS+1e-10
        for j,axis in enumerate(VARIABLES):
            for k,endpoint in ((0,"core"),(1,"outside_core")):
                q=delta[:,k*len(VARIABLES)+j]
                out["transition_balance"].append(dict(**meta,axis=axis,endpoint=endpoint,
                    common_transitions=len(l),native_mean=float(dn["x"][pnn[k],j].mean()),
                    introduced_mean=float(di["x"][pii[k],j].mean()),
                    signed_difference_native_sd=float(q.mean()),
                    mean_abs_difference_native_sd=float(abs(q).mean())))
    if ok:
        bn2=fit_contrast(out,dict(**meta,range="native",comparison="shared_transition"),dn,pnn,center,scale,qq)
        bi2=fit_contrast(out,dict(**meta,range="introduced",comparison="shared_transition"),di,pii,center,scale,qq)
    else:
        bn2=bi2=None
        for side in ("native","introduced"):
            for trait in NAMES:
                out["coefficients"].append(dict(**meta,range=side,comparison="shared_transition",trait=trait,
                     pairs=len(l),regions=nn if side=="native" else ni2,status="limited_support"))
    for comparison,a,b in (("geographic_pairs",bn,bi),("shared_transition",bn2,bi2)):
        for j,trait in enumerate(NAMES):
            q=dict(**meta,comparison=comparison,trait=trait,status="descriptive" if a is not None and b is not None else "limited_support")
            if a is not None and b is not None:
                q.update(introduced_minus_native_linear=float(b[j]-a[j]),
                         introduced_minus_native_quadratic=float(b[j+3]-a[j+3]))
            out["curve_contrast"].append(q)

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
    out={k:[] for k in ("cohort","climate_profiles","phenotype_profiles","marginal_density","coefficients","curves","transition_support","transition_balance","curve_contrast")}
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
