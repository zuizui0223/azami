"""Fixed exploration of mechanisms underlying B's relative-enrichment curves.

Same five climate axes, native reference, matched cells and 125 support bins.
Separate within-cell, among-cell within-region and among-region variation.
Ask whether climate-dependent directional/curvature terms improve
conditional allocation of phenotype mass in held-out geographic blocks.
No selection gradients, no causal mediation, no new threshold search.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit
from analysis.v5.phenotype_water_interaction import (
    ROOT,SPECIES,DEFINITIONS,VARIABLES,RAIN,SEASON,NAMES,COLS,EXPECTED,
    EXPECTED_BETA,joint_records,eligible,load_three_axis,make_cells,match_cells,
    head_cells,background_mesh,_background_species_audit,_pca_env_2d,
    _probability_content,apply_original_five_CHELSA,authority_source,
    angular_cells,X,K,native_cuts,rank_codes,cell_hist,fit,
    climate_cells,common_scale,
)
from analysis.v5.phenotype_three_axis_dispersion import (
    reference_scale,cell_moments,
)
OUT=ROOT/"phenotype_curve_mechanism"
FOCAL=np.array([NAMES.index("chroma"),NAMES.index("angle"),
                NAMES.index("chroma")+3,NAMES.index("angle")+3])
OTHER=tuple(j for j in range(5) if j not in (RAIN,SEASON))
ALL=tuple(range(5))
# Fixed moderation penalty: 0.5*4*sum(gamma^2) in pair-weighted
# profiled log loss, equivalent to a .5-scale weak ridge on the rank terms.
# It is not tuned on the data; the baseline coefficients are unpenalized.
RIDGE=4.
SEED=20261127
MODEL_AXES={"constant":(),"other_climate":OTHER,"all_climate":ALL}
# Match the existing aggregate raw dispersion (native-SD squared).
EXPECTED_DELTA={
    (SPECIES[0],"native_core95"):[[-.24955437,.01888015],[.07759093,-.08181781]],
    (SPECIES[0],"density_surplus"):[[-.22806389,.13116686],[.03134796,-.02795196]],
    (SPECIES[1],"native_core95"):[[.01914196,-.15119085],[.17463344,.01740958]],
    (SPECIES[1],"density_surplus"):[[.01525498,.12668516],[-.05467250,.07820684]],
}
COMPONENTS=("total","between_region","between_cells_within_region","within_cell")

def region_moments(m,blocks):
    ids=np.unique(blocks);n=np.array([(blocks==g).sum() for g in ids])
    # m contains cell means, second moments and squared cell means.
    s=np.stack([m[blocks==g].mean(0) for g in ids])
    return s,n

def components(regions,n,weights):
    # Repeated blocks in a bootstrap carry their original sample sizes.
    # weights: (resamples, blocks), integer multiplicities.
    a,q,b=np.split(regions,3,axis=1)
    w=np.asarray(weights,float)*n[None,:]
    w=w/w.sum(1,keepdims=True)
    mean=w@a
    within=w@(q-b)
    mid=w@(b-a*a)
    between=w@(a*a)-mean*mean
    total=w@q-mean*mean
    result=np.stack([total,between,mid,within],axis=1)
    np.testing.assert_allclose(result[:,0],result[:,1:].sum(1),atol=1e-9)
    assert np.min(result)>-1e-8
    return result

def variance_diagnostic(left,right,blocks):
    a,n=region_moments(left,blocks);b,n2=region_moments(right,blocks)
    np.testing.assert_array_equal(n,n2)
    g=len(n)
    rng=np.random.default_rng(SEED)
    w=np.vstack([np.ones(g),rng.multinomial(g,np.full(g,1/g),size=3999)])
    ca=components(a,n,w);cb=components(b,n,w)
    diff=cb-ca
    rows=[]
    for i,component in enumerate(COMPONENTS):
        for j,name in enumerate(NAMES):
            lo,hi=np.quantile(diff[1:,i,j],[.025,.975])
            rows.append(dict(trait=name,component=component,
                             analogous=float(ca[0,i,j]),novel=float(cb[0,i,j]),
                             difference=float(diff[0,i,j]),lo=float(lo),hi=float(hi)))
    # Absolute region-level mean shifts vs pooled signed shift; equal-cell
    # weighted. No locations or regional IDs exported. Signs are descriptive.
    d=b[:,:3]-a[:,:3]
    ww=n/n.sum()
    shifts=[]
    for j,name in enumerate(NAMES):
        v=d[:,j];avg=float(ww@v);absolute=float(ww@np.abs(v))
        shifts.append(dict(trait=name,mean_shift=avg,
                           mean_abs_region_shift=absolute,
                           cancellation_fraction=(1-abs(avg)/absolute if absolute>1e-12 else 0.),
                           positive_regions=int((v>0).sum()),negative_regions=int((v<0).sum()),
                           zero_regions=int((v==0).sum()),positive_pair_weight=float(ww[v>0].sum())))
    return rows,shifts

def eta_profile(a,b,beta,z,axes):
    # beta = six original rank coefficients, followed by four
    # focal-term modifiers for each selected climate difference.
    m=a+b
    eta=np.broadcast_to(X@beta[:6],a.shape).copy()
    if len(axes):
        gamma=beta[6:].reshape(len(axes),len(FOCAL))
        eta+=(z[:,axes]@gamma)@X[:,FOCAL].T
    lower=-eta.max(1)-40.;upper=-eta.min(1)+40.
    for _ in range(48):
        mid=(lower+upper)/2
        p=expit(eta+mid[:,None])
        hi=(m*p).sum(1)>1.
        upper=np.where(hi,mid,upper)
        lower=np.where(hi,lower,mid)
    eta=eta+((lower+upper)/2)[:,None]
    p=expit(eta)
    np.testing.assert_allclose((m*p).sum(1),1.,atol=1e-10,rtol=0)
    loss=(m*np.logaddexp(0.,eta)-b*eta).sum(1)
    score=(b-m*p)@X
    grad=-score.sum(0)
    if len(axes):
        grad=np.r_[grad,-(z[:,axes].T@score[:,FOCAL]).ravel()]
    return loss,grad,p

def objective(beta,a,b,z,axes):
    loss,grad,_=eta_profile(a,b,beta,z,axes)
    penalty=.5*RIDGE*np.sum(beta[6:]**2)
    grad[6:]+=RIDGE*beta[6:]
    return float(loss.sum()+penalty),grad

def fit_model(a,b,z,axes):
    start=fit(a,b,X,np.arange(6))[0]
    if not len(axes):
        return start
    start=np.r_[start,np.zeros(len(axes)*len(FOCAL))]
    opt=minimize(objective,start,args=(a,b,z,axes),jac=True,method="L-BFGS-B",
                 options=dict(maxiter=500,ftol=1e-13,gtol=1e-7,maxls=40))
    val,grad=objective(opt.x,a,b,z,axes)
    if (not opt.success or not np.isfinite(val)
            or np.max(np.abs(grad))>1e-4 or np.max(np.abs(opt.x))>12):
        raise ValueError("Unqualified climate moderation fit")
    return opt.x

def marginals(h):
    return np.stack([h.reshape(-1,5,5,5).sum(axis=tuple(1+k for k in range(3) if k!=j))
                     for j in range(3)],axis=1)

def climate_models(a,b,z,blocks):
    ids=np.unique(blocks)
    labels=("null",*MODEL_AXES)
    loss={k:np.full(len(a),np.nan) for k in labels}
    pred={k:(np.full_like(a,np.nan),np.full_like(b,np.nan)) for k in labels}
    completed={k:0 for k in labels}
    coefs=[]
    for model,axes in MODEL_AXES.items():
        try:
            beta=fit_model(a,b,z,axes)
        except (ValueError,np.linalg.LinAlgError):
            coefs.append(dict(model=model,status="fit_failed"))
            continue
        for j,name in enumerate(NAMES):
            for q,term in enumerate(("directional","quadratic")):
                k=j+3*q
                coefs.append(dict(model=model,status="estimated",trait=name,term=term,
                     predictor="constant",coefficient=float(beta[k])))
        if len(axes):
            gamma=beta[6:].reshape(len(axes),len(FOCAL))
            for j,axis in enumerate(axes):
                for l,k in enumerate(FOCAL):
                    coefs.append(dict(model=model,status="estimated",trait=NAMES[k%3],
                        term="directional" if k<3 else "quadratic",
                        predictor=VARIABLES[axis],coefficient=float(gamma[j,l])))
    # Conditional discrimination, not prediction of absolute phenotype at an
    # unseen site: each test pair supplies its pooled phenotype histogram.
    # Only its intercept is profiled to balance total mass, with no test-label
    # fitting of slopes. Leave out the entire 20-degree block.
    for g in ids:
        train=blocks!=g;test=~train
        if len(np.unique(blocks[train]))<7 or train.sum()<20:
            continue
        m=a[test]+b[test]
        loss["null"][test]=2*np.log(2)
        pred["null"][0][test]=m/2;pred["null"][1][test]=m/2
        completed["null"]+=1
        for model,axes in MODEL_AXES.items():
            try:
                beta=fit_model(a[train],b[train],z[train],axes)
                ls,_,p=eta_profile(a[test],b[test],beta,z[test],axes)
            except (ValueError,np.linalg.LinAlgError):
                continue
            loss[model][test]=ls
            pred[model][0][test]=m*(1-p);pred[model][1][test]=m*p
            completed[model]+=1
    rows=[];profiles=[]
    obsa,obsb=marginals(a).mean(0),marginals(b).mean(0)
    for model in labels:
        row=dict(model=model,planned_blocks=len(ids),valid_blocks=completed[model],
                 status="estimated" if completed[model]==len(ids) else "incomplete_folds")
        if row["status"]=="estimated":
            ls=loss[model]
            null=loss["null"]
            byblock=lambda v:np.array([v[blocks==g].mean() for g in ids])
            row.update(mean_loss=float(ls.mean()),
                       equal_block_loss=float(byblock(ls).mean()),
                       skill_vs_null_pct=float(100*(null.mean()-ls.mean())/null.mean()))
            for comparator in ("constant","other_climate"):
                if completed[comparator]==len(ids):
                    ref=loss[comparator]
                    improvement=byblock(ref-ls)
                    row.update({f"skill_vs_{comparator}_pct":float(100*(ref.mean()-ls.mean())/ref.mean()),
                                f"equal_block_gain_vs_{comparator}":float(improvement.mean()),
                                f"blocks_better_than_{comparator}":int((improvement>1e-12).sum())})
            aa,bb=pred[model]
            np.testing.assert_allclose(aa.sum(1),1.,atol=1e-10)
            np.testing.assert_allclose(bb.sum(1),1.,atol=1e-10)
            pa,pb=marginals(aa).mean(0),marginals(bb).mean(0)
            for j,name in enumerate(NAMES):
                for k in range(5):
                    profiles.append(dict(model=model,trait=name,native_quintile=k+1,
                      observed_analogous=float(obsa[j,k]),observed_novel=float(obsb[j,k]),
                      predicted_analogous=float(pa[j,k]),predicted_novel=float(pb[j,k]),
                      observed_difference=float(obsb[j,k]-obsa[j,k]),
                      predicted_difference=float(pb[j,k]-pa[j,k])))
        rows.append(row)
    return rows,profiles,coefs

def analyze(sp,definition,n,c,pairs,fn,fi,out):
    nv,ni,nw=joint_records(fn,n,COLS,False)
    iv,ii,_=joint_records(fi,c,COLS,False)
    assert len(nv)==EXPECTED[sp][4]
    cuts,lo,hi=native_cuts(nv,nw)
    rc,valid=rank_codes(nv,cuts,lo,hi)
    assert valid.all()
    support=np.bincount(rc,minlength=K)>0
    assert support.sum()==125
    h,tot,kept=cell_hist(iv,ii,len(c),cuts,lo,hi,support)
    left=np.concatenate([l for l,r in pairs]);right=np.concatenate([r for l,r in pairs])
    blocks=c.iloc[left].block.to_numpy()
    np.testing.assert_array_equal(blocks,c.iloc[right].block.to_numpy())
    assert np.all(kept[np.r_[left,right]]>0)
    beta,_,_=fit(h[left],h[right],X,np.arange(6))
    np.testing.assert_allclose(beta,EXPECTED_BETA[(sp,definition)],atol=1e-6,rtol=0)
    mu,sd=reference_scale(nv,ni,len(n))
    moments,counts=cell_moments((iv-mu)/sd,ii,len(c))
    assert np.all(counts[np.r_[left,right]]>0)
    # Reproduce the old dispersion decomposition before subdividing the
    # between-cell component. Compare chroma and angle in fixed order.
    def old_components(index):
        v=moments[index].mean(0);a,q,b=np.split(v,3)
        return np.column_stack([b-a*a,q-b])[[NAMES.index("chroma"),NAMES.index("angle")]]
    np.testing.assert_allclose(old_components(right)-old_components(left),
        EXPECTED_DELTA[(sp,definition)],atol=1e-6,rtol=0)
    climate=climate_cells(fi,c)
    out[(sp,definition)]=dict(native_climate=climate_cells(fn,n),
        left=climate[left],right=climate[right],blocks=blocks,
        hist_left=h[left],hist_right=h[right],
        moment_left=moments[left],moment_right=moments[right])

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
    out={}
    for sp,rows in wild.groupby("taxon_name"):
        env,_=_background_species_audit(auth,land.dropna(subset=list(VARIABLES)),sp,land)
        complete=rows.dropna(subset=list(VARIABLES))
        sn=complete.loc[complete.range_state=="native"];si=complete.loc[complete.range_state=="introduced"]
        transform,edges,nat,inv=_pca_env_2d(env[env._range=="native"],env[env._range=="introduced"],sn,si,return_grid=True)
        assert abs(1-abs(nat-inv).sum()/2-frozen["species"][sp]["PCA_env_background_corrected_Schoener_D"])<1e-5
        run_species(sp,sn,si,transform,edges,nat,inv,out)
    variances=[];shifts=[];cv=[];predictions=[];coefficients=[];cohort=[]
    for definition in DEFINITIONS:
        ca=out[(SPECIES[0],definition)];cvd=out[(SPECIES[1],definition)]
        center,scale=common_scale(ca["native_climate"],cvd["native_climate"])
        for sp,d in zip(SPECIES,(ca,cvd)):
            z=(d["right"]-d["left"])/scale
            meta=dict(species=sp,definition=definition,pairs=len(z),
                      blocks=len(np.unique(d["blocks"])))
            cohort.append(dict(**meta,**{str(axis)+"_mean_delta":float(z[:,j].mean()) for j,axis in enumerate(VARIABLES)}))
            vv,ss=variance_diagnostic(d["moment_left"],d["moment_right"],d["blocks"])
            variances.extend(dict(**meta,**r) for r in vv)
            shifts.extend(dict(**meta,**r) for r in ss)
            fitrows,pred,coef=climate_models(d["hist_left"],d["hist_right"],z,d["blocks"])
            cv.extend(dict(**meta,**r) for r in fitrows)
            predictions.extend(dict(**meta,**r) for r in pred)
            coefficients.extend(dict(**meta,**r) for r in coef)
    pd.DataFrame(variances).to_csv(OUT/"variance_components.csv",index=False)
    pd.DataFrame(shifts).to_csv(OUT/"regional_shift_diagnostics.csv",index=False)
    pd.DataFrame(cv).to_csv(OUT/"heldout_climate_models.csv",index=False)
    pd.DataFrame(predictions).to_csv(OUT/"heldout_quintile_profiles.csv",index=False)
    pd.DataFrame(coefficients).to_csv(OUT/"shape_model_coefficients.csv",index=False)
    pd.DataFrame(cohort).to_csv(OUT/"cohort.csv",index=False)
    print("Complete; aggregate outputs only.",flush=True)

if __name__=="__main__":
    main()
