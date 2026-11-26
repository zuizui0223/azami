"""Within-species directional and nonlinear enrichment of standing phenotypes.

A matched case-control density-ratio model of phenotype observations in
introduced novel vs analogous climates. These are relative enrichment
coefficients, NOT fitness, selection gradients, or evolutionary rates.
No predictions transferred from interspecific Result 2, no time axis.
"""
from pathlib import Path
from itertools import product
import json
import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.stats import t as student_t
from statsmodels.stats.multitest import multipletests
from analysis.v5.phenotype_three_axis_dispersion import (
    ROOT, eligible, load_three_axis, NAMES, COLS, EXPECTED, joint_records, T,
    make_cells, match_cells, head_cells, background_mesh, _background_species_audit,
    _pca_env_2d, _probability_content, apply_original_five_CHELSA,
    authority_source, angular_cells, VARIABLES,
)
from analysis.v5.phenotype_redistribution import weighted_quantile
OUT=ROOT/"phenotype_within_species_enrichment"
RANKS=np.array(list(product(range(-2,3),repeat=3)),float)
# Other two phenotype ranks held fixed. Curvature centered for numerical stability.
X=np.column_stack([RANKS,RANKS**2-2])
K=len(X)

def native_cuts(v,w):
    cuts=np.column_stack([weighted_quantile(v[:,j],w,[.2,.4,.6,.8]) for j in range(3)]).T
    if not np.all(np.diff(cuts,axis=1)>0):
        raise ValueError("Degenerate native quintiles")
    return cuts,v.min(0),v.max(0)

def rank_codes(v,cuts,lo,hi):
    v=np.asarray(v,float)
    rank=np.stack([np.searchsorted(cuts[j],v[:,j],side="right") for j in range(3)],axis=1)
    code=rank@np.array([25,5,1])
    valid=np.isfinite(v).all(1)&np.all((v>=lo)&(v<=hi),axis=1)
    return code,valid

def cell_hist(v,ix,n,cuts,lo,hi,support):
    code,valid=rank_codes(v,cuts,lo,hi)
    valid&=support[code]
    total=np.bincount(ix,minlength=n)
    kept=np.bincount(ix[valid],minlength=n)
    h=np.bincount(ix[valid]*K+code[valid],minlength=n*K).reshape(n,K).astype(float)
    h=np.divide(h,kept[:,None],out=np.full_like(h,np.nan),where=kept[:,None]>0)
    return h,total,kept

def profile_state(a,b,x,beta):
    """Profile one intercept for each matched pair, balanced cell weights.

    In a cell each retained observation has weight 1/n_cell. Both cells in
    every pair have total weight 1. The pair-specific intercept is solved
    to make predicted novel mass exactly 1. A common log density-ratio
    surface across pairs is then estimated from their within-pair contrasts.
    """
    m=a+b
    eta=x@beta
    lower=np.full(len(a),-np.max(eta)-40.)
    upper=np.full(len(a),-np.min(eta)+40.)
    for _ in range(48):
        mid=(lower+upper)/2
        p=expit(eta[None,:]+mid[:,None])
        bigger=(m*p).sum(1)>1.
        upper=np.where(bigger,mid,upper)
        lower=np.where(bigger,lower,mid)
    z=eta[None,:]+((lower+upper)/2)[:,None]
    p=expit(z)
    np.testing.assert_allclose((m*p).sum(1),1,atol=1e-10)
    loss=float(np.sum(m*np.logaddexp(0,z)-b*z))
    score=(b-m*p)@x
    w=m*p*(1-p)
    wx=w@x
    # Schur complement removes pair-specific intercepts from information.
    h=(x.T*w.sum(0))@x-(wx.T/w.sum(1))@wx
    return loss,score,h

def fit(a,b,x,active):
    active=np.asarray(active,int)
    beta=np.zeros(x.shape[1])
    for _ in range(50):
        loss,s,h=profile_state(a,b,x,beta)
        g=s.sum(0)[active]
        if not len(active) or np.max(np.abs(g))<1e-7:
            return beta,s,h
        hh=h[np.ix_(active,active)]
        if np.linalg.eigvalsh(hh).min()<1e-9:
            raise ValueError("Singular profile information")
        step=np.linalg.solve(hh,g)
        accepted=False
        for k in range(24):
            new=beta.copy()
            new[active]+=step/2**k
            if np.max(np.abs(new))>12:
                continue
            f,ss,hh2=profile_state(a,b,x,new)
            if f<=loss+1e-10:
                beta=new
                accepted=True
                break
        if not accepted:
            raise ValueError("Profile Newton line search failed")
    raise ValueError("Profile fit did not converge")

def group_sums(s,blocks):
    return np.stack([s[blocks==g].sum(0) for g in np.unique(blocks)])

def sign_patterns(blocks):
    b=len(np.unique(blocks))
    if not 8<=b<=15:
        raise ValueError("Frozen regional support changed")
    return np.array(list(product([1.,-1.],repeat=b)))

def score_omnibus(a,b,x,blocks):
    # Exact whole-block label swap at the no-enrichment null, with a
    # sign-invariant quadratic score statistic. This tests distributional
    # association captured by these fixed linear and quadratic features.
    scores=group_sums(.5*(b-a)@x,blocks)
    h=scores.T@scores
    rank=int(np.linalg.matrix_rank(h))
    if rank==0:
        return dict(statistic=0.,p=1.,score_rank=0,patterns=len(sign_patterns(blocks)))
    s=sign_patterns(blocks)@scores
    value=np.einsum("ij,jk,ik->i",s,np.linalg.pinv(h),s)
    return dict(statistic=float(value[0]),p=float(np.mean(value>=value[0]-1e-12)),
                score_rank=rank,patterns=len(value))

def fit_coefficients(a,b,x,blocks):
    np.testing.assert_allclose(a.sum(1),1)
    np.testing.assert_allclose(b.sum(1),1)
    p=x.shape[1]
    beta,s,h=fit(a,b,x,np.arange(p))
    nblocks=len(np.unique(blocks))
    bs=group_sums(s,blocks)
    bread=np.linalg.inv(h)
    cov=bread@(bs.T@bs)@bread*nblocks/(nblocks-1)
    se=np.sqrt(np.maximum(np.diag(cov),0))
    crit=student_t.ppf(.975,nblocks-1)
    flips=sign_patterns(blocks)
    rows=[]
    for j in range(p):
        nuisance=np.delete(np.arange(p),j)
        _,ss,hh=fit(a,b,x,nuisance)
        # Approximate conditional score-sign test with nuisance fitted under
        # the one-coefficient null. Orthogonalize nuisance score contributions
        # using the profiled information. Not an exact randomized experiment.
        q=ss[:,j]-ss[:,nuisance]@np.linalg.solve(hh[np.ix_(nuisance,nuisance)],hh[nuisance,j])
        q=group_sums(q[:,None],blocks).ravel()
        null=flips@q
        pv=float(np.mean(np.abs(null)>=abs(null[0])-1e-12))
        lo,hi=beta[j]-crit*se[j],beta[j]+crit*se[j]
        rows.append(dict(trait=NAMES[j%3],term="directional" if j<3 else "quadratic",
            coefficient=float(beta[j]),cluster_se=float(se[j]),lo=float(lo),hi=float(hi),
            p_score=pv,contrast_ratio=float(np.exp(4*beta[j])),
            contrast_ratio_lo=float(np.exp(4*lo)),contrast_ratio_hi=float(np.exp(4*hi)),
            positive_score_blocks=int((q>0).sum()),score_patterns=len(null),status="estimated"))
    return beta,rows

def analyze(sp,definition,n,c,pairs,fn,fi,out):
    nv,ni,nw=joint_records(fn,n,COLS,False)
    iv,ii,_=joint_records(fi,c,COLS,False)
    assert len(nv)==EXPECTED[sp][4]
    cuts,lo,hi=native_cuts(nv,nw)
    rc,valid=rank_codes(nv,cuts,lo,hi)
    assert valid.all()
    support=np.bincount(rc,minlength=K)>0
    # All 125 interior combinations were occupied in the preceding full
    # native reference. Re-check rather than assume.
    assert support.sum()==125
    nh,nt,nk=cell_hist(nv,ni,len(n),cuts,lo,hi,support)
    h,tot,kept=cell_hist(iv,ii,len(c),cuts,lo,hi,support)
    left=np.concatenate([l for l,r in pairs]);right=np.concatenate([r for l,r in pairs])
    blocks=c.iloc[left].block.to_numpy()
    np.testing.assert_array_equal(blocks,c.iloc[right].block.to_numpy())
    assert np.all(kept[np.r_[left,right]]>0)
    assert np.all(nk==nt) and np.isfinite(nh).all()
    a,b=h[left],h[right]
    meta=dict(species=sp,definition=definition,pairs=len(left),blocks=len(pairs),
              native_cells=len(n),native_records=len(nv),native_observed_combinations=int(support.sum()),
              analogous_records=int(tot[left].sum()),novel_records=int(tot[right].sum()),
              analogous_excluded_fraction=float(np.mean(1-kept[left]/tot[left])),
              novel_excluded_fraction=float(np.mean(1-kept[right]/tot[right])))
    out["omnibus"].append(dict(**meta,**score_omnibus(a,b,X,blocks)))
    # Do not tune a penalty, use a different trait scale, or select a model
    # if numerical identification fails. Preserve a row for each planned test.
    try:
        beta,rows=fit_coefficients(a,b,X,blocks)
    except (ValueError,np.linalg.LinAlgError) as exc:
        beta=np.full(6,np.nan)
        rows=[dict(trait=t,term=k,status="fit_failed",failure_reason=str(exc),p_score=1.)
              for k in ("directional","quadratic") for t in NAMES]
    out["coefficients"].extend(dict(**meta,**r) for r in rows)
    for j,trait in enumerate(NAMES):
        for rank in range(-2,3):
            idx=RANKS[:,j]==rank
            # Conditional fitted enrichment relative to this trait's middle
            # bin, holding the other two trait ranks fixed.
            z=beta[j]*rank+beta[3+j]*rank**2
            out["profiles"].append(dict(**meta,trait=trait,native_quintile=rank+3,
                    native_frequency=float(nh[:,idx].sum(1).mean()),
                    analogous_frequency=float(a[:,idx].sum(1).mean()),
                    novel_frequency=float(b[:,idx].sum(1).mean()),
                    fitted_relative_enrichment_vs_middle=float(np.exp(z)) if np.isfinite(z) else np.nan,
                    native_bin_lower=float(lo[j] if rank==-2 else cuts[j,rank+1]),
                    native_bin_upper=float(hi[j] if rank==2 else cuts[j,rank+2])))

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
    out={"coefficients":[],"omnibus":[],"profiles":[]}
    for sp,rows in wild.groupby("taxon_name"):
        env,_=_background_species_audit(auth,land.dropna(subset=list(VARIABLES)),sp,land)
        complete=rows.dropna(subset=list(VARIABLES))
        sn=complete.loc[complete.range_state=="native"];si=complete.loc[complete.range_state=="introduced"]
        transform,edges,nat,inv=_pca_env_2d(env[env._range=="native"],env[env._range=="introduced"],sn,si,return_grid=True)
        assert abs(1-abs(nat-inv).sum()/2-frozen["species"][sp]["PCA_env_background_corrected_Schoener_D"])<1e-5
        run_species(sp,sn,si,transform,edges,nat,inv,out)
    d=pd.DataFrame(out["coefficients"])
    assert len(d)==24
    d["holm_all24"]=multipletests(d.p_score.fillna(1),method="holm")[1]
    d["holm_family12"]=d.groupby("definition").p_score.transform(lambda p:multipletests(p.fillna(1),method="holm")[1])
    d.to_csv(OUT/"enrichment_coefficients.csv",index=False)
    g=pd.DataFrame(out["omnibus"])
    assert len(g)==4
    g["holm_all4"]=multipletests(g.p.fillna(1),method="holm")[1]
    g["holm_family2"]=g.groupby("definition").p.transform(lambda p:multipletests(p.fillna(1),method="holm")[1])
    g.to_csv(OUT/"enrichment_omnibus.csv",index=False)
    q=pd.DataFrame(out["profiles"])
    assert len(q)==60
    q.to_csv(OUT/"enrichment_profiles.csv",index=False)
    print("Complete; aggregate outputs only.",flush=True)

if __name__=="__main__":
    main()
