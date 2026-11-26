"""Geographic transfer of the frozen within-species enrichment surface.

The only new question: does the fitted relative enrichment surface generalize
to held-out introduced geographic blocks, and are coefficient signs stable
when one block is deleted? No new response, phenotype or climate variable.
The blocks are NOT known independent introduction events or populations.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from analysis.v5.phenotype_within_species_enrichment import (
    ROOT, NAMES, COLS, EXPECTED, joint_records, T, eligible, load_three_axis,
    make_cells, match_cells, head_cells, background_mesh, _background_species_audit,
    _pca_env_2d, _probability_content, apply_original_five_CHELSA,
    authority_source, angular_cells, VARIABLES, X, K,
    native_cuts, rank_codes, cell_hist, fit, profile_state,
)
OUT=ROOT/"phenotype_enrichment_transfer"
EXPECTED_BETA={
    ("Cirsium arvense","native_core95"):[.078932,-.207146,.097324,-.195635,.167441,-.032402],
    ("Cirsium arvense","density_surplus"):[.010634,-.013189,-.049352,-.002273,.041548,-.002733],
    ("Cirsium vulgare","native_core95"):[.269631,.042021,-.054898,-.115440,.126021,-.006817],
    ("Cirsium vulgare","density_surplus"):[.085949,-.079453,-.010954,.009071,.032332,-.011042],
}
# Four fixed evaluations of the same model. No selection of a best subset.
MODELS={"full":np.arange(6),"linear_only":np.arange(3),
        "drop_chroma":np.array([1,2,4,5]),"drop_angle":np.array([0,2,3,5])}

def heldout_loss(a,b,beta):
    # Balanced within-pair discrimination: condition on pooled phenotype mass
    # (a+b) and equal analogue/novel weights. Only the intercept is profiled in
    # the held-out pair; no held-out label-specific phenotype slope is fitted.
    f,_,_=profile_state(a,b,X,beta)
    return f

def cross_validate(a,b,blocks):
    ids=np.unique(blocks)
    base,_,_=fit(a,b,X,np.arange(6))
    betas=[];scores={k:[] for k in MODELS}
    for g in ids:
        train=blocks!=g
        test=~train
        # Freeze model; all folds retained. If a fit fails, mark that
        # specification unsupported, not a favorable subset of folds.
        full_beta=np.full(6,np.nan)
        for name,active in MODELS.items():
            try:
                beta,_,_=fit(a[train],b[train],X,active)
                if name=="full":full_beta=beta.copy()
                loss=heldout_loss(a[test],b[test],beta)
                null_loss=2*int(test.sum())*np.log(2.)
                scores[name].append((int(test.sum()),loss,null_loss))
            except (ValueError,np.linalg.LinAlgError):
                scores[name].append((int(test.sum()),np.nan,2*int(test.sum())*np.log(2.)))
        betas.append(full_beta)
    return base,np.stack(betas),scores

def summarize_transfer(scores):
    result={}
    full=np.array(scores["full"],float)
    for name,rows in scores.items():
        m=np.array(rows,float)
        n=len(m);valid=np.isfinite(m[:,1])
        rec=dict(model=name,blocks=n,valid_blocks=int(valid.sum()),
                 pairs=int(m[:,0].sum()),status="estimated" if valid.all() else "fit_failed")
        if valid.all():
            gain=m[:,2]-m[:,1]
            rec.update(null_loss=float(m[:,2].sum()),heldout_loss=float(m[:,1].sum()),
                      loss_reduction_pct=float(100*gain.sum()/m[:,2].sum()),
                      equal_block_reduction_pct=float(100*np.mean(gain/m[:,2])),
                      positive_blocks=int((gain>1e-12).sum()),
                      worst_block_reduction_pct=float(100*np.min(gain/m[:,2])),
                      best_block_reduction_pct=float(100*np.max(gain/m[:,2])))
            if np.isfinite(full[:,1]).all():
                delta=m[:,1]-full[:,1]
                rec.update(full_advantage_pct_of_null=float(100*delta.sum()/m[:,2].sum()),
                           blocks_full_better=int((delta>1e-12).sum()))
        result[name]=rec
    return result

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
    a,b=h[left],h[right]
    base,betas,scores=cross_validate(a,b,blocks)
    np.testing.assert_allclose(base,EXPECTED_BETA[(sp,definition)],rtol=0,atol=1e-6)
    meta=dict(species=sp,definition=definition,pairs=len(left),blocks=len(pairs))
    valid=np.isfinite(betas).all(1)
    for j,beta in enumerate(base):
        v=betas[valid,j]
        rec=dict(**meta,trait=NAMES[j%3],term="directional" if j<3 else "quadratic",
                 full_coefficient=float(beta),valid_folds=int(valid.sum()),
                 status="estimated" if valid.all() else "fit_failed")
        if valid.all():
            rec.update(loo_min=float(v.min()),loo_max=float(v.max()),
                       same_sign_folds=int((v*beta>0).sum()),
                       opposite_sign_folds=int((v*beta<0).sum()),
                       max_abs_change=float(np.max(np.abs(v-beta))),
                       loo_ratio_min=float(np.exp(4*v.min())),
                       loo_ratio_max=float(np.exp(4*v.max())))
        out["coefficients"].append(rec)
    for rec in summarize_transfer(scores).values():
        # rec already supplies pairs and blocks.
        out["transfer"].append(dict(species=sp,definition=definition,**rec))

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
    out={"coefficients":[],"transfer":[]}
    for sp,rows in wild.groupby("taxon_name"):
        env,_=_background_species_audit(auth,land.dropna(subset=list(VARIABLES)),sp,land)
        complete=rows.dropna(subset=list(VARIABLES))
        sn=complete.loc[complete.range_state=="native"];si=complete.loc[complete.range_state=="introduced"]
        transform,edges,nat,inv=_pca_env_2d(env[env._range=="native"],env[env._range=="introduced"],sn,si,return_grid=True)
        assert abs(1-abs(nat-inv).sum()/2-frozen["species"][sp]["PCA_env_background_corrected_Schoener_D"])<1e-5
        run_species(sp,sn,si,transform,edges,nat,inv,out)
    d=pd.DataFrame(out["coefficients"])
    assert len(d)==24
    d.to_csv(OUT/"leave_region_out_coefficients.csv",index=False)
    g=pd.DataFrame(out["transfer"])
    assert len(g)==16
    g.to_csv(OUT/"heldout_region_transfer.csv",index=False)
    print("Complete; aggregate outputs only.",flush=True)

if __name__=="__main__":
    main()
