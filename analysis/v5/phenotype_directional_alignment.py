"""Two fixed native climate-phenotype direction hypotheses in introduced ranges.

Cross-sectional; use only the three original traits. No new dates, models of
ancestry, inferred fitness, or genomic selection. Directional hypotheses come
from existing result 2, but this is a focused exploratory reanalysis, not
independent confirmation. Keep all old tests and their corrections.
"""
from pathlib import Path
from itertools import product
import json
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from analysis.v5.phenotype_three_axis_dispersion import (
    ROOT, SEED, eligible, load_three_axis, NAMES, COLS, EXPECTED,
    joint_records, T, make_cells, match_cells, head_cells, background_mesh,
    _background_species_audit, _pca_env_2d, _probability_content,
    apply_original_five_CHELSA, authority_source, angular_cells, VARIABLES,
)
from analysis.v5.phenotype_redistribution import cellkeys, weighted_quantile
OUT=ROOT/"phenotype_directional_alignment"
PREDICTIONS={"angle":"chelsa_bio12","chroma":"chelsa_bio15"}

def existing_tails(v,ix,n,reference,w):
    """Same full native range + weighted q20/q80 thresholds for all groups.

    Only values within observed native extrema can contribute to existing tails.
    Between-threshold values and out-of-range values do not enter the tail score.
    All observations have equal weights within each geographic cell.
    """
    v=np.asarray(v,float);reference=np.asarray(reference,float)
    low,high=weighted_quantile(reference,w,[.2,.8])
    lo,hi=reference.min(),reference.max()
    if not lo<=low<high<=hi:
        raise ValueError("Degenerate native thresholds")
    flags=np.column_stack([(v>=lo)&(v<low),(v>=high)&(v<=hi),(v<lo)|(v>hi)])
    count=np.bincount(ix,minlength=n)
    h=np.stack([np.bincount(ix,weights=flags[:,j],minlength=n) for j in range(3)],axis=1)
    h=np.divide(h,count[:,None],out=np.full_like(h,np.nan),where=count[:,None]>0)
    assert np.nanmax(h.sum(1))<=1+1e-10
    return h,(float(lo),float(low),float(high),float(hi))

def exact_block_test(num,den):
    """Whole-block symmetric score sign null, conditional on fixed climate.

    For alignment this flips phenotype contrasts, not climate differences.
    For the within-block diagnostic it flips within-region covariance scores.
    Not a randomized experiment or a selection null.
    """
    num=np.asarray(num,float);den=np.asarray(den,float)
    b=len(num)
    if b<8 or b>15 or not np.isfinite(num).all() or not np.isfinite(den).all() or den.sum()<=1e-12:
        return dict(status="insufficient_support",effect=np.nan,p=1.,blocks=b,
                    patterns=0,lo=np.nan,hi=np.nan,positive_blocks=int((num>0).sum()))
    signs=np.array(list(product([1.,-1.],repeat=b)))
    effect=num.sum()/den.sum()
    null=signs@num/den.sum()
    p=float(np.mean(np.abs(null)>=abs(effect)-1e-12))
    rng=np.random.default_rng(SEED)
    boots=rng.multinomial(b,np.full(b,1/b),size=3999)
    bd=boots@den
    vals=(boots@num)[bd>1e-12]/bd[bd>1e-12]
    lo,hi=np.quantile(vals,[.025,.975])
    return dict(status="estimated",effect=float(effect),p=p,blocks=b,
                patterns=len(signs),lo=float(lo),hi=float(hi),positive_blocks=int((num>0).sum()))

def alignment(dx,dy,block):
    """Climate-direction weighted difference of predicted vs opposite tails.

    dy[:,0]=change in low native tail, dy[:,1]=change in high native tail.
    For wetter or more seasonal differences predict the high tail; reverse the
    prediction for negative climate differences. Weight by |climate difference|
    so negligible climate differences are not counted as strong direction tests.
    """
    dx=np.asarray(dx,float);dy=np.asarray(dy,float);block=np.asarray(block)
    assert dy.shape==(len(dx),2) and block.shape==dx.shape
    assert np.isfinite(dx).all() and np.isfinite(dy).all()
    w=np.abs(dx)
    # Climate differences in native climate SD units. A predicted tail is a
    # phenotype associated with the climate direction, not a measured fitness benefit.
    pred=np.where(dx>0,dy[:,1],dy[:,0])
    opp=np.where(dx>0,dy[:,0],dy[:,1])
    score=w*(pred-opp)
    np.testing.assert_allclose(score,dx*(dy[:,1]-dy[:,0]),atol=1e-12)
    ids=np.unique(block)
    active=[g for g in ids if w[block==g].sum()>1e-12]
    a=exact_block_test([score[block==g].sum() for g in active],
                       [w[block==g].sum() for g in active])
    a["predicted_tail_change"]=float(w@pred/w.sum()) if w.sum()>0 else np.nan
    a["opposite_tail_change"]=float(w@opp/w.sum()) if w.sum()>0 else np.nan
    # Guard against a global shift of all novel cells and geographic group
    # structure. Does variation among pairs within the SAME block follow the
    # gradient? A block with only one pair supplies no within-block slope.
    num=[];den=[];informative_pairs=0
    y=dy[:,1]-dy[:,0]
    for g in ids:
        k=block==g
        x=dx[k]-dx[k].mean()
        z=y[k]-y[k].mean()
        s=x@x
        if s<=1e-12:
            continue
        num.append(x@z);den.append(s)
        informative_pairs+=int(k.sum())
    b=exact_block_test(num,den)
    b["informative_pairs"]=informative_pairs
    return a,b

def climate(frame,cells):
    d=frame.groupby("cellkey")[list(VARIABLES)].median().reindex(cellkeys(cells))
    if not np.isfinite(d.to_numpy()).all():
        raise ValueError("Missing climate values in frozen cells")
    return d

def analyze(sp,definition,n,c,pairs,fn,fi,out):
    nv,ni,nw=joint_records(fn,n,COLS,False)
    iv,ii,_=joint_records(fi,c,COLS,False)
    assert len(nv)==EXPECTED[sp][4]
    left=np.concatenate([l for l,r in pairs])
    right=np.concatenate([r for l,r in pairs])
    b=c.iloc[left].block.to_numpy()
    np.testing.assert_array_equal(b,c.iloc[right].block.to_numpy())
    cn,ci=climate(fn,n),climate(fi,c)
    for j,trait in enumerate(NAMES):
        ref,th=existing_tails(nv[:,j],ni,len(n),nv[:,j],nw)
        h,th2=existing_tails(iv[:,j],ii,len(c),nv[:,j],nw)
        assert th==th2 and np.isfinite(ref).all() and np.isfinite(h[np.r_[left,right]]).all()
        a,z,r=h[left].mean(0),h[right].mean(0),ref.mean(0)
        row=dict(species=sp,definition=definition,trait=trait,native_cells=len(n),
                 pairs=len(left),blocks=len(pairs),native_min=th[0],native_q20=th[1],
                 native_q80=th[2],native_max=th[3])
        for label,arr in [("native",r),("analogous",a),("novel",z),("difference",z-a)]:
            for i,name in enumerate(["low","high","outside"]):
                row[label+"_"+name]=float(arr[i])
        out["frequencies"].append(row)
        # The outline index has no established signed climate prediction.
        # Report its descriptive frequencies; do not select a climate axis by p.
        if trait not in PREDICTIONS:
            continue
        axis=PREDICTIONS[trait]
        sd=float(cn[axis].std(ddof=0))
        assert sd>0
        raw=ci[axis].to_numpy()[right]-ci[axis].to_numpy()[left]
        dx=raw/sd
        aa,bb=alignment(dx,h[right,:2]-h[left,:2],b)
        info=dict(species=sp,definition=definition,trait=trait,climate_axis=axis,
                  pairs=len(left),blocks=len(pairs),native_climate_sd=sd,
                  mean_climate_difference=float(raw.mean()),
                  median_climate_difference=float(np.median(raw)),
                  positive_climate_pairs=int((dx>0).sum()),negative_climate_pairs=int((dx<0).sum()),
                  zero_climate_pairs=int((dx==0).sum()))
        for suffix,d in [("alignment",aa),("within_block",bb)]:
            info.update({k+"_"+suffix:v for k,v in d.items()})
        out["tests"].append(info)

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
    out={"frequencies":[],"tests":[]}
    for sp,rows in wild.groupby("taxon_name"):
        env,_=_background_species_audit(auth,land.dropna(subset=list(VARIABLES)),sp,land)
        complete=rows.dropna(subset=list(VARIABLES))
        sn=complete.loc[complete.range_state=="native"];si=complete.loc[complete.range_state=="introduced"]
        transform,edges,nat,inv=_pca_env_2d(env[env._range=="native"],env[env._range=="introduced"],sn,si,return_grid=True)
        assert abs(1-abs(nat-inv).sum()/2-frozen["species"][sp]["PCA_env_background_corrected_Schoener_D"])<1e-5
        run_species(sp,sn,si,transform,edges,nat,inv,out)
    d=pd.DataFrame(out["tests"])
    assert len(d)==8
    for p in ["p_alignment","p_within_block"]:
        d["holm8_"+p]=multipletests(d[p].fillna(1),method="holm")[1]
        d["holm4_"+p]=d.groupby("definition")[p].transform(lambda v:multipletests(v.fillna(1),method="holm")[1])
    d.to_csv(OUT/"directional_alignment.csv",index=False)
    g=pd.DataFrame(out["frequencies"])
    assert len(g)==12
    g.to_csv(OUT/"existing_tail_frequencies.csv",index=False)
    print("Complete; aggregate outputs only.",flush=True)

if __name__=="__main__":
    main()
