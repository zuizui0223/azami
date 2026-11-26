"""Direct existing-bin gains versus native-unobserved-bin gains.

Cross-sectional only. Native reference fixed independently of introduced
phenotypes. Full reference and equal-size geographic native holdouts.
No dates, prediction models, inferred genes, ancestry or fitness.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from analysis.v5.phenotype_redistribution import ROOT, SEED, eligible, weighted_quantile
from analysis.v5.phenotype_endpoint_scan import load_measures, permute
from analysis.v5.phenotype_native_calibration import native_splits
from analysis.v5.phenotype_reservoir import records
from analysis.v5.novel_head_overlap import T, make_cells, match_cells
from analysis.v5.invasion_partition import head_cells
from analysis.v5.native_pca_env import background_mesh, _background_species_audit, _pca_env_2d, _probability_content
from analysis.v5.native_pca_env_fullcoverage import apply_original_five_CHELSA
from analysis.v4.native_introduced_opportunity import authority_source, angular_cells
from analysis.v4.native_cached_climate import VARIABLES
from analysis.v4.native_integration import FEATURES

OUT=ROOT/"phenotype_existing_vs_novel"
SPACES={**FEATURES,
        "colour_joint":FEATURES["chroma"]+FEATURES["lightness"],
        "colour_angle_joint":FEATURES["chroma"]+FEATURES["lightness"]+FEATURES["angle"]}
REFS=["full_native","rarefied_native"]
DEFS=["native_core95","density_surplus"]
N_REF,N_CAL,N_DRAWS=256,128,32

def codes(ref,w,x,circular=False):
    """Five native quantile intervals plus explicit below/above observed bounds.
    For hue use 12 fixed circular sectors; no artificial linear endpoints.
    The same grid is used for reference and evaluation records.
    """
    if circular:
        a=np.mod(np.arctan2(x[:,0],x[:,1]),2*np.pi)
        return np.minimum(np.floor(a/(2*np.pi)*12).astype(int),11),12
    c=np.zeros(len(x),int)
    for j in range(ref.shape[1]):
        q=weighted_quantile(ref[:,j],w,[.2,.4,.6,.8])
        lo,hi=ref[:,j].min(),ref[:,j].max()
        b=1+np.searchsorted(q,x[:,j],side="right")
        b[x[:,j]<lo]=0
        b[x[:,j]>hi]=6
        c=c*7+b
    return c,7**ref.shape[1]

def hist_by_cell(code,idx,n,k):
    h=np.bincount(idx*k+code,minlength=n*k).reshape(n,k).astype(float)
    counts=h.sum(1)
    result=np.full_like(h,np.nan)
    np.divide(h,counts[:,None],out=result,where=counts[:,None]>0)
    return result

def gains(a,b,support):
    positive=np.maximum(b-a,0)
    old,new=positive[support].sum(),positive[~support].sum()
    assert np.isclose(old+new,.5*np.abs(b-a).sum(),atol=1e-10)
    return float(old),float(new)

def joint_records(frame,c,cols,circular):
    v,ix,w=records(frame,c,cols)
    if circular:
        ok=np.linalg.norm(v,axis=1)>1e-8
        v,ix=v[ok],ix[ok]
        cnt=np.bincount(ix,minlength=len(c));w=1/cnt[ix]
    return v,ix,w

def one_per_cell(v,idx,cells,rng,lookup):
    return v[np.array([rng.choice(lookup[int(c)]) for c in cells],dtype=int)]

def evaluate(ref,rw,cal,iv,ii,n,left,right,circular):
    rc,k=codes(ref,rw,ref,circular)
    nc,_=codes(ref,rw,iv,circular)
    ref_freq=np.bincount(rc,weights=rw,minlength=k)/rw.sum()
    support=ref_freq>0
    hist=hist_by_cell(nc,ii,n,k)
    assert np.isfinite(hist[np.r_[left,right]]).all()
    u,v=hist[left].mean(0),hist[right].mean(0)
    old,new=gains(u,v,support)
    old_nat,new_nat=gains(ref_freq,v,support)
    # Partition both total positive redistribution and signed outside-bin mass.
    novel_score=hist[:,~support].sum(1)
    novel_score[~np.isfinite(hist).all(1)]=np.nan
    cal_mass=np.nan
    if cal is not None:
        cc,_=codes(ref,rw,cal,circular)
        cal_mass=float((~support[cc]).mean())
    return dict(native_calibration_unobserved=cal_mass,
        analogue_unobserved=float(u[~support].sum()),novel_unobserved=float(v[~support].sum()),
        observed_bins=int(support.sum()),possible_bins=k,existing_gain=old,new_gain=new,
        total_redistribution=old+new,existing_gain_from_native=old_nat,new_gain_from_native=new_nat,
        ref_records=len(ref)),novel_score

def analyze(sp,definition,space,n,c,pairs,fn,fi,out):
    cols=["scan_"+x for x in SPACES[space]];circular=space=="hue"
    nv,ni,nw=joint_records(fn,n,cols,circular)
    iv,ii,iw=joint_records(fi,c,cols,circular)
    good_native=np.unique(ni);good_intro=np.bincount(ii,minlength=len(c))>0
    total=sum(len(l) for l,r in pairs)
    matched=[(l[good_intro[l]&good_intro[r]],r[good_intro[l]&good_intro[r]]) for l,r in pairs]
    matched=[(l,r) for l,r in matched if len(l)]
    count=sum(len(l) for l,r in matched)
    meta=dict(species=sp,definition=definition,space=space,native_cells=len(good_native),
              native_records=len(nv),fixed_pairs=total,pairs=count,blocks=len(matched))
    if len(good_native)<N_REF+N_CAL or count<20 or count<.8*total or len(matched)<8:
        for ref in REFS:out.append(dict(**meta,reference=ref,status="insufficient_support",p=1.))
        return
    left=np.concatenate([l for l,r in matched]);right=np.concatenate([r for l,r in matched])
    ns=np.array([len(l) for l,r in matched])
    # Equal reference count across both species and all endpoints, one archived
    # observation per cell; geographically distinct reference/calibration sets.
    splits=native_splits(n.block.to_numpy()[good_native],nref=N_REF,ncal=N_CAL,draws=N_DRAWS)
    lookup={int(k):np.flatnonzero(ni==k) for k in good_native}
    rng=np.random.default_rng(SEED+1)
    for mode in REFS:
        summaries=[];scores=[]
        if mode=="full_native":
            m,y=evaluate(nv,nw,None,iv,ii,len(c),left,right,circular)
            summaries.append(m);scores.append(y)
        else:
            for tr,ca in splits:
                ref=one_per_cell(nv,ni,good_native[tr],rng,lookup)
                cal=one_per_cell(nv,ni,good_native[ca],rng,lookup)
                m,y=evaluate(ref,np.ones(len(ref)),cal,iv,ii,len(c),left,right,circular)
                summaries.append(m);scores.append(y)
        avg=pd.DataFrame(summaries).mean().to_dict()
        # Only paired cells enter inference. Avoid averaging NaNs for other cells.
        v=np.asarray(scores)
        block_diff=np.array([(v[:,r]-v[:,l]).mean(0).sum() for l,r in matched])
        delta=float(block_diff.sum()/count)
        np.testing.assert_allclose(delta,avg["novel_unobserved"]-avg["analogue_unobserved"],atol=1e-10)
        _,p,exact,nperm=permute(block_diff[:,None],ns)
        rng_ci=np.random.default_rng(SEED)
        boot=rng_ci.multinomial(len(ns),np.full(len(ns),1/len(ns)),size=3999)
        bd=(boot@block_diff)/(boot@ns)
        lo,hi=np.quantile(bd,[.025,.975])
        gain=avg["total_redistribution"]
        out.append(dict(**meta,reference=mode,status="estimated",draws=len(summaries),**avg,
            share_existing_gain=avg["existing_gain"]/gain if gain>0 else np.nan,
            novelty_difference=delta,lo=float(lo),hi=float(hi),p=p,exact=exact,permutations=nperm))

def run_species(sp,sn,si,transform,edges,nat,inv,out):
    fn,fi=eligible(sn),eligible(si)
    n=make_cells(sn,True);full=make_cells(si,True)
    expected={"Cirsium arvense":(1366,860,30,183),"Cirsium vulgare":(1113,1139,146,316)}[sp]
    assert (len(n),len(full))==expected[:2]
    coords=transform(si)
    uv=np.column_stack([np.searchsorted(e,coords[:,k],side="right")-1 for k,e in enumerate(edges)])
    for k,e in enumerate(edges):uv[coords[:,k]==e[-1],k]=len(e)-2
    inside=np.all((uv>=0)&(uv<np.array([len(e)-1 for e in edges])),axis=1)
    labeled=si.loc[inside].copy()
    labeled["novel"]=(~_probability_content(nat,.95)[uv[inside,0],uv[inside,1]]).astype(int)
    c=make_cells(labeled);c=c.loc[c.novel!=.5].reset_index(drop=True)
    pairs=match_cells(c);assert sum(len(l) for l,r in pairs)==expected[2]
    for space in SPACES:analyze(sp,"native_core95",space,n,c,pairs,fn,fi,out)
    hc=head_cells(si,transform,edges)[["i","j","u","v"]]
    c=full.merge(hc,on=["i","j"],validate="one_to_one")
    u,v=c.u.to_numpy(),c.v.to_numpy()
    c=c.loc[inv[u,v]>0].reset_index(drop=True);u,v=c.u.to_numpy(),c.v.to_numpy()
    c["novel"]=(inv[u,v]>nat[u,v]).astype(int)
    pairs=match_cells(c);assert sum(len(l) for l,r in pairs)==expected[3]
    for space in SPACES:analyze(sp,"density_surplus",space,n,c,pairs,fn,fi,out)

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    auth=authority_source((ROOT/"private/authority350.zip").read_bytes())
    wild=pd.read_csv(ROOT/"private/focal_wild.csv.gz",dtype={"obs_id":str},low_memory=False)
    assert len(wild)==297622
    wild["_cell"]=angular_cells(wild.latitude,wild.longitude)
    wild=wild.sort_values("obs_id").drop_duplicates(["taxon_name","range_state","_cell","user_id"]).copy()
    land=background_mesh(auth)
    wild,land,_=apply_original_five_CHELSA(wild,land,source_meta=Path("analysis/v4/environment_source_metadata.json"),work=ROOT/"private/rasters")
    measured,mapping=load_measures(ROOT/"private/original-measured.csv.gz")
    wild=wild.drop(columns=[c for c in measured if c!="obs_id"],errors="ignore").merge(measured,on="obs_id",how="left",validate="one_to_one")
    frozen=json.loads(Path("analysis/v4/results/20261021_unified_abcde/S2_TWO_SPECIES_FULL_WCVP_NICHE_VERIFIED_20261021.json").read_text())
    out=[]
    for sp,rows in wild.groupby("taxon_name"):
        env,_=_background_species_audit(auth,land.dropna(subset=list(VARIABLES)),sp,land)
        complete=rows.dropna(subset=list(VARIABLES))
        sn=complete.loc[complete.range_state=="native"];si=complete.loc[complete.range_state=="introduced"]
        transform,edges,nat,inv=_pca_env_2d(env[env._range=="native"],env[env._range=="introduced"],sn,si,return_grid=True)
        assert abs(1-abs(nat-inv).sum()/2-frozen["species"][sp]["PCA_env_background_corrected_Schoener_D"])<1e-5
        run_species(sp,sn,si,transform,edges,nat,inv,out)
    d=pd.DataFrame(out)
    assert len(d)==2*len(SPACES)*len(REFS)*len(DEFS)==80
    d["holm_all80"]=multipletests(d.p.fillna(1),method="holm")[1]
    d["holm_family20"]=d.groupby(["reference","definition"]).p.transform(lambda p:multipletests(p.fillna(1),method="holm")[1])
    d.to_csv(OUT/"existing_vs_novel.csv",index=False)
    print("Complete; aggregate outputs only.",flush=True)

if __name__=="__main__":
    main()
