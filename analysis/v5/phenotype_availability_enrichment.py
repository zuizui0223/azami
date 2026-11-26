"""Fixed opportunity-aware within-species enrichment comparison.
Use the original native-99%-background and native-95%-occupancy grids.
Newly occupied common background vs common core; and common-background
novel vs outside-native-background AND outside-native-core. No selection
coefficients, fitness, ancestry, time axis or interspecific predictions.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from analysis.v5.phenotype_within_species_enrichment import (
    ROOT,NAMES,COLS,EXPECTED,joint_records,T,eligible,load_three_axis,
    make_cells,match_cells,head_cells,background_mesh,_background_species_audit,
    _pca_env_2d,_probability_content,apply_original_five_CHELSA,authority_source,
    angular_cells,VARIABLES,X,K,RANKS,native_cuts,rank_codes,cell_hist,
    fit_coefficients,score_omnibus,
)
from analysis.v5.invasion_partition import background
OUT=ROOT/"phenotype_availability_enrichment"
GROUPS=("common_core","common_expansion","outside_background_novel",
        "outside_background_core","outside_own_background")
CONTRASTS={"common_expansion_vs_core":("common_core","common_expansion"),
           "common_vs_outside_background_novel":("outside_background_novel","common_expansion")}

def partition(an,ai,on):
    groups=[ai&an&on,ai&an&~on,ai&~an&~on,ai&~an&on,~ai]
    assert np.all(np.sum(groups,axis=0)==1)
    return dict(zip(GROUPS,groups))

def group_profile(v,ix,cells,indices,cuts,lo,hi):
    """Unconditional observed-record frequency; do not hide outside-native mass."""
    selected=np.isin(ix,indices)
    v,ix=v[selected],ix[selected]
    if not len(indices):
        return []
    count=np.bincount(ix,minlength=cells)
    assert np.all(count[indices]>0)
    w=1/count[ix]/len(indices)
    np.testing.assert_allclose(w.sum(),1)
    rows=[]
    for j,trait in enumerate(NAMES):
        low=(v[:,j]>=lo[j])&(v[:,j]<cuts[j,0])
        high=(v[:,j]>=cuts[j,-1])&(v[:,j]<=hi[j])
        outside=(v[:,j]<lo[j])|(v[:,j]>hi[j])
        rows.append(dict(trait=trait,cells=len(indices),records=len(ix),
            mean=float(w@v[:,j]),existing_low=float(w@low),
            existing_high=float(w@high),outside=float(w@outside)))
    return rows

def emit_unavailable(out,meta,status):
    out["omnibus"].append(dict(**meta,status=status,p=1.))
    out["coefficients"].extend(dict(**meta,trait=trait,term=term,status=status,p_score=1.)
        for term in ("directional","quadratic") for trait in NAMES)

def run_species(sp,sn,si,env,transform,edges,nat,inv,out):
    fn,fi=eligible(sn),eligible(si)
    n=make_cells(sn,True)
    full=make_cells(si,True)
    assert (len(n),len(full))==EXPECTED[sp][:2]
    nv,ni,nw=joint_records(fn,n,COLS,False)
    assert len(nv)==EXPECTED[sp][4]
    cuts,lo,hi=native_cuts(nv,nw)
    rc,valid=rank_codes(nv,cuts,lo,hi)
    assert valid.all()
    support=np.bincount(rc,minlength=K)>0
    assert support.sum()==125
    bn,bi=[background(env[env._range==r],transform,edges) for r in ("native","introduced")]
    an,ai=[_probability_content(v,.99) for v in (bn,bi)]
    on=_probability_content(nat,.95)
    groups=partition(an,ai,on)
    c=full.merge(head_cells(si,transform,edges)[["i","j","u","v"]],on=["i","j"],validate="one_to_one")
    u,v=c.u.to_numpy(),c.v.to_numpy()
    c["availability"]="unclassified"
    for name,mask in groups.items():
        c.loc[mask[u,v],"availability"]=name
    assert c.availability.ne("unclassified").all()
    observed=c.availability.value_counts()
    expected={"Cirsium arvense":(822,28,7),"Cirsium vulgare":(877,146,110)}[sp]
    assert (int(observed.get("common_core",0)),int(observed.get("common_expansion",0)),
            int(observed.get("outside_background_novel",0)+observed.get("outside_background_core",0)))==expected
    iv,ii,_=joint_records(fi,c,COLS,False)
    h,tot,kept=cell_hist(iv,ii,len(c),cuts,lo,hi,support)
    for g,idx in [("native",np.arange(len(n)))]+[(g,np.flatnonzero(c.availability==g)) for g in GROUPS]:
        if g=="native":
            profiles=group_profile(nv,ni,len(n),idx,cuts,lo,hi)
        else:
            out["cohort"].append(dict(species=sp,group=g,cells=len(idx),
                eligible_records=int(tot[idx].sum()),outside_joint_support=int((tot-kept)[idx].sum()),
                blocks=int(c.iloc[idx].block.nunique()),unclassified_cells=len(full)-len(c)))
            if len(idx):
                assert np.all(tot[idx]>0)
            profiles=group_profile(iv,ii,len(c),idx,cuts,lo,hi)
        out["frequencies"].extend(dict(species=sp,contrast="all_groups",group=g,**r) for r in profiles)
    for contrast,(left_group,right_group) in CONTRASTS.items():
        sub=c[c.availability.isin([left_group,right_group])].copy()
        sub["original_index"]=sub.index
        sub=sub.reset_index(drop=True)
        sub["novel"]=(sub.availability==right_group).astype(int)
        pairs=match_cells(sub)
        left=np.concatenate([l for l,r in pairs]) if pairs else np.array([],int)
        right=np.concatenate([r for l,r in pairs]) if pairs else np.array([],int)
        li=sub.original_index.to_numpy()[left];ri=sub.original_index.to_numpy()[right]
        blocks=sub.iloc[left].block.to_numpy()
        np.testing.assert_array_equal(blocks,sub.iloc[right].block.to_numpy())
        meta=dict(species=sp,contrast=contrast,reference_group=left_group,target_group=right_group,
            pairs=len(li),blocks=len(pairs),native_cells=len(n),native_records=len(nv))
        for g,idx in [(left_group,li),(right_group,ri)]:
            out["frequencies"].extend(dict(species=sp,contrast=contrast,group=g,**r)
                for r in group_profile(iv,ii,len(c),idx,cuts,lo,hi))
        if len(li)<20 or not 8<=len(pairs)<=15:
            emit_unavailable(out,meta,"insufficient_support")
            continue
        if np.any(kept[np.r_[li,ri]]==0):
            emit_unavailable(out,meta,"missing_supported_observations")
            continue
        meta.update(reference_records=int(tot[li].sum()),target_records=int(tot[ri].sum()),
            reference_excluded_fraction=float(np.mean(1-kept[li]/tot[li])),
            target_excluded_fraction=float(np.mean(1-kept[ri]/tot[ri])))
        a,b=h[li],h[ri]
        out["omnibus"].append(dict(**meta,status="estimated",**score_omnibus(a,b,X,blocks)))
        try:
            _,rows=fit_coefficients(a,b,X,blocks)
        except (ValueError,np.linalg.LinAlgError) as e:
            rows=[dict(trait=t,term=k,status="fit_failed",failure_reason=str(e),p_score=1.)
                    for k in ("directional","quadratic") for t in NAMES]
        out["coefficients"].extend(dict(**meta,**r) for r in rows)

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
    out={"cohort":[],"frequencies":[],"coefficients":[],"omnibus":[]}
    for sp,rows in wild.groupby("taxon_name"):
        env,_=_background_species_audit(auth,land.dropna(subset=list(VARIABLES)),sp,land)
        complete=rows.dropna(subset=list(VARIABLES))
        sn=complete.loc[complete.range_state=="native"];si=complete.loc[complete.range_state=="introduced"]
        transform,edges,nat,inv=_pca_env_2d(env[env._range=="native"],env[env._range=="introduced"],sn,si,return_grid=True)
        assert abs(1-abs(nat-inv).sum()/2-frozen["species"][sp]["PCA_env_background_corrected_Schoener_D"])<1e-5
        run_species(sp,sn,si,env,transform,edges,nat,inv,out)
    d=pd.DataFrame(out["coefficients"])
    assert len(d)==24
    d["holm_all24"]=multipletests(d.p_score.fillna(1),method="holm")[1]
    d.to_csv(OUT/"availability_coefficients.csv",index=False)
    g=pd.DataFrame(out["omnibus"])
    assert len(g)==4
    g["holm_all4"]=multipletests(g.p.fillna(1),method="holm")[1]
    g.to_csv(OUT/"availability_omnibus.csv",index=False)
    for name in ["cohort","frequencies"]:
        pd.DataFrame(out[name]).to_csv(OUT/(name+".csv"),index=False)
    print("Complete; aggregate outputs only.",flush=True)

if __name__=="__main__":
    main()
