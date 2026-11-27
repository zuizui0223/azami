"""Direct empirical climate occupancy of the five native phenotype ranks.

No model-based reconstruction, new fitness proxy, climate-effect regression,
or new significance tests. The species' sampled photographic cells are the
reference, not all available habitat. All locations remain private.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from analysis.v5.phenotype_curve_mechanism import (
    ROOT,SPECIES,DEFINITIONS,VARIABLES,RAIN,SEASON,NAMES,COLS,EXPECTED,
    authority_source,angular_cells,background_mesh,apply_original_five_CHELSA,
    load_three_axis,_background_species_audit,_pca_env_2d,run_species,
    X,K,eligible,make_cells,joint_records,native_cuts,rank_codes,cell_hist,
    climate_cells,marginals,
)
OUT=ROOT/"phenotype_water_occupancy"
WATER=(RAIN,SEASON)
SEED=20261127

def wquant(x,w,probs):
    x=np.asarray(x,float);w=np.asarray(w,float)
    good=w>0
    x,w=x[good],w[good]
    if not len(x) or not np.isfinite(x).all() or not np.isfinite(w).all():
        raise ValueError("No finite positive phenotype mass")
    ix=np.argsort(x,kind="stable")
    x,w=x[ix],w[ix]
    c=np.cumsum(w)/w.sum()
    return x[np.minimum(np.searchsorted(c,probs,side="left"),len(x)-1)]

def full_cohort(sp,sn,si):
    fn,fi=eligible(sn),eligible(si)
    n=make_cells(sn,True);c=make_cells(si,True)
    nv,ni,nw=joint_records(fn,n,COLS,False)
    iv,ii,_=joint_records(fi,c,COLS,False)
    assert (len(n),len(c),len(nv))==(EXPECTED[sp][0],EXPECTED[sp][1],EXPECTED[sp][4])
    cuts,lo,hi=native_cuts(nv,nw)
    rc,valid=rank_codes(nv,cuts,lo,hi)
    assert valid.all()
    support=np.bincount(rc,minlength=K)>0
    assert support.sum()==125
    out={}
    for label,v,ix,cells,frame in (
        ("native",nv,ni,n,fn),("introduced",iv,ii,c,fi)
    ):
        h,tot,kept=cell_hist(v,ix,len(cells),cuts,lo,hi,support)
        good=kept>0
        if label=="native":
            assert good.all() and np.array_equal(tot,kept)
        # Explicitly record any cells without retained interior records.
        assert good.any()
        np.testing.assert_allclose(h[good].sum(1),1.,atol=1e-10)
        x=climate_cells(frame,cells)
        assert np.isfinite(x[good]).all()
        out[label]=dict(h=marginals(h[good]),x=x[good],
                       g=cells.block.to_numpy()[good],
                       source_cells=len(cells),source_records=int(tot.sum()),
                       retained_records=int(kept.sum()))
    return out

def common_cuts(full,axis):
    native=[full[s]["native"]["x"][:,axis] for s in SPECIES]
    # Same physical cutpoints in both species and ranges, from climate
    # alone, with equal total native-cell weight per species.
    x=np.concatenate(native)
    w=np.concatenate([np.full(len(v),.5/len(v)) for v in native])
    q=wquant(x,w,[.2,.4,.6,.8])
    if not (np.diff(q)>0).all():
        raise ValueError("Non-distinct climate bands; do not tune")
    return q

def environment_quantiles(d,axis):
    x,h,g=d["x"][:,axis],d["h"],d["g"]
    n=len(x)
    base_q=wquant(x,np.ones(n),[.1,.5,.9])
    rows=[]
    for j,trait in enumerate(NAMES):
        for k in range(5):
            w=h[:,j,k];pos=w>0
            if not pos.any():
                rows.append(dict(axis=VARIABLES[axis],trait=trait,native_quintile=k+1,
                                 status="absent",cells=n,phenotype_cells=0))
                continue
            q=wquant(x,w,[.1,.5,.9])
            rows.append(dict(axis=VARIABLES[axis],trait=trait,native_quintile=k+1,
                cells=n,phenotype_cells=int(pos.sum()),phenotype_regions=len(np.unique(g[pos])),
                phenotype_frequency=float(w.mean()),effective_cells=float(w.sum()**2/(w@w)),
                climate_q10=float(q[0]),climate_median=float(q[1]),climate_q90=float(q[2]),
                climate_mean=float(w@x/w.sum()),reference_q10=float(base_q[0]),
                reference_median=float(base_q[1]),reference_q90=float(base_q[2]),
                reference_mean=float(x.mean()),status="descriptive"))
    return rows

def climate_band_profile(d,axis,cuts):
    x,h,g=d["x"][:,axis],d["h"],d["g"]
    band=np.searchsorted(cuts,x,side="right")
    rows=[]
    for b in range(5):
        mask=band==b;n=int(mask.sum())
        ids=np.unique(g[mask])
        for j,trait in enumerate(NAMES):
            for k in range(5):
                row=dict(axis=VARIABLES[axis],climate_band=b+1,trait=trait,native_quintile=k+1,
                    cells=n,regions=len(ids),status=("empty" if n==0 else
                    "sparse" if n<20 or len(ids)<3 else "descriptive"))
                if n:
                    row.update(climate_median=float(np.median(x[mask])),
                        frequency=float(h[mask,j,k].mean()),
                        equal_region_frequency=float(np.mean([h[mask&(g==z),j,k].mean() for z in ids])),
                        overall_frequency=float(h[:,j,k].mean()))
                rows.append(row)
    return rows

def signed_profiles(d,axis):
    a,b=marginals(d["hist_left"]),marginals(d["hist_right"])
    left,right=d["left"][:,axis],d["right"][:,axis]
    diff=right-left;g=d["blocks"];n=len(diff)
    rows=[];pieces=[]
    # B's original contrast in each water-direction stratum, then its
    # orientation reversal. Zero differences are retained in accounting.
    for label,mask in (("new_wetter_or_more_seasonal",diff>0),
                       ("new_drier_or_less_seasonal",diff<0),
                       ("equal",diff==0)):
        m=int(mask.sum());ids=np.unique(g[mask])
        if m:
            delta=(b[mask]-a[mask]).mean(0)
            pieces.append(m/n*delta)
        for j,trait in enumerate(NAMES):
            for k in range(5):
                row=dict(axis=VARIABLES[axis],direction=label,trait=trait,native_quintile=k+1,
                    pairs=m,regions=len(ids),all_pairs=n,status=("empty" if not m else
                        "sparse" if m<10 or len(ids)<3 else "descriptive"))
                if m:
                    val=float(delta[j,k])
                    row.update(analogous_frequency=float(a[mask,j,k].mean()),
                      novel_frequency=float(b[mask,j,k].mean()),
                      novel_minus_analogous=val,contribution_to_total=m/n*val,
                      toward_higher_climate=(val if label.startswith("new_wetter") else -val if label.startswith("new_drier") else np.nan),
                      analogous_climate_median=float(np.median(left[mask])),
                      novel_climate_median=float(np.median(right[mask])),
                      climate_difference_median=float(np.median(diff[mask])))
                rows.append(row)
    np.testing.assert_allclose(np.sum(pieces,axis=0),(b-a).mean(0),atol=1e-10)
    nonzero=diff!=0
    if nonzero.any():
        ids=np.unique(g[nonzero]);m=int(nonzero.sum())
        aligned=(b-a)[nonzero]*np.sign(diff[nonzero])[:,None,None]
        left_h=np.where(diff[nonzero,None,None]>0,a[nonzero],b[nonzero])
        right_h=np.where(diff[nonzero,None,None]>0,b[nonzero],a[nonzero])
        np.testing.assert_allclose(aligned,right_h-left_h,atol=1e-10)
        for j,trait in enumerate(NAMES):
            for k in range(5):
                y=aligned[:,j,k];gg=g[nonzero]
                rows.append(dict(axis=VARIABLES[axis],direction="aligned_lower_to_higher",
                    trait=trait,native_quintile=k+1,pairs=m,regions=len(ids),all_pairs=n,
                    status="descriptive" if m>=10 and len(ids)>=3 else "sparse",
                    lower_climate_frequency=float(left_h[:,j,k].mean()),
                    higher_climate_frequency=float(right_h[:,j,k].mean()),
                    toward_higher_climate=float(y.mean()),
                    equal_region_toward_higher=float(np.mean([y[gg==z].mean() for z in ids])),
                    lower_climate_median=float(np.median(np.minimum(left,right)[nonzero])),
                    higher_climate_median=float(np.median(np.maximum(left,right)[nonzero])),
                    climate_difference_median=float(np.median(abs(diff[nonzero])))))
    return rows

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
    out={};full={}
    for sp,rows in wild.groupby("taxon_name"):
        env,_=_background_species_audit(auth,land.dropna(subset=list(VARIABLES)),sp,land)
        complete=rows.dropna(subset=list(VARIABLES))
        sn=complete.loc[complete.range_state=="native"];si=complete.loc[complete.range_state=="introduced"]
        transform,edges,nat,inv=_pca_env_2d(env[env._range=="native"],env[env._range=="introduced"],sn,si,return_grid=True)
        assert abs(1-abs(nat-inv).sum()/2-frozen["species"][sp]["PCA_env_background_corrected_Schoener_D"])<1e-5
        full[sp]=full_cohort(sp,sn,si)
        run_species(sp,sn,si,transform,edges,nat,inv,out)
    quantiles=[];bands=[];directions=[];coverage=[];cuts_rows=[]
    for axis in WATER:
        cuts=common_cuts(full,axis)
        cuts_rows.append(dict(axis=VARIABLES[axis],**{f"q{k*20}":float(c) for k,c in enumerate(cuts,1)}))
        for sp in SPECIES:
            for range_name,d in full[sp].items():
                quantiles.extend(dict(species=sp,range=range_name,**r)
                    for r in environment_quantiles(d,axis))
                bands.extend(dict(species=sp,range=range_name,**r)
                    for r in climate_band_profile(d,axis,cuts))
            for definition in DEFINITIONS:
                d=out[(sp,definition)]
                directions.extend(dict(species=sp,definition=definition,**r)
                    for r in signed_profiles(d,axis))
    for sp in SPECIES:
        for range_name,d in full[sp].items():
            coverage.append(dict(species=sp,range=range_name,
                source_cells=d["source_cells"],retained_cells=len(d["h"]),
                regions=len(np.unique(d["g"])),source_records=d["source_records"],
                retained_records=d["retained_records"]))
    pd.DataFrame(quantiles).to_csv(OUT/"phenotype_water_distributions.csv",index=False)
    pd.DataFrame(bands).to_csv(OUT/"observed_water_band_composition.csv",index=False)
    pd.DataFrame(directions).to_csv(OUT/"matched_signed_water_composition.csv",index=False)
    pd.DataFrame(coverage).to_csv(OUT/"cohort.csv",index=False)
    pd.DataFrame(cuts_rows).to_csv(OUT/"shared_water_bands.csv",index=False)
    print("Complete; aggregate outputs only.",flush=True)

if __name__=="__main__":
    main()
