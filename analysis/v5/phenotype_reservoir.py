"""Cross-sectional native phenotype reservoir and source-mixture diagnostic.

No observation dates used. Fixed climate-only assignments and original pairs.
The mixture is a descriptive feasibility model, not inferred introduction
ancestry, natural selection, fitness or genetic variant frequencies.
"""
from pathlib import Path
from itertools import product
import json
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from analysis.v5.phenotype_redistribution import ROOT, eligible, cellkeys, weighted_quantile
from analysis.v5.phenotype_endpoint_scan import load_measures, cell_average_matrix, permute
from analysis.v5.novel_head_overlap import T, make_cells, match_cells
from analysis.v5.invasion_partition import head_cells
from analysis.v5.native_pca_env import background_mesh, _background_species_audit, _pca_env_2d, _probability_content
from analysis.v5.native_pca_env_fullcoverage import apply_original_five_CHELSA
from analysis.v4.native_introduced_opportunity import authority_source, angular_cells
from analysis.v4.native_cached_climate import VARIABLES
from statsmodels.stats.multitest import multipletests

OUT=ROOT/"phenotype_reservoir"
COLS={"chroma":"scan_corolla_lab_chroma","lightness":"scan_corolla_lab_lightness",
      "angle":"scan_orientation_image_vertical_angle"}
SPACES={"colour":("chroma","lightness"),"colour_angle":("chroma","lightness","angle")}
LABELS=("low","mid","high")
MIN_REGION=20

def records(f,c,cols):
    ix=pd.Index(cellkeys(c)).get_indexer(f.cellkey)
    v=f[cols].to_numpy(float)
    ok=(ix>=0)&np.isfinite(v).all(1)
    ix,v=ix[ok],v[ok]
    counts=np.bincount(ix,minlength=len(c))
    return v,ix,1/counts[ix]

def hist(v,ix,n,cuts):
    code=np.zeros(len(ix),int)
    for j,q in enumerate(cuts):
        code=3*code+np.searchsorted(q,v[:,j],side="right")
    k=3**len(cuts)
    h=cell_average_matrix(np.eye(k)[code],ix,n)
    good=np.isfinite(h).all(1)
    np.testing.assert_allclose(h[good].sum(1),1)
    return h

def tv(a,b):
    return .5*np.abs(np.asarray(a)-np.asarray(b)).sum()

def convex_fit(regions,target,initial):
    regions=np.asarray(regions,float)
    target=np.asarray(target,float)
    initial=np.asarray(initial,float)
    np.testing.assert_allclose(regions.sum(1),1)
    if len(regions)==1:return regions[0].copy()
    def obj(w):
        e=w@regions-target
        return float(e@e)
    def jac(w):
        return 2*regions@(w@regions-target)
    r=minimize(obj,initial,jac=jac,method="SLSQP",bounds=[(0,1)]*len(regions),
        constraints=[{"type":"eq","fun":lambda w:w.sum()-1,"jac":lambda w:np.ones_like(w)}],
        options={"ftol":1e-12,"maxiter":1000})
    if not r.success or abs(r.x.sum()-1)>1e-6:
        raise ValueError("Native source mixture optimization did not converge")
    return r.x@regions

def supported_reweight(native,target):
    # The target defines empirical frequency weights on already observed bins.
    # Unsupported target mass is explicitly counted, not invented in native.
    supported=native>0
    q=target*supported
    return (q/q.sum() if q.sum()>0 else native.copy()),float(target[~supported].sum())

def analyze(sp,definition,space,n,c,pairs,fn,fi,out):
    names=SPACES[space];cols=[COLS[k] for k in names]
    nv,ni,nw=records(fn,n,cols);iv,ii,iw=records(fi,c,cols)
    if not len(nv) or not len(iv):
        out["cohort"].append(dict(species=sp,definition=definition,space=space,status="missing_joint_measurements"))
        return
    cuts=[weighted_quantile(nv[:,j],nw,[1/3,2/3]) for j in range(len(cols))]
    a=hist(nv,ni,len(n),cuts);b=hist(iv,ii,len(c),cuts)
    finite=np.isfinite(b).all(1)
    pairs2=[(l[finite[l]&finite[r]],r[finite[l]&finite[r]]) for l,r in pairs]
    pairs2=[(l,r) for l,r in pairs2 if len(l)]
    total=sum(len(l) for l,r in pairs)
    count=sum(len(l) for l,r in pairs2)
    good=np.isfinite(a).all(1)
    region_ids=n.block.to_numpy()
    region_counts=pd.Series(region_ids[good]).value_counts()
    keep=region_counts[region_counts>=MIN_REGION].sort_index()
    if not count or good.sum()<50:
        out["cohort"].append(dict(species=sp,definition=definition,space=space,status="insufficient_joint_measurements",fixed_pairs=total,eligible_pairs=count))
        return
    ref=a[good].mean(0)
    source=np.stack([a[good&(region_ids==g)].mean(0) for g in keep.index]) if len(keep) else np.empty((0,len(ref)))
    initial=keep.to_numpy(float)/keep.sum() if len(keep) else np.array([])
    left=np.concatenate([l for l,r in pairs2]);right=np.concatenate([r for l,r in pairs2])
    u,v=b[left].mean(0),b[right].mean(0)
    ready=count>=20 and count>=.8*total and len(pairs2)>=8 and len(keep)>=3
    out["cohort"].append(dict(species=sp,definition=definition,space=space,
       status="estimated" if ready else "insufficient_support",
       native_cells=int(good.sum()),native_records=len(nv),native_regions=len(keep),
       native_cells_in_regions=int(keep.sum()),introduced_records=len(iv),
       fixed_pairs=total,eligible_pairs=count,blocks=len(pairs2),
       native_observed_bins=int((ref>0).sum()),bins=len(ref),
       analogue_unobserved_native_mass=float(u[ref==0].sum()),
       novel_unobserved_native_mass=float(v[ref==0].sum()),
       tv_native_analogue=tv(ref,u),tv_native_novel=tv(ref,v)))
    labels=["/".join(x) for x in product(LABELS,repeat=len(names))]
    # One fixed biological prediction from preceding observations: joint high
    # chroma and low lightness, in the 2D colour space. All bins remain visible.
    candidate=np.array([x.split("/")[:2]==["high","low"] for x in labels])
    for j,label in enumerate(labels):
        out["bins"].append(dict(species=sp,definition=definition,space=space,bin=label,
            native=ref[j],analogue=u[j],novel=v[j],difference=v[j]-u[j],
            candidate_colour=bool(candidate[j]),
            analogue_over_native=(u[j]/ref[j] if ref[j]>0 else np.nan),
            novel_over_native=(v[j]/ref[j] if ref[j]>0 else np.nan)))
    if not ready:return
    # Distribution-level assessment on held-out INTRODUCED regions. Train on
    # other paired regions only. Same pairs/folds for both sides; no date use.
    source_pool=initial@source
    accum={}
    for side in [0,1]:
        for k,(l,r) in enumerate(pairs2):
            te=(l,r)[side]
            tr=np.concatenate([x[side] for j,x in enumerate(pairs2) if j!=k])
            target=b[tr].mean(0);truth=b[te]
            mixture=convex_fit(source,target,initial)
            rw,unsupported=supported_reweight(ref,target)
            predictions={"native_pool":source_pool,"native_region_mixture":mixture,
                         "native_bin_reweight":rw,"introduced_empirical":target}
            for model,p in predictions.items():
                e=np.square(truth-p).sum(1)
                key=(side,model)
                accum.setdefault(key,[]).append((len(te),float(e.sum()),tv(truth.mean(0),p),unsupported if model=="native_bin_reweight" else np.nan))
    for (side,model),rows in accum.items():
        z=np.array(rows)
        baseline=np.array(accum[(side,"native_pool")])
        # Equal-cell squared error and equal-region mean TV; improvement is
        # descriptive cross-validation, not genetic/source selection evidence.
        out["heldout"].append(dict(species=sp,definition=definition,space=space,
            group=["analogue","novel"][side],model=model,blocks=len(rows),cells=int(z[:,0].sum()),
            mean_cell_squared_error=z[:,1].sum()/z[:,0].sum(),mean_block_tv=z[:,2].mean(),
            skill_vs_native_pool=1-z[:,1].sum()/baseline[:,1].sum() if baseline[:,1].sum()>0 else np.nan,
            max_unsupported_target_mass=float(np.nanmax(z[:,3])) if model=="native_bin_reweight" else np.nan))
    # Candidate contrast is explicitly post-hoc, not a confirmatory test.
    # A separate omnibus uses all 2D/3D bins with the same block flips.
    counts=np.array([len(l) for l,r in pairs2])
    ds=np.stack([(b[r]-b[l]).sum(0) for l,r in pairs2])
    effect,p,exact,nperm=permute(ds,counts)
    out["tests"].append(dict(species=sp,definition=definition,space=space,test="joint_distribution",
        pairs=count,blocks=len(pairs2),effect=effect,p=p,exact=exact,permutations=nperm))
    if space=="colour":
        d=ds[:,candidate].sum(1)
        _,p,exact,nperm=permute(d[:,None],counts)
        out["tests"].append(dict(species=sp,definition=definition,space=space,test="high_chroma_low_lightness",
            pairs=count,blocks=len(pairs2),effect=float(d.sum()/count),p=p,exact=exact,permutations=nperm))

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
    assert all(mapping[key][0] and mapping[key][1] for key in ["orientation_image_vertical_angle","corolla_lab_chroma","corolla_lab_lightness"])
    wild=wild.drop(columns=[c for c in measured if c!="obs_id"],errors="ignore").merge(measured,on="obs_id",how="left",validate="one_to_one")
    frozen=json.loads(Path("analysis/v4/results/20261021_unified_abcde/S2_TWO_SPECIES_FULL_WCVP_NICHE_VERIFIED_20261021.json").read_text())
    out={k:[] for k in ["cohort","bins","heldout","tests"]}
    for sp,rows in wild.groupby("taxon_name"):
        env,_=_background_species_audit(auth,land.dropna(subset=list(VARIABLES)),sp,land)
        complete=rows.dropna(subset=list(VARIABLES))
        sn=complete.loc[complete.range_state=="native"];si=complete.loc[complete.range_state=="introduced"]
        transform,edges,nat,inv=_pca_env_2d(env[env._range=="native"],env[env._range=="introduced"],sn,si,return_grid=True)
        target=frozen["species"][sp]["PCA_env_background_corrected_Schoener_D"]
        assert abs(1-abs(nat-inv).sum()/2-target)<1e-5
        run_species(sp,sn,si,transform,edges,nat,inv,out)
    # Preserve all 12 planned tests, including any unsupported tests as p=1.
    planned=pd.DataFrame([(sp,d,s,t) for sp in ["Cirsium arvense","Cirsium vulgare"]
        for d in ["native_core95","density_surplus"] for s in SPACES
        for t in (["joint_distribution","high_chroma_low_lightness"] if s=="colour" else ["joint_distribution"])],
        columns=["species","definition","space","test"])
    tested=pd.DataFrame(out["tests"])
    d=planned.merge(tested,on=["species","definition","space","test"],how="left",validate="one_to_one") if len(tested) else planned.assign(p=np.nan)
    d["status"]=np.where(d.p.notna(),"estimated","not_eligible")
    d["p_planned"]=d.p.fillna(1.)
    d["holm_all12"]=multipletests(d.p_planned,method="holm")[1]
    d.to_csv(OUT/"tests.csv",index=False)
    for name in ["cohort","bins","heldout"]:
        pd.DataFrame(out[name]).to_csv(OUT/(name+".csv"),index=False)
    print("Complete; aggregate outputs only.",flush=True)

if __name__=="__main__":
    main()
