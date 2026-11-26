"""Rarefied native geographic holdout calibration of phenotypic support.

Contemporary native observations only. Fixed climatic pairs; no fitted selection
coefficients, no inference about ancestral allele frequencies, no private rows.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from statsmodels.stats.multitest import multipletests
from analysis.v5.phenotype_redistribution import ROOT, SEED, eligible, cellkeys
from analysis.v5.phenotype_endpoint_scan import permute
from analysis.v5.novel_head_overlap import T, make_cells, match_cells
from analysis.v5.invasion_partition import head_cells
from analysis.v5.native_pca_env import background_mesh, _background_species_audit, _pca_env_2d, _probability_content
from analysis.v5.native_pca_env_fullcoverage import apply_original_five_CHELSA
from analysis.v4.native_introduced_opportunity import authority_source, angular_cells
from analysis.v4.native_cached_climate import VARIABLES

OUT=ROOT/"phenotype_native_calibration"
N_REF,N_CAL,N_DRAWS=256,128,32
LEVELS=np.array([.90,.95,.99])

def native_splits(blocks,nref=N_REF,ncal=N_CAL,draws=N_DRAWS,seed=SEED):
    blocks=np.asarray(blocks)
    unique=np.unique(blocks)
    if len(unique)<4:raise ValueError("Too few native blocks")
    rng=np.random.default_rng(seed)
    result=[]
    for attempt in range(10000):
        ordered=rng.permutation(unique)
        held=ordered[:max(2,len(unique)//4)]
        cal=np.flatnonzero(np.isin(blocks,held))
        train=np.flatnonzero(~np.isin(blocks,held))
        if len(cal)<ncal or len(train)<nref:continue
        a=rng.choice(train,nref,replace=False);b=rng.choice(cal,ncal,replace=False)
        assert not set(blocks[a])&set(blocks[b])
        result.append((a,b))
        if len(result)==draws:return result
    raise ValueError("Insufficient geographically distinct native support")

def distances(reference,values,k=5):
    mu=reference.mean(0);sd=reference.std(0)
    if (sd<=0).any():raise ValueError("Constant native feature")
    ref=(reference-mu)/sd
    z=(values-mu)/sd
    if not 1<=k<=len(ref):raise ValueError("Invalid nearest-neighbor count")
    # Sequence k preserves the neighbor axis even for k=1 (scipy squeezes scalar k=1).
    return cKDTree(ref).query(z,k=list(range(1,k+1)))[0].mean(axis=1)

def calibration(n,c,cols,splits):
    a=n[cols].to_numpy(float);b=c[cols].to_numpy(float)
    continuous=[];tails=[];heldout=[];thresholds=[]
    for tr,ca in splits:
        h=distances(a[tr],a[ca])
        z=distances(a[tr],b)
        q=np.quantile(h,LEVELS)
        if q[1]<=0:raise ValueError("Degenerate support calibration")
        thresholds.append(q)
        heldout.append((h[:,None]>q).mean(0))
        tails.append(z[:,None]>q)
        continuous.append(z/q[1])
    return (np.mean(continuous,axis=0),np.mean(tails,axis=0),
            np.mean(heldout,axis=0),np.mean(thresholds,axis=0))

def compare(sp,definition,n,c,pairs,splits,out):
    total=sum(len(a) for a,b in pairs)
    for name,cols in [("angle_chroma",T[:2]),("three_traits",T)]:
        distance,tail,baseline,threshold=calibration(n,c,cols,splits)
        values=np.column_stack([distance,tail])
        counts=np.array([len(a) for a,b in pairs])
        left=np.concatenate([a for a,b in pairs]);right=np.concatenate([b for a,b in pairs])
        for group,ix in [("analogue_all",np.flatnonzero(c.novel==0)),("novel_all",np.flatnonzero(c.novel==1)),
                         ("analogue_matched",left),("novel_matched",right)]:
            for j,label in enumerate(["normalized_knn_distance","above_native_q90","above_native_q95","above_native_q99"]):
                out["profiles"].append(dict(species=sp,definition=definition,trait_space=name,group=group,endpoint=label,
                    cells=len(ix),pairs=total,blocks=len(pairs),value=values[ix,j].mean(),
                    native_calibration=baseline[j-1] if j else np.nan,
                    mean_distance_threshold=threshold[j-1] if j else np.nan))
        ds=np.stack([(values[b]-values[a]).sum(0) for a,b in pairs])
        rng=np.random.default_rng(SEED)
        bs=rng.multinomial(len(counts),np.full(len(counts),1/len(counts)),size=3999)
        bd=(bs@ds)/(bs@counts)[:,None]
        for j,label in enumerate(["normalized_knn_distance","above_native_q90","above_native_q95","above_native_q99"]):
            x=ds[:,j:j+1];_,p,exact,nperm=permute(x,counts)
            delta=float(x.sum()/total)
            lo,hi=np.quantile(bd[:,j],[.025,.975])
            loo=[float((x.sum()-x[k,0])/(total-counts[k])) for k in range(len(counts)) if total>counts[k]]
            out["tests"].append(dict(species=sp,definition=definition,trait_space=name,endpoint=label,pairs=total,
              blocks=len(pairs),difference=delta,lo=lo,hi=hi,p=p,exact=exact,permutations=nperm,
              loo_min=min(loo),loo_max=max(loo)))

def run_species(sp,sn,si,transform,edges,nat,inv,out):
    n=make_cells(sn,True);full=make_cells(si,True)
    expected={"Cirsium arvense":(1366,860,30,183),"Cirsium vulgare":(1113,1139,146,316)}[sp]
    assert (len(n),len(full))==expected[:2]
    splits=native_splits(n.block.to_numpy())
    out["native_design"].append(dict(species=sp,n_native=len(n),native_blocks=n.block.nunique(),reference_cells=N_REF,
       calibration_cells=N_CAL,draws=len(splits),minimum_train_blocks=min(n.iloc[a].block.nunique() for a,b in splits),
       minimum_heldout_blocks=min(n.iloc[b].block.nunique() for a,b in splits)))
    coords=transform(si)
    uv=np.column_stack([np.searchsorted(e,coords[:,k],side="right")-1 for k,e in enumerate(edges)])
    for k,e in enumerate(edges):uv[coords[:,k]==e[-1],k]=len(e)-2
    inside=np.all((uv>=0)&(uv<np.array([len(e)-1 for e in edges])),axis=1)
    labeled=si.loc[inside].copy();on=_probability_content(nat,.95)
    labeled["novel"]=(~on[uv[inside,0],uv[inside,1]]).astype(int)
    c=make_cells(labeled);c=c.loc[c.novel!=.5].reset_index(drop=True)
    pairs=match_cells(c);assert sum(len(a) for a,b in pairs)==expected[2]
    compare(sp,"native_core95",n,c,pairs,splits,out)
    hc=head_cells(si,transform,edges)[["i","j","u","v"]]
    c=full.merge(hc,on=["i","j"],validate="one_to_one")
    u,v=c.u.to_numpy(),c.v.to_numpy()
    c=c.loc[inv[u,v]>0].reset_index(drop=True);u,v=c.u.to_numpy(),c.v.to_numpy()
    c["novel"]=(inv[u,v]>nat[u,v]).astype(int)
    pairs=match_cells(c);assert sum(len(a) for a,b in pairs)==expected[3]
    compare(sp,"density_surplus",n,c,pairs,splits,out)

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    auth=authority_source((ROOT/"private/authority350.zip").read_bytes())
    wild=pd.read_csv(ROOT/"private/focal_wild.csv.gz",dtype={"obs_id":str},low_memory=False)
    assert len(wild)==297622
    wild["_cell"]=angular_cells(wild.latitude,wild.longitude)
    wild=wild.sort_values("obs_id").drop_duplicates(["taxon_name","range_state","_cell","user_id"]).copy()
    land=background_mesh(auth)
    wild,land,_=apply_original_five_CHELSA(wild,land,
        source_meta=Path("analysis/v4/environment_source_metadata.json"),work=ROOT/"private/rasters")
    measured=pd.read_csv(ROOT/"private/original-measured.csv.gz",usecols=["obs_id",*T],dtype={"obs_id":str})
    wild=wild.drop(columns=T,errors="ignore").merge(measured,on="obs_id",how="left",validate="one_to_one")
    frozen=json.loads(Path("analysis/v4/results/20261021_unified_abcde/S2_TWO_SPECIES_FULL_WCVP_NICHE_VERIFIED_20261021.json").read_text())
    out={k:[] for k in ["native_design","profiles","tests"]}
    for sp,rows in wild.groupby("taxon_name"):
        env,_=_background_species_audit(auth,land.dropna(subset=list(VARIABLES)),sp,land)
        complete=rows.dropna(subset=list(VARIABLES))
        sn=complete.loc[complete.range_state=="native"];si=complete.loc[complete.range_state=="introduced"]
        transform,edges,nat,inv=_pca_env_2d(env[env._range=="native"],env[env._range=="introduced"],sn,si,return_grid=True)
        assert abs(1-abs(nat-inv).sum()/2-frozen["species"][sp]["PCA_env_background_corrected_Schoener_D"])<1e-5
        run_species(sp,sn,si,transform,edges,nat,inv,out)
    for name,rows in out.items():
        d=pd.DataFrame(rows)
        if name=="tests":
            assert len(d)==32
            d["holm_all32"]=multipletests(d.p,method="holm")[1]
        d.to_csv(OUT/(name+".csv"),index=False)
    print("Complete; aggregate outputs only.",flush=True)

if __name__=="__main__":
    main()
