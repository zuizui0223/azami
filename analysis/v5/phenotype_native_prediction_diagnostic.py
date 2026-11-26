"""Final native-only model diagnostic after failed rarefied kNN transport.

All available native cells in leave-one-geographic-block-out validation.
Four fixed climate models, no choice using introduced outcomes. Aggregate only.
"""
from pathlib import Path
from itertools import product
import json
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from statsmodels.stats.multitest import multipletests
from sklearn.linear_model import Ridge
from analysis.v5.phenotype_redistribution import ROOT, SEED, eligible, cellkeys
from analysis.v5.phenotype_endpoint_scan import load_measures, values, design_matrices
from analysis.v5.novel_head_overlap import T, make_cells
from analysis.v5.native_pca_env import background_mesh, _background_species_audit, _pca_env_2d
from analysis.v5.native_pca_env_fullcoverage import apply_original_five_CHELSA
from analysis.v4.native_introduced_opportunity import authority_source, angular_cells
from analysis.v4.native_cached_climate import VARIABLES
from analysis.v4.native_integration import FEATURES

OUT=ROOT/"phenotype_native_prediction_diagnostic"
MODELS=["knn20","knn100","ridge_0.1n","ridge_1n"]
MIN_TEST=8
MIN_TRAIN=100

def predictions(train_x,train_y,test_x):
    mu=train_x.mean(0);sd=train_x.std(0)
    if not np.isfinite(sd).all() or (sd<=0).any() or len(train_x)<100:
        raise ValueError("Insufficient native climate reference")
    x=(train_x-mu)/sd;z=(test_x-mu)/sd
    _,ix=cKDTree(x).query(z,k=list(range(1,101)))
    # Ridge is fixed a priori; outputs clipped only for squared-error scoring
    # in valid CDF range in the caller, identical treatment in each fold.
    return {"knn20":train_y[ix[:,:20]].mean(1),
            "knn100":train_y[ix].mean(1),
            "ridge_0.1n":Ridge(alpha=.1*len(x)).fit(x,train_y).predict(z),
            "ridge_1n":Ridge(alpha=len(x)).fit(x,train_y).predict(z)}

def fold_scan(x,blocks,v,ids,w,circ):
    n=len(x);rows=[]
    for b in np.unique(blocks):
        train=np.flatnonzero(blocks!=b);test=np.flatnonzero(blocks==b)
        if len(train)<MIN_TRAIN or len(test)<MIN_TEST:continue
        mask=np.isin(ids,train)
        ref=v[mask];rw=w[mask]
        if not len(ref):continue
        f,_,_=design_matrices(v,ids,n,ref,rw,circ)
        train=train[np.isfinite(f[train]).all(1)];test=test[np.isfinite(f[test]).all(1)]
        if len(train)<MIN_TRAIN or len(test)<MIN_TEST:continue
        y=f[train];truth=f[test];base=y.mean(0)
        baseline=np.square(truth-base).sum()
        if baseline<=1e-12:continue
        pred=predictions(x[train],y,x[test])
        for name,p in pred.items():
            if not circ:p=np.clip(p,0,1)
            err=np.square(truth-p).sum()
            rows.append(dict(model=name,test_cells=len(test),train_cells=len(train),
                             baseline_sse=baseline,sse=err))
    return pd.DataFrame(rows)

def aggregate(g):
    e=g.sse.to_numpy();b=g.baseline_sse.to_numpy()
    delta=(b-e)/b
    k=len(g)
    if not k:return {}
    # Block-equal skill is not dominated by large geographic blocks.
    # Sign flips are an approximate descriptive uncertainty diagnostic, not
    # independent experimental replicates (LOBO training sets overlap).
    exact=k<=15
    rng=np.random.default_rng(SEED)
    flips=np.array(list(product([-1.,1.],repeat=k))) if exact else rng.choice([-1.,1.],(9999,k))
    val=float(delta.mean());null=(flips@delta)/k
    p=(np.count_nonzero(abs(null)>=abs(val)-1e-12)+int(not exact))/(len(null)+int(not exact))
    return dict(blocks=k,test_cells=int(g.test_cells.sum()),min_train_cells=int(g.train_cells.min()),
                max_train_cells=int(g.train_cells.max()),q2_cell_weighted=1-e.sum()/b.sum(),
                q2_block_equal=val,positive_blocks=int((delta>0).sum()),
                min_block_q2=float(delta.min()),max_block_q2=float(delta.max()),
                p_block_sign=p,sign_patterns=len(null))

def run_species(sp,sn,mapping,out):
    n=make_cells(sn,True)
    expected={"Cirsium arvense":1366,"Cirsium vulgare":1113}[sp]
    assert len(n)==expected
    fn=eligible(sn)
    x=fn.groupby("cellkey")[list(VARIABLES)].median().reindex(cellkeys(n)).to_numpy(float)
    assert np.isfinite(x).all()
    blocks=n.block.to_numpy()
    for name,components in FEATURES.items():
        circ=name=="hue"
        cols=["scan_"+a for a in components]
        v,ids,w,_=values(fn,n,cols,"record_equalcell",circ)
        meta=dict(species=sp,trait=name,total_native_cells=len(n),total_native_blocks=len(np.unique(blocks)),
                  source_columns=";".join(str(mapping[a][0] or "unavailable") for a in components))
        folds=fold_scan(x,blocks,v,ids,w,circ) if len(v) else pd.DataFrame()
        for model in MODELS:
            g=folds.loc[folds.model==model] if not folds.empty else folds
            stats=aggregate(g)
            if not stats or stats["blocks"]<8:
                out.append(dict(**meta,model=model,status="insufficient_native_blocks",**stats,p_planned=1.))
            else:
                out.append(dict(**meta,model=model,status="estimated",**stats,p_planned=stats["p_block_sign"]))

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
        run_species(sp,sn,mapping,out)
    d=pd.DataFrame(out)
    assert len(d)==64
    d["holm_all64"]=multipletests(d.p_planned,method="holm")[1]
    d["diagnostic_pass"]=(d.status.eq("estimated") & (d.q2_cell_weighted>0) & (d.q2_block_equal>0) &
                          (d.positive_blocks>=.5*d.blocks) & (d.holm_all64<.05))
    d.to_csv(OUT/"native_prediction_diagnostic.csv",index=False)
    print("Complete; aggregate outputs only.",flush=True)

if __name__=="__main__":
    main()
