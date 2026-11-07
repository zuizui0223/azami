"""V4 H3c: test outcome-blind selected climate variables in original S1.

BLIND variable list is frozen in the X-only selection JSON before observing
H3c traits. Seven selected climate variables are converted to 3 unsupervised
PCs (no flower measurements enter transform). Compare against the original
BIO12/BIO01/RSDS representation using identical source-QC taxon/cell support.
Climate-only quantile calibration matches expected pairwise-analogue frequency.
The established H3c source run and its published P remain untouched.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
from analysis.v4.species_environment_existing_inputs import load_existing_inputs
from analysis.v4 import s1_h3_climate_analogue_alignment as h3
from analysis.v4.s1_climate_trait_hypervolume import PHENOTYPE,CLIMATE

SEED=20261107
DRAWS=24
PERMS=199
def aggregate(d,climate_cols,degree):
    y=list(PHENOTYPE)
    env=list(dict.fromkeys([*CLIMATE,*climate_cols]))
    needed={"obs_id","accepted_key","latitude","longitude",*y,*env,
            *(v+"__eligible" for v in y)}
    if not needed<=set(d):raise ValueError("H3_ORIGINAL_SOURCE_OR_CLIMATE_MISSING")
    valid=(d.accepted_key.notna()&d.latitude.between(-90,90)&d.longitude.between(-180,180))
    valid=valid.to_numpy(bool)
    values=d[[*y,*env]].apply(pd.to_numeric,errors="coerce").to_numpy(float)
    valid &= np.isfinite(values).all(axis=1)
    for v in y:valid &= d[v+"__eligible"].eq(True).fillna(False).to_numpy(bool)
    valid &= d[y[0]].between(0,180).fillna(False).to_numpy(bool)
    selected=d.loc[valid,["accepted_key","latitude","longitude",*env,*y]].copy()
    selected["cell_lat"]=np.floor((selected.latitude+90)/degree).astype(int)
    selected["cell_lon"]=np.floor((selected.longitude+180)/degree).astype(int)
    groups=selected.groupby(["accepted_key","cell_lat","cell_lon"],observed=True,sort=True)
    cols=[*env,*y,"latitude","longitude"]
    cells=groups[cols].median().reset_index()
    cells["n_photos"]=groups.size().to_numpy(int)
    counts=cells.groupby("accepted_key",sort=True,observed=True).agg(
        photos=("n_photos","sum"),cell_n=("cell_lat","size"))
    qualified=counts.index[(counts.photos>=40)&(counts.cell_n>=8)]
    cells=cells.loc[cells.accepted_key.isin(qualified)].copy()
    if cells.accepted_key.nunique()<25:raise ValueError("H3_NEW_CLIMATE_SUPPORT_TOO_SMALL")
    return cells,{"source_joint_QC_photos_complete_X":int(len(selected)),
                  "source_joint_QC_taxa":int(selected.accepted_key.nunique()),
                  "eligible_8cells_40photos_taxa":int(len(qualified)),
                  "eligible_taxon_cells":int(len(cells))}

def standardized(cells,names,w):
    X=cells[list(names)].to_numpy(float)
    mu=np.average(X,axis=0,weights=w)
    sd=np.sqrt(np.average((X-mu)**2,axis=0,weights=w))
    if not np.isfinite(sd).all() or (sd<1e-10).any():
        raise ValueError("H3_XONLY_CLIMATE_OR_PHENOTYPE_VARIANCE_INVALID")
    return (X-mu)/sd

def climate_pcs(cells,selected,w):
    X=standardized(cells,selected,w)
    _,singular,vt=np.linalg.svd(X*np.sqrt(w[:,None]),full_matrices=False)
    raw=X@vt[:3].T
    scales=np.sqrt(np.average(raw**2,weights=w,axis=0))
    if (scales<1e-10).any():raise ValueError("H3_CLIMATE_PCA_DEGENERATE")
    return raw/scales,{
        "selected_source_climate_columns":list(selected),
        "first_three_PC_variance_fraction":float((singular[:3]**2).sum()/(singular**2).sum()),
        "PC_weights":vt[:3].tolist(),
        "PCA_outcome_blind":True}

def _prelim_pair_distances(x,groups,sector,seed):
    rng=np.random.default_rng(seed)
    picks=[rng.choice(len(arr),size=8,replace=False) for arr in groups]
    result=[]
    for i in range(len(groups)):
        for j in range(i+1,len(groups)):
            if sector[i]==sector[j]:continue
            a=x[groups[i][picks[i]]];b=x[groups[j][picks[j]]]
            result.extend(np.linalg.norm(a[:,None,:]-b[None,:,:],axis=2).ravel())
    return np.asarray(result,float)

def calibrate_threshold(cells,x_original,x_selected,seed):
    ix=[i.to_numpy(int) for _,i in cells.groupby("accepted_key",sort=True).groups.items()]
    # groupby indices are original row indices; caller resets row index.
    sectors=[]
    for _,g in cells.groupby("accepted_key",sort=True):
        lon=g.longitude.to_numpy(float)
        codes=h3.geographic_sector(lon)
        sectors.append(int(np.bincount(codes,minlength=3).argmax()))
    sec=np.asarray(sectors,int)
    d0=_prelim_pair_distances(x_original,ix,sec,seed)
    d1=_prelim_pair_distances(x_selected,ix,sec,seed)
    q=float(np.mean(d0<=h3.CLIMATE_ANALOG_Z))
    q=np.clip(q,.002,.998)
    threshold=float(np.quantile(d1,q))
    return threshold,{"baseline_pairwise_analogue_fraction":float(q),
                      "outcome_blind_X_distance_quantile_calibration":True,
                      "selected_distance_threshold":threshold}

def matched_pairs_raw_nd(x:np.ndarray,y:np.ndarray,sector:np.ndarray,threshold:float):
    """Same H3c pair/support rules, allowing 7-dimensional CLIMATE and 3D Y."""
    if (x.ndim!=3 or y.ndim!=3 or x.shape[:2]!=y.shape[:2] or
        x.shape[1]!=8 or y.shape[2]!=3 or x.shape[2]<3 or
        sector.shape!=(len(x),) or threshold<=0):
        raise ValueError("H3_DIMENSIONAL_PAIR_SHAPE_UNSUPPORTED")
    ii=[];jj=[];distances=[];masks=[];counts=[];sector_types=[]
    for i in range(len(x)):
        for j in range(i+1,len(x)):
            if sector[i]==sector[j]:continue
            d_env=np.linalg.norm(x[i,:,None,:]-x[j,None,:,:],axis=2)
            matched=d_env<=threshold
            if (int(matched.sum())<h3.MIN_ANALOG_CELL_PAIRS or
                np.flatnonzero(matched.any(axis=1)).size<h3.MIN_DISTINCT_CELLS_PER_SIDE or
                np.flatnonzero(matched.any(axis=0)).size<h3.MIN_DISTINCT_CELLS_PER_SIDE):
                continue
            d_y=np.linalg.norm(y[i,:,None,:]-y[j,None,:,:],axis=2)
            ii.append(i);jj.append(j)
            distances.append(d_y);masks.append(matched)
            counts.append(int(matched.sum()))
            sector_types.append(tuple(sorted((int(sector[i]),int(sector[j])))))
    if not ii:
        return None,{"status":"no_analogous_pair_support","n_pairs":0,
                     "n_sector_pair_types":0}
    n_sector_types=len(set(sector_types))
    meta={"status":"qualified" if
          len(ii)>=h3.MIN_ELIGIBLE_TAXON_PAIRS and
          n_sector_types>=h3.MIN_SECTOR_PAIR_TYPES
          else "insufficient_analogue_pairs",
          "n_pairs":len(ii),"n_sector_pair_types":n_sector_types,
          "n_analog_cell_pairs":int(sum(counts)),
          "mean_analog_cell_pairs_per_taxon_pair":float(np.mean(counts)),
          "n_taxa_represented":len(set(ii)|set(jj))}
    pair={"i":np.asarray(ii,int),"j":np.asarray(jj,int),
          "d":np.stack(distances),"mask":np.stack(masks),
          "analog_n":np.asarray(counts,float)}
    return pair,meta


def _compare(cells,climate_columns,degree):
    cells=cells.reset_index(drop=True)
    counts=cells.groupby("accepted_key",sort=True,observed=True).size()
    w=cells.accepted_key.map(1/counts).to_numpy(float)
    w=w/w.sum()
    xorig=standardized(cells,CLIMATE,w)
    xnew,pca=climate_pcs(cells,climate_columns,w)
    xraw=standardized(cells,climate_columns,w)
    Y=standardized(cells,PHENOTYPE,w)
    threshold,calibrate=calibrate_threshold(cells,xorig,xnew,SEED+degree)
    raw_threshold,raw_calibration=calibrate_threshold(cells,xorig,xraw,SEED+degree)
    taxa_groups=[g.index.to_numpy(int) for _,g in cells.groupby("accepted_key",sort=True,observed=True)]
    sectors=np.asarray([int(np.bincount(h3.geographic_sector(
                        cells.loc[ix,"longitude"].to_numpy(float)),minlength=3).argmax())
                        for ix in taxa_groups],int)
    if len(taxa_groups)<25:raise ValueError("H3_PAIR_TAXA_INSUFFICIENT")
    # Same rarefied source cells in both climate representations within draw.
    rng=np.random.default_rng(SEED+20+degree)
    metrics={name:{"adv":[],"null":[],"support":[],"overlap":[],"overlap_null":[]}
             for name in ("original3","blind_selected_PCA3","blind_selected_raw7")}
    for rep in range(DRAWS):
        take=[rng.choice(ix,8,replace=False) for ix in taxa_groups]
        y=np.stack([Y[ids] for ids in take])
        for name,X,cutoff in (("original3",xorig,h3.CLIMATE_ANALOG_Z),
                              ("blind_selected_PCA3",xnew,threshold),
                              ("blind_selected_raw7",xraw,raw_threshold)):
            xx=np.stack([X[ids] for ids in take])
            if name=="blind_selected_raw7":
                pairs,support=matched_pairs_raw_nd(xx,y,sectors,cutoff)
            else:
                previous=h3.CLIMATE_ANALOG_Z
                h3.CLIMATE_ANALOG_Z=cutoff
                try:
                    pairs,support=h3.matched_pairs(xx,y,sectors)
                finally:
                    h3.CLIMATE_ANALOG_Z=previous
            metrics[name]["support"].append(support)
            if support["status"]!="qualified":continue
            nullrng=np.random.default_rng(SEED+degree*100000+rep*3)
            observed,null=h3.observed_and_permutations(pairs,len(taxa_groups),
                                                      nullrng,reps=PERMS)
            # KDE uses same y and matched masks, with a separately seeded null.
            kdenullrng=np.random.default_rng(SEED+degree*200000+rep*5)
            overlap,overlap_null=h3.kde_overlap_and_permutations(pairs,y,
                                            kdenullrng,reps=PERMS)
            metrics[name]["adv"].append(observed)
            metrics[name]["null"].append(null)
            metrics[name]["overlap"].append(overlap["conditional_minus_unconditional_mean"])
            metrics[name]["overlap_null"].append(overlap_null)
    result={}
    for name,data in metrics.items():
        support=[s for s in data["support"] if s["status"]=="qualified"]
        row={"n_eligible_source_taxa":int(len(taxa_groups)),
             "sector_counts":{str(k):int((sectors==k).sum()) for k in range(3)},
             "successful_rarefied_draws":int(len(support)),
             "n_requested_draws":DRAWS,
             "mean_eligible_taxon_pairs":float(np.mean([s["n_pairs"] for s in support])) if support else None,
             "mean_eligible_analog_cell_pairs":float(np.mean([s["n_analog_cell_pairs"] for s in support])) if support else None}
        if len(support)<12:
            row["status"]="insufficient_matched_support"
            result[name]=row;continue
        obs=float(np.mean(data["adv"]))
        null=np.mean(np.stack(data["null"]),axis=0)
        nexp=float(null.mean())
        p=float((1+np.sum(np.abs(null-nexp)>=abs(obs-nexp)))/(PERMS+1))
        overlap=float(np.mean(data["overlap"]))
        onull=np.mean(np.stack(data["overlap_null"]),axis=0)
        onexp=float(onull.mean())
        po=float((1+np.sum(np.abs(onull-onexp)>=abs(overlap-onexp)))/(PERMS+1))
        adj=multipletests([p,po],method="holm")[1]
        row.update(status="estimated",
                   trait_distance_advantage=float(obs),
                   null_expected_advantage=float(nexp),
                   two_sided_exploratory_permutation_p=p,
                   two_metric_holm_p=float(adj[0]),
                   KDE_conditional_minus_all=float(overlap),
                   KDE_null_expected=float(onexp),
                   KDE_two_sided_exploratory_p=po,
                   KDE_two_metric_holm_p=float(adj[1]),
                   evidence_is_observational_not_adaptation=True)
        result[name]=row
    return {"degree":degree,"matched_source_support":int(len(taxa_groups)),
            "source_climate_axes_metadata":pca,
            "climate_distance_support_calibration":{
                "selected_PCA3":calibrate,"selected_raw7":raw_calibration},
            "tested_original_selected_seven_dimensions_without_PCA":True,
            "results":result}

def run_existing(native,decoded_receipt,seven,nine,nine_manifest,source_receipt,
                 source_sha,selected,output):
    frame,rec=load_existing_inputs(native,decoded_receipt,seven,nine,
                                   source_receipt,nine_manifest,
                                   source_receipt_sha256=source_sha,expected_n=313350)
    geo=pd.read_csv(native,usecols=["obs_id","latitude","longitude"],
                    dtype={"obs_id":str},low_memory=False)
    if len(geo)!=len(frame) or not geo.obs_id.is_unique:raise ValueError("H3_GEO_SOURCE_INVALID")
    original_ids=frame.obs_id.copy()
    frame=frame.merge(geo,on="obs_id",validate="one_to_one",how="left",sort=False)
    if not frame.obs_id.equals(original_ids):raise ValueError("H3_COHORT_JOIN_INVALID")
    blind=json.loads(Path(selected).read_text())
    if blind.get("status")!="original_S1_authenticated_X_only_blind_selection" or \
        blind.get("climate_source_receipt_sha256")!=source_sha:
        raise ValueError("H3_VARIABLE_SET_NOT_BLIND_OR_AUTHENTICATED")
    results={}
    for degree in (1,5):
        newvars=blind["section33"]["by_resolution"][str(degree)+"degree"]["selected_climate_X"]["selected_variables"]
        cells,support=aggregate(frame,newvars,degree)
        details=_compare(cells,newvars,degree)
        results[str(degree)+"degree"]={"source":support,**details}
    out={"schema":"azami_v4_s1_blind_X_climate_PCA3_H3c_conditional_alignment_v1",
         "status":"original_S1_Xonly_selection_post_selection_H3c_exploration",
         "original_source_observation_n":313350,
         "source_SHA256":rec["input_digests"]["native"],
         "receipt_sha256":source_sha,
         "H3_climate_variable_selection_outcome_blind":True,
         "new_climate_distance_PCA_uses_X_only":True,
         "new_analogue_cutoff_climate_distance_only_calibrated":True,
         "existing_H3c_original_Holm_result_untouched":True,
         "n_permutations_per_draw":PERMS,
         "results":results}
    banned={"obs_id","accepted_key","taxon","taxon_name","latitude","longitude",
            "user_id","photographer_id","photo_id","individual_prediction"}
    def privacy(x):
        if isinstance(x,dict):
            for k,v in x.items():
                if k in banned:raise ValueError("ORIGINAL_H3_PRIVATE_KEY")
                privacy(v)
        elif isinstance(x,list):
            for v in x:privacy(v)
    privacy(out)
    Path(output).parent.mkdir(parents=True,exist_ok=True)
    Path(output).write_text(json.dumps(out,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    return out

def main():
    a=argparse.ArgumentParser()
    for arg in ("native","decoded-receipt","seven","nine","nine-manifest",
                "source-receipt","source-sha","selected","output"):
        a.add_argument("--"+arg,required=True)
    p=a.parse_args()
    x=run_existing(p.native,p.decoded_receipt,p.seven,p.nine,
                   p.nine_manifest,p.source_receipt,p.source_sha,p.selected,p.output)
    print(json.dumps({"status":x["status"],"results":{
        d:{k:v["status"] for k,v in a["results"].items()}
        for d,a in x["results"].items()}}),flush=True)
if __name__=="__main__":main()
