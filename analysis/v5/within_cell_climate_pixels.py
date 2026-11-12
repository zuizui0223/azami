"""V5 S1 independent ~1-km climate-pixel contrasts WITHIN the same taxon×cell.

Existing V5 analyses use one median per taxon×1°/5° geographic cell.
This new source-authenticated supplemental analysis retains the SAME
157,822 complete native photographs, but compares independently sampled
CHELSA raster pixels within their respective taxon×cell groups.

No data-dependent choice of signs, strata or climatic axes. No photograph,
taxon identity, user ID, exact coordinate, pixel ID or group ID is exported.
All 5 annual climate predictors are original CHELSA v2.1 pixel values.
Photo within-pixel replication is collapsed to one median phenotype;
independent contrast support comes from different raster pixels.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from analysis.v4.s1_climate_trait_hypervolume import (
    AXES, PHENOTYPE, geo_aggregate, load_source
)

EXTRA=("chelsa_bio04","chelsa_bio15")
CLIMATE=(
    "chelsa_bio12","chelsa_bio01","chelsa_rsds_mean",
    "chelsa_bio04","chelsa_bio15"
)
OUTCOME=(
    "orientation_image_vertical_angle","corolla_lab_chroma"
)
FOCAL=(
    ("rain_vs_angle","chelsa_bio12","orientation_image_vertical_angle"),
    ("rain_seasonality_vs_chroma","chelsa_bio15","corolla_lab_chroma")
)
# Exact original CHELSA V2.1 1981-2010 affine transform in
# analysis/v4/environment_source_metadata.json (EPSG:4326).
X_LEFT=-180.00013888885002
Y_TOP=83.99986041515001
PIX_DEG=0.0083333333
WIDTH=43200
HEIGHT=20880
PHOTO_N=157822
N_BOOT=399
MIN_PIXELS_PER_TAXON_CELL=5
MIN_PHOTOS_PER_TAXON_CELL=10
MIN_TAXON_CELLS=1
AXES_REPORTED=("BIO12","BIO01","RSDS","BIO04","BIO15")
MAX_CONDITION=1e8
ROOT=Path("work/v5_within_cell")
OUT=ROOT/"public"
GROUPCOLS=("_taxon","_latcell","_loncell")
PIXCOLS=(*GROUPCOLS,"_pixel_row","_pixel_col")


def _climate_pixel_index(frame):
    lon=frame.longitude.to_numpy(dtype=float)
    lat=frame.latitude.to_numpy(dtype=float)
    col=np.floor((lon-X_LEFT)/PIX_DEG).astype(np.int64)
    row=np.floor((Y_TOP-lat)/PIX_DEG).astype(np.int64)
    if ((col<0)|(col>=WIDTH)|(row<0)|(row>=HEIGHT)).any():
        raise ValueError("V5_WITHIN_CELL_ORIGINAL_CHELSA_PIXEL_INDEX_OOB")
    return row,col


def _read_frozen_photo_population(source,receipt):
    photo,proof=load_source(source,receipt,EXTRA)
    if len(photo)!=PHOTO_N:
        raise ValueError("V5_WITHIN_CELL_SOURCE_PHOTO_COHORT_CHANGED")
    if not np.isfinite(photo[list(CLIMATE)+list(OUTCOME)]).all().all():
        raise ValueError("V5_WITHIN_CELL_NONFINITE_FROZEN_CLIMATE_OR_TRAIT")
    return photo,proof


def build_pixel_cells(photo, degrees, *, min_pixels=MIN_PIXELS_PER_TAXON_CELL,
                      min_photos=MIN_PHOTOS_PER_TAXON_CELL):
    """Within-cell independent climate locations; no photo-count weights."""
    if degrees not in (1,5):raise ValueError("BAD_S1_GRID_SIZE")
    cells,original_support=geo_aggregate(photo,degrees,EXTRA)
    # Keep only taxon×cell groups which were eligible in the original V5
    # multivariate analysis; do not change the original study population.
    g=photo.rename(columns={"accepted_key":"_taxon"}).copy()
    g["_latcell"]=np.floor((g.latitude.to_numpy(float)+90)/degrees).astype(int)
    g["_loncell"]=np.floor((g.longitude.to_numpy(float)+180)/degrees).astype(int)
    pixels=(_climate_pixel_index(g))
    g["_pixel_row"],g["_pixel_col"]=pixels
    keys=cells.rename(columns={"accepted_key":"_taxon"})[list(GROUPCOLS)].drop_duplicates()
    g=g.merge(keys,on=list(GROUPCOLS),how="inner",validate="many_to_one")
    if g.empty:raise ValueError("NO_ORIGINAL_S1_ELIGIBLE_GEOGRAPHIC_CELLS")
    group=g.groupby(list(PIXCOLS),sort=True,observed=True)
    pix=group.agg(
        **{v:(v,"median") for v in (*CLIMATE,*OUTCOME)},
        n_photos=("latitude","size"),
        lat=("latitude","median"),
        lon=("longitude","median")
    ).reset_index()
    # An original climate raster pixel must carry a single physical climate
    # vector across all photographs assigned to that raster pixel.
    for v in CLIMATE:
        spread=group[v].max().to_numpy()-group[v].min().to_numpy()
        if np.nanmax(np.abs(spread))>1e-6:
            raise ValueError("INCONSISTENT_ORIGINAL_RASTER_PIXEL_CLIMATE_"+v)
    sg=pix.groupby(list(GROUPCOLS),sort=True,observed=True)
    support=sg.agg(n_pixels=("n_photos","size"),
                   n_photos=("n_photos","sum")).reset_index()
    yes=support[(support.n_pixels>=min_pixels)&(support.n_photos>=min_photos)]
    px=pix.merge(yes,on=list(GROUPCOLS),how="inner",validate="many_to_one",
                 suffixes=("","_cell"))
    if px.empty:
        raise ValueError("NO_TAXON_CELLS_WITH_REPLICATE_CLIMATE_PIXELS")
    # Every taxon×cell is weighted equally within its taxon, and each
    # taxon contributes equal overall weight; pixels are equal within cell.
    groupcount=px[list(GROUPCOLS)].drop_duplicates().groupby("_taxon").size()
    px["_n_taxon_cells"]=px["_taxon"].map(groupcount)
    px["_weight"]=1/(
        px["_n_taxon_cells"] * px["n_pixels"] * len(groupcount))
    assert np.isclose(px["_weight"].sum(),1)
    # Each 1° or 5° cell is assigned a region block via its centroid.
    px["_blocklat"]=np.floor((px.lat.to_numpy(float)+90)/20).astype(int)
    px["_blocklon"]=np.floor((px.lon.to_numpy(float)+180)/20).astype(int)
    px["_block"]=px["_blocklat"].astype(str)+"_"+px["_blocklon"].astype(str)
    return px,{
        "degrees":degrees,
        "source_photo_n":int(len(photo)),
        "source_taxa_original_supported":int(cells.accepted_key.nunique()),
        "source_taxon_cells_original_supported":int(len(cells)),
        "eligible_source_photos_in_original_supported_cells":int(len(g)),
        "photos_supported_by_independent_pixels":int(px.n_photos.sum()),
        "supported_taxon_cells":int(len(yes)),
        "supported_taxa":int(px._taxon.nunique()),
        "supported_unique_taxon_cell_pixels":int(len(px)),
        "supported_geographic_blocks":int(px._block.nunique()),
        "median_pixels_per_taxon_cell":float(yes.n_pixels.median()),
        "median_photos_per_taxon_cell":float(yes.n_photos.median()),
        "fraction_of_original_taxon_cells_with_pixel_support":
            float(len(yes)/len(cells)),
        "min_different_original_pixels":int(min_pixels),
        "min_photos_per_group":int(min_photos),
        "original_photo_eligibility_unchanged":True,
        "same_original_pixel_each_photo_is_not_independent":True
    }


def _group_demean(px):
    """Weighted within taxon×geographical-cell means, pixel as observation."""
    names=[*CLIMATE,*OUTCOME]
    # Within each group the contribution of each original raster pixel
    # is exactly equal. Repeated photos cannot inflate the weight.
    means=px.groupby(list(GROUPCOLS),sort=True,observed=True)[names].transform("mean")
    residual=px[names].to_numpy(dtype=float)-means.to_numpy(dtype=float)
    wt=px._weight.to_numpy(dtype=float)
    sd=np.sqrt((wt[:,None]*residual**2).sum(axis=0))
    if not np.isfinite(sd).all():
        raise ValueError("NONFINITE_CLIMATE_WITHIN_CELL_SPREAD")
    if (sd<=1e-12).any():
        raise ValueError("ONE_CLIMATE_OR_TRAIT_HAS_NO_WITHIN_CELL_SUPPORT")
    z=residual/sd
    x=z[:,:len(CLIMATE)]
    y=z[:,len(CLIMATE):]
    return x,y,wt,sd


def _cluster_moments(x,y,w,labels):
    codes,unique=pd.factorize(labels,sort=True)
    n=len(unique)
    xx=np.zeros((n,x.shape[1],x.shape[1]),float)
    xy=np.zeros((n,x.shape[1],y.shape[1]),float)
    np.add.at(xx,codes,w[:,None,None]*x[:,:,None]*x[:,None,:])
    np.add.at(xy,codes,w[:,None,None]*x[:,:,None]*y[:,None,:])
    uni=np.zeros((n,x.shape[1],x.shape[1]),float)
    uny=np.zeros((n,x.shape[1],y.shape[1]),float)
    for k in range(x.shape[1]):
        uni[:,k,k]=xx[:,k,k]
        uny[:,k,:]=xy[:,k,:]
    return xx,xy,unique


def _solve_boot(xx,xy,cidx,tidx):
    """Coefficient vector for two focal endpoints; no phenotype-driven selection."""
    totalx=np.sum(xx,axis=0)
    totaly=np.sum(xy,axis=0)
    if not np.isfinite(totalx).all():return None
    cond=np.linalg.cond(totalx)
    try:
        b=np.linalg.solve(totalx,totaly) if cond<MAX_CONDITION else np.full((5,2),np.nan)
    except np.linalg.LinAlgError:
        b=np.full((5,2),np.nan)
    focal_ind=np.array([totaly[cidx[0],tidx[0]]/totalx[cidx[0],cidx[0]],
                        totaly[cidx[1],tidx[1]]/totalx[cidx[1],cidx[1]]])
    focal_adj=np.array([b[cidx[0],tidx[0]],b[cidx[1],tidx[1]]])
    return focal_ind,focal_adj,cond


def _bootstrap_group(xx,xy,labels_seed,nboot=N_BOOT):
    rng=np.random.default_rng(labels_seed)
    n=len(xx)
    if n<10:return np.empty((0,2)),np.empty((0,2))
    blockcounts=rng.multinomial(n,np.full(n,1/n),size=nboot)
    cidx=(CLIMATE.index(FOCAL[0][1]),CLIMATE.index(FOCAL[1][1]))
    tidx=(OUTCOME.index(FOCAL[0][2]),OUTCOME.index(FOCAL[1][2]))
    out_one=[];out_multi=[]
    for counts in blockcounts:
        momentsx=np.tensordot(counts,xx,axes=(0,0))
        momentsy=np.tensordot(counts,xy,axes=(0,0))
        got=_solve_boot(momentsx[None,:,:],momentsy[None,:,:],cidx,tidx)
        if got is None:continue
        out_one.append(got[0]);out_multi.append(got[1])
    return np.asarray(out_one,float),np.asarray(out_multi,float)


def _summarise_ci(bootstrap,index):
    a=bootstrap[:,index] if bootstrap.ndim==2 and bootstrap.shape[1]>index else np.array([])
    a=a[np.isfinite(a)]
    if len(a)<.8*N_BOOT:return {"status":"NOT_ENOUGH_BOOTSTRAPS","valid_draws":int(len(a))}
    a95=np.quantile(a,[.025,.975])
    a975=np.quantile(a,[.0125,.9875])
    return {"status":"ESTIMATED","valid_draws":int(len(a)),
            "ci95":[float(a95[0]),float(a95[1])],
            "ci97_5_two_endpoint_bonferroni":[float(a975[0]),float(a975[1])]}


def fit_within_cells(px,*,nboot=N_BOOT):
    x,y,w,sd=_group_demean(px)
    gram=x.T@(w[:,None]*x)
    cross=x.T@(w[:,None]*y)
    cidx=(CLIMATE.index(FOCAL[0][1]),CLIMATE.index(FOCAL[1][1]))
    tidx=(OUTCOME.index(FOCAL[0][2]),OUTCOME.index(FOCAL[1][2]))
    estimates=_solve_boot(gram[None,:,:],cross[None,:,:],cidx,tidx)
    condition=float(estimates[2])
    condition_warning=condition>1000
    cov=gram/np.sqrt(np.outer(np.diag(gram),np.diag(gram)))
    try:vifs=np.diag(np.linalg.inv(cov))
    except np.linalg.LinAlgError:vifs=np.full(5,np.nan)
    xt,yt,taxa=_cluster_moments(x,y,w,px._taxon)
    xb,yb,blocks=_cluster_moments(x,y,w,px._block)
    # Distinct cluster systems; confidence is only called robust when BOTH
    # taxon- and geographic-block bootstraps exclude zero after the
    # family-wise 2-endpoint Bonferroni interval adjustment.
    nboot_before=N_BOOT
    bt1,bt5=_bootstrap_group(xt,yt,20261112,nboot)
    bb1,bb5=_bootstrap_group(xb,yb,20261113,nboot)
    results=[]
    for j,(name,predictor,response) in enumerate(FOCAL):
        taxon_ci=_summarise_ci(bt1,j)
        block_ci=_summarise_ci(bb1,j)
        taxon_partial=_summarise_ci(bt5,j)
        block_partial=_summarise_ci(bb5,j)
        def robust(a,b):
            if a["status"]!="ESTIMATED" or b["status"]!="ESTIMATED":
                return "INSUFFICIENT_CLUSTER_SUPPORT"
            def sign(z):
                lo,hi=z["ci97_5_two_endpoint_bonferroni"]
                if lo>0:return "positive"
                if hi<0:return "negative"
                return "includes_zero"
            a_s,b_s=sign(a),sign(b)
            return a_s if a_s==b_s and a_s!="includes_zero" else "NOT_TWO_CLUSTER_SUPPORTED"
        results.append({
            "hypothesis":name,"climate":predictor,"trait":response,
            "univariate_within_cell_std_beta":float(estimates[0][j]),
            "five_axis_partial_within_cell_std_beta":
                float(estimates[1][j]) if np.isfinite(estimates[1][j]) else None,
            "univariate_taxon_cluster":taxon_ci,
            "univariate_geographic_block_cluster":block_ci,
            "five_axis_taxon_cluster":taxon_partial,
            "five_axis_geographic_block_cluster":block_partial,
            "two_cluster_support_univariate":robust(taxon_ci,block_ci),
            "two_cluster_support_five_axis":robust(taxon_partial,block_partial),
            "within_cell_SD_climate_source_units":float(sd[CLIMATE.index(predictor)]),
            "within_cell_SD_trait_source_units":float(sd[len(CLIMATE)+OUTCOME.index(response)]),
            "family_size_for_bonferroni_ci":2,
            "not_a_causal_or_heritable_response":True
        })
    return {
        "results":results,
        "axes":list(AXES_REPORTED),
        "condition_number_of_five_axis_within_cell_gram":condition,
        "condition_number_over_1000":condition_warning,
        "axis_VIF_after_within_cell_demeaning":
            {k:float(z) if np.isfinite(z) else None for k,z in zip(AXES_REPORTED,vifs)},
        "n_taxon_clusters":int(len(taxa)),
        "n_20degree_geographic_clusters":int(len(blocks)),
        "n_cluster_bootstrap_draws":int(nboot_before),
        "cell_pixel_within_group_demeaning":True,
        "equal_weight_per_original_climate_pixel_within_taxon_cell":True,
        "equal_total_weight_per_taxon":True,
        "no_raw_site_or_taxon_identifiers_exported":True
    }


def run(table:Path,receipt:Path,out:Path,*,boot=N_BOOT):
    photo,proof=_read_frozen_photo_population(table,receipt)
    report={
        "schema":"V5_NATIVE_S1_WITHIN_GEOGRAPHIC_CELL_ORIGINAL_CLIMATE_PIXELS_20261112_V1",
        "source_authentication":"decoded_S1_SHA256_AND_OBSERVATION_COUNT_VERIFIED",
        "source_photo_count":len(photo),
        "source_decoded_sha256":proof["decoded_sha256"],
        "observed_climate_resolution":"original_CHELSA_30_arcsecond_pixel",
        "climate_normal_period":"1981-2010",
        "adopted_five_axes":list(AXES_REPORTED),
        "method":"taxon_by_angular_cell_fixed_effects_and_pixel_medians",
        "spatial_units":{},
        "outcomes_not_fitness_or_genetics":True,
        "posthoc_supplement_to_existing_V5":True,
        "no_original_photo_coordinates_taxon_ids_exported":True,
    }
    for deg in (1,5):
        px,support=build_pixel_cells(photo,deg)
        fit=fit_within_cells(px,nboot=boot)
        report["spatial_units"][str(deg)+"degree"]={
            "support":support,"inference":fit}
        print({"degrees":deg,
            "supported_groups":support["supported_taxon_cells"],
            "supported_taxa":support["supported_taxa"],
            "independent_pixel_groups":support["supported_unique_taxon_cell_pixels"],
            "conditions":fit["condition_number_of_five_axis_within_cell_gram"],
            "endpoints":[{k:r[k] for k in (
                "hypothesis","univariate_within_cell_std_beta",
                "five_axis_partial_within_cell_std_beta",
                "two_cluster_support_univariate",
                "two_cluster_support_five_axis")}
                for r in fit["results"]]},flush=True)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,ensure_ascii=False,allow_nan=False,indent=2)+"\n")
    return report


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--receipt",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    run(args.source,args.receipt,args.out)


if __name__=="__main__":main()
