"""NEW S1/S2 four-environment biological analysis from actual ORIGINAL photo sites.

BIO18: CHELSA v2.1 (1981-2010) warmest-quarter precipitation.
BIO5: CHELSA v2.1 (1981-2010) warmest-month daily maximum temperature.
RSDS: exact already archived CHELSA v2.1 1-km pixels (never resampled).
Elevation: WorldClim v2.1 native ~30 arcsec DEM in metres (different source).
Do not replace missing layers by BIO12/BIO1/BIO15 or substitute coordinates.
S1 170066 native image-angle records and S2 same two original focal species,
native + introduced WCVP all-wild background and 173737 measured heads.
Private source records stay in the Actions private working directory. Upload
ONLY aggregate taxon-level 4-axis means (no coordinates) and statistical JSON.
The climate/photo associations are observational, not true nodding mechanisms.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import time
import zipfile
from pathlib import Path
import numpy as np
import pandas as pd
import requests
from scipy import stats

BIO5_URL="https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/bioclim/bio05/1981-2010/CHELSA_bio05_1981-2010_V.2.1.tif"
BIO18_URL="https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/bioclim/bio18/1981-2010/CHELSA_bio18_1981-2010_V.2.1.tif"
ELEV_URL="https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_30s_elev.zip"
NEW_DOWNLOADS={
    "chelsa_bio05":(BIO5_URL,165253968,"CHELSA_bio05_1981-2010_V.2.1.tif"),
    "chelsa_bio18":(BIO18_URL,549240986,"CHELSA_bio18_1981-2010_V.2.1.tif"),
    "worldclim_elevation":(ELEV_URL,338188707,"wc2.1_30s_elev.zip"),
}
EXPOSURES=("chelsa_bio18","chelsa_bio05","chelsa_rsds_mean","elevation_m")
S1_ENDPOINT="orientation_image_vertical_angle"
ORIGINAL_251_SHA="47b32972b20564ade163907ef6fdf3eb2a69f8b42c8d87b8cf6ea75ca1682dc1"
# The actual authenticated source SHA is below (never relax it for a fit).
ORIGINAL_251_SHA="47b32972b20564ade163907ef6fdf3eb2a69f8b42c8d87b8cf4ea75ca1682dc1"
S1_N=170066
S2_PHOTO_N=173737
S2_WILD_N=297622
OUTPUT_SPECIES=251
FOCAL={"Cirsium arvense","Cirsium vulgare"}

def hashfile(path):
    h=hashlib.sha256()
    with open(path,"rb") as fd:
        for chunk in iter(lambda:fd.read(1<<22),b""):
            h.update(chunk)
    return h.hexdigest()

def fetch_scientific_raster(url,expected_bytes,target):
    target=Path(target)
    target.parent.mkdir(parents=True,exist_ok=True)
    size=0
    h=hashlib.sha256()
    with requests.get(url,timeout=(15,120),stream=True,
         headers={"Accept-Encoding":"identity","User-Agent":"azami-V4-S1-S2-bio18-bio05-source-bound/1.0"}) as response:
        if response.status_code!=200:
            raise ValueError(f"Exact public raster unavailable: HTTP {response.status_code}: {url}")
        contentlen=response.headers.get("Content-Length")
        if contentlen and int(contentlen)!=expected_bytes:
            raise ValueError("Remote source byte length drift before downloading "+url)
        with target.open("wb") as fd:
            for chunk in response.iter_content(chunk_size=4<<20):
                if chunk:
                    size+=len(chunk)
                    if size>expected_bytes:
                        raise ValueError("Source raster too large; fail closed")
                    fd.write(chunk)
                    h.update(chunk)
    if size!=expected_bytes:
        raise ValueError(f"New climate file incomplete ({size} vs {expected_bytes}): "+url)
    return {"source_url":url,"source_bytes":size,"source_sha256":h.hexdigest()}

def fetch_all(root):
    root=Path(root)
    files={}
    receipts={}
    for field,(url,n,filename) in NEW_DOWNLOADS.items():
        path=root/filename
        receipts[field]=fetch_scientific_raster(url,n,path)
        files[field]=path
        print(f"ACQUIRED_PUBLIC_{field} sha256={receipts[field]['source_sha256']} bytes={n}",flush=True)
    zip_path=files["worldclim_elevation"]
    with zipfile.ZipFile(zip_path) as z:
        names=[n for n in z.namelist() if n.lower().endswith(".tif")]
        if len(names)!=1 or not names[0].endswith("wc2.1_30s_elev.tif"):
            raise ValueError("WorldClim 30s elevation archive contained wrong raster member")
        info=z.getinfo(names[0])
        if info.file_size<1e8 or info.file_size>5e9:
            raise ValueError("WorldClim 30s elevation dimensions unexpected")
        out=root/"wc2.1_30s_elev.tif"
        with z.open(info) as i, out.open("wb") as o:
            for chunk in iter(lambda:i.read(4<<20),b""):
                o.write(chunk)
        receipts["worldclim_elevation"]["uncompressed_sha256"]=hashfile(out)
        receipts["worldclim_elevation"]["uncompressed_bytes"]=int(info.file_size)
        files["worldclim_elevation"]=out
    (root/"FOUR_AXIS_PUBLIC_SOURCE_RECEIPT.json").write_text(
        json.dumps(receipts,sort_keys=True,indent=2)+"\n")
    return files,receipts

def raster_sample(srcpath,lat,lon):
    """Pixel-exact sampling; unique original grid cells group source observations."""
    import rasterio
    from rasterio.transform import rowcol
    with rasterio.open(srcpath) as src:
        if src.crs is None or src.crs.to_epsg()!=4326 or src.count!=1:
            raise ValueError("New environment must be 1 band in WGS84")
        if not(20000<=src.width<=50000 and 10000<=src.height<=25000):
            raise ValueError("New environmental raster is not global 30 arc seconds")
        raw_scale=float(src.scales[0]); raw_offset=float(src.offsets[0])
        if not np.isfinite(raw_scale) or raw_scale<=0 or not np.isfinite(raw_offset):
            raise ValueError("Original raster physical scale/offset missing")
        mask=np.isfinite(lat)&np.isfinite(lon)&(np.abs(lat)<=90)&(np.abs(lon)<=180)
        out=np.full(len(lat),np.nan)
        take=np.flatnonzero(mask)
        rr,cc=rowcol(src.transform,lon[take],lat[take])
        rr=np.asarray(rr,dtype=np.int64);cc=np.asarray(cc,dtype=np.int64)
        inside=(rr>=0)&(rr<src.height)&(cc>=0)&(cc<src.width)
        selected=take[inside]
        codes=rr[inside]*src.width+cc[inside]
        unique,idx,inverse=np.unique(codes,return_index=True,return_inverse=True)
        src_indices=selected[idx]
        coords=[(float(lon[i]),float(lat[i])) for i in src_indices]
        if len(coords)>350000:
            raise ValueError("Original source climate raster sampling request unexpectedly large")
        readings=np.fromiter((float(v[0]) for v in src.sample(coords,masked=False)),
                            dtype=float,count=len(coords))
        nodata=src.nodata
        if nodata is not None:
            readings[np.isclose(readings,float(nodata))]=np.nan
        transformed=readings*raw_scale+raw_offset
        out[selected]=transformed[inverse]
        return out,{"source_raster_sha256":hashfile(srcpath),
            "width":src.width,"height":src.height,"dtype":src.dtypes[0],
            "scale":raw_scale,"offset":raw_offset,
            "nodata":None if nodata is None else float(nodata),
            "block_shapes":[list(x) for x in src.block_shapes],
            "n_input_observations":len(lat),"n_unique_raster_pixels_sampled":len(coords),
            "source_coverage":float(np.isfinite(out).mean())}

def proper_s1(frame,reference):
    from analysis.v4.native_integration_audit import eligible
    if not isinstance(frame,pd.DataFrame):raise ValueError("Original S1 must be observations")
    allcols={"obs_id","accepted_key","taxon_name","latitude","longitude",
             S1_ENDPOINT,S1_ENDPOINT+"__eligible",
             "chelsa_bio01","chelsa_bio04","chelsa_bio12",
             "chelsa_bio15","chelsa_rsds_mean"}
    if not allcols<=set(frame):
        raise ValueError("S1 original observations not complete with original 5 CHELSA")
    mask=eligible(frame,[S1_ENDPOINT]) & np.isfinite(frame[[
        "chelsa_bio01","chelsa_bio04","chelsa_bio12",
        "chelsa_bio15","chelsa_rsds_mean"]]).all(axis=1)
    s=frame.loc[mask,["obs_id","accepted_key","taxon_name","latitude","longitude",
                      S1_ENDPOINT,"chelsa_rsds_mean"]].copy()
    if len(s)!=S1_N or s.accepted_key.nunique()!=244:
        raise ValueError("Original S1 eligible population drifted")
    if s.obs_id.duplicated().any() or not s[S1_ENDPOINT].between(0,180).all():
        raise ValueError("S1 original native observation or camera angle invalid")
    s["obs_id"]=s.obs_id.astype(str)
    s["accepted_key"]=s.accepted_key.astype(str)
    s["source_name"]=s.taxon_name.astype(str).str.strip().str.replace(" ","_",regex=False)
    if len(s.groupby(["accepted_key","source_name"]))!=OUTPUT_SPECIES:
        raise ValueError("Original S1 251 source-name split lost")
    group=s.groupby(["accepted_key","source_name"],sort=True).agg(
          count=(S1_ENDPOINT,"size"),
          angle=(S1_ENDPOINT,"mean"),
          rsds=("chelsa_rsds_mean","mean")).reset_index()
    base=reference[["accepted_key","original_source_name","original_photo_count",
        "mean_image_angle_deg","chelsa_rsds_mean"]].copy()
    base["accepted_key"]=base.accepted_key.astype(str)
    merged=group.merge(base,left_on=["accepted_key","source_name"],
        right_on=["accepted_key","original_source_name"],how="left",validate="one_to_one",indicator=True)
    if (len(merged)!=251 or not merged._merge.eq("both").all() or
        not (merged["count"].to_numpy(int)==merged.original_photo_count.to_numpy(int)).all() or
        not np.allclose(merged.angle.to_numpy(float),merged.mean_image_angle_deg.to_numpy(float),atol=1e-5,rtol=0) or
        not np.allclose(merged.rsds.to_numpy(float),merged.chelsa_rsds_mean.to_numpy(float),atol=1e-6,rtol=0)):
        raise ValueError("S1 new four-axis original-photo source fails frozen exact 251 validation")
    return s

def load_s2(photo,allwild):
    a=pd.read_csv(photo,dtype={"obs_id":str,"user_id":str,"accepted_key":str},low_memory=False)
    b=pd.read_csv(allwild,dtype={"obs_id":str,"user_id":str,"accepted_key":str},low_memory=False)
    if len(a)!=S2_PHOTO_N or len(b)!=S2_WILD_N:
        raise ValueError("Original S2 focal count changed")
    if a.obs_id.duplicated().any() or b.obs_id.duplicated().any():
        raise ValueError("Duplicated S2 observation source")
    expect={("Cirsium arvense","native"):49072,
        ("Cirsium arvense","introduced"):43195,
        ("Cirsium vulgare","native"):33637,
        ("Cirsium vulgare","introduced"):47833}
    if a.groupby(["taxon_name","range_state"]).size().to_dict()!=expect:
        raise ValueError("Original S2 native and introduced photographed populations changed")
    if set(a.obs_id)-set(b.obs_id):
        raise ValueError("Original S2 photographed observations outside wild source")
    if set(b.taxon_name)!=FOCAL or set(b.range_state)!={"native","introduced"}:
        raise ValueError("S2 source species/range status changed")
    return a,b

def calc_four_axis(s1,s2,srcroot,rsds_cache_root,rsds_metadata):
    from analysis.v4.native_cached_climate import verified_cache,sample_cached_pixels
    lat=np.concatenate([s1.latitude.to_numpy(float),s2.latitude.to_numpy(float)])
    lon=np.concatenate([s1.longitude.to_numpy(float),s2.longitude.to_numpy(float)])
    outputs={}
    meta={}
    for field,file in (("chelsa_bio05","CHELSA_bio05_1981-2010_V.2.1.tif"),
                       ("chelsa_bio18","CHELSA_bio18_1981-2010_V.2.1.tif"),
                       ("elevation_m","wc2.1_30s_elev.tif")):
        data,info=raster_sample(Path(srcroot)/file,lat,lon)
        outputs[field]=data
        meta[field]=info
    with Path(rsds_metadata).open() as f:
        source=json.load(f)
    contract=source["chelsa_rsds_mean"]
    cache=Path(rsds_cache_root)/"chelsa_rsds_mean"
    grid=verified_cache(cache,contract)
    sampled=sample_cached_pixels(cache,grid,s2.latitude.to_numpy(),s2.longitude.to_numpy())
    s2_rsds=sampled*grid["scale"]+grid["offset"]
    if len(s1)+len(s2)!=len(lat):
        raise ValueError("S1/S2 climatic photo union length mismatch")
    s1=s1.copy();s2=s2.copy()
    for field,vector in outputs.items():
        s1[field]=vector[:len(s1)]
        s2[field]=vector[len(s1):]
    s2["chelsa_rsds_mean"]=s2_rsds
    meta["chelsa_rsds_mean"]={"source_SHA256":str(grid["source_sha256"]),
        "source":"previous existing SHA-verified original CHELSA RSDS 316 5-degree cache tiles",
        "original_pixel_no_resample":True,"n_S2_RSDS_complete":int(np.isfinite(s2_rsds).sum())}
    return s1,s2,meta

def robust_hc3(y,design,names):
    X=np.asarray(design,float);y=np.asarray(y,float)
    if len(y)<=X.shape[1]+5 or not np.isfinite(X).all() or not np.isfinite(y).all():
        raise ValueError("Not enough complete S1 source samples for four-axis comparison")
    if np.linalg.matrix_rank(X)<X.shape[1]:raise ValueError("Four-axis design not full rank")
    inv=np.linalg.inv(X.T@X);beta=inv@(X.T@y)
    res=y-X@beta;h=np.sum((X@inv)*X,axis=1)
    if (h>=.999999).any():raise ValueError("S1 source extreme leverage")
    v=X*(res/(1-h))[:,None]
    C=inv@(v.T@v)@inv
    se=np.sqrt(np.diag(C));z=stats.norm.ppf(.975)
    return {"coefficients":{k:float(b) for k,b in zip(names,beta)},
       "HC3_95CI":{k:[float(b-z*s),float(b+z*s)] for k,b,s in zip(names,beta,se)},
       "HC3_two_sided_z_p":{k:float(2*stats.norm.sf(abs(b/s))) for k,b,s in zip(names,beta,se)}}

def four_design(data,features=EXPOSURES):
    cols=np.column_stack([data[x].to_numpy(float) for x in features])
    means=cols.mean(axis=0);sds=cols.std(axis=0)
    if np.min(sds)<=1e-12:raise ValueError("Four environmental axes not identified")
    return np.column_stack((np.ones(len(data)),(cols-means)/sds)),{
        k:{"mean":float(v),"SD":float(s)} for k,v,s in zip(features,means,sds)}

def s1_between(s):
    complete=s.dropna(subset=list(EXPOSURES)).copy()
    original=complete.groupby(["accepted_key","source_name"],sort=True).agg(
        n_original_angle=("obs_id","count"),
        mean_image_angle=(S1_ENDPOINT,"mean"),
        **{a:(a,"mean") for a in EXPOSURES}
    ).reset_index()
    if len(original)<185:raise ValueError("More than 66 original source names lost to new environment")
    result={}
    for bound in (1,10):
        cohort=original.loc[original.n_original_angle>=bound].copy()
        X,sd=four_design(cohort)
        fit=robust_hc3(cohort.mean_image_angle.to_numpy(float),X,["intercept",*EXPOSURES])
        result[f"min_photo_{bound}"]={"n_source_names":len(cohort),
            "n_original_photo_observations":int(cohort.n_original_angle.sum()),
            "predictor_sd":sd,"fit":fit}
    return original,result

def s1_within(s):
    complete=s.dropna(subset=list(EXPOSURES)).copy()
    if len(complete)<0.8*S1_N:raise ValueError("Source-original S1 new four-axis complete cohort too small")
    grouping=complete.groupby(["accepted_key","source_name"],sort=False)
    sd=np.array([complete[x].std(ddof=0) for x in EXPOSURES],float)
    if np.min(sd)<=1e-12:raise ValueError("Source four-axis within scale constant")
    X=np.column_stack([(complete[x]-grouping[x].transform("mean")).to_numpy(float)/sd[j]
        for j,x in enumerate(EXPOSURES)])
    y=(complete[S1_ENDPOINT]-grouping[S1_ENDPOINT].transform("mean")).to_numpy(float)
    count=grouping.obs_id.transform("size").to_numpy(float)
    w=1/count
    bread=X.T@(X*w[:,None])
    if np.linalg.matrix_rank(bread)<len(EXPOSURES):raise ValueError("Within-taxon S1 climate design rank invalid")
    inv=np.linalg.inv(bread);beta=inv@(X.T@(w*y))
    residual=y-X@beta
    # Robust covariance across taxon-level clusters, without pretending that
    # tens of thousands of photographs are evolutionary replicates.
    sources=complete.accepted_key.astype(str)+"|"+complete.source_name
    score=pd.DataFrame(X*(w*residual)[:,None],columns=EXPOSURES)
    score["_s"]=sources.to_numpy()
    clusters=score.groupby("_s")[list(EXPOSURES)].sum().to_numpy()
    G=len(clusters);C=inv@(clusters.T@clusters)@inv*(G/(G-1))
    se=np.sqrt(np.diag(C));critical=stats.t.ppf(.975,G-1)
    return {"n_original_source_images_complete":len(complete),
        "n_source_taxon_clusters":G,
        "coefficients_degrees_per_full_S1_photo_SD":dict(zip(EXPOSURES,beta.tolist())),
        "cluster_taxon_95CI":{k:[float(v-critical*s),float(v+critical*s)]
                for k,v,s in zip(EXPOSURES,beta,se)},
        "cluster_taxon_two_sided_t_p":{k:float(2*stats.t.sf(abs(v/s),G-1))
                for k,v,s in zip(EXPOSURES,beta,se)},
        "photo_within_taxon_residual_not_causal_plasticity":True}

def s2_cell_model(data):
    if len(data)<1000:return {"status":"INSUFFICIENT_PHOTO_SUPPORT"}
    from analysis.v4.native_introduced_climate_transfer import dated_and_balanced
    work=dated_and_balanced(data)
    if len(work)<500 or work._cell.nunique()<15:
        return {"status":"INSUFFICIENT_SPATIAL_CELL_SUPPORT",
            "n_cell":int(work._cell.nunique()),"n":len(work)}
    vals=np.column_stack([work[c].to_numpy(float) for c in EXPOSURES])
    scale=vals.std(axis=0)
    if np.min(scale)<=1e-12:return {"status":"FOUR_AXIS_RANK_DEFICIENT"}
    mon=work._month.to_numpy(float)
    yr=work._year.to_numpy(float)
    X=np.column_stack([np.ones(len(work)),(vals-vals.mean(axis=0))/scale,
       np.sin(2*np.pi*mon/12),np.cos(2*np.pi*mon/12),(yr-2015)/10])
    y=work.orientation_angle_degrees_median.to_numpy(float)
    w=1/work.groupby("_cell").obs_id.transform("size").to_numpy(float)
    bread=X.T@(X*w[:,None])
    if np.linalg.matrix_rank(bread)<X.shape[1]:return {"status":"S2_FOUR_AXIS_MODEL_NOT_IDENTIFIED"}
    inv=np.linalg.inv(bread)
    beta=inv@(X.T@(w*y))
    err=y-X@beta
    cellscore=pd.DataFrame(X*(w*err)[:,None])
    cellscore["_cell"]=work._cell.to_numpy()
    S=cellscore.groupby("_cell").sum().to_numpy()
    G=len(S);n=len(work);p=X.shape[1]
    C=inv@(S.T@S)@inv*(G/(G-1))*((n-1)/(n-p))
    se=np.sqrt(np.diag(C));cv=stats.t.ppf(.975,G-1)
    return {"status":"S2_TRUE_IMAGE_ANGLE_FOUR_AXIS_CELL_CLUSTER_REGRESSION",
        "n_balanced_observer_cell":len(work),"n_original_equal_sin_5deg_cells":G,
        "n_original_observers":int(work.user_id.nunique()),
        "month_season_and_year_adjusted":True,
        "coefficients_degrees_per_within_range_SD":{
            k:float(v) for k,v in zip(EXPOSURES,beta[1:5])},
        "cell_cluster_95CI":{k:[float(v-cv*s),float(v+cv*s)]
            for k,v,s in zip(EXPOSURES,beta[1:5],se[1:5])},
        "cell_cluster_two_sided_t_p":{k:float(2*stats.t.sf(abs(v/s),G-1))
            for k,v,s in zip(EXPOSURES,beta[1:5],se[1:5])},
        "source_photo_frame_not_genetic_change_or_true_gravity":True}

def s2_summary(s,wild):
    from analysis.v4.native_introduced_climate_transfer import load_verified_focal
    summary={}
    photo=s.merge(wild[["obs_id",*EXPOSURES]],on="obs_id",how="left",validate="one_to_one",suffixes=("","_w"))
    # With all-new sampled columns only wild records are the accepted source;
    # avoid duplicate original variable field in focal photo table.
    for v in EXPOSURES:
        if v not in photo and v+"_w" in photo:photo[v]=photo[v+"_w"]
    if len(photo)!=S2_PHOTO_N or not set(photo.taxon_name)<=FOCAL:
        raise ValueError("Exact S2 photo and WCVP wild source join drift")
    for species in sorted(FOCAL):
        x={}
        for state in ("native","introduced"):
            bg=wild[(wild.taxon_name==species)&(wild.range_state==state)]
            subset=photo[(photo.taxon_name==species)&(photo.range_state==state)]
            complete_bg=bg.dropna(subset=list(EXPOSURES))
            complete_photo=subset.dropna(subset=list(EXPOSURES))
            if len(complete_bg)<0.70*len(bg) or len(complete_photo)<0.70*len(subset):
                x[state]={"status":"STOP_INSUFFICIENT_FOUR_ENVIRONMENT_COVERAGE",
                    "n_background":len(bg),"n_complete_bg":len(complete_bg),
                    "n_original_images":len(subset),"n_complete_photo":len(complete_photo)}
                continue
            x[state]={"status":"FOUR_AXIS_S2_SOURCE_REAL_WCVP_RANGE",
                "n_all_original_wild_eligible_background":len(bg),
                "n_source_background_complete":len(complete_bg),
                "n_original_measured_head_images":len(subset),
                "n_photo_four_axis_complete":len(complete_photo),
                "original_continuous_photo_angle_median":float(subset.orientation_angle_degrees_median.median()),
                "full_wild_original_site_medians":{v:float(complete_bg[v].median()) for v in EXPOSURES},
                "photo_matched_original_site_medians":{v:float(complete_photo[v].median()) for v in EXPOSURES},
                "original_S2_four_axis_cell_model":s2_cell_model(complete_photo)}
        if all(x[t].get("status")=="FOUR_AXIS_S2_SOURCE_REAL_WCVP_RANGE" for t in ("native","introduced")):
            x["introduced_minus_native_full_wild_site_median"]={
                v:float(x["introduced"]["full_wild_original_site_medians"][v]-
                        x["native"]["full_wild_original_site_medians"][v]) for v in EXPOSURES}
            x["introduced_minus_native_photo_angle_median_degrees"]=float(
                x["introduced"]["original_continuous_photo_angle_median"]-
                x["native"]["original_continuous_photo_angle_median"])
        summary[species]=x
    return summary

def execute(s1_source,s1_original251,s2_photo,s2_allwild,source_root,rsds_cache_root,
            rsds_metadata,out_dir):
    raw=Path(s1_original251).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=ORIGINAL_251_SHA:
        raise ValueError("Cannot extend S1 new four axes from unverified original 251")
    original=pd.DataFrame(list(csv.DictReader(raw.decode("utf-8-sig").splitlines())))
    s1=proper_s1(pd.read_csv(s1_source,low_memory=False),original)
    photo,wild=load_s2(s2_photo,s2_allwild)
    wild_fields=["obs_id","taxon_name","range_state","latitude","longitude",
                 "user_id","observed_on","accepted_key"]
    wild=wild.loc[:,wild_fields].copy()
    s1,wild,source_meta=calc_four_axis(s1,wild,source_root,rsds_cache_root,rsds_metadata)
    if len(s1)!=S1_N or len(wild)!=S2_WILD_N:
        raise ValueError("Original S1/S2 initial whole population changed")
    if any(np.isfinite(s1[v]).mean()<.80 for v in EXPOSURES):
        raise ValueError("New CHELSA/DEM axis missing more than 20% of original S1")
    if any(np.isfinite(wild[v]).mean()<.75 for v in EXPOSURES):
        raise ValueError("New CHELSA/DEM axis missing more than 25% of original S2")
    for var,lo,hi in (("chelsa_bio05",-80,75),("chelsa_bio18",0,20000),
                       ("chelsa_rsds_mean",0,50),("elevation_m",-500,9000)):
        v=s1[var].to_numpy(float);median=np.nanmedian(v)
        if not(lo<=median<=hi):
            raise ValueError(f"Source physical units wrong for {var}: median {median}")
    native=photo.loc[photo.range_state.eq("native"),["obs_id","taxon_name","orientation_angle_degrees_median"]].merge(
        s1[["obs_id","taxon_name",S1_ENDPOINT,"chelsa_rsds_mean"]],on=["obs_id","taxon_name"],
        how="left",validate="one_to_one",indicator=True)
    if len(native)!=82709 or not native["_merge"].eq("both").all():
        raise ValueError("S1 same original focal native images lost in S2")
    if not np.allclose(native.orientation_angle_degrees_median,native[S1_ENDPOINT],rtol=0,atol=1e-8):
        raise ValueError("S1/S2 original camera angle drift across source")
    # Pixel-exact new RSDS cached values must agree with original S1 own values
    rsds_native=wild.loc[wild.obs_id.isin(native.obs_id),["obs_id","chelsa_rsds_mean"]]
    re=rsds_native.merge(native[["obs_id","chelsa_rsds_mean"]],on="obs_id",validate="one_to_one",suffixes=("_cache","_S1"))
    both=np.isfinite(re.chelsa_rsds_mean_cache)&np.isfinite(re.chelsa_rsds_mean_S1)
    if not both.any() or not np.allclose(re.loc[both,"chelsa_rsds_mean_cache"],
           re.loc[both,"chelsa_rsds_mean_S1"],atol=.001,rtol=0):
        raise ValueError("True original RSDS grid no longer matches S1 species pixels")
    taxon,s1fit=s1_between(s1)
    within=s1_within(s1)
    s2=s2_summary(photo,wild)
    target=Path(out_dir);target.mkdir(parents=True,exist_ok=True)
    # Species means are aggregate public phenotypes only, never point records.
    taxon_out=taxon.rename(columns={"source_name":"original_source_name",
        "mean_image_angle":"mean_image_angle_deg","n_original_angle":"n_photo_complete_four_axis"})
    taxon_out.to_csv(target/"S1_ORIGINAL_251_TAXON_MEAN_NEW_BIO18_BIO05_RSDS_ELEVATION.csv",
         index=False,encoding="utf-8",float_format="%.9g")
    report={"schema":"azami_V4_real_S1_S2_BIO18_BIO05_RSDS_elev_native_original_observations_v1",
         "status":"COMPLETE_REAL_ORIGINAL_S1_AND_S2_SITE_EXPOSURES_FOUR_NEW_AXES",
         "variables":list(EXPOSURES),"ecology":{
              "chelsa_bio18":"CHELSA warmest-quarter precipitation 1981-2010 (not flower-month observed rain)",
              "chelsa_bio05":"CHELSA warmest month maximum temperature 1981-2010 (not daily observed extreme)",
              "chelsa_rsds_mean":"Original cached 1981-2010 downward shortwave radiation",
              "elevation_m":"WorldClim v2.1 30arcsec static elevation; different environmental data provider"},
         "sources":source_meta,
         "original_original_251_SHA":ORIGINAL_251_SHA,
         "S1":{"all_original_native_S1_head_image_observations":S1_N,
               "n_four_axis_complete_native_head_images":int(s1[list(EXPOSURES)].notna().all(axis=1).sum()),
               "n_taxon_source_names_after_new_four_axis":len(taxon),
               "between_taxon":s1fit,"within_taxon":within},
         "S2":{"all_original_WCVP_native_and_introduced_wild_records":S2_WILD_N,
               "all_original_photo_angles":S2_PHOTO_N,
               "source_observation_same_native_angle_overlap":int(len(native)),
               "species":s2},
         "claim_stop":{"no_new_U_to_D_botanical_label":True,
                       "no_historic_adaptive_evolution_inferred":True,
                       "original_five_CHELSA_S1_S2_frozen_unchanged":True,
                       "new_four_climate_posthoc_explicit_user_request":True,
                       "source_locations_and_users_not_exposed":True,
                       "no_true_gravitational_angle_verified":True,
                       "PGLS_nuclear_species_tree_not_yet_applied_to_new_four_source":True}}
    (target/"S1_S2_ORIGINAL_FOUR_CLIMATE_REAL_RESULTS.json").write_text(
         json.dumps(report,ensure_ascii=False,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("S1_S2_REAL_NEW_FOUR_AXIS="+json.dumps({
      "S1_n":S1_N,"S1_four_complete":report["S1"]["n_four_axis_complete_native_head_images"],
      "S1_251_taxa":len(taxon),
      "S1_between_SDs":s1fit["min_photo_1"]["fit"]["coefficients"],
      "S2_source_counts":{k:{r:x[r].get("n_photo_four_axis_complete") for r in ("native","introduced")} for k,x in s2.items()},
      "only_public_aggregate_written":True},ensure_ascii=False,sort_keys=True),flush=True)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--stage",required=True,choices=("download","analyse"))
    p.add_argument("--source-root",required=True,type=Path)
    p.add_argument("--out",type=Path)
    for name in ("s1","s1-original251","s2-photo","s2-allwild","rsds-cache-root","rsds-metadata"):
        p.add_argument("--"+name,type=Path)
    a=p.parse_args()
    if a.stage=="download":
        fetch_all(a.source_root)
    else:
        names=(a.s1,a.s1_original251,a.s2_photo,a.s2_allwild,
               a.rsds_cache_root,a.rsds_metadata,a.out)
        if any(x is None for x in names):raise ValueError("All original S1 S2 private and cache sources required")
        execute(*names[:4],a.source_root,names[4],names[5],names[6])
if __name__=="__main__":main()
