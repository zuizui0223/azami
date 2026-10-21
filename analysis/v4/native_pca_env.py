"""Availability-corrected PCA-env-style native vs introduced climate niche pilot.

Strictly preserve distinction between the original 665115 OBSERVATION universe,
the 350-name WCVP authority map and the raster climate BACKGROUND available
within each TDWG3 native or explicitly introduced region. No head angle or
other image phenotype enters the niche estimation.

The historical five CHELSA v2.1 source TIFFs were extracted into an exact-pixel
cache for 316 5-degree storage tiles, initially focused on native occurrences.
For reproducible background comparisons we audit all *potential* WCVP polygon
mesh points and report full/partial cache coverage. This is PCA-env STYLE, not
a bitwise ecospat.grid.clim.dyn reproduction or a calibrated MaxEnt model.

PRIVATE: original sighting metadata, photo IDs, human users, exact occurrence
coordinates, TDWG geom JSON and all raster pixel arrays remain in the runner.
PUBLIC: named species and continent counts, PCA/environmental summary metrics.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
from shapely import points
from shapely.geometry import shape
from shapely.strtree import STRtree

from analysis.v4.native_cached_climate import (
    VARIABLES, sample_cached_pixels, verified_cache
)
from analysis.v4.native_introduced_opportunity import (
    authority_source, CONTINENTS
)
from analysis.v4.native_range_status import distribution_status

FROZEN_SOURCE_N=297622
SPECIES_KEYS={
    "Cirsium arvense":"206745077",
    "Cirsium vulgare":"206739917"
}
SOURCE_COUNTS={
    ("Cirsium arvense","native"):78646,
    ("Cirsium arvense","introduced"):70964,
    ("Cirsium vulgare","native"):62368,
    ("Cirsium vulgare","introduced"):85644
}
STEP_DEGREES=0.5
SEED=20261021
BINS=80
KDE_SIGMA_BINS=1.5
MIN_BACKGROUND_POINTS=200
MIN_OCCURRENCE_POINTS=100
FULL_BACKGROUND_COVERAGE_REQUIRED=0.8
INTRO_CONTINENT_MIN_OCC=80


def _native_and_introduced_codes(authority:dict,accepted_key:str)->dict:
    rows=authority["distributions"].get(accepted_key)
    if not isinstance(rows,list) or not rows:
        raise ValueError("No WCVP distribution for pinned accepted species "+accepted_key)
    by_code={}
    for rec in rows:
        key=str(rec.get("locationId","")).removeprefix("TDWG:").upper()
        status=distribution_status(rec)
        if not key:
            continue
        by_code.setdefault(key,set()).add(status)
    allowed={k:next(iter(s)) for k,s in by_code.items()
             if len(s)==1 and next(iter(s)) in ("native","introduced")}
    if not any(v=="native" for v in allowed.values()):
        raise ValueError("Species has no unambiguous WCVP-native available geography")
    if not any(v=="introduced" for v in allowed.values()):
        raise ValueError("Species has no unambiguous WCVP-introduced geography")
    return allowed


def background_mesh(authority:dict,*,step:float=STEP_DEGREES)->pd.DataFrame:
    """Global fixed lat/lon mesh, WCVP TDWG3 polygon query; no occurrence bias.

    Every potential cell centroid is accounted for, not merely 316 cached
    storage tiles. Ties to two TDWG3 codes are dropped (conservative).
    """
    if not .1<=step<=1.:
        raise ValueError("Geographical mesh step must be 0.1-1 degree")
    features=authority["geometry"]["features"]
    polys=[shape(feature["geometry"]) for feature in features]
    codes=sorted({str(f["properties"]["LEVEL3_COD"]).upper() for f in features})
    dictionary={v:i for i,v in enumerate(codes)}
    geometry_to_code=np.array([
        dictionary[str(f["properties"]["LEVEL3_COD"]).upper()]
        for f in features
    ],dtype=np.int32)
    geo_continent={}
    for f in features:
        k=str(f["properties"]["LEVEL3_COD"]).upper()
        continent=int(f["properties"]["LEVEL1_COD"])
        if continent not in CONTINENTS:
            raise ValueError("Unexpected TDWG Level1 background continent")
        if k in geo_continent and geo_continent[k]!=continent:
            raise ValueError("One TDWG region belongs to two continents")
        geo_continent[k]=continent

    # 0.5-degree world grid: 250,560 land-or-ocean points; world extent is the
    # latitude coverage of the original CHELSA 43200x20880 raster.
    lats=np.arange(-89.75,83.999999,step)
    lons=np.arange(-179.75,179.999999,step)
    lon,lat=np.meshgrid(lons,lats)
    la=lat.ravel()
    lo=lon.ravel()
    pts=points(lo,la)
    hits=STRtree(polys).query(pts,predicate="within")
    source_idx=hits[0].astype(int)
    poly_idx=hits[1].astype(int)
    region_id=np.full(len(la),-1,np.int32)
    if len(source_idx):
        order=np.argsort(source_idx,kind="mergesort")
        idx=source_idx[order]
        coded=geometry_to_code[poly_idx[order]]
        unique,offset=np.unique(idx,return_index=True)
        lower=np.minimum.reduceat(coded,offset)
        upper=np.maximum.reduceat(coded,offset)
        keep=lower==upper
        region_id[unique[keep]]=lower[keep]
    valid=region_id>=0
    land=pd.DataFrame({
        "latitude":la[valid],
        "longitude":lo[valid],
        "TDWG_code":np.asarray(codes,dtype=object)[region_id[valid]],
    })
    land["TDWG_level1"]=land.TDWG_code.map(geo_continent).astype(int)
    land["_cell_key"]=(
        ((land.latitude.to_numpy()+90)//5).astype(int).astype(str)+"_"+
        ((land.longitude.to_numpy()+180)//5).astype(int).astype(str)
    )
    land["_area_weight"]=np.cos(np.deg2rad(land.latitude.to_numpy(float)))
    if land.empty or (land._area_weight<=0).any():
        raise ValueError("No WCVP land climatology mesh")
    return land


def _occurrence_frame(path:Path)->pd.DataFrame:
    data=pd.read_csv(path,dtype={"obs_id":str,"user_id":str,"accepted_key":str},
                     low_memory=False)
    if len(data)!=FROZEN_SOURCE_N or data.obs_id.isna().any() or not data.obs_id.is_unique:
        raise ValueError("Original two-species native plus introduced wild roster changed")
    if set(data.taxon_name)!=set(SPECIES_KEYS):
        raise ValueError("Original source taxon names differ from frozen all-species audit")
    if set(data.range_state)!={"native","introduced"} or data.user_id.isna().any():
        raise ValueError("WCVP native/introduced status or original photographer missing")
    for (name,status),n in SOURCE_COUNTS.items():
        block=data.loc[data.taxon_name.eq(name)&data.range_state.eq(status)]
        if len(block)!=n:
            raise ValueError("Original source population changed for "+name+" "+status)
        if set(block.accepted_key)!={SPECIES_KEYS[name]}:
            raise ValueError("Original source WCVP accepted key changed")
    from analysis.v4.native_introduced_opportunity import angular_cells
    data["_cell"]=angular_cells(data.latitude,data.longitude)
    result=data.sort_values("obs_id").drop_duplicates(
        ["taxon_name","range_state","_cell","user_id"],keep="first"
    ).copy()
    if result.obs_id.duplicated().any():
        raise ValueError("Observer/cell balancing did not preserve unique original records")
    return result


def source_background_climate(
        authority:dict,cache_root:Path,contract_path:Path
    )->tuple[pd.DataFrame,dict]:
    landscape=background_mesh(authority)
    contracts=json.loads(contract_path.read_text(encoding="utf-8"))
    configs={}
    grid_sets=[]
    for var in VARIABLES:
        grid=verified_cache(cache_root/var,contracts[var])
        configs[var]=grid
        grid_sets.append(set(grid["tiles"]))
    cache_tiles=set.intersection(*grid_sets)
    if len(cache_tiles)!=316 or any(s!=cache_tiles for s in grid_sets):
        raise ValueError("Five frozen original-pixel rasters do not share identical 316 tiles")
    landscape["_cache_candidate"]=landscape._cell_key.isin(cache_tiles)
    sample=landscape.loc[landscape._cache_candidate].copy()
    for var in VARIABLES:
        grid=configs[var]
        raw=sample_cached_pixels(cache_root/var,grid,
                                 sample.latitude.to_numpy(float),
                                 sample.longitude.to_numpy(float))
        sample[var]=raw*float(grid["scale"])+float(grid["offset"])
    complete=np.isfinite(sample[list(VARIABLES)].to_numpy(float)).all(axis=1)
    sample=sample.loc[complete].copy()
    return sample,{
        "potential_background_mesh_points":len(landscape),
        "cached_candidate_background_mesh_points":int(landscape._cache_candidate.sum()),
        "complete_original_5variable_climate_background_points":len(sample),
        "original_source_raster_cache_tile_count":316,
        "background_mesh_spacing_degrees":STEP_DEGREES,
        "background_metadata":"Private WCVP TDWG3 polygon membership plus exact original raster pixels; 5-degree tiles are storage only",
    }


def _background_species_audit(authority:dict,background:pd.DataFrame,
                              species:str,land:pd.DataFrame)->tuple[pd.DataFrame,dict]:
    source_status=_native_and_introduced_codes(authority,SPECIES_KEYS[species])
    full=land.loc[land.TDWG_code.isin(source_status)].copy()
    full["_range"]=full.TDWG_code.map(source_status)
    usable=background.loc[background.TDWG_code.isin(source_status)].copy()
    usable["_range"]=usable.TDWG_code.map(source_status)
    audit={}
    for st in ("native","introduced"):
        all_=full.loc[full._range.eq(st)]
        ok=usable.loc[usable._range.eq(st)]
        if not len(all_):
            raise ValueError("TDWG region has no available raster-grid mesh points")
        full_weight=float(all_._area_weight.sum())
        covered=float(ok._area_weight.sum())
        audit[st]={
            "potential_area_weighted_mesh":full_weight,
            "cached_complete_area_weighted_mesh":covered,
            "complete_pixel_background_fraction_of_all_WCVP_available_area":covered/full_weight,
            "potential_mesh_points":len(all_),
            "complete_raster_cache_points":len(ok),
            "fully_available_environment_sampled":bool(covered/full_weight>=FULL_BACKGROUND_COVERAGE_REQUIRED),
        }
        for cont,g in all_.groupby("TDWG_level1",sort=True):
            q=ok.loc[ok.TDWG_level1.eq(cont)]
            audit[f"{st}:{CONTINENTS[int(cont)]}"]={
                "potential_mesh_points":len(g),
                "complete_cached_background_mesh_points":len(q),
                "fraction_original_background_tile_coverage_area_weighted":float(q._area_weight.sum()/g._area_weight.sum())
            }
    return usable,audit


def _sample_occupied_climate(sightings:pd.DataFrame,cache_root:Path,
                             contracts:dict)->pd.DataFrame:
    result=sightings.copy()
    for var in VARIABLES:
        g=verified_cache(cache_root/var,contracts[var])
        raw=sample_cached_pixels(cache_root/var,g,
                                 result.latitude.to_numpy(float),
                                 result.longitude.to_numpy(float))
        result[var]=raw*float(g["scale"])+float(g["offset"])
    return result


def _weights_per_5degree_cell(frame:pd.DataFrame)->np.ndarray:
    if frame.empty or frame._cell.isna().any():
        raise ValueError("No original observer×coarse-cell sightings")
    n=frame.groupby("_cell").obs_id.transform("size").to_numpy(float)
    if np.any(n<=0):
        raise ValueError("Original photo/cell weighting invalid")
    return 1/n


def _weighted_quantile(values,weights,probabilities=(.05,.5,.95)):
    x=np.asarray(values,float);w=np.asarray(weights,float)
    if len(x)!=len(w) or not len(x) or not np.isfinite(x).all() or not np.isfinite(w).all() or (w<=0).any():
        raise ValueError("Invalid original climate values/weights")
    ix=np.argsort(x,kind="stable")
    z=np.cumsum(w[ix])-0.5*w[ix]
    z/=w.sum()
    return [float(np.interp(q,z,x[ix],left=x[ix][0],right=x[ix][-1]))
            for q in probabilities]


def _weighted_projected_pca(values:np.ndarray,weights:np.ndarray):
    if values.ndim!=2 or values.shape[1]!=len(VARIABLES) or len(values)!=len(weights):
        raise ValueError("Five identical original climate predictors are required")
    if not np.isfinite(values).all() or not np.isfinite(weights).all() or (weights<=0).any():
        raise ValueError("Source climate PCA received non-finite values")
    w=weights/weights.sum()
    center=np.average(values,axis=0,weights=w)
    scale=np.sqrt(np.average((values-center)**2,axis=0,weights=w))
    if np.any(scale<=1e-10):
        raise ValueError("Original CHELSA pooled climate range lacks a variable")
    z=(values-center)/scale
    cov=(z*w[:,None]).T@z
    eig,evec=np.linalg.eigh(cov)
    order=np.argsort(eig)[::-1]
    eig=eig[order]
    basis=evec[:,order[:2]].copy()
    for k in range(2):
        j=int(np.argmax(np.abs(basis[:,k])))
        if basis[j,k]<0:
            basis[:,k]*=-1
    return center,scale,basis,[float(x) for x in eig[:2]/eig.sum()]


def _probability_content(z:np.ndarray,coverage:float=.95)->np.ndarray:
    a=np.asarray(z,float)
    if a.ndim!=2 or abs(a.sum()-1)>1e-5 or not 0<coverage<1:
        raise ValueError("PCA-env occupancy density not normalized")
    ix=np.argsort(a.ravel())[::-1]
    ordered=a.ravel()[ix]
    threshold=np.searchsorted(np.cumsum(ordered),coverage,side="left")
    out=np.zeros(a.size,bool)
    out[ix[:min(threshold+1,len(ix))]]=True
    return out.reshape(a.shape)


def _pca_env_2d(env_native:pd.DataFrame,env_intro:pd.DataFrame,
                occ_native:pd.DataFrame,occ_intro:pd.DataFrame)->dict:
    if min(len(env_native),len(env_intro))<MIN_BACKGROUND_POINTS:
        return {"status":"STOP_BACKGROUND_INSUFFICIENT_OR_RASTER_COVERAGE"}
    if min(len(occ_native),len(occ_intro))<MIN_OCCURRENCE_POINTS:
        return {"status":"STOP_OCCURRENCE_CLIMATE_SUPPORT_TOO_SMALL"}
    background=pd.concat([env_native,env_intro],ignore_index=True)
    weights=np.r_[
        env_native._area_weight.to_numpy(float)/env_native._area_weight.sum()/2,
        env_intro._area_weight.to_numpy(float)/env_intro._area_weight.sum()/2
    ]
    center,scale,load,variance=_weighted_projected_pca(
        background[list(VARIABLES)].to_numpy(float),weights)
    def transform(data):
        return ((data[list(VARIABLES)].to_numpy(float)-center)/scale)@load
    eb=transform(background)
    ob0=transform(occ_native)
    ob1=transform(occ_intro)
    low=np.quantile(eb,[.005],axis=0)[0]
    high=np.quantile(eb,[.995],axis=0)[0]
    margin=np.maximum((high-low)*.1,.1)
    edges=[np.linspace(low[j]-margin[j],high[j]+margin[j],BINS+1) for j in range(2)]
    def distribution(env,occ):
        ebasis=transform(env)
        oxyz=transform(occ)
        ew=env._area_weight.to_numpy(float)
        ew=ew/ew.sum()
        ow=_weights_per_5degree_cell(occ)
        ow=ow/ow.sum()
        hist_e=np.histogram2d(ebasis[:,0],ebasis[:,1],bins=edges,weights=ew)[0]
        hist_o=np.histogram2d(oxyz[:,0],oxyz[:,1],bins=edges,weights=ow)[0]
        if hist_e.sum()<.95 or hist_o.sum()<.80:
            raise ValueError("Most native/introduced observations outside common PCA reference bounds")
        smoothed_background=gaussian_filter(hist_e,KDE_SIGMA_BINS,mode="constant")
        smoothed_occurrence=gaussian_filter(hist_o,KDE_SIGMA_BINS,mode="constant")
        floor=max(float(smoothed_background.max())*.0001,1e-12)
        corrected=smoothed_occurrence/np.maximum(smoothed_background,floor)
        corrected[smoothed_background<floor]=0
        if corrected.sum()<=0 or smoothed_occurrence.sum()<=0:
            raise ValueError("No joint ecological occupancy support in background")
        corrected/=corrected.sum()
        raw=smoothed_occurrence/smoothed_occurrence.sum()
        return corrected,raw,{
            "photographed_occurrences_projected_inside_common_grid":float(hist_o.sum()),
            "available_environment_area_projected_inside_common_grid":float(hist_e.sum()),
            "background_density_floor_fraction_of_mode":.0001
        }
    z0,raw0,area0=distribution(env_native,occ_native)
    z1,raw1,area1=distribution(env_intro,occ_intro)
    support0=_probability_content(z0)
    support1=_probability_content(z1)
    overlap=1-float(np.abs(z0-z1).sum())/2
    raw_overlap=1-float(np.abs(raw0-raw1).sum())/2
    return {
        "status":"EXECUTED_PCA_ENV_STYLE_CACHED_BACKGROUND_NOT_EXACT_R_ECOSPAT",
        "pca_axes_explained_background_variance":variance,
        "shared_background_standardization_center":{
            v:float(x) for v,x in zip(VARIABLES,center)},
        "shared_background_standardization_sd":{
            v:float(x) for v,x in zip(VARIABLES,scale)},
        "two_climate_PC_loadings":{
            v:[float(x) for x in row] for v,row in zip(VARIABLES,load)},
        "common_density_grid_bins_per_axis":BINS,
        "gaussian_smoothing_bins":KDE_SIGMA_BINS,
        "Schoener_D_available_environment_corrected_density":float(np.clip(overlap,0,1)),
        "Schoener_D_uncorrected_photo_observation_density":float(np.clip(raw_overlap,0,1)),
        "introduced_probability_outside_native_95pct_probability_content":float(z1[~support0].sum()),
        "native_probability_outside_introduced_95pct_probability_content":float(z0[~support1].sum()),
        "reference_definition":"PCA pooled WCVP climate environments with area weighting and equal weight per range; occupancy density corrected by within-range smoothed climate availability",
        "range_native_projection":area0,
        "range_introduced_projection":area1,
        "not_ecospat_exact_implementation":True,
        "nonoverlap_probabilities_not_significance_tests":True,
    }


def climatic_shift_summary(native:pd.DataFrame,introduced:pd.DataFrame)->dict:
    if min(len(native),len(introduced))<MIN_OCCURRENCE_POINTS:
        return {"status":"STOP_OCCUPANCY_TOO_SMALL"}
    out={}
    weights={
        "native":_weights_per_5degree_cell(native),
        "introduced":_weights_per_5degree_cell(introduced),
    }
    for var in VARIABLES:
        a=_weighted_quantile(native[var].to_numpy(float),weights["native"])
        b=_weighted_quantile(introduced[var].to_numpy(float),weights["introduced"])
        iqr=np.diff(_weighted_quantile(native[var].to_numpy(float),weights["native"],(.25,.75)))[0]
        out[var]={
            "native_cell_equal_q05_q50_q95":a,
            "introduced_cell_equal_q05_q50_q95":b,
            "introduced_minus_native_median_original_units":b[1]-a[1],
            "introduced_minus_native_median_native_IQR_units":(b[1]-a[1])/iqr if iqr>0 else None,
            "CHELSA_units":"climatological, physical source units per original raster metadata"
        }
    return out


def run_pca_env(occurrences:pd.DataFrame,source:dict,
                cache_root:Path,source_metadata:Path)->dict:
    source_background,all_bg_audit=source_background_climate(source,cache_root,source_metadata)
    land=background_mesh(source)
    contracts=json.loads(source_metadata.read_text(encoding="utf-8"))
    occ_frame=_sample_occupied_climate(occurrences,cache_root,contracts)
    out={}
    for species in SPECIES_KEYS:
        env,audit=_background_species_audit(source,source_background,species,land)
        data=occ_frame.loc[occ_frame.taxon_name.eq(species)].copy()
        data["_complete"]=np.isfinite(data[list(VARIABLES)].to_numpy(float)).all(axis=1)
        support={}
        for st in ("native","introduced"):
            sub=data.loc[data.range_state.eq(st)]
            supported=sub.loc[sub._complete]
            support[st]={
                "original_observer_cell_balanced_sightings":len(sub),
                "original_five_CHELSA_pixels_complete":len(supported),
                "original_distinct_5degree_cells_with_complete_climate":int(
                    supported._cell.nunique()),
                "pixels_missing_not_imputed":len(sub)-len(supported),
            }
        n=data.loc[data.range_state.eq("native")&data._complete].copy()
        i=data.loc[data.range_state.eq("introduced")&data._complete].copy()
        en=env.loc[env._range.eq("native")]
        ei=env.loc[env._range.eq("introduced")]
        label="PARTIAL_CLIMATE_BACKGROUND_HOLDOUT_NOT_FULL_TDWG_PCA_ENV"
        if audit["native"]["fully_available_environment_sampled"] and audit["introduced"]["fully_available_environment_sampled"]:
            label="PCA_ENV_BACKGROUND_COVERAGE_AT_LEAST_80PCT"
        result={
            "status":label,
            "source_sighting_climate_coverage":support,
            "available_climatic_background_audit":audit,
            "five_original_climate_variable_shift":climatic_shift_summary(n,i),
            "pca_env_style":_pca_env_2d(en,ei,n,i),
            "introduced_continents":{}
        }
        for group,part in i.groupby("tdwg1",sort=True):
            continent=CONTINENTS[int(group)]
            eb=ei.loc[ei.TDWG_level1.eq(int(group))]
            denom=audit.get(f"introduced:{continent}",{})
            regional={"introduction_photo_occurrences_original_pixel_complete":len(part),
                      "background_potential_mesh_points":denom.get("potential_mesh_points",0),
                      "background_climate_complete_mesh_points":len(eb),
                      "available_background_geographic_area_fraction":denom.get(
                        "fraction_original_background_tile_coverage_area_weighted",0.)}
            if len(part)>=INTRO_CONTINENT_MIN_OCC and len(eb)>=MIN_BACKGROUND_POINTS:
                regional["pca_env_style"]=_pca_env_2d(en,eb,n,part)
                regional["five_original_climate_variable_shift"]=climatic_shift_summary(n,part)
                regional["status"]="EXECUTED_REGIONAL_CLIMATE_PCA_ENV_CACHED_SUPPORT_ONLY"
            else:
                regional["pca_env_style"]={"status":"STOP_BACKGROUND_OR_ORIGINAL_OCCURRENCE_NOT_AVAILABLE"}
                regional["status"]="STOP_SOURCE_CHELSA_COVERAGE"
            result["introduced_continents"][continent]=regional
        out[species]=result
    return {
        "schema_version":"v4_two_Cirsium_native_introduced_PCA_env_available_climate_20261021_v1",
        "status":"EXECUTED_GEOREFERENCED_AVAILABLE_BACKGROUND_AUDIT_AND_PCA_ENV_PILOT",
        "taxon_results":out,
        "source_background":all_bg_audit,
        "species_original_sightings_not_conditioned_on_photo_trait_eligibility":True,
        "original_sampling_universe":"native/explicit WCVP introduced wild species-rank original sightings; one per original observer×5degree cell per species and status",
        "raster_source":"same five CHELSA v2.1 frozen exact-pixel 1981-2010 source TIFF tile caches",
        "limitation":[
            "Accessible environment E is restricted to WCVP-TDWG3 classified land and original cached tiles. Reported area coverage is required before any full-niche claim.",
            "This is PCA-env-style smoothed and availability-corrected density, not byte-exact ecospat R PCA-env or MaxEnt.",
            "Environmental availability is approximated by deterministic 0.5-degree polygon-intersected mesh with latitude cosine area weights, not full 1km background.",
            "Original human-source observer×5degree-cell deduplication does not remove detector and geospatial taxonomic sampling bias.",
            "PCA-env density overlap and exclusive 95% probability-content mass are descriptive, not ecology-wide confidence intervals.",
            "Data availability in formerly native-targeted caches is spatially biased and may have poor introduced background completeness, including Africa and South America.",
            "Neither niche-shift patterns nor same-species geographic contrasts identify heritable head-angle changes, plasticity or fast evolutionary adaptation.",
        ],
        "private_individual_observations_or_user_coordinates_in_export":False
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("occurrences","authority","cache_root","metadata","output"):
        p.add_argument("--"+name.replace("_","-"),type=Path,required=True)
    args=p.parse_args()
    sample=_occurrence_frame(args.occurrences)
    a=authority_source(args.authority.read_bytes())
    result=run_pca_env(sample,a,args.cache_root,args.metadata)
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/"S2_two_species_native_introduced_pca_env.json").write_text(
        json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n",
        encoding="utf-8")
    print(json.dumps({
        "status":result["status"],
        "species":{name:{
            "availability_n":z["available_climatic_background_audit"]["native"]["complete_raster_cache_points"],
            "availability_i":z["available_climatic_background_audit"]["introduced"]["complete_raster_cache_points"],
            "background_coverage_i":z["available_climatic_background_audit"]["introduced"]["complete_pixel_background_fraction_of_all_WCVP_available_area"],
            "Schoener_D_cached_background":z["pca_env_style"].get("Schoener_D_available_environment_corrected_density")
        } for name,z in result["taxon_results"].items()}
    },ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
