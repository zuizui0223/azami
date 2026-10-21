"""S2 whole rainfall versus image-angle joint observational landscape.

The main S2 comparison remains native vs ALL introduced range, per species.
Joint original CHELSA BIO12 / photo-image-angle data test whether separated
marginal ECDFs concealed responses associated with precipitation at one site.
Climate plotting bins arise from the original wild source, independent of angle.
No ecological response metric, width or extreme threshold selected. Regional
TDWG1 data reconstruct the ALL-introduced joint distribution exactly. Only
aggregate output escapes the source-authenticated GitHub Actions runner.
"""
from __future__ import annotations

import json
from pathlib import Path
import numpy as np
import pandas as pd

from analysis.v4.native_cached_climate import VARIABLES
from analysis.v4.native_climate_trait_coupling import (
    original_cell_weights, weighted_quantiles
)
from analysis.v4.native_introduced_climate_transfer import dated_and_balanced
from analysis.v4.native_intro_continent_orientation import CONTINENTS

SPECIES=("Cirsium arvense","Cirsium vulgare")
RAIN="chelsa_bio12"
ANGLE="orientation_angle_degrees_median"
RAIN_DISPLAY_KNOTS=41
ANGLE_EDGES=np.linspace(0,180,37)
MIN_DRAW_ROWS=30
MIN_DRAW_CELLS=5
EXPECTED_SOURCE=144178
EXPECTED_PHOTO=173737


def _check(wild:pd.DataFrame,photo:pd.DataFrame)->None:
    for data in (wild,photo):
        required={"obs_id","user_id","taxon_name","range_state","tdwg1",
                  "_cell",RAIN}
        if not required<=set(data):
            raise ValueError("Original WCVP observation and climate metadata unavailable")
        if data.obs_id.duplicated().any() or data.user_id.isna().any():
            raise ValueError("Original source observation IDs or photographers invalid")
        if set(data.taxon_name)!=set(SPECIES) or set(data.range_state)!={"native","introduced"}:
            raise ValueError("Original WCVP native/introduced source classes changed")
    if not set(VARIABLES)<=set(wild) or ANGLE not in photo:
        raise ValueError("Original full CHELSA climates/photo-angle source missing")
    if len(wild)!=EXPECTED_SOURCE or len(photo)!=EXPECTED_PHOTO:
        raise ValueError("Original two-species source and angle archives changed")


def _climate_source(wild:pd.DataFrame)->pd.DataFrame:
    mask=np.isfinite(wild[list(VARIABLES)].to_numpy(float)).all(axis=1)
    return wild.loc[mask].copy()


def _measured(photo:pd.DataFrame)->pd.DataFrame:
    v=photo[ANGLE].to_numpy(float)
    rain=photo[RAIN].to_numpy(float)
    keep=np.isfinite(v)&(v>=0)&(v<=180)&np.isfinite(rain)
    data=dated_and_balanced(photo.loc[keep])
    if data.empty:
        raise ValueError("No authenticated, dated photographed capitula")
    return data


def _rain_grid(wild:pd.DataFrame,photo:pd.DataFrame)->np.ndarray:
    """Plotting axis derives from phenotype-unselected climate source."""
    groups=[]
    for label in ("native","introduced"):
        data=wild.loc[wild.range_state.eq(label)]
        weight=original_cell_weights(data)
        groups.append((data[RAIN].to_numpy(float),weight/weight.sum()/2))
    climate=np.concatenate([a for a,w in groups])
    weights=np.concatenate([w for a,w in groups])
    rank_cut=weighted_quantiles(
        climate,weights,np.linspace(0,1,RAIN_DISPLAY_KNOTS))
    p=photo[RAIN].to_numpy(float)
    low=min(float(np.min(climate)),float(np.min(p)))
    high=max(float(np.max(climate)),float(np.max(p)))
    edges=np.unique(np.r_[low,rank_cut,high])
    if len(edges)<12:
        raise ValueError("Insufficient climate display resolution")
    return edges


def _source_climate_probabilities(frame:pd.DataFrame,edges:np.ndarray)->list[float]:
    if frame.empty:
        raise ValueError("No original non-photo-selected climate reference")
    hist=np.histogram(frame[RAIN],bins=edges,weights=original_cell_weights(frame))[0]
    if hist.sum()<=0:
        raise ValueError("Climate background bins have no mass")
    return (hist/hist.sum()).tolist()


def _raw_histogram(frame:pd.DataFrame,edges:np.ndarray,
                   weights:np.ndarray)->tuple[np.ndarray,np.ndarray,np.ndarray]:
    """Same-observation CHELSA rain versus image-angle histogram."""
    rain=frame[RAIN].to_numpy(float)
    angle=frame[ANGLE].to_numpy(float)
    if (not len(frame) or len(weights)!=len(frame) or
        not np.isfinite(weights).all() or (weights<0).any() or
        not np.isfinite(rain).all() or not np.isfinite(angle).all() or
        ((angle<0)|(angle>180)).any()):
        raise ValueError("Invalid original rainfall × head angle photo sample")
    hist=np.histogram2d(rain,angle,bins=(edges,ANGLE_EDGES),
                        weights=weights)[0]
    if not np.isclose(hist.sum(),weights.sum(),atol=1e-8):
        raise ValueError("Original pixels outside global display axes")
    bin_ids=np.clip(np.searchsorted(edges,rain,side="right")-1,0,len(edges)-2)
    obs=np.bincount(bin_ids,minlength=len(edges)-1)
    cells=pd.DataFrame({"c":frame._cell.to_numpy(int),"b":bin_ids})
    counts=cells.groupby("b").c.nunique()
    independent=np.zeros(len(obs),int)
    independent[counts.index.to_numpy(int)]=counts.to_numpy(int)
    return hist,obs,independent


def _panel(frame:pd.DataFrame,source:pd.DataFrame,edges:np.ndarray)->dict:
    weights=original_cell_weights(frame)
    raw,obs,cells=_raw_histogram(frame,edges,weights)
    mass=raw/raw.sum()
    by_climate=mass.sum(axis=1)
    conditional=np.divide(mass,by_climate[:,None],
                          out=np.zeros_like(mass),
                          where=by_climate[:,None]>0)
    valid=(obs>=MIN_DRAW_ROWS)&(cells>=MIN_DRAW_CELLS)
    # Missing conditional support is NULL, not a biological zero.
    conditional[~valid,:]=np.nan
    return {
        "status":"ORIGINAL_SAME_OBSERVATION_JOINT_IMAGE_ANGLE_AND_BIO12",
        "original_photo_observations":int(len(frame)),
        "original_geographic_cells":int(frame._cell.nunique()),
        "original_photographers":int(frame.user_id.nunique()),
        "original_wild_source_rain_bin_mass":
            _source_climate_probabilities(source,edges),
        "original_photo_rain_bin_mass":by_climate.tolist(),
        "original_photo_joint_mass":mass.tolist(),
        "photo_angle_conditional_on_rain_or_null":[
            [float(v) if np.isfinite(v) else None for v in row]
            for row in conditional],
        "photo_n_per_display_rain_bin":obs.tolist(),
        "photo_cells_per_display_rain_bin":cells.tolist(),
        "masked_conditionals_only_if": {
            "at_least_original_photos":MIN_DRAW_ROWS,
            "at_least_geographic_cells":MIN_DRAW_CELLS
        },
        "rain_and_angle_measured_on_same_photo_observation":True
    }


def _region_mix(source:pd.DataFrame,photo:pd.DataFrame,
                edges:np.ndarray)->dict:
    """Exact all-introduced mixture, even for straddling 5-degree cells."""
    photo_names=photo.tdwg1.map(CONTINENTS).fillna("other_unmapped").to_numpy(str)
    source_names=source.tdwg1.map(CONTINENTS).fillna("other_unmapped").to_numpy(str)
    global_weights=original_cell_weights(photo)
    all_raw=_raw_histogram(photo,edges,global_weights)[0]
    available=[]
    suppressed=[]
    for name in sorted(set(photo_names)):
        subset=photo.loc[photo_names==name]
        if len(subset)>=MIN_DRAW_ROWS and subset._cell.nunique()>=MIN_DRAW_CELLS:
            available.append((name,photo_names==name))
        else:
            suppressed.append(name)
    if suppressed:
        available.append(("other_low_support",np.isin(photo_names,suppressed)))
    recon=np.zeros_like(all_raw)
    groups={}
    for name,mask in available:
        data=photo.loc[mask]
        w=global_weights[mask]
        hist,obs,cells=_raw_histogram(data,edges,w)
        recon+=hist
        p=hist/hist.sum()
        climate_mask=(np.isin(source_names,suppressed)
                      if name=="other_low_support" else source_names==name)
        groups[name]={
            "status":"GROUPED_REGIONAL_DISPLAY" if name!="other_low_support"
                    else "POOLED_LOW_SUPPORT_REGIONS",
            "source_wild_n":int(climate_mask.sum()),
            "photo_n":int(len(data)),
            "photo_cells":int(data._cell.nunique()),
            "fraction_global_cell_weight":float(w.sum()/global_weights.sum()),
            "joint_mass_within_region":p.tolist(),
            "rain_mass_within_region":p.sum(axis=1).tolist(),
            "not_independent_invasion_events":True,
        }
    if not np.allclose(recon,all_raw,rtol=0,atol=1e-8):
        raise ValueError("Original globally weighted regional photo mixture mismatched")
    return {
        "status":"ALL_INTRODUCED_REGIONAL_MIXTURE_EXACTLY_RECONSTRUCTED",
        "groups":groups,
        "low_support_region_names_suppressed":True,
        "same_globally_cell_equal_weights_used_before_region_split":True,
        "all_introduced_joint_mass_recovered_exactly":True
    }


def evaluate(wild:pd.DataFrame,photo:pd.DataFrame,
             *,source_sha:dict|None=None)->dict:
    _check(wild,photo)
    source=_climate_source(wild)
    qualified=_measured(photo)
    results={}
    for species in SPECIES:
        w=source.loc[source.taxon_name.eq(species)]
        p=qualified.loc[qualified.taxon_name.eq(species)]
        grid=_rain_grid(w,p)
        native_w=w.loc[w.range_state.eq("native")]
        intro_w=w.loc[w.range_state.eq("introduced")]
        native=p.loc[p.range_state.eq("native")]
        intro=p.loc[p.range_state.eq("introduced")]
        if min(len(native),len(intro))<1000 or min(native._cell.nunique(),intro._cell.nunique())<15:
            raise ValueError("Original native/introduced photographic support insufficient")
        results[species]={
            "display_climate_BIO12_edges_mm_year":grid.tolist(),
            "display_angle_edges_image_degrees":ANGLE_EDGES.tolist(),
            "native":_panel(native,native_w,grid),
            "all_introduced":_panel(intro,intro_w,grid),
            "introduced_continents_supplementary":
                _region_mix(intro_w,intro,grid),
        }
    return {
        "schema":"azami_v4_S2_original_same_photo_rain_orientation_joint_v1",
        "status":"ORIGINAL_S2_JOINT_RAIN_AND_IMAGE_ANGLE_ALL_INTRODUCED",
        "two_species_original_wild_source":EXPECTED_SOURCE,
        "two_species_original_angle_eligible":EXPECTED_PHOTO,
        "balanced_photo_observations":int(len(qualified)),
        "primary_comparison":"native versus all WCVP explicitly introduced, for each species",
        "visualization_only_no_trait_width_or_tail_indices_selected":True,
        "global_rain_display_source_angle_blind":True,
        "source_raster_SHA256":source_sha or {},
        "privacy":{"obs_ids":False,"photographers":False,"precise_locations":False,
                   "photos":False,"individual_image_angles":False},
        "limitations":[
            "Joint photograph BIO12 and image angle link does not imply adaptive evolution",
            "BIO12 is 1981–2010 annual rainfall normal, not date-specific precipitation",
            "Original EXIF-corrected image-vertical reference is not true gravity",
            "Camera orientation, observer selection and flowering stage remain unknown",
            "Regional source categories do not represent independent invasion events",
            "The display grid and masking are technical support controls, not ecological response indices",
        ],
        "species":results,
    }


def save(wild:pd.DataFrame,photo:pd.DataFrame,output:Path,
         *,source_sha:dict|None=None)->dict:
    out=evaluate(wild,photo,source_sha=source_sha)
    output.mkdir(parents=True,exist_ok=True)
    (output/"S2_full_joint_rain_angle_landscape.json").write_text(
        json.dumps(out,indent=2,ensure_ascii=False,allow_nan=False)+"\n",
        encoding="utf-8")
    return out
