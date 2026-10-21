"""S2: How does capitulum PHOTO orientation vary across the wetter invaded climate?

Simple, evidence-first two-panel summary: (1) climate-occupancy redistribution
in original WCVP-native / explicitly introduced wild Cirsium records independent
of trait; (2) original orientation distribution in the SAME native-defined
rainfall bands. No biomean inference, rare q95-only selection, or causal claim.

One original observation / photographer x 5-degree cell is maintained. The
default introduced focus is North America (both source-supported species);
other continents and all introduced ranges are descriptive sensitivities.
Returned outputs contain aggregated source and photo counts only.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.v4.native_cached_climate import VARIABLES
from analysis.v4.native_climate_trait_coupling import (
    original_cell_weights, weighted_quantiles,
)
from analysis.v4.native_introduced_climate_transfer import dated_and_balanced
from analysis.v4.native_intro_continent_orientation import CONTINENTS

SPECIES = ("Cirsium arvense", "Cirsium vulgare")
RAIN = "chelsa_bio12"
ANGLE = "orientation_angle_degrees_median"
Q = (0.50, 0.75, 0.95)
BANDS = (
    "below_native_median",
    "native_q50_to_q75",
    "native_q75_to_q95",
    "above_native_q95",
)
MIN_ANGLE_N = 30
MIN_ANGLE_CELLS = 5
MIN_ANGLE_OBSERVERS = 10
MIN_COMPARE_CELLS = 15
BOOTSTRAPS = 249
SEED = 20261021


def _verify(table: pd.DataFrame, *, photo: bool) -> None:
    required = {"obs_id","taxon_name","range_state","user_id","_cell",
                "tdwg1","observed_on",RAIN,*VARIABLES}
    if photo:
        required.add(ANGLE)
    if not required <= set(table):
        raise ValueError("Original wild/climate/observation provenance incomplete")
    if table.obs_id.duplicated().any() or table.user_id.isna().any():
        raise ValueError("Original observation/photographer identity not unique")
    if not set(table.taxon_name) <= set(SPECIES):
        raise ValueError("Unexpected species in original source")
    if not set(table.range_state) <= {"native","introduced"}:
        raise ValueError("Unverified WCVP range category")


def _finite(table: pd.DataFrame) -> pd.DataFrame:
    keep = np.isfinite(table[list(VARIABLES)].to_numpy(float)).all(axis=1)
    return table.loc[keep].copy()


def _band(rain: np.ndarray, edges: tuple[float,float,float]) -> np.ndarray:
    """Native-defined rain bands; no photo angles enter bin selection."""
    return np.asarray(BANDS,object)[np.searchsorted(
        np.asarray(edges,float),rain,side="right"
    )]


def _source_rain_bands(source_native: pd.DataFrame) -> tuple[float,float,float]:
    if len(source_native)<200 or source_native._cell.nunique()<20:
        raise ValueError("Original native climate source too sparse")
    q=weighted_quantiles(
        source_native[RAIN].to_numpy(float),
        original_cell_weights(source_native),
        quantiles=Q,
    )
    if not (q[0]<q[1]<q[2]):
        raise ValueError("Native rain bands are degenerate")
    return tuple(float(v) for v in q)


def _source_summary(source: pd.DataFrame, edges: tuple[float,float,float],
                    eligible_ids: set[str]) -> dict:
    if not len(source):
        return {"status":"NO_SOURCE_OBSERVATIONS"}
    cats=_band(source[RAIN].to_numpy(float),edges)
    base_weights=original_cell_weights(source)
    output={
        "n_original_source_observations":int(len(source)),
        "n_5degree_cells":int(source._cell.nunique()),
        "n_original_observers":int(source.user_id.nunique()),
        "BIO12_median_mm_cell_balanced":float(weighted_quantiles(
            source[RAIN].to_numpy(float),base_weights,quantiles=(.5,))[0]),
        "bands":{},
    }
    for name in BANDS:
        group=source.loc[cats==name]
        w=base_weights[cats==name]
        output["bands"][name]={
            "source_n":int(len(group)),
            "source_5degree_cells":int(group._cell.nunique()),
            "source_observers":int(group.user_id.nunique()),
            "source_fraction_observations":float(len(group)/len(source)),
            "source_fraction_cell_equal_weight":float(w.sum()/base_weights.sum()),
            "with_archived_angle_obs_id":int(group.obs_id.isin(eligible_ids).sum()),
        }
    return output


def _image_profile(frame: pd.DataFrame) -> dict:
    if (len(frame)<MIN_ANGLE_N or frame._cell.nunique()<MIN_ANGLE_CELLS
        or frame.user_id.nunique()<MIN_ANGLE_OBSERVERS):
        return {
            "status":"INSUFFICIENT_IMAGE_ANGLE_OBSERVATIONS",
            "angle_n":int(len(frame)),
            "independent_cells":int(frame._cell.nunique()),
        }
    y=frame[ANGLE].to_numpy(float)
    if not np.isfinite(y).all() or ((y<0)|(y>180)).any():
        raise ValueError("Invalid EXIF image vertical orientation")
    w=original_cell_weights(frame)
    w/=w.sum()
    q=weighted_quantiles(y,w,quantiles=(.1,.5,.9))
    return {
        "status":"DESCRIPTIVE_CAMERA_RELATIVE_ANGLE",
        "angle_n":int(len(frame)),
        "independent_cells":int(frame._cell.nunique()),
        "original_observers":int(frame.user_id.nunique()),
        "image_angle_Q10_Q50_Q90_deg":q,
        "fraction_at_or_below_image_horizontal_90_deg":
            float(w[y<=90].sum()),
        "fraction_at_or_above_image_horizontal_90_deg":
            float(w[y>=90].sum()),
    }


def _profiles(frame: pd.DataFrame, edges:tuple[float,float,float]) -> dict:
    if frame.empty:
        return {name:{"status":"NO_PHOTO_SAMPLES"} for name in BANDS}
    label=_band(frame[RAIN].to_numpy(float),edges)
    return {name:_image_profile(frame.loc[label==name]) for name in BANDS}


def _cell_bootstrap_contrast(a:pd.DataFrame,b:pd.DataFrame,
                             *,seed:int) -> dict:
    """Simple within-band introduced-minus-native median/fraction, cell CI.

    Independent 5-degree cells are resampled separately across ranges.
    Not independent colonisation events; no camera/selection-error inference.
    """
    if min(a._cell.nunique(),b._cell.nunique())<MIN_COMPARE_CELLS:
        return {"status":"NOT_IDENTIFIED_INDEPENDENT_CELLS"}
    aa=a[ANGLE].to_numpy(float)
    bb=b[ANGLE].to_numpy(float)
    wa=original_cell_weights(a)
    wb=original_cell_weights(b)
    ca,ai=np.unique(a._cell.to_numpy(),return_inverse=True)
    cb,bi=np.unique(b._cell.to_numpy(),return_inverse=True)
    rng=np.random.default_rng(seed)
    vals=[]
    for _ in range(BOOTSTRAPS):
        wa2=wa*np.bincount(rng.integers(len(ca),size=len(ca)),
                           minlength=len(ca))[ai]
        wb2=wb*np.bincount(rng.integers(len(cb),size=len(cb)),
                           minlength=len(cb))[bi]
        if wa2.sum()==0 or wb2.sum()==0:
            continue
        za=wa2>0
        zb=wb2>0
        w1=wa2[za]/wa2[za].sum()
        w0=wb2[zb]/wb2[zb].sum()
        median_delta=(
            weighted_quantiles(aa[za],w1,(.5,))[0]
            -weighted_quantiles(bb[zb],w0,(.5,))[0]
        )
        down_delta=float(w1[aa[za]>=90].sum()-w0[bb[zb]>=90].sum())
        vals.append((median_delta,down_delta))
    if len(vals)<.90*BOOTSTRAPS:
        return {"status":"NOT_IDENTIFIED_BOOTSTRAP_DRAWS"}
    v=np.asarray(vals,float)
    return {
        "status":"EXPLORATORY_CELL_BOOTSTRAP_NOT_CAMERA_ADJUSTED",
        "original_cells_introduced":int(len(ca)),
        "original_cells_native":int(len(cb)),
        "replicates":len(v),
        "median_image_angle_introduced_minus_native_95pct":
            [float(x) for x in np.quantile(v[:,0],[.025,.5,.975])],
        "fraction_image_angle_ge90_introduced_minus_native_95pct":
            [float(x) for x in np.quantile(v[:,1],[.025,.5,.975])],
        "not_independent_introductions_or_adaptation":True,
    }


def assess(source: pd.DataFrame, photo: pd.DataFrame,
           *,climate_shas:dict|None=None) -> dict:
    _verify(source,photo=False)
    _verify(photo,photo=True)
    if len(source)<1000 or len(photo)<1000:
        raise ValueError("Original 2-species study source missing")
    wild=_finite(source)
    photographed=_finite(photo)
    photographed=photographed.loc[
        photographed[ANGLE].between(0,180) &
        np.isfinite(photographed[ANGLE].to_numpy(float))
    ].copy()
    photographed=dated_and_balanced(photographed)
    photo_ids=set(photo.obs_id.astype(str))
    focal={}
    for name in SPECIES:
        src=wild.loc[wild.taxon_name.eq(name)]
        ref=src.loc[src.range_state.eq("native")]
        intro=src.loc[src.range_state.eq("introduced")]
        edges=_source_rain_bands(ref)
        pi=photographed.loc[photographed.taxon_name.eq(name)]
        pn=pi.loc[pi.range_state.eq("native")]
        introduced=pi.loc[pi.range_state.eq("introduced")]
        result={
            "native_rain_reference_BIO12_q50_q75_q95_mm":list(edges),
            "climate_occupancy":{
                "native":_source_summary(ref,edges,photo_ids),
                "introduced_all":_source_summary(intro,edges,photo_ids),
            },
            "orientation_by_native_rain_band":{
                "native":_profiles(pn,edges),
                "introduced_all":_profiles(introduced,edges),
            },
            "introduced_minus_native_median_and_down_fraction_by_band":{},
            "introduced_by_continent":{},
        }
        for k,band in enumerate(BANDS):
            aa=introduced.loc[_band(introduced[RAIN].to_numpy(float),edges)==band]
            bb=pn.loc[_band(pn[RAIN].to_numpy(float),edges)==band]
            if (_image_profile(aa)["status"]=="DESCRIPTIVE_CAMERA_RELATIVE_ANGLE"
                and _image_profile(bb)["status"]=="DESCRIPTIVE_CAMERA_RELATIVE_ANGLE"):
                result["introduced_minus_native_median_and_down_fraction_by_band"][band]=(
                    _cell_bootstrap_contrast(aa,bb,seed=SEED+k+sum(map(ord,name)))
                )
            else:
                result["introduced_minus_native_median_and_down_fraction_by_band"][band]={
                    "status":"NOT_IDENTIFIED_ANGLE_SUPPORT"}
        for tdwg,srcpart in intro.groupby("tdwg1",sort=True):
            if pd.isna(tdwg) or int(tdwg) not in CONTINENTS:
                raise ValueError("Unrecognised source TDWG continent")
            continent=CONTINENTS[int(tdwg)]
            observed=introduced.loc[introduced.tdwg1.eq(tdwg)]
            result["introduced_by_continent"][continent]={
                "climate_occupancy":_source_summary(srcpart,edges,photo_ids),
                "orientation_by_native_rain_band":_profiles(observed,edges)
            }
        focal[name]=result
    return {
        "schema":"azami_S2_rain_occupancy_and_image_orientation_by_native_rain_band_v1",
        "status":"SOURCE_VERIFIED_RAIN_OCCUPANCY_AND_ANGLE_BANDS_EXPLORATORY",
        "total_original_wild_two_species_source":297622,
        "total_original_source_cell_balanced":int(len(source)),
        "total_original_photo_eligible":173737,
        "photo_balanced_climate_complete":int(len(photographed)),
        "rainfall_bands":"Native wild (not photo-selected) original WCVP source, q50/q75/q95, cell balanced; same edges for both native and introduced",
        "trait_response":"EXIF image-vertical orientation, 0 up 90 horizontal 180 down, not gravity",
        "first_scientific_question":"Have native and introduced populations occupied different rainfall bands?",
        "second_scientific_question":"Across those same rainfall bands, is the distribution of original photographic head angles different?",
        "focal":focal,
        "CHELSA_original_TIFF_SHA256":climate_shas or {},
        "privacy":{"individual_observation_ids":False,"original_photographer_ids":False,
                   "coordinates":False,"source_photos":False},
        "limits":[
            "Source records represent sampled geographic occupancy, not plant population frequencies or true realized range",
            "Photographic eligibility differs across rainfall bands and status",
            "Climatic normals are not rainfall on flowering days",
            "Differences by rainfall band are observational and do not establish plasticity, genetic adaptation or rain shelter",
            "Only 5degree-cell conditional intervals; photographic pose and genotype uncalibrated",
        ]
    }


def save(source:pd.DataFrame,photo:pd.DataFrame,
         output:Path,climate_shas:dict|None=None) -> dict:
    result=assess(source,photo,climate_shas=climate_shas)
    output.mkdir(parents=True,exist_ok=True)
    (output/"S2_rain_occupancy_vs_head_angle_bands.json").write_text(
        json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n",
        encoding="utf-8"
    )
    return result
