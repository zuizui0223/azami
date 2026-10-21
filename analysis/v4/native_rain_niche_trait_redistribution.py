"""S2: distinguish rainfall-niche redistribution from conditional head-photo phenotypes.

Source rain classes are determined from the original observation/observer
balanced WCVP-native roster before any photographic phenotype selection.
The original CHELSA pixel values also define area-weighted available climate.
A transport-standardisation identity contrasts rainfall composition and
within-rain-band image-angle distribution: F_intro - F_native =
(F_counterfactual - F_native) + (F_intro - F_counterfactual).
This is descriptive, NOT a causal invasion, plasticity, or evolution estimate.
All original observations, photographers, precise geography remain private.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.v4.native_cached_climate import VARIABLES
from analysis.v4.native_climate_trait_coupling import weighted_quantiles
from analysis.v4.native_introduced_climate_transfer import dated_and_balanced
from analysis.v4.native_pca_env import (
    SPECIES_KEYS, CONTINENTS, _native_and_introduced_codes
)

RAIN="chelsa_bio12"
ANGLE="orientation_angle_degrees_median"
OTHER=tuple(v for v in VARIABLES if v!=RAIN)
BANDS=("below_native_q25","native_q25_to_q50","native_q50_to_q75",
       "native_q75_to_q90","native_q90_to_q95","above_native_q95")
N_BINS=len(BANDS)
MIN_CELLS_PER_BAND=10
MIN_PHOTO_PER_BAND=30
MIN_SHARED_SOURCE_MASS=.65
MIN_VALID_BOOTSTRAPS=.80
BOOTSTRAPS=199
ANGLE_SUPPORT=np.arange(181,dtype=float)
RAIN_QUANTILES=(.25,.50,.75,.90,.95)
SEED=20261021

def _source_cells(frame):
    if (not len(frame) or frame.obs_id.duplicated().any()
        or frame._cell.isna().any() or frame.user_id.isna().any()):
        raise ValueError("Untrusted full-source individual/observer/cell keys")
    x=frame.copy()
    finite=np.isfinite(x[list(VARIABLES)].to_numpy(float)).all(axis=1)
    return x.loc[finite].copy()

def _cell_weight(frame):
    n=frame.groupby("_cell").obs_id.transform("size").to_numpy(float)
    if not len(n) or not np.isfinite(n).all() or (n<1).any():
        raise ValueError("Invalid observer-balanced spatial cells")
    return 1/n

def _rain_classes(frame,cutpoints):
    v=frame[RAIN].to_numpy(float)
    if not np.isfinite(v).all():
        raise ValueError("Missing source rain: do not assign artificial dry zero")
    return np.searchsorted(np.asarray(cutpoints,float),v,side="right").astype(int)

def _histogram(frame,cutpoints,*,angle=False):
    """Per independent 5-degree cell × native-defined rain stratum counts."""
    vals=_rain_classes(frame,cutpoints)
    loc,idx=np.unique(frame._cell.to_numpy(),return_inverse=True)
    w=_cell_weight(frame)
    counts=np.zeros((len(loc),N_BINS),float)
    raw=np.zeros((len(loc),N_BINS),int)
    np.add.at(counts,(idx,vals),w)
    np.add.at(raw,(idx,vals),1)
    result={"cells":loc,"band_cell_mass":counts,
            "band_cell_n":raw,
            "band_n":raw.sum(axis=0),
            "band_cells":np.array([np.unique(loc[idx[vals==k]]).size
                                   for k in range(N_BINS)],int),
            "rows":int(len(frame))}
    if angle:
        y=frame[ANGLE].to_numpy(float)
        if not np.isfinite(y).all() or ((y<0)|(y>180)).any():
            raise ValueError("Uncalibrated source image angle missing or invalid")
        # 1-degree resolution for CDF integration, not a replacement for the
        # continuous original measurement in the older within-species model.
        pos=np.clip(np.rint(y).astype(int),0,180)
        angle_mass=np.zeros((len(loc),N_BINS,181),float)
        np.add.at(angle_mass,(idx,vals,pos),w)
        result["angle_cell_band"]=angle_mass
    return result

def _background_profile(bg,cutpoints):
    if not len(bg) or not np.isfinite(bg[RAIN]).all():
        return None
    if "_area_weight" not in bg or not np.isfinite(bg._area_weight).all():
        raise ValueError("No original WCVP area-weighted background available")
    weights=bg._area_weight.to_numpy(float)
    b=np.bincount(_rain_classes(bg,cutpoints),weights=weights,minlength=N_BINS)
    if b.sum()<=0:
        return None
    return b/b.sum()

def _counts(h):
    b=h["band_cell_mass"].sum(axis=0)
    return b/b.sum()

def _area_and_source(source,available):
    occ=_counts(source)
    if available is None:
        return {"status":"NO_AVAILABLE_ORIGINAL_CLIMATE_BACKGROUND"}
    return {
        "status":"DESCRIPTIVE_SOURCE_VS_AVAILABLE_CLIMATE",
        "source_occurrence_probabilities":occ.tolist(),
        "available_area_probabilities":available.tolist(),
        "relative_occupancy_density_to_available_area":
            [float(a/b) if b>1e-10 else None
             for a,b in zip(occ,available)],
        "not_organism_independent_or_causal_environmental_selection":True,
    }

def _profile(photo):
    cell=photo["angle_cell_band"].sum(axis=0)
    weights=cell.sum(axis=1)
    conditional=np.divide(cell,weights[:,None],
                          out=np.zeros_like(cell),
                          where=weights[:,None]>0)
    return conditional

def _summary_prob(prob):
    p=np.asarray(prob,float)
    p=p/p.sum()
    c=np.cumsum(p)
    idx=[int(np.searchsorted(c,t,side="left"))
         for t in (.10,.50,.90)]
    return {
        "Q10_Q50_Q90_image_degrees":idx,
        "P_image_horizontal_or_down":float(p[90:].sum()),
        "P_image_at_least_120_degrees":float(p[120:].sum())
    }

def _composition(native_source,intro_source,native_photo,intro_photo,
                 *,common=None):
    pn=_counts(native_source)
    pi=_counts(intro_source)
    ng=_profile(native_photo)
    ig=_profile(intro_photo)
    if common is None:
        common=(np.minimum(native_photo["band_cells"],
                           intro_photo["band_cells"])>=MIN_CELLS_PER_BAND)
        common&=(native_photo["band_n"]>=MIN_PHOTO_PER_BAND)
        common&=(intro_photo["band_n"]>=MIN_PHOTO_PER_BAND)
    if common.sum()<2 or pn[common].sum()<MIN_SHARED_SOURCE_MASS or pi[common].sum()<MIN_SHARED_SOURCE_MASS:
        return {"status":"STOP_NO_RAIN_BAND_SHARED_PHOTOGRAPHIC_SUPPORT",
                "supported_bands":[BANDS[i] for i in np.flatnonzero(common)],
                "native_source_mass_supported":float(pn[common].sum()),
                "introduced_source_mass_supported":float(pi[common].sum())}
    pna=pn[common]/pn[common].sum()
    pia=pi[common]/pi[common].sum()
    fn=pna@ng[common] # What was observed in supported native rain environments
    fcf=pia@ng[common] # Composition-only counterfactual, native conditional phenotype
    fi=pia@ig[common] # Observed introduced within-band phenotypes
    if any(abs(x.sum()-1)>1e-8 for x in (fn,fcf,fi)):
        raise ValueError("Climate-standardised phenotype density not normalized")
    cf_delta=fcf.cumsum()-fn.cumsum()
    conditional_delta=fi.cumsum()-fcf.cumsum()
    observed_delta=fi.cumsum()-fn.cumsum()
    if not np.allclose(cf_delta+conditional_delta,observed_delta,atol=1e-10):
        raise ValueError("Conditional/occupancy mixture identity violated")
    p_down=lambda x:float(x[90:].sum())
    return {
        "status":"SUPPORTED_DESCRIPTIVE_RAIN_MIXTURE_DECOMPOSITION",
        "supported_bands":[BANDS[i] for i in np.flatnonzero(common)],
        "native_source_mass_supported":float(pn[common].sum()),
        "introduced_source_mass_supported":float(pi[common].sum()),
        "native_observed":_summary_prob(fn),
        "introduced_observed":_summary_prob(fi),
        "introduced_rain_distribution_with_native_angle_response":_summary_prob(fcf),
        "image_down_fraction_observed_introduced_minus_native":
            p_down(fi)-p_down(fn),
        "image_down_fraction_rain_composition_part":
            p_down(fcf)-p_down(fn),
        "image_down_fraction_within_rain_band_part":
            p_down(fi)-p_down(fcf),
        "cdf_difference_at_60_90_120":{
            str(t):{
                "observed":float(observed_delta[t]),
                "rain_occupancy_part":float(cf_delta[t]),
                "within_rain_band_part":float(conditional_delta[t])
            } for t in (60,90,120)
        },
        "angle_distribution_W1_native_vs_introduced_degrees":
            float(np.abs(observed_delta[:-1]).sum()),
        "one_degree_histogram_approximation":True,
        "not_causal_niche_shift_or_plasticity":True,
        "not_adjusted_for_four_other_climate_variables_or_camera":True,
    }

def _resampled_hist(hist,multipliers):
    m=np.asarray(multipliers,int)
    if len(m)!=len(hist["cells"]) or (m<0).any():
        raise ValueError("Same-region geographic block multipliers misaligned")
    out=dict(hist)
    out["band_cell_mass"]=hist["band_cell_mass"]*m[:,None]
    out["band_cell_n"]=hist["band_cell_n"]*m[:,None]
    out["band_n"]=out["band_cell_n"].sum(axis=0)
    if "angle_cell_band" in hist:
        out["angle_cell_band"]=hist["angle_cell_band"]*m[:,None,None]
        out["band_cells"]=np.count_nonzero(
            out["band_cell_n"],axis=0)
    return out

def _draw_source_photo_joint(source,photo,rng):
    """Draw the same source 5° blocks for the species × establishment region."""
    union=np.union1d(source["cells"],photo["cells"])
    multiplicities=np.bincount(rng.integers(len(union),size=len(union)),
                               minlength=len(union))
    srcmult=multiplicities[np.searchsorted(union,source["cells"])]
    phmult=multiplicities[np.searchsorted(union,photo["cells"])]
    return (_resampled_hist(source,srcmult),
            _resampled_hist(photo,phmult))

def _bootstrap(native_source,intro_source,native_photo,intro_photo,*,label):
    """Exploratory 5°-cell intervals; rain edges and eligible bins frozen."""
    seed=int.from_bytes(hashlib.sha256(label.encode()).digest()[:8],"little") ^ SEED
    rng=np.random.default_rng(seed)
    original=_composition(native_source,intro_source,native_photo,intro_photo)
    if original["status"]!="SUPPORTED_DESCRIPTIVE_RAIN_MIXTURE_DECOMPOSITION":
        return {"status":"STOP_NOT_IDENTIFIED_ORIGINAL_RAIN_BAND_COHORT"}
    # Preserve the pre-phenotype-analysed band support in the sensitivity
    # estimand. Resampling must not discover a more favorable set of bands.
    common=np.array([b in original["supported_bands"] for b in BANDS],bool)
    diffs=[]
    blocked=0
    for _ in range(BOOTSTRAPS):
        a,c=_draw_source_photo_joint(native_source,native_photo,rng)
        b,d=_draw_source_photo_joint(intro_source,intro_photo,rng)
        q=_composition(a,b,c,d,common=common)
        if q["status"]!="SUPPORTED_DESCRIPTIVE_RAIN_MIXTURE_DECOMPOSITION":
            blocked+=1
            continue
        diffs.append([
            q["image_down_fraction_observed_introduced_minus_native"],
            q["image_down_fraction_rain_composition_part"],
            q["image_down_fraction_within_rain_band_part"]
        ])
    if len(diffs)<BOOTSTRAPS*MIN_VALID_BOOTSTRAPS:
        return {"status":"NOT_ESTIMABLE_TOO_FEW_VALID_GEO_BOOTSTRAP_DRAWS",
                "valid":len(diffs),"requested":BOOTSTRAPS}
    v=np.array(diffs,float)
    return {
        "status":"EXPLORATORY_FIXED_BAND_REGION_CELL_RESAMPLING",
        "n_valid":len(diffs),"n_blocked":blocked,
        "image_down_probability_delta_intervals_2p5_50_97p5":{
            k:[float(z) for z in np.quantile(v[:,i],[.025,.5,.975])]
            for i,k in enumerate(("observed","rain_composition","within_band"))
        },
        "not_independent_plant_confidence_intervals":True,
        "native_rain_bin_edges_fixed_before_resampling":True,
        "geographic_resampling_joint_for_source_photo_same_region":True,
        "no_camera_genetic_or_climate_confounder_inference":True,
    }

def evaluate(wild:pd.DataFrame,land:pd.DataFrame,authority:dict,
             photo:pd.DataFrame,*,bootstrap:bool=True)->dict:
    """Conservative climate-availability × head-angle analysis on two species."""
    expected_source=144178
    expected_photo=173737
    if len(wild)!=expected_source or len(photo)!=expected_photo:
        raise ValueError("Full original 2-species source identity changed")
    if not set(VARIABLES)<=(set(wild)&set(land)&set(photo)):
        raise ValueError("Five original pixel CHELSA variables missing")
    if not {"taxon_name","range_state","tdwg1","obs_id","user_id","_cell"}<=set(wild):
        raise ValueError("Original source WCVP/observer scope changed")
    original_photo=photo.loc[np.isfinite(photo[ANGLE].to_numpy(float))&
        photo[ANGLE].between(0,180)].copy()
    balance=dated_and_balanced(original_photo)
    if len(balance)<30000:
        raise ValueError("No valid observer/site photo response sample")
    source=_source_cells(wild)
    output={}
    for species in SPECIES_KEYS:
        native=source.loc[source.taxon_name.eq(species)&source.range_state.eq("native")]
        introduced=source.loc[source.taxon_name.eq(species)&source.range_state.eq("introduced")]
        pw=balance.loc[balance.taxon_name.eq(species)]
        pn=pw.loc[pw.range_state.eq("native")]
        pi=pw.loc[pw.range_state.eq("introduced")]
        if min(len(native),len(introduced),len(pn),len(pi))<200:
            raise ValueError("Insufficient independently verified native and introduced cases")
        percentiles=weighted_quantiles(native[RAIN].to_numpy(float),
                                      _cell_weight(native),
                                      quantiles=RAIN_QUANTILES)
        if not all(a<b for a,b in zip(percentiles[:-1],percentiles[1:])):
            raise ValueError("Native source rainfall cutoffs degenerate")
        source_status=_native_and_introduced_codes(authority,SPECIES_KEYS[species])
        region=land.loc[land.TDWG_code.isin(source_status)].copy()
        region["_range"]=region.TDWG_code.map(source_status)
        nbg=region.loc[region._range.eq("native")]
        ibg=region.loc[region._range.eq("introduced")]
        if min(len(nbg),len(ibg))<200:
            raise ValueError("Original TDWG available background absent")
        nsh=_histogram(native,percentiles)
        nph=_histogram(pn,percentiles,angle=True)
        summary={
            "native_rain_source_cutoffs_q25_q50_q75_q90_q95_mm":percentiles,
            "rain_classes":list(BANDS),
            "native_source_n":len(native),
            "native_photo_n":len(pn),
            "native_source_rain":_area_and_source(nsh,_background_profile(nbg,percentiles)),
            "native_photo_n_by_rain_band":nph["band_n"].tolist(),
            "native_photo_cells_by_rain_band":nph["band_cells"].tolist(),
            "introduced":{}
        }
        for name,scope,back in [
            ("all_introduced",introduced,ibg),
            *[(CONTINENTS[int(c)],s,
               ibg.loc[ibg.TDWG_level1.eq(int(c))])
              for c,s in introduced.groupby("tdwg1",sort=True)]
        ]:
            iphotos=pi if name=="all_introduced" else pi.loc[pi.tdwg1.eq(
                next(k for k,v in CONTINENTS.items() if v==name)
            )]
            if not len(scope):
                continue
            ish=_histogram(scope,percentiles)
            entry={
                "source_n":len(scope),"photo_n":len(iphotos),
                "source_rain":_area_and_source(
                    ish,_background_profile(back,percentiles)),
                "source_n_by_band":ish["band_n"].tolist(),
                "source_cells_by_band":ish["band_cells"].tolist(),
            }
            if len(iphotos)>=200 and iphotos._cell.nunique()>=15:
                iph=_histogram(iphotos,percentiles,angle=True)
                entry["photo_n_by_band"]=iph["band_n"].tolist()
                entry["photo_cells_by_band"]=iph["band_cells"].tolist()
                entry["angle_distribution_rain_redistribution"]=_composition(
                    nsh,ish,nph,iph)
                if bootstrap and entry["angle_distribution_rain_redistribution"]["status"]=="SUPPORTED_DESCRIPTIVE_RAIN_MIXTURE_DECOMPOSITION":
                    entry["cell_resampling"]=_bootstrap(
                        nsh,ish,nph,iph,label=species+"|"+name)
            else:
                entry["angle_distribution_rain_redistribution"]={
                    "status":"STOP_INSUFFICIENT_INTRODUCED_PHOTO_CELLS"}
            summary["introduced"][name]=entry
        output[species]=summary
    return {
        "schema_version":"S2_native_intro_rain_niche_phenotype_decomposition_20261021_v1",
        "status":"EXECUTED_ORIGINAL_SOURCE_RAIN_OCCUPANCY_AND_IMAGE_TRAIT_MIXTURE",
        "original_source_two_species_n":297622,
        "observer_cell_balanced_source_n":len(wild),
        "original_angle_eligible_observations":len(photo),
        "date_and_original_observer_balanced_angle_n":len(balance),
        "background_mesh_point_n":len(land),
        "sources_are_angle_blind":True,
        "response":"Uncalibrated image-relative angle; not true gravity",
        "source_rainfall_1981_2010":"CHELSA v2.1 BIO12 climatic mean; not flowering-day rainfall",
        "species":output,
        "privacy":{"individual_ids_exported":False,"observer_ids_exported":False,
                   "coordinates_exported":False,"photos_exported":False},
        "limits":[
            "Native-introduced rainfall increase can follow a shift in available habitat alone",
            "Cell/photographer sampling does not fully remove observation effort selection",
            "Area-weighted WCVP availability is not true colonisable or surveyed area",
            "Composition versus conditional photo-angle decomposition is an arithmetic identity, not causal attribution",
            "Within rain-band differences retain four-climate/season/geography/camera confounding",
            "Rainfall bands outside overlap remain explicitly unsupported, never extrapolated",
            "Image horizontal-or-down fraction is not field-calibrated physical nodding",
            "Two species and multiple continents are not independent evolutionary experiments"
        ]
    }

def save(wild:pd.DataFrame,land:pd.DataFrame,authority:dict,photo:pd.DataFrame,
         output:Path)->dict:
    result=evaluate(wild,land,authority,photo)
    output.mkdir(parents=True,exist_ok=True)
    (output/"S2_rain_niche_phenotype_decomposition.json").write_text(
        json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n",
        encoding="utf-8")
    return result
