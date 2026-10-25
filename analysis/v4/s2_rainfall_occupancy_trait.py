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



def _dated_metadata(frame:pd.DataFrame)->pd.DataFrame:
    """Observation timing, NOT a measured flowering phenophase."""
    out=frame.copy()
    date=pd.to_datetime(out.observed_on,errors="coerce",format="mixed")
    out["_audit_month"]=date.dt.month.fillna(0).astype(int)
    out["_audit_grade"]=out.quality_grade.fillna("unknown").astype(str)
    out["_audit_key"]=out["_audit_month"].astype(str).str.zfill(2)+":"+out["_audit_grade"]
    return out


def _weighted_subset_fraction(frame:pd.DataFrame,mask:np.ndarray)->float:
    if frame.empty:
        return 0.0
    w=original_cell_weights(frame)
    return float(w[np.asarray(mask,bool)].sum()/w.sum())


def _source_eligibility_technical_profile(frame:pd.DataFrame,
                                           eligible_obs:set[str])->dict:
    """Only SOURCE-selected rows: do not call record grade camera image quality."""
    if frame.empty:
        return {"status":"NO_ORIGINAL_SOURCE_RECORDS"}
    d=_dated_metadata(frame)
    month=d["_audit_month"].to_numpy(int)
    nphoto=pd.to_numeric(d.source_photo_count,errors="coerce").to_numpy(float)
    if not np.isfinite(nphoto).all() or np.any(nphoto<1):
        raise ValueError("Authenticated original linked photograph count invalid")
    count=d.user_id.value_counts()
    top_n=max(1,int(np.ceil(.01*len(count))))
    top_observers=set(count.head(top_n).index)
    return {
        "status":"SOURCE_REPRESENTATIVE_OBSERVATION_TECHNICAL_PROFILE",
        "original_source_observations":int(len(d)),
        "original_5degree_cells":int(d._cell.nunique()),
        "original_reporters":int(d.user_id.nunique()),
        "same_selected_source_observation_has_eligible_archived_angle":
            _weighted_subset_fraction(d,d.obs_id.isin(eligible_obs).to_numpy()),
        "dated_observation_fraction":
            _weighted_subset_fraction(d,month>0),
        "summer_June_to_August_observation_fraction":
            _weighted_subset_fraction(d,(month>=6)&(month<=8)),
        "original_quality_grade_research_fraction_NOT_photo_pose":
            _weighted_subset_fraction(d,d["_audit_grade"].eq("research").to_numpy()),
        "original_multiple_linked_photographs_fraction_NOT_same_head":
            _weighted_subset_fraction(d,nphoto>=2),
        "top_1pct_original_reporters_fraction_of_source_observations":
            float(d.user_id.isin(top_observers).sum()/len(d)),
        "not_a_measured_camera_roll_pitch_or_gravity_angle":True,
    }


def _source_stratum_probs(frame:pd.DataFrame)->dict:
    d=_dated_metadata(frame)
    d=d.loc[d["_audit_month"].between(1,12)].copy()
    if d.empty:
        return {}
    w=original_cell_weights(d)
    labels=d["_audit_key"].to_numpy(str)
    result={}
    for k in np.unique(labels):
        result[str(k)]=float(w[labels==k].sum()/w.sum())
    return result


def _eligible_photo_strata(frame:pd.DataFrame)->dict:
    d=_dated_metadata(frame)
    d=d.loc[d["_audit_month"].between(1,12)].copy()
    result={}
    for name,g in d.groupby("_audit_key",sort=True):
        if (len(g)>=20 and g._cell.nunique()>=5 and
            g.user_id.nunique()>=10):
            result[str(name)]=g
    return result


def _normalised_stratum_composition(a:pd.DataFrame,b:pd.DataFrame,
                                     target:dict)->tuple[pd.DataFrame,pd.DataFrame]:
    def reweight(data):
        arr=[]
        for key,prob in target.items():
            sub=data.loc[data["_audit_key"].eq(key)].copy()
            w=original_cell_weights(sub)
            sub["_audit_weight"]=w/w.sum()*prob
            arr.append(sub)
        return pd.concat(arr,ignore_index=True)
    return reweight(a),reweight(b)


def _technical_standardised_angles(source_n:pd.DataFrame,source_i:pd.DataFrame,
                                   measured_n:pd.DataFrame,measured_i:pd.DataFrame,
                                   *,seed:int)->dict:
    """Same rainfall band; compare angle after shared month×record-grade composition.

    Common target is computed from original non-angle-selected occurrences.
    This does NOT adjust photographic camera orientation, biological flowering
    stage, spatial/country history, temperature or actual genetic provenance.
    """
    pn=_source_stratum_probs(source_n)
    pi=_source_stratum_probs(source_i)
    a=_eligible_photo_strata(measured_n)
    b=_eligible_photo_strata(measured_i)
    overlap={k:min(pn.get(k,0.),pi.get(k,0.)) for k in set(pn)&set(pi)}
    eligible={k:v for k,v in overlap.items()
              if k in a and k in b and v>0}
    total=float(sum(eligible.values()))
    meta={"common_source_month_grade_overlap_fraction":float(sum(overlap.values())),
          "admitted_common_source_month_grade_fraction":total,
          "supported_month_grade_strata":int(len(eligible)),
          "record_quality_grade_is_not_photo_pose":True,
          "month_is_observation_time_not_verified_anthesis":True}
    if total<.40 or not eligible:
        return dict(meta,status="NOT_IDENTIFIABLE_MONTH_GRADE_OVERLAP")
    target={k:eligible[k]/total for k in sorted(eligible)}
    src_a=_dated_metadata(measured_n)
    src_b=_dated_metadata(measured_i)
    x,y=_normalised_stratum_composition(src_a,src_b,target)
    if (min(len(x),len(y))<30 or min(x._cell.nunique(),y._cell.nunique())<15
        or min(x.user_id.nunique(),y.user_id.nunique())<10):
        return dict(meta,status="NOT_IDENTIFIABLE_MONTH_GRADE_SPATIAL_SUPPORT")
    def summary(frame,weights):
        angle=frame[ANGLE].to_numpy(float)
        if not np.isfinite(angle).all() or np.any((angle<0)|(angle>180)):
            raise ValueError("Archived photo-angle response invalid during QC")
        ww=np.asarray(weights,float)
        if ww.sum()<=0:
            raise ValueError("Technical standardisation lost observation weights")
        ww/=ww.sum()
        return (weighted_quantiles(angle,ww,(.1,.5,.9)),
                float(ww[angle>=90].sum()))
    wx=x["_audit_weight"].to_numpy(float)
    wy=y["_audit_weight"].to_numpy(float)
    qx,fx=summary(x,wx.copy())
    qy,fy=summary(y,wy.copy())
    qdiff=[float(i-j) for i,j in zip(qy,qx)]
    out=dict(meta,status="EXPLORATORY_COMMON_MONTH_GRADE_STANDARDISATION",
             source_common_target_fixed_without_photo_angle=True,
             standardised_native_n=int(len(x)),
             standardised_introduced_n=int(len(y)),
             standardised_native_cells=int(x._cell.nunique()),
             standardised_introduced_cells=int(y._cell.nunique()),
             standardised_native_effective_n=float(wx.sum()**2/np.square(wx).sum()),
             standardised_introduced_effective_n=float(wy.sum()**2/np.square(wy).sum()),
             standardised_native_Q10_Q50_Q90=qx,
             standardised_introduced_Q10_Q50_Q90=qy,
             introduced_minus_native_Q10_Q50_Q90=qdiff,
             introduced_minus_native_fraction_ge90=float(fy-fx))
    if min(out["standardised_native_effective_n"],
           out["standardised_introduced_effective_n"])<30:
        out["status"]="NOT_IDENTIFIABLE_LOW_EFFECTIVE_SAMPLE"
        return out
    # An explicit spatial sensitivity, with the prespecified common
    # source-month/record-grade mix kept constant in each resample.
    rng=np.random.default_rng(seed)
    prep=[]
    for frame in (x,y):
        cells,inverse=np.unique(frame._cell.to_numpy(),return_inverse=True)
        strata=frame["_audit_key"].astype(str).to_numpy()
        strata_ix={k:np.flatnonzero(strata==k) for k in target}
        prep.append((frame,cells,inverse,strata_ix,
                     frame["_audit_weight"].to_numpy(float)))
    estimates=[]
    for _ in range(199):
        outcome=[]
        for frame,cells,inv,by_stratum,base in prep:
            counts=np.bincount(rng.integers(len(cells),size=len(cells)),
                               minlength=len(cells))
            w=base*counts[inv]
            good=True
            for key,idx in by_stratum.items():
                z=float(w[idx].sum())
                if z<=0:
                    good=False
                    break
                w[idx]*=target[key]/z
            if not good:
                outcome=[]
                break
            outcome.append(summary(frame,w))
        if len(outcome)==2:
            qnative,pnative=outcome[0]
            qintro,pintro=outcome[1]
            estimates.append([*(float(a-b) for a,b in zip(qintro,qnative)),
                              float(pintro-pnative)])
    if len(estimates)<180:
        out["spatial_interval_status"]="INSUFFICIENT_VALID_MONTH_GRADE_SPATIAL_RESAMPLES"
        out["spatial_resamples_valid"]=len(estimates)
        return out
    values=np.asarray(estimates,float)
    out["spatial_interval_status"]="EXPLORATORY_FIXED_SEASON_GRADE_CELL_BOOTSTRAP"
    out["spatial_resamples_valid"]=len(estimates)
    out["introduced_minus_native_Q10_Q50_Q90_95pct_sensitivity"]=[
        [float(q) for q in np.quantile(values[:,i],(.025,.975))]
        for i in range(3)]
    out["introduced_minus_native_fraction_ge90_95pct_sensitivity"]=[
        float(q) for q in np.quantile(values[:,3],(.025,.975))]
    out["camera_pose_or_head_identity_measured"]=False
    return out


def _top_reporter_exclusion(source_n:pd.DataFrame,source_i:pd.DataFrame,
                            measured_n:pd.DataFrame,measured_i:pd.DataFrame)->dict:
    """Directional reporter-concentration sensitivity, not independent cameras."""
    result={}
    medians=[]
    for label,source,angles in (("native",source_n,measured_n),
                                ("introduced",source_i,measured_i)):
        names=source.user_id.value_counts()
        k=max(1,int(np.ceil(len(names)*.01)))
        exclude=set(names.head(k).index)
        subset=angles.loc[~angles.user_id.isin(exclude)].copy()
        n=len(subset)
        result[label]={"observers_excluded":k,
                       "observations_removed":int(len(angles)-n),
                       "retained_cells":int(subset._cell.nunique())}
        if (n<30 or subset._cell.nunique()<15
            or subset.user_id.nunique()<10):
            result["status"]="NOT_IDENTIFIABLE_AFTER_ORIGINAL_REPORTER_REMOVAL"
            return result
        w=original_cell_weights(subset)
        medians.append(weighted_quantiles(
            subset[ANGLE].to_numpy(float),w,(.5,))[0])
    result["status"]="EXPLORATORY_REMOVE_TOP_1PCT_ORIGINAL_REPORTERS"
    result["introduced_minus_native_median_image_degrees"]=float(medians[1]-medians[0])
    result["not_a_camera_pose_or_identity_validation"]=True
    return result


def rainfall_wet_band_technical_audit(
    source:pd.DataFrame,measured:pd.DataFrame,photo_ids:set[str],
    edges:tuple[float,float,float],*,seed:int
)->dict:
    """Only the already-used q75–q95 band and original photographic metadata."""
    required={"quality_grade","source_photo_count"}
    if not required<=set(source) or not required<=set(measured):
        return {"status":"MISSING_ORIGINAL_GRADE_OR_PHOTO_LINK_METADATA",
                "no_camera_pose_inferred":True}
    low,high=edges[1],edges[2]
    s=source.loc[source[RAIN].between(low,high)].copy()
    m=measured.loc[measured[RAIN].between(low,high)].copy()
    sn=s.loc[s.range_state.eq("native")]
    si=s.loc[s.range_state.eq("introduced")]
    mn=m.loc[m.range_state.eq("native")]
    mi=m.loc[m.range_state.eq("introduced")]
    if min(len(sn),len(si),len(mn),len(mi))<30:
        return {"status":"INSUFFICIENT_COMMON_WET_BAND_SOURCE"}
    return {
        "status":"ORIGINAL_SOURCE_WET_BAND_TECHNICAL_BIAS_AUDIT",
        "native_source_metadata":_source_eligibility_technical_profile(sn,photo_ids),
        "introduced_source_metadata":_source_eligibility_technical_profile(si,photo_ids),
        "same_wet_band_month_x_record_grade_standardised_angle":
            _technical_standardised_angles(sn,si,mn,mi,seed=seed),
        "same_wet_band_excluding_top_1pct_source_observers":
            _top_reporter_exclusion(sn,si,mn,mi),
        "verified_camera_roll_pitch_or_actual_gravity":False,
        "verified_flowering_stage":False,
        "individual_head_identifiers_linked":False,
        "original_linked_photo_count_is_not_same_head_repeat":True,
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
        result["wet_band_observation_bias_audit"]=rainfall_wet_band_technical_audit(
            src,pi,photo_ids,edges,seed=SEED+sum(map(ord,name)))
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
