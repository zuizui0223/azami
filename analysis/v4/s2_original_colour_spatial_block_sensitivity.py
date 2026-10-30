"""S2 original frozen photographic-colour gradients: spatial *block* sensitivity.

Do not add latitude or longitude as regression covariates.  Use source
coordinates only to predefine 20-degree spatial resampling blocks, retaining
the original CHELSA BIO18/BIO05/RSDS/elevation, date nuisance and observer
× 5-degree-cell reduction, inverse-5-degree-cell WLS weighting.

For each of 3 original Holm-passing introduced colour associations and their
native counterparts (6 panels), require exact regeneration of the original
coefficient before any spatial result is labelled source-verified.  Then run
a 20-degree spatial block bootstrap (399 draws) and five source-blind
highest-contribution-block holdouts. Neither bootstrap CI nor block jackknife
is part of the original S2 128-Holm family, nor a test of evolved differences.

Original individuals, observations, photographer IDs, and locations remain in
the authenticated GitHub Actions runner. Export aggregate estimates only.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm

from analysis.v4.s1_s2_extreme_lego import (
    ALL_ENV, S2_RENAMING
)
from analysis.v4.s1_s2_extreme_range_models import _get_dated_balanced, _formula
from analysis.v4.s1_s2_vif5_si10_alltraits import (
    FOUR_WARM, standardized_y
)

SOURCE_N = 297622
RANDOM_SEED = 20261029
BOOTSTRAPS = 399
DEGREES_BLOCK = 20
MAX_LARGEST_HOLDOUTS = 5
EXPECTED = (
 ("Cirsium arvense","native","corolla_lab_lightness",
  "topo_elevation",22936,0.05626673194073065),
 ("Cirsium arvense","introduced","corolla_lab_lightness",
  "topo_elevation",26171,-0.17373187423223405),
 ("Cirsium vulgare","native","corolla_lab_chroma",
  "chelsa_rsds_mean",17335,0.1368829067685122),
 ("Cirsium vulgare","introduced","corolla_lab_chroma",
  "chelsa_rsds_mean",27200,0.1345037874025455),
 ("Cirsium vulgare","native","corolla_lab_lightness",
  "chelsa_rsds_mean",17335,-0.04880796672586844),
 ("Cirsium vulgare","introduced","corolla_lab_lightness",
  "chelsa_rsds_mean",27200,-0.14280073287481068),
)


def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for blob in iter(lambda:f.read(1<<20),b""):
            h.update(blob)
    return h.hexdigest()


def source_cohort(focal_path: Path, extra_path: Path, measured_path: Path):
    focal=pd.read_csv(focal_path,dtype={"obs_id":str},low_memory=False)
    extra=pd.read_csv(extra_path,dtype={"obs_id":str},low_memory=False)
    measurements=pd.read_csv(measured_path,dtype={"obs_id":str},low_memory=False)
    if (len(focal)!=SOURCE_N or not focal.obs_id.is_unique or
        not extra.obs_id.is_unique or
        not measurements.obs_id.is_unique or
        set(ALL_ENV[:5])-set(focal) or
        set(ALL_ENV[5:])-set(extra) or
        set(["taxon_name","range_state","latitude","longitude","user_id",
             "observed_on"])-set(focal)):
        raise ValueError("S2_AUTHENTIC_SOURCE_JOIN_OR_COVARIATE_FAILURE")
    if set(focal.range_state)!={"native","introduced"} or set(focal.taxon_name)!={
        "Cirsium arvense","Cirsium vulgare"
    }:
        raise ValueError("ORIGINAL_S2_TWO_SPECIES_FOUR_STRATA_CHANGED")
    # Some original S2 measured files contain original canonical names, others
    # archived median aliases. Use the exact original S2 alias-preference order.
    used={}
    for canonical in ("corolla_lab_lightness","corolla_lab_chroma"):
        available=[name for name in (canonical,*S2_RENAMING[canonical])
                   if name in measurements]
        if not available:
            raise ValueError("ORIGINAL_S2_COLOUR_MEASUREMENT_NOT_FOUND")
        used[canonical]=available[0]
    measured=measurements[["obs_id",*used.values()]].rename(
        columns={v:k for k,v in used.items()})
    sub=extra[["obs_id",*[e for e in ALL_ENV if e not in focal]]]
    combined=focal.merge(sub,on="obs_id",how="left",validate="one_to_one")
    combined=combined.merge(measured,on="obs_id",how="left",validate="one_to_one")
    if len(combined)!=SOURCE_N:
        raise ValueError("S2_IDENTITY_CHANGED_DURING_MERGE")
    return combined,{
        "original_focal_wild_records":len(focal),
        "original_focal_source_sha256":sha(focal_path),
        "original_new_climate_source_sha256":sha(extra_path),
        "original_colour_measurements_sha256":sha(measured_path),
        "measurement_source_columns":used,
        "no_original_ids_or_locations_exported":True
    }


def original_sample(focal: pd.DataFrame, species: str, state: str,
                    trait: str):
    selection=(focal.taxon_name==species)&(focal.range_state==state)
    stratum=focal.loc[selection].copy()
    if stratum.empty:
        raise ValueError("ORIGINAL_S2_STRATUM_MISSING")
    subset=np.isfinite(stratum[list(ALL_ENV)].to_numpy(float)).all(axis=1)
    part=stratum.loc[subset].copy()
    part=part.loc[np.isfinite(part[trait])].copy()
    balanced=_get_dated_balanced(part)
    if (balanced._cell.nunique()<15 or
        balanced.user_id.nunique()<10 or len(balanced)<100):
        raise ValueError("UNEXPECTED_S2_ELIGIBLE_COLOUR_COHORT_SUPPORT")
    return balanced,{
        "wild_source_observations":len(stratum),
        "all12_complete_source":int(subset.sum()),
        "measured_colour_all12_complete":len(part),
        "balanced_eligible_photographs":len(balanced),
        "balanced_5degree_cells":int(balanced._cell.nunique()),
        "balanced_source_photographers":int(balanced.user_id.nunique()),
    }


def spatial_block_labels(frame: pd.DataFrame, degrees=DEGREES_BLOCK):
    """Only a resampling cluster label, NOT an explanatory variable."""
    lat,lon=frame.latitude.to_numpy(float),frame.longitude.to_numpy(float)
    if (not np.isfinite(lat).all() or not np.isfinite(lon).all() or
        np.any(np.abs(lat)>90) or np.any(np.abs(lon)>180)):
        raise ValueError("INVALID_COORDS_FOR_SPATIAL_BOOTSTRAP")
    if degrees!=20:
        raise ValueError("UNREGISTERED_BLOCK_SCALE")
    lb=np.floor((lat+90)/degrees).astype(int)
    ob=np.floor((lon+180)/degrees).astype(int)
    return lb*int(np.ceil(360/degrees))+ob


def _original_fit(frame: pd.DataFrame, trait: str, exposure: str):
    X,_,_,month,omitted=_formula(frame,list(FOUR_WARM))
    y=standardized_y(frame[trait].to_numpy(float))
    cells=frame._cell.to_numpy()
    counts=frame.groupby("_cell").size()
    weights=frame._cell.map(1/counts).to_numpy(float)
    fit=sm.WLS(y,X,weights=weights).fit(
        cov_type="cluster",cov_kwds={"groups":cells},use_t=False)
    index=list(FOUR_WARM).index(exposure)+1
    result={"beta":float(fit.params[index]),"cluster_robust_p_uncorrected":float(fit.pvalues[index]),
            "pointwise_95CI":[float(v) for v in np.asarray(fit.conf_int())[index]],
            "source_date_nuisance_cols":month,
            "source_date_nuisance_omitted":omitted,
            "n_original_5degree_clusters":int(len(counts))}
    return X,y,weights,index,result


def _sufficient(X, y, wt, labels):
    codes,unique=pd.factorize(labels,sort=True)
    nblock=len(unique);p=X.shape[1]
    gram=np.zeros((nblock,p,p))
    cross=np.zeros((nblock,p))
    contrib=np.zeros(nblock)
    # Each ORIGINAL 5deg cell has equal weight in the baseline; the spatial
    # block bootstrap resamples entire 20deg blocks of preweighted photo rows.
    for k in range(nblock):
        ix=(codes==k)
        x=X[ix];w=wt[ix];z=y[ix]
        gram[k]=x.T@(w[:,None]*x)
        cross[k]=x.T@(w*z)
        contrib[k]=w.sum()
    return gram,cross,contrib


def bootstrap_blocks(X,y,weights,index,block,*,reps=BOOTSTRAPS,seed=RANDOM_SEED):
    gram,cross,contrib=_sufficient(X,y,weights,block)
    p=gram.shape[1]
    if len(gram)<8:
        return {"status":"insufficient_spatial_blocks",
                "n_blocks":len(gram),"min_required":8}
    full=np.linalg.solve(gram.sum(axis=0),cross.sum(axis=0))[index]
    rng=np.random.default_rng(seed)
    coefficient=[]
    for _ in range(reps):
        indices=rng.integers(0,len(gram),size=len(gram))
        mult=np.bincount(indices,minlength=len(gram))
        G=np.einsum("b,bij->ij",mult,gram)
        h=np.einsum("b,bj->j",mult,cross)
        if np.linalg.matrix_rank(G)==p and np.linalg.cond(G)<1e10:
            coefficient.append(float(np.linalg.solve(G,h)[index]))
    nvalid=len(coefficient)
    if nvalid<max(50,int(reps*.8)):
        return {"status":"insufficient_bootstrap_identifiability",
                "n_blocks":len(gram),"n_valid_bootstraps":nvalid}
    coefs=np.asarray(coefficient)
    ci=[float(x) for x in np.quantile(coefs,[.025,.975])]
    # Source-blind influence: drop 5 blocks with greatest ORIGINAL total
    # aggregate cell weights, not largest fitted effects or smallest P.
    order=np.argsort(-contrib,kind="stable")[:min(MAX_LARGEST_HOLDOUTS,len(gram))]
    holds=[]
    for j in order:
        G=gram.sum(axis=0)-gram[j]
        h=cross.sum(axis=0)-cross[j]
        if np.linalg.matrix_rank(G)<p:
            holds.append({"status":"not_estimable"})
        else:
            coef=float(np.linalg.solve(G,h)[index])
            holds.append({"status":"estimated","coefficient":coef,
                          "source_weight_fraction_omitted":float(contrib[j]/contrib.sum()),
                          "sign_retained":bool(np.sign(coef)==np.sign(full))})
    return {
        "status":"estimated",
        "n_original_source_geo_blocks":len(gram),
        "grid_cell_resolution_degrees":20,
        "block_partition_used_for_resampling_not_a_regression_predictor":True,
        "n_valid_bootstraps":nvalid,
        "n_requested_bootstraps":reps,
        "source_block_bootstrap_95_percentile_interval":ci,
        "fraction_bootstraps_same_sign_as_original":float(np.mean(np.sign(coefs)==np.sign(full))),
        "block_bootstrap_center":float(np.mean(coefs)),
        "n_largest_source_weight_blocks_heldout":len(order),
        "largest_source_weight_block_holdout_results":holds,
        "n_largest_holdouts_retaining_sign":sum(h.get("sign_retained",False) for h in holds),
        "not_original_S2_Holm_family_or_causal_spatial_partition":True
    }


def evaluate(focal,*,bootstrap_reps=BOOTSTRAPS):
    original=[];summaries={}
    for i,(species,state,trait,exposure,n,expected) in enumerate(EXPECTED):
        cohort,counts=original_sample(focal,species,state,trait)
        if len(cohort)!=n:
            raise ValueError(
                f"S2_ORIGINAL_SOURCE_COHORT_DRIFT_IN_{i}: {len(cohort)} != {n}"
            )
        X,y,w,idx,fit=_original_fit(cohort,trait,exposure)
        if not np.isclose(fit["beta"],expected,rtol=0,atol=2e-7):
            raise ValueError(
                f"S2_ORIGINAL_HOLM_POSITIVE_REFERENCE_COEFFICIENT_DRIFT_IN_{i}: "
                f"{fit['beta']} != {expected}"
            )
        blocks=spatial_block_labels(cohort)
        sens=bootstrap_blocks(X,y,w,idx,blocks,
                              reps=bootstrap_reps,seed=RANDOM_SEED+i)
        original.append({
            "species":species,"range_status":state,"photographic_trait":trait,
            "environment":exposure,"n":len(cohort),
            "source_verification":{"exact_frozen_beta":expected,
                                   "recovered_beta":fit["beta"],
                                   "source_5degree_cluster_p_uncorrected":fit["cluster_robust_p_uncorrected"],
                                   "same_original_data_and_nuisance_design":True},
            "source_support":counts,
            "spatial_block_sensitivity":sens,
        })
    return {"status":"SOURCE_EXACT_S2_COLOUR_SPATIAL_BLOCK_BOOTSTRAP",
            "registered_panels":len(EXPECTED),
            "original_wild_source_observations":SOURCE_N,
            "original_holm_family_not_refitted":True,
            "spatial_degrees_only_used_to_define_resampling_blocks":True,
            "model_covariates":list(FOUR_WARM),
            "original_sample_reweighted_equal_5degree_cells":True,
            "panel_rows":original,
            "no_photo_observer_taxon_or_precise_geographic_ids_in_output":True}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--focal",type=Path,required=True)
    ap.add_argument("--extra",type=Path,required=True)
    ap.add_argument("--measured",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()
    source,provenance=source_cohort(a.focal,a.extra,a.measured)
    receipt=evaluate(source)
    receipt["source_integrity"]=provenance
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(receipt,indent=2,ensure_ascii=False,allow_nan=False)+"\n",
                     encoding="utf-8")
    print(json.dumps({"status":receipt["status"],
                      "source":provenance,
                      "panels":[{
                          "species":x["species"],"range_status":x["range_status"],
                          "trait":x["photographic_trait"],
                          "exposure":x["environment"],
                          "n":x["n"],
                          "beta":x["source_verification"]["recovered_beta"],
                          "spatial":x["spatial_block_sensitivity"],
                      } for x in receipt["panel_rows"]]},ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
