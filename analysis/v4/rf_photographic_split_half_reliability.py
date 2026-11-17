"""Protected original-S1 split-half spatial-reliability audit for eight photo traits.

Two disjoint, outcome-blind halves of >=4 photographs in each taxon x 1/5deg
cell. Held-out half-B cell medians are predicted from half-A cell medians, and
vice versa, versus training-half *taxon* median. No original rows, IDs, raster,
individual cell scores or species names may leave the protected job.

Neither replicate split variation nor within-cell photo variation by itself
identifies camera error, genotype, local selection, or ecology.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.v4.native_integration import FEATURES


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open('rb') as source:
        for block in iter(lambda: source.read(8 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def _normalise(a: np.ndarray) -> np.ndarray:
    arr = np.asarray(a, dtype=float)
    length = np.linalg.norm(arr, axis=-1, keepdims=True)
    return np.divide(arr, length, out=np.full_like(arr, np.nan), where=length > 1e-8)


def _digest_distribution(values) -> dict:
    data = np.asarray(values, dtype=float)
    data = data[np.isfinite(data)]
    if data.size == 0:
        return {'n': 0, 'q05': None, 'q50': None, 'q95': None}
    a = np.quantile(data, [0.05, 0.5, 0.95])
    return {'n': int(data.size), 'q05': float(a[0]), 'q50': float(a[1]), 'q95': float(a[2])}

SOURCE_N = 313350
SALT = "azami-original-S1-photo-spatial-split-half-20261031-v1"
REPLICATES = 8
GRAINS = (1, 5)
GROUP = ["taxon", "lat_cell", "lon_cell"]
MIN_PHOTOS_CELL = 4
MIN_CELLS_TAXON = 8
MIN_PHOTOS_TAXON = 40
HUE_MIN_VECTOR_NORM = 0.20


def _source_mask(frame: pd.DataFrame, trait: str) -> np.ndarray:
    if trait not in FEATURES:
        raise ValueError("SPLIT_HALF_UNREGISTERED_CONSTRUCT")
    names = FEATURES[trait]
    required = {"obs_id", "accepted_key", "latitude", "longitude", *names,
                *(c + "__eligible" for c in names)}
    if not required.issubset(frame.columns):
        raise ValueError("SPLIT_HALF_REQUIRED_ORIGINAL_PHOTO_COLUMNS_MISSING")
    if frame.obs_id.isna().any() or not frame.obs_id.is_unique:
        raise ValueError("SPLIT_HALF_DUPLICATE_OR_MISSING_PHOTO_ID")
    lat = pd.to_numeric(frame.latitude, errors="coerce").to_numpy(float)
    lon = pd.to_numeric(frame.longitude, errors="coerce").to_numpy(float)
    y = frame[names].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    mask = (frame.accepted_key.notna().to_numpy(bool) & np.isfinite(y).all(axis=1)
            & np.isfinite(lat) & np.isfinite(lon)
            & (np.abs(lat) <= 90) & (np.abs(lon) <= 180))
    for c in names:
        mask &= frame[c + "__eligible"].eq(True).to_numpy(bool)
    if trait == "hue":
        mask &= np.linalg.norm(y, axis=1) >= HUE_MIN_VECTOR_NORM
    return mask


def _spatial_photo_cohort(frame: pd.DataFrame, trait: str, degrees: int,
                          *, min_photos_cell: int = MIN_PHOTOS_CELL,
                          min_cells_taxon: int = MIN_CELLS_TAXON,
                          min_photos_taxon: int = MIN_PHOTOS_TAXON) -> tuple[pd.DataFrame, dict]:
    if degrees not in GRAINS or min_photos_cell != 4 or min_cells_taxon != 8 or min_photos_taxon != 40:
        raise ValueError("SPLIT_HALF_UNREGISTERED_SUPPORT_SCOPE")
    mask = _source_mask(frame, trait)
    names = FEATURES[trait]
    good = frame.loc[mask, ["obs_id", "accepted_key", "latitude", "longitude", *names]].copy()
    good["taxon"] = good.accepted_key.astype(str)
    good["lat_cell"] = np.floor((pd.to_numeric(good.latitude).to_numpy(float) + 90) / degrees).astype(int)
    good["lon_cell"] = np.floor((pd.to_numeric(good.longitude).to_numpy(float) + 180) / degrees).astype(int)
    if trait == "hue":
        uv = _normalise(good[names].to_numpy(float))
        good.loc[:, names] = uv
    # Identical original photographs in both splits are never allowed.
    g = good.groupby(GROUP, sort=False, observed=True)
    count = g.size().rename("photos").reset_index()
    qualifying = count.loc[count.photos >= min_photos_cell, GROUP + ["photos"]]
    if trait == "hue":
        m = g[names].median().reset_index()
        qualifying = qualifying.merge(m, on=GROUP, how="left", validate="one_to_one")
        qualifying = qualifying.loc[np.linalg.norm(qualifying[names].to_numpy(float),axis=1)
                                    >= HUE_MIN_VECTOR_NORM, GROUP + ["photos"]]
    tax = qualifying.groupby("taxon").agg(cells=("photos", "size"), photos=("photos", "sum"))
    keep = tax.index[(tax.cells >= min_cells_taxon) & (tax.photos >= min_photos_taxon)]
    qualifying = qualifying.loc[qualifying.taxon.isin(keep)]
    selected = good.merge(qualifying[GROUP], on=GROUP, how="inner",validate="many_to_one")
    # Drop private coordinates + taxon source identifier from the return; only anonymised taxon label remains privately.
    selected = selected[["obs_id", *GROUP, *names]].copy()
    summary = {
        "n_source_eligible_photo_QC": int(mask.sum()),
        "n_qualified_taxa": int(len(keep)),
        "n_qualified_cells": int(len(qualifying)),
        "n_original_photos_in_supported_cells": int(len(selected)),
        "photos_per_supported_cell": _digest_distribution(qualifying.photos.to_numpy(int)),
        "cells_per_supported_taxon": _digest_distribution(tax.loc[keep, "cells"].to_numpy(int)),
        "min_source_photos_per_cell": 4,
        "min_cells_per_taxon": 8,
        "min_source_photos_per_taxon": 40,
        "source_rows_per_split_disjoint": True,
    }
    return selected, summary


def _hash_keys(obs_id: pd.Series, rep: int) -> np.ndarray:
    # Fixed opaque photo-ID-derived hash, not a sort by outcome, environment, or geography.
    h = (hashlib.blake2b(f"{SALT}|{rep}|{v}".encode(), digest_size=8).digest()
         for v in obs_id.astype(str))
    return np.fromiter((int.from_bytes(x, "big") for x in h), dtype=np.uint64, count=len(obs_id))


def _paired_half_cells(photo: pd.DataFrame, trait: str, *, rep: int) -> tuple[pd.DataFrame, dict]:
    names = FEATURES[trait]
    frame = photo.copy()
    frame["hash_rank"] = _hash_keys(frame.obs_id, rep)
    frame = frame.sort_values(GROUP + ["hash_rank", "obs_id"], kind="stable")
    frame["half"] = frame.groupby(GROUP, sort=False).cumcount().mod(2).astype("int8")
    # At least two unique original photos per half (full-cell gate >=4).
    group_count = frame.groupby(GROUP + ["half"], sort=False).size()
    if group_count.lt(2).any() or frame.obs_id.duplicated().any():
        raise ValueError("SPLIT_HALF_PHOTO_OVERLAP_OR_TOO_FEW_PHOTOS_PER_HALF")
    cells = frame.groupby(GROUP + ["half"], sort=False, observed=True)[names].median().unstack("half")
    flat = {}
    for side in (0, 1):
        for col in names:
            flat[f"{col}__half{side}"] = cells[(col, side)].to_numpy(float)
    pairs = cells.reset_index()[GROUP].copy()
    for key, values in flat.items():
        pairs[key] = values
    before = len(pairs)
    if trait == "hue":
        keep = np.ones(len(pairs), dtype=bool)
        for half in (0, 1):
            cols = [f"{c}__half{half}" for c in names]
            vec = pairs[cols].to_numpy(float)
            keep &= np.linalg.norm(vec,axis=1) >= HUE_MIN_VECTOR_NORM
            pairs.loc[:,cols] = _normalise(vec)
        pairs = pairs.loc[keep].copy()
    # Reapply support after hue-ambiguous cells; do NOT relax based on fitted skill.
    cnt = pairs.groupby("taxon").size()
    names_keep = cnt.index[cnt >= MIN_CELLS_TAXON]
    pairs = pairs.loc[pairs.taxon.isin(names_keep)].copy()
    return pairs, {"cells_before_half_hue_QC": int(before),
                   "cells_after_half_hue_QC": int(len(pairs)),
                   "taxa_after_half_hue_QC": int(pairs.taxon.nunique()),
                   "no_same_observation_in_both_halves": True}


def _centroid(values: np.ndarray, trait: str) -> np.ndarray:
    v = np.nanmedian(values, axis=0)
    if trait == "hue":
        v = _normalise(v[None, :])[0]
    return v


def _taxon_half_stats(cell_pairs: pd.DataFrame, trait: str, *, shuffle_seed: int) -> tuple[list[dict], float, float, float, float]:
    names = FEATURES[trait]
    rng=np.random.default_rng(shuffle_seed)
    diagnostic=[]
    total_loss=total_baseline=total_null_loss=0.
    for _, tax in cell_pairs.groupby("taxon",sort=False):
        a = tax[[f"{c}__half0" for c in names]].to_numpy(float)
        b = tax[[f"{c}__half1" for c in names]].to_numpy(float)
        n = len(a)
        if n < MIN_CELLS_TAXON:
            continue
        ca, cb = _centroid(a, trait), _centroid(b, trait)
        if not np.isfinite(np.r_[ca,cb]).all():
            continue
        # Symmetric cross-validation: each half's cell median is predicted ONLY
        # from other disjoint half; taxon baseline comes from training half.
        sse = float(np.square(a-b).sum() * 2)
        baseline = float(np.square(b-ca).sum() + np.square(a-cb).sum())
        if baseline <= 1e-12:
            continue
        order=rng.permutation(n)
        null=float(2*np.square(a-b[order]).sum())
        # Within-taxon centred cell-rank agreement (Pearson for scalars;
        # vector cosine correlation for circular 2D hue).
        ac=a-a.mean(axis=0)
        bc=b-b.mean(axis=0)
        den=float(np.sqrt(np.square(ac).sum()*np.square(bc).sum()))
        corr=float(np.sum(ac*bc)/den) if den>1e-12 else np.nan
        # Equal-cell loss PER taxon, avoiding domination by prolific taxa.
        total_loss += sse/n
        total_baseline += baseline/n
        total_null_loss += null/n
        diagnostic.append({"correlation":corr,
                           "cross_half_skill":1-sse/baseline,
                           "skill_minus_permuted_cell_label_null":(null-sse)/baseline,
                           "n_cells":n})
    if not diagnostic or total_baseline <= 1e-12:
        return [], float("nan"), float("nan"), float("nan"), float("nan")
    macro_skill = float(np.mean([x["cross_half_skill"] for x in diagnostic]))
    pooled_skill = float(1-total_loss/total_baseline)
    pooled_null_skill=float(1-total_null_loss/total_baseline)
    corr=np.array([x["correlation"] for x in diagnostic])
    median_corr=float(np.nanmedian(corr)) if np.isfinite(corr).any() else float("nan")
    return diagnostic, pooled_skill, pooled_null_skill, macro_skill, median_corr


def audit_split_half(frame: pd.DataFrame, trait: str, *, degrees: int,
                     reps: int=REPLICATES) -> dict:
    if trait not in FEATURES or degrees not in GRAINS or reps != REPLICATES:
        raise ValueError("SPLIT_HALF_ONLY_FIXED_EIGHT_REPLICATES_AND_GRAINS_ALLOWED")
    photo, support = _spatial_photo_cohort(frame, trait, degrees)
    out={"construct":trait,"degrees":int(degrees),"support":support,
         "status":"source_QC_only", "registered_seed_count":int(reps),
         "photo_halves_are_disjoint_but_not_biologically_independent":True,
         "trained_taxon_mean_only_no_test_outcome_leakage":True}
    if support["n_qualified_taxa"]<5:
        out["status"]="insufficient_repeated_taxa_for_split_half"
        return out
    metrics=[]
    all_correlations=[]
    all_skills=[]
    qualified=[]
    for rep in range(reps):
        pair,diag=_paired_half_cells(photo,trait,rep=rep)
        tax_stats, skill, null_skill, macro, corr = _taxon_half_stats(
            pair,trait,shuffle_seed=20261031+rep)
        if not tax_stats or not np.isfinite([skill,null_skill,macro,corr]).all():
            metrics.append({"status":"insufficient_cross_half_cells"})
            continue
        all_correlations.append(corr)
        all_skills.append(skill)
        qualified.append((len(tax_stats),len(pair)))
        metrics.append({"status":"estimated","n_taxa_eligible_this_rep":int(len(tax_stats)),
                        "n_taxon_cells_this_rep":int(len(pair)),
                        "equal_taxon_weighted_two_way_skill":float(skill),
                        "permuted_cell_alignment_skill":float(null_skill),
                        "gain_over_permuted_cell_alignment":float(skill-null_skill),
                        "taxon_mean_cross_half_skill":float(macro),
                        "median_within_taxon_cell_agreement_correlation":float(corr),
                        "fraction_taxa_positive_correlation":float(np.mean([
                            x["correlation"]>0 for x in tax_stats if np.isfinite(x["correlation"])]))})
    if not len(all_skills):
        out["status"]="insufficient_cross_half_spatial_signal_support"
        return out
    out.update({"status":"estimated_source_split_half_QC_not_RF",
                "n_successful_predeclared_splits":int(len(all_skills)),
                "n_taxa_per_split":_digest_distribution([x[0] for x in qualified]),
                "n_cells_per_split":_digest_distribution([x[1] for x in qualified]),
                "two_way_cross_half_skill_vs_train_taxon_median":_digest_distribution(all_skills),
                "within_taxon_cross_half_cell_agreement_correlation":_digest_distribution(all_correlations),
                "n_splits_with_positive_spatial_cell_skill":int(np.sum(np.asarray(all_skills)>0)),
                "n_splits_with_positive_median_taxon_correlation":int(np.sum(np.asarray(all_correlations)>0)),
                "gain_over_shuffled_cell_correspondence":_digest_distribution([
                    v["gain_over_permuted_cell_alignment"] for v in metrics if v["status"]=="estimated"]),
                "random_split_quantiles_are_not_CIs_or_p_values":True,
                "randomization_null_is_not_spatially_autocorrelation_preserving":True,
                "individual_taxon_cell_photo_rows_and_split_ids_not_exported":True})
    return out


def run(native:Path,receipt_path:Path,output:Path,*,expected_n:int=SOURCE_N,
        constructs:tuple[str,...]|None=None, grains:tuple[int,...]=GRAINS) -> dict:
    receipt=json.loads(Path(receipt_path).read_text(encoding="utf-8"))
    if (sha256(native)!=receipt.get("prepared_table_sha256") or
        receipt.get("native_observations")!=expected_n or
        receipt.get("environment_decoding_status")!=
           "physical_units_gsp_invalid_encoding_masked_v1"):
        raise ValueError("S1_ORIGINAL_SPLIT_HALF_SOURCE_NOT_AUTHENTICATED")
    chosen=tuple(FEATURES) if constructs is None else tuple(constructs)
    if (not chosen or len(chosen)!=len(set(chosen)) or not set(chosen)<=set(FEATURES)
        or not grains or len(grains)!=len(set(grains)) or not set(grains)<=set(GRAINS)):
        raise ValueError("SPLIT_HALF_UNREGISTERED_CONSTRUCT_OR_GRAIN_SCOPE")
    names=[c for x in chosen for c in FEATURES[x]]
    cols=["obs_id","accepted_key","latitude","longitude",*names,
          *(x+"__eligible" for x in names)]
    raw=pd.read_csv(native,usecols=list(dict.fromkeys(cols)),dtype={"obs_id":str},low_memory=False)
    if len(raw)!=expected_n or raw.obs_id.isna().any() or not raw.obs_id.is_unique:
        raise ValueError("S1_SPLIT_HALF_UNEXPECTED_ORIGINAL_PHOTO_GRAIN")
    results=[audit_split_half(raw,c,degrees=d) for c in chosen for d in grains]
    report={"schema":"azami_original_s1_photo_split_half_geographic_signal_v1",
            "status":"authenticated_original_S1_exploratory_split_half_photographic_QC",
            "source_original_s1_observations":int(expected_n),
            "source_original_decoded_sha256":sha256(native),
            "design":{"constructs":list(chosen),"degrees":list(grains),
                      "disjoint_photo_halves":True,"seeds":REPLICATES,
                      "min_photos_per_cell":4,"min_cells_taxon":8,"min_photos_taxon":40},
            "not_the_same_RF28_or_RF34_complete_case_population":True,
            "not_independent_cameras_individuals_observers_or_floral_development":True,
            "not_a_direct_measurement_error_or_adaptive_ecological_test":True,
            "new_P_values_or_selection_by_favourable_random_split":False,
            "original_GEB_S1_S2_Holm_H3c_and_RF_untouched":True,
            "private_taxon_photo_coordinate_observer_or_row_predictions_exported":False,
            "results":results}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("native","receipt","output"):
        p.add_argument("--"+key,type=Path,required=True)
    args=p.parse_args()
    report=run(args.native,args.receipt,args.output)
    print(json.dumps({"s1_original_n":report["source_original_s1_observations"],
                      "n_diagnostics":len(report["results"]),
                      "n_estimated":sum(x["status"]=="estimated_source_split_half_QC_not_RF"
                                        for x in report["results"])}),flush=True)

if __name__=="__main__":main()