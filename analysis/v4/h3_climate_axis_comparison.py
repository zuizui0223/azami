"""Source-verified H3c analogue comparison across three elevation-free climate definitions.

Original: BIO12/BIO01/RSDS.
Warm: BIO18/BIO05/RSDS.
Seasonal: BIO18/BIO01/RSDS/BIO04/BIO15.

Matched S1 observations, taxa, taxon-cells, 3D photographic phenotypes and RNG
are shared across scenarios. Climate-only quantile calibration equalises each
scenario's proportion of potential climate-analogue cell pairs *per draw*.
Only source-blind aggregate results leave the private GitHub Actions runner.

The H3c 2026-10-29 numerical result and original S1/S2 Holm families are NOT
recomputed or overwritten by this exploratory three-way sensitivity.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from analysis.v4.s1_climate_trait_hypervolume import (
    SOURCE_N, PHENOTYPE, MIN_CELLS, MIN_PHOTOS, RAREFIED_CELLS,
    _required_cols,
)
from analysis.v4.s1_h3_climate_analogue_alignment import (
    SEED, RESOLUTIONS, CLIMATE_ANALOG_Z, N_DRAWS, N_PERMUTATIONS,
    MIN_ANALOG_CELL_PAIRS, MIN_DISTINCT_CELLS_PER_SIDE,
    MIN_ELIGIBLE_TAXON_PAIRS, MIN_SUCCESSFUL_DRAWS,
    MIN_SECTOR_PAIR_TYPES, geographic_sector,
    observed_and_permutations,
)

CLIMATE_SETS = {
    "annual3_original": ("chelsa_bio12", "chelsa_bio01", "chelsa_rsds_mean"),
    "warm3_bio18_bio05": ("chelsa_bio18", "chelsa_bio05", "chelsa_rsds_mean"),
    "seasonal5_bio18_bio01": (
        "chelsa_bio18", "chelsa_bio01", "chelsa_rsds_mean",
        "chelsa_bio04", "chelsa_bio15",
    ),
}
BASE = "annual3_original"
CLIMATE_UNION = tuple(dict.fromkeys(k for group in CLIMATE_SETS.values() for k in group))
EXTRA = ("chelsa_bio18", "chelsa_bio05")
EXPECTED_UNION = (
    "chelsa_bio12", "chelsa_bio01", "chelsa_rsds_mean",
    "chelsa_bio18", "chelsa_bio05", "chelsa_bio04", "chelsa_bio15",
)
SOURCE_SHA = "342d78e9a8a960e8d9bdcf0bc2cdacd13bee48214a8b962419d5a71d5aef46ed"
EXPECTED_ORIGINAL_ELIGIBLE = 157822
EXPECTED_ORIGINAL_TAXA = 242
MIN_COMMON_CROSS_TAXON_PAIRS = 40
PRIVATE_FIELDS = {
    "obs_id", "user_id", "photo_id", "photographer_id", "accepted_key",
    "taxon_name", "latitude", "longitude", "lat", "lon", "source_record",
}
N_SCENARIOS = 3


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def load_common(native: Path, receipt: Path, extra: Path):
    source_receipt = json.loads(receipt.read_text(encoding="utf-8"))
    source_sha = sha256(native)
    if (source_sha != SOURCE_SHA or
        source_receipt.get("prepared_table_sha256") != source_sha or
        source_receipt.get("native_observations") != SOURCE_N or
        source_receipt.get("environment_decoding_status")
            != "physical_units_gsp_invalid_encoding_masked_v1"):
        raise ValueError("WRONG_ORIGINAL_S1_SOURCE_OR_RECEIPT")
    original = pd.read_csv(native, usecols=[*_required_cols(), "chelsa_bio04", "chelsa_bio15"], dtype={"obs_id": str},
                           low_memory=False)
    if (len(original) != SOURCE_N or original.obs_id.isna().any()
            or original.obs_id.duplicated().any()):
        raise ValueError("ORIGINAL_S1_SOURCE_GRAIN_DRIFT")
    base_axes = [*CLIMATE_SETS[BASE], *PHENOTYPE]
    eligible = np.isfinite(original[base_axes + ["latitude", "longitude"]].to_numpy(float)).all(axis=1)
    for trait in PHENOTYPE:
        eligible &= original[trait + "__eligible"].fillna(False).eq(True).to_numpy(bool)
    eligible &= original[PHENOTYPE[0]].between(0, 180).to_numpy(bool)
    original = original.loc[eligible].copy()
    if len(original) != EXPECTED_ORIGINAL_ELIGIBLE or original.accepted_key.nunique() != EXPECTED_ORIGINAL_TAXA:
        raise ValueError("ORIGINAL_H3C_QC_COHORT_NOT_REPRODUCED")

    expansion = pd.read_csv(extra, usecols=["obs_id", *EXTRA],
                            dtype={"obs_id": str}, low_memory=False)
    if expansion.obs_id.isna().any() or expansion.obs_id.duplicated().any():
        raise ValueError("INVALID_EXTRA_RASTER_JOIN_GRAIN")
    joined = original.merge(expansion, on="obs_id", how="left", validate="one_to_one")
    keep = np.isfinite(joined[list(CLIMATE_UNION)].to_numpy(float)).all(axis=1)
    shared = joined.loc[keep, ["accepted_key", "latitude", "longitude", *CLIMATE_UNION, *PHENOTYPE]].copy()
    if (len(shared) > len(original) or len(shared) < 100000 or
        shared.accepted_key.nunique() < 100 or
        shared[["latitude", "longitude"]].isna().any().any()):
        raise ValueError("NEW_CLIMATE_DIMENSIONS_INSUFFICIENT_SOURCE_SUPPORT")
    return shared, {
        "original_source_sha256": source_sha,
        "original_H3c_qualifying_photographs": int(len(original)),
        "original_H3c_source_taxa": EXPECTED_ORIGINAL_TAXA,
        "new_all_climate_variables_complete_photographs": int(len(shared)),
        "shared_photographs_lost_only_to_new_climate_layers": int(len(original)-len(shared)),
        "shared_taxa_before_8cell_40photo_filter": int(shared.accepted_key.nunique()),
        "original_source_observations": int(SOURCE_N),
        "no_elevation_in_any_matched_source_selection": True,
    }


def aggregate_cells(shared: pd.DataFrame, degrees: int):
    if degrees not in RESOLUTIONS:
        raise ValueError("UNREGISTERED_RESOLUTION")
    d = shared.copy()
    d["_latcell"] = np.floor((d.latitude + 90)/degrees).astype(int)
    d["_loncell"] = np.floor((d.longitude + 180)/degrees).astype(int)
    d = d.groupby(["accepted_key", "_latcell", "_loncell"], sort=True).agg(
        **{name: (name, "median") for name in [*CLIMATE_UNION, *PHENOTYPE]},
        n_photos=(PHENOTYPE[0], "size"),
        longitude=("longitude", "median"),
    ).reset_index()
    group = d.groupby("accepted_key", sort=True).agg(
        n_cells=("n_photos", "size"), n_photos=("n_photos", "sum")
    )
    eligible = group.index[(group.n_cells >= MIN_CELLS) & (group.n_photos >= MIN_PHOTOS)]
    d = d.loc[d.accepted_key.isin(eligible)].copy()
    if d.empty or len(eligible) < 25:
        raise ValueError("INSUFFICIENT_COMMON_CLIMATE_TAXON_CELL_SUPPORT")
    return d, {
        "degrees": degrees, "common_eligible_taxa": int(len(eligible)),
        "common_taxon_cell_values": int(len(d)),
        "common_source_photographic_records": int(d.n_photos.sum()),
        "minimum_taxon_cells": MIN_CELLS, "minimum_taxon_photos": MIN_PHOTOS,
        "same_taxon_cell_cohort_for_every_climate_scenario": True,
    }


def scaled_taxa(cells: pd.DataFrame):
    # All scenarios share each taxon's observations and the exact same global
    # photographic trait SD; climate SD is estimated outcome-blind per axis.
    names = [*CLIMATE_UNION, *PHENOTYPE]
    values = cells[names].to_numpy(float)
    counts = cells.groupby("accepted_key").size()
    weight = cells.accepted_key.map(1/counts).to_numpy(float)
    weight /= weight.sum()
    mu = np.average(values, axis=0, weights=weight)
    scale = np.sqrt(np.average((values-mu)**2, axis=0, weights=weight))
    if not np.isfinite(scale).all() or (scale < 1e-10).any():
        raise ValueError("CLIMATE_OR_PHENOTYPE_AXES_DEGENERATE")
    z = (values-mu)/scale
    for k, name in enumerate(names):
        cells["_common_z" + str(k)] = z[:, k]
    mapping = {name: "_common_z"+str(i) for i, name in enumerate(names)}
    taxa = []
    for taxon, part in cells.groupby("accepted_key", sort=True):
        sector = geographic_sector(part.longitude.to_numpy(float))
        assigned = int(np.argmax(np.bincount(sector, minlength=3)))
        taxa.append({
            "taxon": str(taxon),  # in-memory only; never exported
            "sector": assigned,
            "x": {key:part[[mapping[v] for v in group]].to_numpy(float)
                  for key, group in CLIMATE_SETS.items()},
            "y":part[[mapping[v] for v in PHENOTYPE]].to_numpy(float),
        })
    return taxa


def source_blind_pairs(x: np.ndarray, y: np.ndarray, sector: np.ndarray, threshold: float):
    n_taxa, cells, ndim = x.shape
    if (len(y) != n_taxa or y.shape != (n_taxa, RAREFIED_CELLS, 3)
        or sector.shape != (n_taxa,) or
        cells != RAREFIED_CELLS or ndim not in (3,5)):
        raise ValueError("INVALID_CLIMATE_AXIS_DIMENSIONS")
    rows=[]
    all_dist=[]
    for i in range(n_taxa):
        for j in range(i+1,n_taxa):
            if sector[i] == sector[j]:
                continue
            climate = np.linalg.norm(x[i, :, None, :] - x[j, None, :, :], axis=2)
            all_dist.append(climate.ravel())
            rows.append((i,j,climate))
    return rows, np.concatenate(all_dist)


def pack_pairs(rows, y: np.ndarray, limit: float):
    selected = {}
    for i,j,d_env in rows:
        mask = d_env <= limit
        if (int(mask.sum()) < MIN_ANALOG_CELL_PAIRS or
            np.flatnonzero(mask.any(axis=1)).size < MIN_DISTINCT_CELLS_PER_SIDE or
            np.flatnonzero(mask.any(axis=0)).size < MIN_DISTINCT_CELLS_PER_SIDE):
            continue
        distance = np.linalg.norm(y[i,:,None,:]-y[j,None,:,:],axis=2)
        selected[(i,j)] = (distance,mask,float(mask.sum()))
    return selected


def build_common_pair(pack: dict, shared_keys: list[tuple[int,int]]):
    return {
        "i": np.array([i for i,j in shared_keys], int),
        "j": np.array([j for i,j in shared_keys], int),
        "d": np.stack([pack[key][0] for key in shared_keys]),
        "mask": np.stack([pack[key][1] for key in shared_keys]),
        "analog_n": np.array([pack[key][2] for key in shared_keys], float),
    }


def analyse_paired(cells: pd.DataFrame, *, draws=N_DRAWS,
                   permutations=N_PERMUTATIONS, seed=SEED):
    taxa = scaled_taxa(cells)
    n_taxa = len(taxa)
    results = {k: [] for k in CLIMATE_SETS}
    permutations_by_scenario = {k: [] for k in CLIMATE_SETS}
    coverage = []
    for rep in range(draws):
        # Seed is independent of which climate model passes the match gate:
        # all models share the same rarefied cells within each source taxon.
        rng = np.random.default_rng(seed + rep*9109)
        selected = [rng.choice(len(t["y"]), RAREFIED_CELLS, replace=False) for t in taxa]
        y = np.stack([t["y"][idx] for t,idx in zip(taxa,selected)])
        sector = np.asarray([t["sector"] for t in taxa], int)
        x = {
            key: np.stack([t["x"][key][idx] for t,idx in zip(taxa,selected)])
            for key in CLIMATE_SETS
        }
        rows, distances = {}, {}
        for key in CLIMATE_SETS:
            rows[key],distances[key] = source_blind_pairs(x[key],y,sector,CLIMATE_ANALOG_Z)
        # Equalise the *overall fraction of potential climate-analogue cell pairs*
        # across 3D and 5D models. Calibration is climate-only and blind to y.
        shared_match_rate = float(np.mean(distances[BASE] <= CLIMATE_ANALOG_Z))
        thresholds = {
            BASE: CLIMATE_ANALOG_Z,
            **{key:float(np.quantile(distances[key], shared_match_rate, method="higher"))
               for key in CLIMATE_SETS if key!=BASE}
        }
        packs = {key:pack_pairs(rows[key],y,thresholds[key]) for key in CLIMATE_SETS}
        qualified = sorted(set.intersection(*(set(q) for q in packs.values())))
        sector_types = {tuple(sorted((int(sector[i]),int(sector[j])))) for i,j in qualified}
        coverage.append({
            "matched_rate_baseline": shared_match_rate,
            "matched_rate_quantiles": thresholds,
            "qualified_taxon_pairs_per_climate": {k:len(v) for k,v in packs.items()},
            "common_cross_taxon_pairs": len(qualified),
            "common_sector_pair_types": len(sector_types),
        })
        if len(qualified)<MIN_COMMON_CROSS_TAXON_PAIRS or len(sector_types)<MIN_SECTOR_PAIR_TYPES:
            continue
        for key in CLIMATE_SETS:
            pair = build_common_pair(packs[key], qualified)
            result,null = observed_and_permutations(
                pair, n_taxa, np.random.default_rng(seed + rep*500003 + 777),
                reps=permutations,
            )
            results[key].append(result)
            permutations_by_scenario[key].append(null)
    complete = len(results[BASE])
    if complete < min(MIN_SUCCESSFUL_DRAWS,draws):
        return {
            "status":"INSUFFICIENT_COMMON_CLIMATE_ANALOGUE_SUPPORT",
            "n_taxa":n_taxa, "requested_draws":draws, "common_draws":complete,
            "comparison_no_conclusions":True,
            "diagnostics":{
                "min_common_cross_taxon_pairs":int(min(x["common_cross_taxon_pairs"] for x in coverage)),
                "max_common_cross_taxon_pairs":int(max(x["common_cross_taxon_pairs"] for x in coverage)),
            }
        }
    scenario_outputs={}
    pvals=[]
    for key in CLIMATE_SETS:
        vals = np.asarray(results[key],float)
        nulls = np.stack(permutations_by_scenario[key])
        measured = float(vals.mean())
        perm_average = nulls.mean(axis=0)
        null_mean = float(perm_average.mean())
        p = float((1+np.sum(abs(perm_average-null_mean) >= abs(measured-null_mean)))/(permutations+1))
        scenario_outputs[key] = {
            "climate_axes":list(CLIMATE_SETS[key]),
            "n_climate_dimensions":len(CLIMATE_SETS[key]),
            "observed_matched_distance_advantage":measured,
            "null_mean":null_mean,
            "observed_minus_null":float(measured-null_mean),
            "p_exploratory_unadjusted":p,
            "draws_successful":complete,
            "rarefaction_sd":float(vals.std()),
            "cell_match_fraction_anchored_to_original_3D_climate":True,
            "shared_source_taxon_pairs_in_each_draw":True,
        }
        pvals.append(p)
    corrected=multipletests(pvals,method="holm")[1]
    for (name,output),pholm in zip(scenario_outputs.items(),corrected):
        output["p_holm3_across_posthoc_climate_definitions_exploratory"]=float(pholm)
    common_draws = [row for row in coverage if row["common_cross_taxon_pairs"]>=40 and row["common_sector_pair_types"]>=2]
    return {
        "status":"PAIRED_MATCHED_CLIMATE_AXIS_COMPARISON_EXPLORATORY",
        "n_taxa":n_taxa,
        "shared_climate_analogue_cell_pair_fraction_mean":float(np.mean([x["matched_rate_baseline"] for x in common_draws])),
        "successful_common_rarefaction_draws":complete,
        "required_minimum_successful_draws":min(MIN_SUCCESSFUL_DRAWS,draws),
        "mean_common_eligible_taxon_pairs":float(np.mean([x["common_cross_taxon_pairs"] for x in common_draws])),
        "min_common_eligible_taxon_pairs":int(min(x["common_cross_taxon_pairs"] for x in common_draws)),
        "source_matched_threshold_z":CLIMATE_ANALOG_Z,
        "distance_measurement_axes":list(PHENOTYPE),
        "models":scenario_outputs,
        "paired_effect_differences_observed_minus_null": {
            key+"_minus_"+BASE:float(scenario_outputs[key]["observed_minus_null"]-scenario_outputs[BASE]["observed_minus_null"])
            for key in CLIMATE_SETS if key!=BASE
        },
        "no_new_independent_replicate":True,
        "climate_dimension_matched_by_outcome_blind_match_rate":True,
        "causal_evolutionary_convergence_not_estimated":True,
    }


def assert_public(result):
    def walk(node):
        if isinstance(node,dict):
            for k,v in node.items():
                if k in PRIVATE_FIELDS:
                    raise ValueError("PRIVATE_SOURCE_KEY_IN_AGGREGATE_JSON")
                walk(v)
        elif isinstance(node,list):
            for item in node:
                walk(item)
    walk(result)
    if result["source"]["original_H3c_qualifying_photographs"] != EXPECTED_ORIGINAL_ELIGIBLE:
        raise ValueError("H3C_PHOTO_GATE_CHANGED")


def run(native: Path, receipt: Path, extra: Path, output: Path):
    shared,provenance = load_common(native,receipt,extra)
    panels = {}
    for degree in RESOLUTIONS:
        cells,cohort = aggregate_cells(shared,degree)
        fit=analyse_paired(cells,seed=SEED+degree)
        panels[f"{degree}degree"]={"source":cohort,"analysis":fit}
    result={
        "schema":"azami_h3c_three_vs_five_elevation_free_climate_axes_v1",
        "status":"SOURCE_VERIFIED_EXPLORATORY_POST_H3C_VARIANT_NOT_FROZEN_DISCOVERY",
        "source":provenance,
        "design":{
            "climate_sets":{k:list(v) for k,v in CLIMATE_SETS.items()},
            "all_models_exclude_elevation":True,
            "same_sample_and_geo_cells_for_each_scenario":True,
            "same_eight_cells_rarefied_per_taxon_per_draw":True,
            "same_joint_flower_traits_and_within_taxon_permutation_null":True,
            "analog_cell_pair_fraction_fixed_climate_blind_per_draw":True,
            "original_20261029_frozen_3D_H3c_not_changed":True,
            "S1_80_S2_128_original_Holm_family_not_changed":True,
        },
        "results":panels,
    }
    assert_public(result)
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    return result


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--native",required=True,type=Path)
    p.add_argument("--receipt",required=True,type=Path)
    p.add_argument("--extra",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    a = p.parse_args()
    result=run(a.native,a.receipt,a.extra,a.output)
    print(json.dumps({
        "status": result["status"],
        "source":result["source"],
        "results":{k:v["analysis"] for k,v in result["results"].items()},
    },ensure_ascii=False,allow_nan=False),flush=True)


if __name__=="__main__":
    main()
