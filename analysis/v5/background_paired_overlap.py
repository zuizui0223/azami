"""Availability-corrected regional climate D paired with photographed head D.

Exploratory V5 branch follow-up: same original WCVP wild/photo records and
20-degree regional resampling as analysis.v5.overlap_coupling, but use the
species- and region-specific available climatic background to correct the
climatic occupancy density. Never export source points or individual images.

The original published 2026-11-09 region-pair rho is an uncorrected *occurrence*
climate-overlap statistic, NOT this new available-environment-corrected one.
A model fit here is observational and is not evidence of evolution/adaptation.
"""
from __future__ import annotations

import io
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
from scipy.spatial.distance import cdist
from statsmodels.stats.multitest import multipletests

from analysis.v4.native_introduced_opportunity import (
    ORIGINAL_AAD, ORIGINAL_SEALED_SHA, ORIGINAL_SHA,
    MEASUREMENT_AAD, MEASUREMENT_SEALED_SHA, MEASUREMENT_SHA,
    angular_cells, authority_source, decrypt_checked, source_photos,
    source_roster,
)
from analysis.v4.native_cached_climate import VARIABLES
from analysis.v5.native_pca_env import (
    SOURCE_COUNTS, SPECIES_KEYS, _native_and_introduced_codes,
    background_mesh,
)
from analysis.v5.native_pca_env_fullcoverage import apply_original_five_CHELSA
from analysis.v5.overlap_coupling import TRAITS, SEED, project, grid, density, infer

ROOT = Path("work/submission")
PRIVATE = ROOT / "private"
OUT = ROOT / "public"
MIN_BG_POINTS = 32
MIN_BG_AREA_COVERAGE = 0.80
MIN_SUPPORTED_NATIVE_REGIONS = 3
MIN_SUPPORTED_INTRO_REGIONS = 3
RNG_DRAWS = 32
PER_REGION_DRAW = 8
PCA_BINS = 64
BACKGROUND_FLOOR = 0.0001
MIN_OCCUPANCY_BG_SUPPORT = 0.80
EXPECTED_PHOTOS = {"Cirsium arvense": 44780, "Cirsium vulgare": 39092}
EXPECTED_REGIONS = {"Cirsium arvense": (17, 8), "Cirsium vulgare": (11, 7)}


def block_key(lat, lon):
    """Identical 20-degree integer partitions to the original paired test."""
    latitude = np.asarray(lat, float)
    longitude = np.asarray(lon, float)
    i = (latitude + 90).clip(max=179.99999).astype(int)
    j = ((longitude + 180) % 360).astype(int)
    return (i // 20) * 18 + (j // 20)


def weighted_project_from_background(bg):
    """One species-wide common PCA, 1/2 area weighting per range (not photo PCA)."""
    if set(bg.range_state) != {"native", "introduced"}:
        raise ValueError("BACKGROUND_BOTH_RANGES_REQUIRED")
    arrays = bg[list(VARIABLES)].to_numpy(float)
    area = bg._area_weight.to_numpy(float)
    w = area.copy()
    ranges = bg.range_state.to_numpy()
    for state in ("native", "introduced"):
        mask = ranges == state
        if len(arrays[mask]) < MIN_BG_POINTS or (area[mask] <= 0).any():
            raise ValueError("INSUFFICIENT_SPECIES_BACKGROUND")
        w[mask] /= 2 * area[mask].sum()
    w /= w.sum()
    center = np.average(arrays, axis=0, weights=w)
    scale = np.sqrt(np.average((arrays - center) ** 2, axis=0, weights=w))
    if (scale <= 1e-10).any():
        raise ValueError("BACKGROUND_ZERO_VARIANCE")
    z = (arrays - center) / scale
    eig, vec = np.linalg.eigh((z * w[:, None]).T @ z)
    order = np.argsort(eig)[::-1]
    eig = eig[order]
    components = vec[:, order[:2]]
    if eig[1] <= 1e-12:
        raise ValueError("BACKGROUND_PCA_TWO_AXIS_UNAVAILABLE")
    scorescale = np.sqrt(eig[:2])

    def transform(values):
        v = np.asarray(values, float)
        if v.ndim != 2 or v.shape[1] != len(VARIABLES):
            raise ValueError("BACKGROUND_PCA_INVALID_CLIMATE")
        return ((v - center) / scale) @ components / scorescale

    return transform, float(eig[:2].sum() / eig.sum())


def regional_density(points, edges, h, weights=None):
    """Normalised KDE on a globally fixed PCA grid, with no point disclosure."""
    pts = np.asarray(points, float)
    if pts.ndim != 2 or pts.shape[1] != 2 or not len(pts):
        raise ValueError("EMPTY_REGIONAL_DENSITY")
    if not np.isfinite(pts).all():
        raise ValueError("NONFINITE_REGIONAL_PCA_COORDINATE")
    if weights is not None:
        weights = np.asarray(weights, float)
        if len(weights) != len(pts) or (weights <= 0).any():
            raise ValueError("INVALID_REGIONAL_WEIGHTS")
    hist = np.histogram2d(pts[:, 0], pts[:, 1], bins=edges, weights=weights)[0]
    if hist.sum() <= 0:
        raise ValueError("PCA_POINTS_OUTSIDE_GRID")
    sigma = [h / np.diff(e)[0] for e in edges]
    z = gaussian_filter(hist, sigma=sigma, mode="constant")
    if z.sum() <= 0:
        raise ValueError("DEGENERATE_KDE")
    return z / z.sum()


def correct_occupancy(occupied, available):
    """PCA-env ratio with explicit unsupported-mass gate."""
    occupied = np.asarray(occupied, float)
    available = np.asarray(available, float)
    if occupied.shape != available.shape or not np.isclose(occupied.sum(), 1, atol=1e-7) or not np.isclose(available.sum(), 1, atol=1e-7):
        raise ValueError("OCCUPIED_AND_BACKGROUND_DENSITY_NOT_NORMALISED")
    floor = max(BACKGROUND_FLOOR * available.max(), 1e-12)
    supported = available >= floor
    fraction = float(occupied[supported].sum())
    if fraction < MIN_OCCUPANCY_BG_SUPPORT:
        raise ValueError("REGIONAL_OCCUPANCY_OUTSIDE_AVAILABLE_CLIMATE")
    corrected = np.zeros_like(occupied)
    corrected[supported] = occupied[supported] / np.maximum(available[supported], floor)
    if corrected.sum() <= 0:
        raise ValueError("NO_CORRECTED_REGIONAL_DENSITY")
    return corrected / corrected.sum(), fraction


def overlap_matrix(a, b):
    """Numerically stable Schoener D for two collections of density vectors."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if a.ndim != 2 or b.ndim != 2 or a.shape[1] != b.shape[1]:
        raise ValueError("OVERLAP_DENSITY_SHAPE_MISMATCH")
    if not np.allclose(a.sum(1), 1, atol=1e-7) or not np.allclose(b.sum(1), 1, atol=1e-7):
        raise ValueError("UNNORMALISED_OVERLAP_DENSITY")
    return np.clip(1 - cdist(a, b, metric="cityblock") / 2, 0, 1)


def _region_background(bg, photos):
    """Gate each region's WCVP-available land and original raster completeness.

    The background mask is a species-specific WCVP-native or explicitly
    introduced TDWG3 subset intersected with the fixed original 20deg block.
    Exclude unsupported regions, NEVER replace missing climates with zero.
    """
    full = bg.copy()
    full["block"] = block_key(full.latitude, full.longitude)
    quality = {}
    for (state, block), cell in photos.groupby(["range_state", "block"], sort=True):
        if len(cell) < PER_REGION_DRAW:
            continue
        region = full[(full.range_state == state) & (full.block == block)]
        complete = region[np.isfinite(region[list(VARIABLES)].to_numpy(float)).all(axis=1)]
        tot = float(region._area_weight.sum())
        area = float(complete._area_weight.sum())
        if tot <= 0:
            continue
        fraction = area / tot
        if len(complete) < MIN_BG_POINTS or fraction < MIN_BG_AREA_COVERAGE:
            continue
        quality[(str(state), int(block))] = {
            "available": complete,
            "background_fraction": fraction,
        }
    return quality


def _one_species(c, land, authority, species):
    """Same matched cells and region random draws for all climate/head versions."""
    states = _native_and_introduced_codes(authority, SPECIES_KEYS[species])
    species_land = land[land.TDWG_code.isin(states)].copy()
    species_land["range_state"] = species_land.TDWG_code.map(states)
    if set(species_land.range_state) != {"native", "introduced"}:
        raise ValueError("NO_WCVP_BOTH_RANGES")
    complete = species_land[np.isfinite(species_land[list(VARIABLES)].to_numpy(float)).all(axis=1)].copy()
    transform, variance = weighted_project_from_background(complete)
    xbg = transform(complete[list(VARIABLES)].to_numpy(float))
    xp = transform(c[list(VARIABLES)].to_numpy(float))
    # Fixed common axes and grid for every region and for the original-photo
    # versus availability-corrected comparison.
    allx = np.vstack([xp, xbg])
    edgebg = [np.linspace(allx[:, k].min() - .5, allx[:, k].max() + .5, PCA_BINS + 1)
              for k in range(2)]
    xold, _ = project(c, list(VARIABLES))
    y, _ = project(c, TRAITS)
    y2, _ = project(c, TRAITS[:2])
    edges = {"head3": grid(y), "head2": grid(y2), "original_photo_climate": grid(xold)}

    cc = c.copy()
    cc["_idx"] = np.arange(len(cc), dtype=int)
    regional = _region_background(species_land, cc)
    keys = {}
    for state in ("native", "introduced"):
        keys[state] = sorted((state, int(block)) for block in cc.loc[cc.range_state == state, "block"].unique()
                             if (state, int(block)) in regional)
    if min(len(keys["native"]), len(keys["introduced"])) < 3:
        return {"status": "STOP_REGIONAL_BACKGROUND_SUPPORT_INSUFFICIENT",
                "n_native_supported": len(keys["native"]),
                "n_intro_supported": len(keys["introduced"]),
                "n_photo_cells": int(len(c))}, []
    bg_density = {}
    all_support = {}
    for key in keys["native"] + keys["introduced"]:
        available = regional[key]["available"]
        arr = transform(available[list(VARIABLES)].to_numpy(float))
        bg_density[key] = regional_density(arr, edgebg, .5,
                                           weights=available._area_weight.to_numpy(float))
        idx = cc.loc[(cc.range_state == key[0]) & (cc.block == key[1]), "_idx"].to_numpy(int)
        density_full = regional_density(xp[idx], edgebg, .5)
        _, frac = correct_occupancy(density_full, bg_density[key])
        all_support[key] = frac

    # No geographic identifiers, original observations, or photo pixel values
    # are included in result outputs.
    output = []
    for h in (.5, .35, .75):
        # Area-weighted WCVP background KDE is independent of photo draw.
        # Compute it once per region and bandwidth, not 32x per draw.
        bg_h = {}
        for key in keys["native"] + keys["introduced"]:
            available = regional[key]["available"]
            bg_h[key] = regional_density(
                transform(available[list(VARIABLES)].to_numpy(float)),
                edgebg, h, weights=available._area_weight.to_numpy(float))
        for head, positions in (("three_trait", y), ("angle_chroma", y2)):
            rain = np.zeros((len(keys["native"]), len(keys["introduced"])))
            raw = np.zeros_like(rain)
            phen = np.zeros_like(rain)
            support_min = 1.
            rng = np.random.default_rng(SEED)  # identical draws, including sensitivities
            for _ in range(RNG_DRAWS):
                picked = []
                for side in ("native", "introduced"):
                    for state, block in keys[side]:
                        z = cc.loc[(cc.range_state == state) & (cc.block == block), "_idx"].to_numpy(int)
                        picked.append((state, block, rng.choice(z, PER_REGION_DRAW, replace=False)))
                bg_parts, original_parts, head_parts = [], [], []
                for state, block, sample in picked:
                    occ = regional_density(xp[sample], edgebg, h)
                    corrected, coverage = correct_occupancy(
                        occ, bg_h[(state, block)])
                    support_min = min(support_min, coverage)
                    bg_parts.append(corrected.ravel())
                    original_parts.append(regional_density(xold[sample], edges["original_photo_climate"], h).ravel())
                    head_parts.append(regional_density(positions[sample], edges["head3" if head == "three_trait" else "head2"], h).ravel())
                n = len(keys["native"])
                rain += overlap_matrix(bg_parts[:n], bg_parts[n:]) / RNG_DRAWS
                raw += overlap_matrix(original_parts[:n], original_parts[n:]) / RNG_DRAWS
                phen += overlap_matrix(head_parts[:n], head_parts[n:]) / RNG_DRAWS
            rho, p = infer(rain, phen)
            rho_raw, p_raw = infer(raw, phen)
            loo = []
            for axis, count in enumerate((len(keys["native"]), len(keys["introduced"]))):
                for j in range(count):
                    a = np.delete(rain, j, axis).ravel()
                    b = np.delete(phen, j, axis).ravel()
                    if len(np.unique(a)) > 1 and len(np.unique(b)) > 1:
                        from scipy.stats import rankdata
                        loo.append(float(np.corrcoef(rankdata(a), rankdata(b))[0, 1]))
            if not loo:
                raise ValueError("LEAVE_ONE_REGION_UNIDENTIFIABLE")
            output.append({
                "species": species, "head_dimensions": head, "h_std_pc": h,
                "status": "SOURCE_SUPPORTED_REGION_PAIR_TEST",
                "background_adjusted_rho": float(rho),
                "background_adjusted_permutation_p": float(p),
                "photo_climate_rho_same_regions": float(rho_raw),
                "photo_climate_permutation_p": float(p_raw),
                "rho_adjustment_change": float(rho - rho_raw),
                "leave_one_region_min_rho": min(loo), "leave_one_region_max_rho": max(loo),
                "n_native_regions": len(keys["native"]),
                "n_intro_regions": len(keys["introduced"]),
                "n_region_pairs": int(rain.size),
                "n_photos_in_photo_cohort": int(c.photos.sum()),
                "n_photo_cells": int(len(c)),
                "species_wcvp_background_pc1_pc2_variance": variance,
                "min_original_20deg_available_area_fraction": float(min(regional[k]["background_fraction"] for k in keys["native"] + keys["introduced"])),
                "min_photo_occupancy_on_available_background_fraction": float(min(min(all_support.values()), support_min)),
                "region_label_permutation": True,
                "same_fixed_photo_draws_for_bg_vs_raw_and_head": True,
                "availability_corrected_by_region_wcvp_background": True,
                "not_genetic_or_causal": True,
            })
    return {"status": "SOURCE_SUPPORTED",
            "n_supported_native_regions": len(keys["native"]),
            "n_supported_introduced_regions": len(keys["introduced"]),
            "n_species_WCVP_background_points": int(len(species_land)),
            "n_species_complete_background_points": int(len(complete))}, output


def _source_original():
    """Return confidential original observations and full WCVP land backgrounds."""
    original = (PRIVATE / "original-metadata.aesgcm").read_bytes()
    content = decrypt_checked(original, os.environ["AZAMI_METADATA_KEY"],
                              ORIGINAL_SEALED_SHA, ORIGINAL_SHA, ORIGINAL_AAD)
    photos = source_photos(content)
    del content
    authority = authority_source((PRIVATE / "authority350.zip").read_bytes())
    wild = source_roster(photos, authority)
    del photos
    wild = wild.loc[wild.taxon_name.isin(SPECIES_KEYS) &
                    wild.eligibility_reason.isin(("include", "range_introduced"))].copy()
    wild["range_state"] = np.where(wild.native_primary.eq(True), "native", "introduced")
    if len(wild) != 297622:
        raise ValueError("ORIGINAL_WCVP_WILD_PHOTO_SOURCE_DRIFT")
    for (sp, state), n in SOURCE_COUNTS.items():
        if len(wild[(wild.taxon_name == sp) & (wild.range_state == state)]) != n:
            raise ValueError("ORIGINAL_WCVP_SPECIES_RANGE_COUNTS_DRIFT")
    land = background_mesh(authority)
    wild, land, receipts = apply_original_five_CHELSA(
        wild, land, source_meta=Path("analysis/v4/environment_source_metadata.json"),
        work=PRIVATE / "region_bg_original_tiffs")
    return wild, land, authority, receipts


def _head_cells(wild, *, require_original_cohort=True):
    """Exactly photo-eligible observer×cell source, with optional diagnostic gate.

    The full global TIFF extraction can have a different missingness mask from
    the original five archived 5-degree caches. Compare the two rather than
    silently enlarging the historical photograph cohort.
    """
    src = (PRIVATE / "original-measurements.aesgcm").read_bytes()
    clear = decrypt_checked(src, os.environ["AZAMI_MEASUREMENTS_KEY"],
                            MEASUREMENT_SEALED_SHA, MEASUREMENT_SHA, MEASUREMENT_AAD)
    measured = pd.read_csv(io.BytesIO(clear), dtype={"obs_id": str}, low_memory=False,
                           usecols=["obs_id", *TRAITS])
    del clear
    if measured.obs_id.duplicated().any():
        raise ValueError("MEASUREMENT_OBS_DUPLICATION")
    d = wild.drop(columns=TRAITS, errors="ignore").merge(
        measured, on="obs_id", how="left", validate="one_to_one")
    d = d[d[TRAITS[0]].between(0, 180)].dropna(subset=[*TRAITS, *VARIABLES]).copy()
    d["coarse"] = angular_cells(d.latitude, d.longitude)
    d = d.sort_values("obs_id").drop_duplicates(
        ["taxon_name", "range_state", "coarse", "user_id"]).copy()
    d["i"] = (d.latitude + 90).clip(upper=179.99999).astype(int)
    d["j"] = ((d.longitude + 180) % 360).astype(int)
    d["block"] = d.i.floordiv(20) * 18 + d.j.floordiv(20)
    result = {}
    for sp, sub in d.groupby("taxon_name", sort=True):
        c = sub.groupby(["range_state", "block", "i", "j"])[[*VARIABLES, *TRAITS]].median().reset_index()
        c["photos"] = sub.groupby(["range_state", "block", "i", "j"]).size().to_numpy()
        # Provenance anchor against the original 2026-11-09 paired overlap run.
        if require_original_cohort and int(c.photos.sum()) != EXPECTED_PHOTOS[sp]:
            raise ValueError(
                "V5_ORIGINAL_PAIRED_PHOTO_COHORT_DRIFT_"
                + sp.replace(" ", "_")
                + "_got_" + str(int(c.photos.sum()))
                + "_expected_" + str(EXPECTED_PHOTOS[sp]))
        old_regions = tuple(
            sum(len(g) >= 8 for _, g in c[c.range_state == state].groupby("block"))
            for state in ("native", "introduced"))
        if require_original_cohort and old_regions != EXPECTED_REGIONS[sp]:
            raise ValueError("V5_ORIGINAL_PAIRED_REGIONS_DRIFT_"
                             + sp.replace(" ", "_")
                             + "_got_" + str(old_regions)
                             + "_expected_" + str(EXPECTED_REGIONS[sp]))
        result[sp] = c
    if set(result) != set(SPECIES_KEYS):
        raise ValueError("V5_ORIGINAL_TWO_SPECIES_NOT_PRESENT")
    return result


def _historical_photo_pixel_cache(wild):
    """Recreate the *original* V5 pairing cache, including authentic missingness.

    This intentionally follows .github/workflows/v5-coupling.yml. The separate
    S2 recover_original_five helper has an unrelated >=95% coverage gate that
    the historical V5 regional pairing did NOT impose. Missing pixels are never
    imputed; the eligible photo cohort is chosen only afterwards.
    """
    import subprocess
    import zipfile
    from analysis.v4.s1_s2_restore_original_five import CACHES, sha256
    from analysis.v4.native_cached_climate import physical_climate_for_focal

    root = PRIVATE / "original_five_photo_cache"
    root.mkdir(parents=True, exist_ok=True)
    receipts = {}
    for name, (artifact_id, expected_hash) in CACHES.items():
        dest = root / name
        dest.mkdir(parents=True, exist_ok=True)
        archive_path = root / (name + ".zip")
        with archive_path.open("wb") as f:
            subprocess.run(
                ["gh", "api",
                 f"repos/zuizui0223/azami/actions/artifacts/{artifact_id}/zip"],
                stdout=f, check=True)
        if sha256(archive_path) != expected_hash:
            raise ValueError("ORIGINAL_V5_PHOTO_CACHE_SHA_MISMATCH_" + name)
        with zipfile.ZipFile(archive_path) as z:
            members = z.namelist()
            if any(Path(member).is_absolute() or ".." in Path(member).parts
                   for member in members):
                raise ValueError("UNSAFE_ORIGINAL_V5_CACHE_ARCHIVE")
            z.extractall(dest)
        archive_path.unlink()
        receipts[name] = {"artifact_id": artifact_id,
                          "zip_sha256": expected_hash,
                          "source_cache_verified": True}
    cached, pixel_audit = physical_climate_for_focal(
        wild, root, Path("analysis/v4/environment_source_metadata.json"))
    photo = wild.copy()
    for name in VARIABLES:
        data = np.asarray(cached[name], float)
        if len(data) != len(wild):
            raise ValueError("HISTORICAL_V5_CACHED_CLIMATE_ROW_COUNT_DRIFT")
        photo[name] = data
        receipts[name]["finite_original_v5_pixels"] = int(np.isfinite(data).sum())
        receipts[name]["missingness_retained_without_imputation"] = True
    return photo, {
        "status": "ORIGINAL_V5_FIVE_CACHED_PIXEL_MASK_REPRODUCED",
        "archives": receipts,
        "pixel_diagnostic": pixel_audit,
        "all_original_pixel_missingness_retained": True,
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    wild, land, authority, receipts = _source_original()

    # Diagnostic source-parity gate. The original 2026-11-09 regional overlap
    # consumed five SHA-verified CACHE pixels, not the newly recovered global
    # TIFF completeness mask. Use exact old cache values for the SAME photo
    # cohort while using full original TIFFs for WCVP background availability.
    # No invented climate values or coalescing of missing pixels.
    photo_wild, photo_receipt = _historical_photo_pixel_cache(wild)
    comparison = {}
    for variable in VARIABLES:
        new = wild[variable].to_numpy(float)
        old = photo_wild[variable].to_numpy(float)
        common = np.isfinite(new) & np.isfinite(old)
        max_abs = float(np.max(np.abs(new[common] - old[common]))) if common.any() else None
        comparison[variable] = {
            "old_cache_complete": int(np.isfinite(old).sum()),
            "full_original_complete": int(np.isfinite(new).sum()),
            "full_only": int((np.isfinite(new) & ~np.isfinite(old)).sum()),
            "cache_only": int((~np.isfinite(new) & np.isfinite(old)).sum()),
            "max_abs_original_pixel_difference_on_common": max_abs,
        }
        if max_abs is None or max_abs > 1e-7:
            raise ValueError("ORIGINAL_CACHE_FULL_PIXEL_MISMATCH_" + variable)
    full_cells = _head_cells(wild, require_original_cohort=False)
    cells = _head_cells(photo_wild, require_original_cohort=True)
    parity = {
        species: {
            "full_original_tiff_photo_eligible": int(full_cells[species].photos.sum()),
            "historical_cached_photo_eligible": int(cells[species].photos.sum()),
            "historical_paired_expected": EXPECTED_PHOTOS[species],
        }
        for species in SPECIES_KEYS
    }
    print({"source_parity_checked": True,
           "photo_eligible_original_vs_full": parity,
           "original_pixel_missingness_diagnostics": comparison}, flush=True)
    rows, audits = [], {}
    for sp in SPECIES_KEYS:
        status, result = _one_species(cells[sp], land, authority, sp)
        audits[sp] = status
        rows.extend(result)
    if not rows:
        raise ValueError("NO_SOURCE_SUPPORTED_BACKGROUND_REGION_PAIR_ANALYSIS")
    t = pd.DataFrame(rows)
    main = (t.head_dimensions == "three_trait") & (t.h_std_pc == .5) & (t.status == "SOURCE_SUPPORTED_REGION_PAIR_TEST")
    if int(main.sum()) != 2:
        raise ValueError("PRIMARY_BOTH_SPECIES_UNSUPPORTED")
    t["primary_holm_p"] = np.nan
    t.loc[main, "primary_holm_p"] = multipletests(
        t.loc[main, "background_adjusted_permutation_p"].to_numpy(float), method="holm")[1]
    # No exact coordinates, names of human photographers, rowwise measurements,
    # region keys, or source pixel values are written to public outputs.
    results = {
        "schema": "azami_v5_background_corrected_regional_niche_head_paired_20261110_v1",
        "status": "SOURCE_VERIFIED_NEW_ESTIMATE",
        "main": json.loads(t.loc[main].to_json(orient="records")),
        "sensitivities": json.loads(t.loc[~main].to_json(orient="records")),
        "support": audits,
        "original_source_wild_n": 297622,
        "source_rasters": {key: {
            "source_SHA256": v["source_SHA256"],
            "original_pixels_only": v["original_pixels_only"]}
            for key, v in receipts.items()},
        "historical_photo_cache_source": photo_receipt["status"],
        "original_pixel_missingness_comparison": comparison,
        "photo_cohort_source_parity": parity,
        "no_individual_records_published": True,
        "old_published_original_rhos_are_not_background_corrected": True,
        "not_genetic_adaptation": True,
    }
    (OUT / "v5_background_paired_overlap_summary.json").write_text(
        json.dumps(results, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8")
    print({"status": results["status"],
           "primary": [{k: q[k] for k in ("species", "background_adjusted_rho",
                "primary_holm_p", "n_region_pairs")} for q in results["main"]]}, flush=True)


if __name__ == "__main__":
    main()
