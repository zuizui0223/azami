"""V4 BIO18 photographic-head-angle source-cohort geographic sensitivity.

This does NOT refit the pre-existing S1 80-test Holm family.  It first reproduces
the exact S1 four-variable BIO18 among-taxon coefficient on the original QC
cohort and then compares predefined longitude-sector and spherical-geography
adjustments, balanced source-cell weighting, and regional holdouts.

Only pooled model coefficients and counts are written. Original photograph
locations, taxon keys, observation IDs, and uploader IDs remain private.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

from analysis.v4.s1_s2_extreme_lego import ALL_ENV, EXTRA
from analysis.v4.s1_s2_vif5_si10_alltraits import (
    FOUR_WARM, selected_source, standardized_y, regression_design,
)

ANGLE = "orientation_image_vertical_angle"
RAIN = "chelsa_bio18"
EXPECTED_ORIGINAL = 313350
EXPECTED_ANGLE = 168907
EXPECTED_TAXA = 178
EXPECTED_BETA = 0.04544684307686103
SECTORS = ("W180_W030", "W030_E060", "E060_E180")
GEO_COLS = ("g_x", "g_y", "g_z", "g_xy", "g_xz", "g_yz",
            "g_x2_minus_y2", "g_3z2_minus1")
GEO_CONTROLS = ("none", "longitude_sector", "spherical_degree2")
WEIGHTINGS = ("photo_weighted", "equal_5degree_cell")
COHORT_MIN = 20


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def geography(lat_deg, lon_deg):
    lat, lon = np.asarray(lat_deg, float), np.asarray(lon_deg, float)
    if lat.shape != lon.shape or not np.isfinite(lat).all() or not np.isfinite(lon).all():
        raise ValueError("Missing source geographic coordinates")
    if np.any(np.abs(lat) > 90) or np.any(np.abs(lon) > 180):
        raise ValueError("Geographic coordinate outside geographic range")
    la, lo = np.deg2rad(lat), np.deg2rad(lon)
    x, y, z = np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)
    features = np.column_stack((x, y, z, x*y, x*z, y*z, x*x-y*y, 3*z*z-1))
    sectors = np.select((lon < -30, lon < 60), (0, 1), default=2).astype(int)
    return features, sectors


def exact_eligible(native_path: Path, extended_path: Path):
    native = pd.read_csv(native_path, dtype={"obs_id": str}, low_memory=False)
    env = pd.read_csv(extended_path, dtype={"obs_id": str}, low_memory=False)
    if (len(native) != EXPECTED_ORIGINAL or not native.obs_id.is_unique
        or not env.obs_id.is_unique or len(env) != EXPECTED_ORIGINAL
        or set(EXTRA) - set(env) or set(ALL_ENV[:5]) - set(native)
        or set(EXTRA) & set(native)):
        raise ValueError("ORIGINAL_SOURCE_OR_RASTER_JOIN_MISMATCH")
    sample = native.merge(env[["obs_id", *EXTRA]], on="obs_id",
                          how="left", validate="one_to_one")
    selected, status = selected_source(sample, [ANGLE], min_taxa=20)
    if selected is None or len(selected) != EXPECTED_ANGLE or selected.accepted_key.nunique() != EXPECTED_TAXA:
        raise ValueError(f"FROZEN_S1_ANGLE_COHORT_CHANGED: {status}")
    if selected.obs_id.duplicated().any() or not np.isfinite(selected[list(FOUR_WARM)].to_numpy(float)).all():
        raise ValueError("DUPLICATED_OR_INVALID_PRIMARY_OBSERVATIONS")
    features, sec = geography(selected.latitude.to_numpy(float),
                               selected.longitude.to_numpy(float))
    return selected, features, sec


def reproduce_frozen_baseline(selected):
    design, groups = regression_design(selected, list(FOUR_WARM))
    y = standardized_y(selected[ANGLE].to_numpy(float))
    counts = selected.groupby("accepted_key").size()
    weights = selected.accepted_key.map(1/counts).to_numpy(float)
    fit = sm.WLS(y, design, weights=weights).fit(
        cov_type="cluster", cov_kwds={"groups": groups}, use_t=False)
    estimated = float(fit.params[1])
    if not np.isclose(estimated, EXPECTED_BETA, rtol=0, atol=3e-7):
        raise ValueError(f"FROZEN_S1_BIO18_SLOPE_NOT_REPRODUCED: {estimated}")
    return {
        "n_eligible": len(selected),
        "n_taxa": int(selected.accepted_key.nunique()),
        "standardised_beta": estimated,
        "cluster_robust_p_raw": float(fit.pvalues[1]),
        "unchanged_original_Holm_P": 0.009206350580895119,
        "original_eight_biological_constructs_80_test_family_not_refitted": True,
    }


def taxon_rosters(selected, geofeatures, sectors):
    """Return only private in-memory taxon means; public output is aggregate."""
    keys = selected.accepted_key.astype(str).to_numpy()
    values = selected[list(FOUR_WARM)].to_numpy(float)
    response = standardized_y(selected[ANGLE].to_numpy(float))
    means, standard_deviations = values.mean(axis=0), values.std(axis=0)
    if np.any(standard_deviations <= 0):
        raise ValueError("CONSTANT_SOURCE_CLIMATE")
    photo = pd.DataFrame({
        "taxon": keys,
        "response": response,
        **{k: (values[:, i]-means[i])/standard_deviations[i]
           for i, k in enumerate(FOUR_WARM)},
        **{g: geofeatures[:, i] for i, g in enumerate(GEO_COLS)},
        **{f"sector_{label}": (sectors == i).astype(float)
           for i, label in enumerate(SECTORS)},
    })
    photo["cell_lat"] = np.floor((selected.latitude.to_numpy(float)+90)/5).astype(int)
    photo["cell_lon"] = np.floor((selected.longitude.to_numpy(float)+180)/5).astype(int)
    quantities = ["response", *FOUR_WARM, *GEO_COLS,
                  *(f"sector_{label}" for label in SECTORS)]
    taxa_photo = photo.groupby("taxon", sort=True)[quantities].mean()
    cells = photo.groupby(["taxon", "cell_lat", "cell_lon"], sort=True)[quantities].mean()
    taxa_equal_cell = cells.groupby(level="taxon", sort=True)[quantities].mean()
    rosters = {"photo_weighted": taxa_photo, "equal_5degree_cell": taxa_equal_cell}
    if set(taxa_photo.index) != set(taxa_equal_cell.index) or len(taxa_photo) != EXPECTED_TAXA:
        raise ValueError("SOURCE_TAXA_LOST_IN_BALANCING")
    # The original equal-taxon observation-level between effects are exactly
    # equivalent to a taxon-mean regression with globally standardised inputs.
    return rosters, {"n_source_cell_taxon_pairs": len(cells),
                     "n_original_photo_observations": len(photo)}


def _fixed_design(data: pd.DataFrame, adjustment: str):
    if adjustment not in GEO_CONTROLS:
        raise ValueError("UNKNOWN_GEO_ADJUSTMENT")
    columns = list(FOUR_WARM)
    if adjustment == "longitude_sector":
        columns += [f"sector_{label}" for label in SECTORS[:2]]
    elif adjustment == "spherical_degree2":
        columns += list(GEO_COLS)
    X = data[columns].to_numpy(float)
    if not np.isfinite(X).all():
        raise ValueError("NONFINITE_TAXON_CLIMATE_OR_GEOGRAPHY")
    X = np.column_stack([np.ones(len(data)), X])
    if np.linalg.matrix_rank(X) != X.shape[1] or np.linalg.cond(X) > 1e7:
        raise ValueError("REGIONAL_TAXON_DESIGN_UNIDENTIFIABLE")
    return X


def fit_taxon(data, adjustment="none"):
    """One original-scale BIO18 coefficient and heteroskedasticity-robust CI."""
    if len(data) < COHORT_MIN:
        return {"status": "too_few_taxa", "n_taxa": len(data)}
    try:
        X = _fixed_design(data, adjustment)
        y = data.response.to_numpy(float)
        fit = sm.OLS(y, X).fit(cov_type="HC3", use_t=False)
        ci = fit.conf_int(alpha=0.05)
        coefficient, interval = float(fit.params[1]), [float(v) for v in ci[1]]
        if not all(np.isfinite([coefficient, *interval, float(fit.pvalues[1])])):
            raise ValueError("NONFINITE_GEO_MODEL_RESULT")
        return {
            "status": "estimated",
            "n_taxa": len(data),
            "n_parameters": X.shape[1],
            "standardised_beta": coefficient,
            "HC3_pointwise_95CI": interval,
            "HC3_p_uncorrected": float(fit.pvalues[1]),
            "geographic_covariates": adjustment,
            "taxa_equal_weight": True,
            "not_phylogenetic_or_causal": True,
            "no_holm_correction_for_new_sensitivity_family": True,
        }
    except (ValueError, np.linalg.LinAlgError) as exc:
        return {"status": "not_identifiable", "n_taxa": len(data),
                "reason": str(exc)[:100]}


def summarise_geography(rosters):
    result = {}
    for weighting in WEIGHTINGS:
        frame = rosters[weighting]
        runs = {adj: fit_taxon(frame, adj) for adj in GEO_CONTROLS}
        if weighting == "photo_weighted" and runs["none"]["status"] == "estimated":
            if not np.isclose(runs["none"]["standardised_beta"], EXPECTED_BETA,
                              rtol=0, atol=3e-7):
                raise ValueError("TAXON_MEAN_BASELINE_DIFFERS_FROM_FROZEN_S1")
        dominant = frame[[f"sector_{label}" for label in SECTORS]].to_numpy(float).argmax(axis=1)
        regional = {}
        holdouts = {}
        for j, name in enumerate(SECTORS):
            part = frame.iloc[np.flatnonzero(dominant == j)]
            rest = frame.iloc[np.flatnonzero(dominant != j)]
            regional[name] = fit_taxon(part, "none")
            holdouts[name] = fit_taxon(rest, "none")
            regional[name]["is_longitude_sector_not_independent_lineage_or_continent"] = True
            holdouts[name]["heldout_n_taxa"] = len(part)
        result[weighting] = {
            "global_geographic_adjustments": runs,
            "within_predeclared_longitude_sectors": regional,
            "leave_one_predeclared_sector_out": holdouts,
            "n_taxa_with_max_sector_fraction_below_0_70": int(
                (frame[[f"sector_{label}" for label in SECTORS]].max(axis=1) < 0.70).sum()),
        }
    return result


def execute(native_path, extended_path, source_receipt, output):
    sample, geo, sectors = exact_eligible(native_path, extended_path)
    baseline = reproduce_frozen_baseline(sample)
    rosters, counts = taxon_rosters(sample, geo, sectors)
    sensitivity = summarise_geography(rosters)
    receipt = json.loads(Path(source_receipt).read_text(encoding="utf-8"))
    if receipt.get("status") != "COMPLETE_SEVEN_ORIGINAL_RASTERS" or len(receipt.get("rasters", [])) != len(EXTRA):
        raise ValueError("INCOMPLETE_CHELSA_OR_ELEVATION_RESTORATION")
    output = {
        "schema": "azami_v4_S1_BIO18_exact_cohort_predeclared_geography_sensitivity_v1",
        "status": "EMPIRICAL_BIO18_SOURCE_VERIFIED_NOT_IN_ORIGINAL_HOLM_FAMILY",
        "primary_model": "CHELSA BIO18 BIO05 RSDS and WorldClim elevation",
        "native_original_rows": EXPECTED_ORIGINAL,
        "native_angle_eligible_records": EXPECTED_ANGLE,
        "native_angle_taxa": EXPECTED_TAXA,
        "frozen_original_baseline_reproduced": baseline,
        "taxon_photo_and_cell_counts": counts,
        "source_sha256": {
            "decoded_native_table": digest(native_path),
            "native_only_extra_raster_join": digest(extended_path),
            "rast_receipt": digest(source_receipt),
            "original_source_archive_locked": True,
        },
        "sensitivity": sensitivity,
        "scientific_limits": {
            "different_regression_estimand_with_new_geographic_covariates": True,
            "new_pointwise_CI_not_familywise_80_test_claim": True,
            "three_longitude_sectors_not_biogeographic_independent_replicates": True,
            "spherical_basis_not_phylogenetic_control": True,
            "image_orientation_not_gravity": True,
            "not_rainfall_shielding_adaptation_test": True,
            "source_taxa_not_all_Cirsium_taxa": True,
            "source_original_coordinates_names_photoIDs_not_published": True,
        },
    }
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(output, indent=2, allow_nan=False, ensure_ascii=False)+"\n",
                    encoding="utf-8")
    print(json.dumps({"status": output["status"],
                      "n": EXPECTED_ANGLE, "taxa": EXPECTED_TAXA,
                      "geographic_fits": sensitivity}, ensure_ascii=False), flush=True)
    return output


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("native", "extended", "source_receipt", "output"):
        p.add_argument("--"+name.replace("_","-"), type=Path, required=True)
    args = p.parse_args()
    execute(args.native, args.extended, args.source_receipt, args.output)


if __name__ == "__main__":
    main()
