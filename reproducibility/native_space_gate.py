"""Native-only capitulum phenotype-space FEASIBILITY audit, not functional adaptation.

Only accepts the frozen already-prepared, decrypted native-observation table.
No WCVP queries, original-image downloads, raster acquisition, model searches,
private row export, or inference about genetic/functional modularity.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2, spearmanr

RAW_UNCORRECTED_SHA256 = "8ceed49e4ac73b31f9a5c726ca84f4198bfc3c51e9ee43303e29524e915a5c29"
EXPECTED_SHA256 = "342d78e9a8a960e8d9bdcf0bc2cdacd13bee48214a8b962419d5a71d5aef46ed"
SOURCE_CORRECTION = "PR143 / Actions run 39422378058; GSP-invalid sentinel masked"
EXPECTED_N = 313350
CORE3 = (
    "orientation_image_vertical_angle",
    "corolla_lab_chroma",
    "capitulum_outline_aspect_ratio",
)
CORE9 = (
    "orientation_image_vertical_angle",
    "corolla_lab_lightness",
    "corolla_lab_chroma",
    "corolla_hue_sin",
    "corolla_hue_cos",
    "capitulum_outline_aspect_ratio",
    "capitulum_outline_circularity",
    "capitulum_outline_solidity",
    "capitulum_width_profile_cv",
)
CLIMATE = ("chelsa_bio12", "chelsa_rsds_mean")



def verify_corrected_native_input(path: Path) -> str:
    """Accept only the affine-decoded, checksum-locked existing PR143 source."""
    h = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1 << 20), b""):
            h.update(block)
    sha = h.hexdigest()
    if sha == RAW_UNCORRECTED_SHA256:
        raise ValueError(
            "Uncorrected environment input: GSP invalid sentinels and missing "
            "CHELSA affine decoding. Supply the already-corrected PR143 output."
        )
    if sha != EXPECTED_SHA256:
        raise ValueError("Unknown native input bytes: expected the corrected PR143 SHA256")
    return sha


def _eligible(frame: pd.DataFrame, endpoints: tuple[str, ...]) -> np.ndarray:
    """Require source-provided endpoint QC flags AND finite numeric measurements."""
    ok = np.ones(len(frame), dtype=bool)
    for name in endpoints:
        flags = frame[name + "__eligible"].astype("string").str.lower().str.strip()
        if not flags.isin(("true", "false")).all():
            raise ValueError("Missing or invalid measurement eligibility: " + name)
        values = pd.to_numeric(frame[name], errors="coerce").to_numpy(float)
        ok &= flags.eq("true").to_numpy(dtype=bool) & np.isfinite(values)
    return ok


def _summary(values: list[float]) -> dict:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if not len(arr):
        return {"n": 0, "median": None, "q25": None, "q75": None}
    return {
        "n": int(len(arr)),
        "median": float(np.median(arr)),
        "q25": float(np.quantile(arr, .25)),
        "q75": float(np.quantile(arr, .75)),
    }



def _correlation_pairs(covariance: np.ndarray) -> dict:
    """Off-diagonal correlations for the same three independently labelled photo axes."""
    diagonal = np.diag(covariance)
    answer = {}
    for i in range(3):
        for j in range(i + 1, 3):
            denom = math.sqrt(max(0.0, diagonal[i]) * max(0.0, diagonal[j]))
            estimate = float(covariance[i, j] / denom) if denom > 1e-12 else None
            answer[CORE3[i] + "__" + CORE3[j]] = estimate
    return answer


def _scale_specific_core3_covariance(
    distinct: pd.DataFrame, qualified, centers: np.ndarray, scales: np.ndarray,
) -> dict:
    """Identity W+B=T for equal-taxon/equal-spatial-cell photo phenotypes.

    This is descriptive covariance in the same matched source observations,
    not genetic covariance, a phylogenetic contrast, or a hypothesis test.
    """
    group_means, group_covariances = [], []
    groups = []
    for taxon in qualified:
        values = (
            distinct.loc[distinct.accepted_key.eq(taxon), list(CORE3)].to_numpy(float)
            - centers
        ) / scales
        if len(values) < 2 or not np.isfinite(values).all():
            raise ValueError("Degenerate within-taxon photo covariance")
        mean = values.mean(axis=0)
        centered = values - mean
        group_means.append(mean)
        group_covariances.append(centered.T @ centered / len(values))
        groups.append(values)
    m = np.asarray(group_means)
    w = np.mean(group_covariances, axis=0)
    mean = m.mean(axis=0)
    b = (m - mean).T @ (m - mean) / len(m)
    total = w + b

    direct_pooled = np.zeros((3, 3), dtype=float)
    for values in groups:
        delta = values - mean
        direct_pooled += (delta.T @ delta) / len(values) / len(groups)
    error = float(np.max(np.abs(direct_pooled - total)))
    if not np.isfinite(error) or error > 1e-8:
        raise ValueError("Pooled covariance identity mismatch")

    within = _correlation_pairs(w)
    among = _correlation_pairs(b)
    contrasts = {
        key: (float(among[key] - within[key])
              if among[key] is not None and within[key] is not None else None)
        for key in within
    }
    fractions = np.diag(b) / np.maximum(np.diag(total), 1e-20)
    return {
        "status": "descriptive_native_photo_covariance_decomposition",
        "n_accepted_taxa": int(len(groups)),
        "axes": list(CORE3),
        "sampling": "one original observation per accepted taxon x 0.5-degree geographic cell",
        "within_taxon_pairwise_correlations": within,
        "among_taxon_pairwise_correlations": among,
        "among_minus_within_pairwise_correlations": contrasts,
        "between_taxon_variance_fraction_by_axis": {
            axis: float(value) for axis, value in zip(CORE3, fractions)
        },
        "identity_max_abs_residual": error,
        "interpretation_limit": "W and B describe photographic phenotypic covariance, not separate genetic or developmental modules. No phylogeny/photographer adjustment, uncertainty intervals, or significance tests.",
    }


def _environment_specific_among_covariance(
    means: pd.DataFrame, environments: pd.DataFrame,
) -> dict:
    """Preselected current BIO12/RSDS coarse halves, with no p or causal claim."""
    result = {}
    for name in CLIMATE:
        aligned = means.join(environments[[name]], how="inner").dropna()
        axis = aligned[name]
        midpoint = float(axis.median()) if len(aligned) else None
        lower = aligned.loc[axis <= midpoint] if midpoint is not None else aligned.iloc[:0]
        upper = aligned.loc[axis > midpoint] if midpoint is not None else aligned.iloc[:0]
        entry = {
            "eligible_taxa": int(len(aligned)),
            "lower_n_taxa": int(len(lower)),
            "higher_n_taxa": int(len(upper)),
            "split": "native taxon median climate, lower <= pooled qualified-taxon median",
            "status": "insufficient_taxa_without_inferential_interpretation",
            "lower_among_pairwise_correlations": None,
            "higher_among_pairwise_correlations": None,
        }
        # Ten distinct accepted taxa per side are a feasibility floor, not power.
        if len(lower) >= 10 and len(upper) >= 10:
            entry["status"] = "descriptive_without_phylogenetic_or_spatial_adjustment"
            entry["lower_among_pairwise_correlations"] = _correlation_pairs(
                np.cov(lower[list(CORE3)].to_numpy(float), rowvar=False, bias=True)
            )
            entry["higher_among_pairwise_correlations"] = _correlation_pairs(
                np.cov(upper[list(CORE3)].to_numpy(float), rowvar=False, bias=True)
            )
        result[name] = entry
    return result



def analyze_native_space(
    data: pd.DataFrame, *, expected_n: int | None = None,
    rarefaction_n: int = 20, replicates: int = 100, seed: int = 20261020,
) -> dict:
    """Return ONLY aggregate diagnostics: no taxa, obs IDs or exact coordinates."""
    if rarefaction_n < 6 or replicates < 2:
        raise ValueError("Need at least six cells and two sampling replicates")
    required = set(CORE9) | {x + "__eligible" for x in CORE9}
    required |= {"obs_id", "accepted_key", "latitude", "longitude", *CLIMATE}
    missing = required - set(data.columns)
    if missing:
        raise ValueError("Prepared native table lacks columns: " + ", ".join(sorted(missing)))
    if expected_n is not None and len(data) != expected_n:
        raise ValueError("Not the fixed native-primary frame")
    if data.obs_id.isna().any() or not data.obs_id.is_unique:
        raise ValueError("Observation IDs must be nonmissing and unique")
    if data.accepted_key.isna().any():
        raise ValueError("Accepted native taxon identities must be complete")

    core3 = _eligible(data, CORE3)
    core9 = _eligible(data, CORE9)
    numeric = data.loc[:, [*CORE3, *CLIMATE, "latitude", "longitude"]].apply(
        pd.to_numeric, errors="coerce"
    )
    env = np.isfinite(numeric[list(CLIMATE)].to_numpy(float)).all(axis=1)
    coord = np.isfinite(numeric[["latitude", "longitude"]].to_numpy(float)).all(axis=1)
    coord &= numeric.latitude.between(-90, 90).to_numpy() & numeric.longitude.between(-180, 180).to_numpy()
    subset = data.loc[core3 & env & coord, ["obs_id", "accepted_key"]].copy()
    for col in [*CORE3, *CLIMATE, "latitude", "longitude"]:
        subset[col] = numeric.loc[core3 & env & coord, col].to_numpy(float)

    report = {
        "status": "native_image_phenotype_space_feasibility_not_functional_validation",
        "input_native_observations": int(len(data)),
        "source_is_native_only": expected_n == EXPECTED_N,
        "core3": list(CORE3),
        "complete_core3_observations": int(core3.sum()),
        "complete_core9_observations": int(core9.sum()),
        "core3_with_bio12_rsds_and_coordinates": int(len(subset)),
        "core3_taxa": int(data.loc[core3, "accepted_key"].nunique()),
        "core9_taxa": int(data.loc[core9, "accepted_key"].nunique()),
        "environment_axes": list(CLIMATE),
        "spatial_thinning": "one source observation per accepted taxon × 0.5-degree cell",
        "rarefaction_n_distinct_cells_per_taxon": int(rarefaction_n),
        "replicates": int(replicates),
        "measurement_semantics": "image-relative angle, uncalibrated JPEG chroma, viewpoint-dependent 2D aspect",
        "inference": "exploratory, no p values or adaptation/functional claims",
    }
    if len(subset) == 0:
        return dict(report, occupancy_status="no_complete_joint_records")

    subset["_latcell"] = np.floor((subset.latitude + 90) * 2).astype(int)
    subset["_loncell"] = np.floor((subset.longitude + 180) * 2).astype(int)
    subset["_stable"] = pd.util.hash_pandas_object(
        subset.obs_id.astype(str), index=False
    ).to_numpy(dtype="uint64")
    distinct = (
        subset.sort_values(["accepted_key", "_latcell", "_loncell", "_stable"])
        .drop_duplicates(["accepted_key", "_latcell", "_loncell"], keep="first")
        .sort_values(["accepted_key", "_latcell", "_loncell"])
    )
    cell_counts = distinct.groupby("accepted_key", sort=True).size()
    qualified = cell_counts[cell_counts >= rarefaction_n].index
    report["core3_taxa_with_enough_distinct_cells"] = int(len(qualified))
    report["core3_distinct_taxon_cells"] = int(len(distinct))
    if len(qualified) < 4:
        return dict(report, occupancy_status="too_few_geographically_supported_taxa")

    # Equal taxon contribution defines z-scales, not observation-count weighting.
    medians = distinct.groupby("accepted_key", sort=True)[list(CORE3)].median()
    centers = medians.median().to_numpy(float)
    deviations = np.abs(medians.to_numpy(float) - centers)
    scales = 1.4826 * np.median(deviations, axis=0)
    fallback = medians.std(axis=0, ddof=1).to_numpy(float)
    scales = np.where(scales > 1e-10, scales, fallback)
    if not np.isfinite(scales).all() or np.any(scales <= 1e-10):
        return dict(report, occupancy_status="degenerate_module_scaling")

    # Primary matched-cohort estimand: covariation at two biological scales.
    # Volume is a secondary summary of this same joint phenotypic geometry.
    report["core3_scale_specific_covariance"] = _scale_specific_core3_covariance(
        distinct, qualified, centers, scales
    )

    rng = np.random.default_rng(seed)
    constant = (4.0 * math.pi / 3.0) * chi2.ppf(.95, df=3) ** 1.5
    volume_by_taxon: dict[object, float] = {}
    rao_by_taxon: dict[object, float] = {}
    for taxon in qualified:
        values = (
            distinct.loc[distinct.accepted_key == taxon, list(CORE3)]
            .to_numpy(float)
            - centers
        ) / scales
        volumes, raos = [], []
        for _ in range(replicates):
            selected = values[rng.choice(len(values), size=rarefaction_n, replace=False)]
            covariance = np.cov(selected, rowvar=False)
            determinant = float(np.linalg.det(covariance))
            if determinant <= 1e-12 or not np.isfinite(determinant):
                continue
            volumes.append(constant * math.sqrt(determinant))
            raos.append(math.sqrt(2.0 * float(np.trace(covariance))))
        if volumes:
            volume_by_taxon[taxon] = float(np.median(volumes))
            rao_by_taxon[taxon] = float(np.median(raos))

    report["core3_taxa_with_nondegenerate_covariance"] = len(volume_by_taxon)
    report["rarefied_gaussian_95pct_ellipsoid_volume_proxy"] = _summary(list(volume_by_taxon.values()))
    report["rarefied_rms_pairwise_distance_proxy"] = _summary(list(rao_by_taxon.values()))
    report["occupancy_status"] = (
        "descriptive_space_proxy_available" if len(volume_by_taxon) >= 4
        else "insufficient_nondegenerate_taxa"
    )
    # Between-taxon climate is computed over the WHOLE native frame,
    # not selectively restricted to observations with image-measurable traits.
    climate = data[["accepted_key", *CLIMATE]].copy()
    for col in CLIMATE:
        climate[col] = pd.to_numeric(climate[col], errors="coerce")
    by_taxon = climate.groupby("accepted_key", sort=True)[list(CLIMATE)].median()
    rows = [(taxon, vol) for taxon, vol in volume_by_taxon.items() if taxon in by_taxon.index]
    # Environment can change the *joint arrangement*, not merely the volume.
    # This conditional comparison is descriptive, not a fitted eco-evolutionary model.
    qualified_means = (
        distinct.loc[distinct.accepted_key.isin(qualified)]
        .groupby("accepted_key", sort=True)[list(CORE3)].mean()
    )
    report["core3_among_covariance_by_climate_halves"] = _environment_specific_among_covariance(
        qualified_means, by_taxon
    )
    report["exploratory_volume_climate_spearman_no_pvalues"] = {}
    for col in CLIMATE:
        pairs = [(math.log(vol), float(by_taxon.at[taxon, col])) for taxon, vol in rows
                 if vol > 0 and np.isfinite(by_taxon.at[taxon, col])]
        if len(pairs) >= 6 and np.std([x[0] for x in pairs]) > 0 and np.std([x[1] for x in pairs]) > 0:
            rho = float(spearmanr([x[0] for x in pairs], [x[1] for x in pairs]).statistic)
        else:
            rho = None
        report["exploratory_volume_climate_spearman_no_pvalues"][col] = {
            "taxa": len(pairs), "rho": rho,
        }
    report["qualifications"] = [
        "Gaussian ellipsoid volume in a three-proxy robustly standardized PHOTO space, NOT a density-based functional hypervolume",
        "Rarefaction quantiles are sampling variation, NOT independent ecological confidence intervals",
        "0.5-degree spatial cells are angular approximations, not equal-area ecological replicates",
        "Source observation is not a validated unique plant, and photographer effects are not identified",
        "Volume-environment Spearman estimates are descriptive, uncorrected and conditioned on image measurability",
        "No kernel hypervolume, genetic covariance, selection, adaptation or geographical transfer is demonstrated",
        "A pooled covariance is partitioned into equal-taxon within and between components; statistical identity is not biological causation",
        "Climate-half trait correlations are descriptive and may be driven by different related taxa in each half",
    ]
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native-input", required=True, type=Path, help="Existing decrypted private native-observation CSV.gz")
    parser.add_argument("--output", required=True, type=Path, help="Aggregate JSON ONLY; no observation export")
    parser.add_argument("--rarefaction-n", type=int, default=20)
    parser.add_argument("--replicates", type=int, default=100)
    args = parser.parse_args()
    verify_corrected_native_input(args.native_input)
    frame = pd.read_csv(args.native_input, dtype={"obs_id": str}, low_memory=False)
    report = analyze_native_space(
        frame, expected_n=EXPECTED_N,
        rarefaction_n=args.rarefaction_n, replicates=args.replicates,
    )
    report["source_prepared_table_sha256"] = EXPECTED_SHA256
    report["source_environment_decoding_status"] = SOURCE_CORRECTION
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": report["occupancy_status"],
        "complete_core3": report["complete_core3_observations"],
        "qualified_taxa": report["core3_taxa_with_enough_distinct_cells"]
    }))


if __name__ == "__main__":
    main()
