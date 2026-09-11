#!/usr/bin/env python3
"""Pilot external validation of image-derived Cirsium head orientation.

This targeted Supplementary validation compares the frozen v2 image-derived
species medians with independent botanical descriptions that explicitly state
capitulum orientation. It then asks whether the same literature state is
aligned with frozen species-level BIO12 medians. It does not rescreen traits or
environmental predictors, alter the primary v2 atlas, or close the independent
measurement-validation gate.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

ORIENTATION_ENDPOINT = "orientation_image_vertical_angle"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--contract", required=True, type=Path)
    p.add_argument("--evidence", required=True, type=Path)
    p.add_argument("--traits-long", required=True, type=Path)
    p.add_argument("--environment", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    return p.parse_args()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def exact_group_permutation(values: np.ndarray, group: np.ndarray, statistic: str) -> tuple[float, float, int]:
    values = np.asarray(values, dtype=float)
    group = np.asarray(group, dtype=int)
    k = int(group.sum())
    if k <= 0 or k >= len(group):
        raise ValueError("Both literature-orientation groups are required")

    def stat(mask: np.ndarray) -> float:
        a = values[mask]
        b = values[~mask]
        if statistic == "median":
            return float(np.median(a) - np.median(b))
        if statistic == "mean":
            return float(np.mean(a) - np.mean(b))
        raise ValueError(statistic)

    observed = stat(group.astype(bool))
    simulated: list[float] = []
    for selected in itertools.combinations(range(len(values)), k):
        mask = np.zeros(len(values), dtype=bool)
        mask[list(selected)] = True
        simulated.append(stat(mask))
    sims = np.asarray(simulated, dtype=float)
    # Direction is fixed by the frozen v2 candidate: nodding > upright.
    p_one_sided = float(np.mean(sims >= observed - 1e-12))
    return observed, p_one_sided, int(len(sims))


def summarize_response(frame: pd.DataFrame, column: str) -> dict[str, object]:
    upright = frame.loc[frame["state_score"].eq(0), column].to_numpy(float)
    nodding = frame.loc[frame["state_score"].eq(1), column].to_numpy(float)
    rho, rho_p = spearmanr(frame["state_score"].to_numpy(float), frame[column].to_numpy(float))
    u, u_p = mannwhitneyu(nodding, upright, alternative="greater", method="exact")
    median_diff, median_p, n_perms = exact_group_permutation(
        frame[column].to_numpy(float), frame["state_score"].to_numpy(int), "median"
    )
    mean_diff, mean_p, _ = exact_group_permutation(
        frame[column].to_numpy(float), frame["state_score"].to_numpy(int), "mean"
    )
    return {
        "n_upright_taxa": int(len(upright)),
        "n_nodding_taxa": int(len(nodding)),
        "upright_median": float(np.median(upright)),
        "nodding_median": float(np.median(nodding)),
        "nodding_minus_upright_median": median_diff,
        "exact_one_sided_permutation_p_median_difference": median_p,
        "upright_mean": float(np.mean(upright)),
        "nodding_mean": float(np.mean(nodding)),
        "nodding_minus_upright_mean": mean_diff,
        "exact_one_sided_permutation_p_mean_difference": mean_p,
        "exact_permutation_assignments": n_perms,
        "spearman_rho_binary_state": float(rho),
        "spearman_two_sided_p": float(rho_p),
        "mann_whitney_u": float(u),
        "mann_whitney_exact_one_sided_p": float(u_p),
        "direction_concordant_with_v2_candidate": bool(np.median(nodding) > np.median(upright)),
    }


def main() -> None:
    args = parse_args()
    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    if contract.get("analysis_id") != "literature_orientation_external_validation_pilot_v1":
        raise ValueError("Unexpected literature-validation contract")
    expected_trait_sha = str(contract["trait_long_sha256"])
    expected_env_sha = str(contract["strict_spatial_chelsa_sha256"])
    if sha256(args.traits_long) != expected_trait_sha:
        raise ValueError("Trait-long input does not match the frozen v2 artifact")
    if sha256(args.environment) != expected_env_sha:
        raise ValueError("CHELSA input does not match the frozen v2 artifact")
    allowed_strict_states = {str(k): int(v) for k, v in contract["strict_states"].items()}
    minimum_orientation_observations = int(contract["minimum_image_orientation_observations_per_taxon"])

    evidence = pd.read_csv(args.evidence, dtype=str, keep_default_na=False)
    required = {
        "taxon_name", "trait_id", "literature_state", "strict_include", "source_id",
        "source_title", "source_url", "source_authority", "source_language",
        "evidence_quote", "source_scope", "notes", "review_status", "source_accessed_date",
    }
    missing = required - set(evidence.columns)
    if missing:
        raise ValueError(f"Evidence table missing columns: {sorted(missing)}")
    if not evidence["review_status"].eq("accepted_manual").all():
        raise ValueError("Every pilot evidence row must be manually accepted")
    if evidence["taxon_name"].duplicated().any():
        raise ValueError("Pilot evidence must contain at most one adjudicated orientation row per taxon")
    if not evidence["trait_id"].eq("capitulum_orientation").all():
        raise ValueError("This pilot accepts capitulum_orientation evidence only")

    traits = pd.read_csv(
        args.traits_long,
        usecols=["obs_id", "taxon_name", "endpoint_id", "value", "measurement_available"],
        low_memory=False,
    )
    traits["value"] = pd.to_numeric(traits["value"], errors="coerce")
    orientation = traits.loc[
        traits["endpoint_id"].eq(ORIENTATION_ENDPOINT)
        & as_bool(traits["measurement_available"])
        & traits["value"].notna()
    ].copy()
    image_summary = (
        orientation.groupby("taxon_name", as_index=False)
        .agg(
            n_image_orientation_observations=("value", "size"),
            image_orientation_median_degrees=("value", "median"),
        )
    )

    environment = pd.read_csv(
        args.environment,
        usecols=["obs_id", "taxon_name", "chelsa_bio12"],
        low_memory=False,
    )
    environment["chelsa_bio12"] = pd.to_numeric(environment["chelsa_bio12"], errors="coerce")
    env_summary = (
        environment.dropna(subset=["chelsa_bio12"])
        .groupby("taxon_name", as_index=False)
        .agg(
            n_environment_observations=("obs_id", "size"),
            bio12_taxon_median=("chelsa_bio12", "median"),
        )
    )

    joined = evidence.merge(image_summary, on="taxon_name", how="left", validate="one_to_one")
    joined = joined.merge(env_summary, on="taxon_name", how="left", validate="one_to_one")
    joined["strict_include_bool"] = as_bool(joined["strict_include"])
    joined["state_score"] = joined["literature_state"].map(allowed_strict_states)
    joined["minimum_image_n_pass"] = (
        pd.to_numeric(joined["n_image_orientation_observations"], errors="coerce")
        .ge(minimum_orientation_observations)
    )
    joined["strict_analysis_eligible"] = (
        joined["strict_include_bool"]
        & joined["literature_state"].isin(allowed_strict_states)
        & joined["minimum_image_n_pass"]
        & pd.to_numeric(joined["image_orientation_median_degrees"], errors="coerce").notna()
        & pd.to_numeric(joined["bio12_taxon_median"], errors="coerce").notna()
    )

    strict = joined.loc[joined["strict_analysis_eligible"]].copy()
    if len(strict) < 8:
        raise ValueError(f"Pilot requires at least 8 strict taxa, observed {len(strict)}")
    if strict["state_score"].nunique() != 2:
        raise ValueError("Pilot requires both upright and nodding taxa")

    image_validation = summarize_response(strict, "image_orientation_median_degrees")
    bio12_validation = summarize_response(strict, "bio12_taxon_median")

    report = {
        "analysis_id": str(contract["analysis_id"]),
        "status": "pilot_external_validation_not_gate_closure",
        "purpose": (
            "Targeted Supplementary check of whether independent botanical orientation descriptions "
            "agree with image-derived taxon medians and the frozen orientation-BIO12 candidate direction."
        ),
        "targeted_candidate": "larger image-referenced orientation angle with higher annual precipitation",
        "evidence_rows_total": int(len(joined)),
        "strict_taxa": int(len(strict)),
        "strict_state_counts": strict["literature_state"].value_counts().to_dict(),
        "minimum_image_orientation_observations_per_taxon": minimum_orientation_observations,
        "image_orientation_validation": image_validation,
        "bio12_directional_replication": bio12_validation,
        "interpretation": (
            "Independent literature state is directionally concordant with both the image-derived orientation median "
            "and BIO12 in this small pilot, but neither targeted comparison reaches conventional significance. "
            "This is therefore feasibility and external-consistency evidence, not confirmation of the v2 candidate."
        ),
        "excluded_or_non_strict_taxa": joined.loc[~joined["strict_analysis_eligible"], [
            "taxon_name", "literature_state", "strict_include", "notes"
        ]].to_dict("records"),
        "input_sha256": {
            "contract": sha256(args.contract),
            "evidence": sha256(args.evidence),
            "traits_long": sha256(args.traits_long),
            "strict_spatial_chelsa": sha256(args.environment),
        },
        "claim_boundary": list(contract["claim_boundary"]),
    }

    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    joined.to_csv(out / "literature_orientation_external_validation_taxa.csv", index=False, lineterminator="\n")
    strict_columns = [
        "taxon_name", "literature_state", "source_authority", "source_title", "source_url",
        "evidence_quote", "n_image_orientation_observations", "image_orientation_median_degrees",
        "n_environment_observations", "bio12_taxon_median", "strict_analysis_eligible",
    ]
    strict.loc[:, strict_columns].to_csv(
        out / "Table_S4_literature_orientation_external_validation.csv", index=False, lineterminator="\n"
    )
    (out / "literature_orientation_external_validation_report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
