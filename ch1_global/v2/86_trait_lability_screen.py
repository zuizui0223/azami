#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

TRAIT_MODULES = {
    "orientation_angle_degrees_median": "orientation",
    "corolla_lab_lightness_median": "colour",
    "corolla_lab_chroma_median": "colour",
    "corolla_hue_sin_median": "colour",
    "corolla_hue_cos_median": "colour",
    "shape_aspect_ratio_median": "shape",
    "shape_circularity_median": "shape",
    "shape_solidity_median": "shape",
    "shape_width_cv_median": "shape",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--environment", required=True)
    p.add_argument("--spde-fixed", required=True)
    p.add_argument("--spde-stability", required=True)
    p.add_argument("--spde-hyper", required=False)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--min-species-n", type=int, default=5)
    p.add_argument("--bootstrap-repeats", type=int, default=200)
    p.add_argument("--random-state", type=int, default=20260714)
    return p.parse_args()


def variance_partition(df: pd.DataFrame, trait: str, min_species_n: int) -> dict:
    z = df[["taxon_name", trait]].copy()
    z[trait] = pd.to_numeric(z[trait], errors="coerce")
    z = z[np.isfinite(z[trait])]
    counts = z.groupby("taxon_name")[trait].size()
    keep = counts[counts >= min_species_n].index
    z = z[z.taxon_name.isin(keep)]
    if len(z) < 100 or z.taxon_name.nunique() < 3:
        return {"n_observations": len(z), "n_species": z.taxon_name.nunique()}

    grand = z[trait].mean()
    total_ss = float(((z[trait] - grand) ** 2).sum())
    means = z.groupby("taxon_name")[trait].transform("mean")
    within_ss = float(((z[trait] - means) ** 2).sum())
    between_ss = max(total_ss - within_ss, 0.0)
    within_fraction = within_ss / total_ss if total_ss > 0 else np.nan
    between_fraction = between_ss / total_ss if total_ss > 0 else np.nan

    species_stats = z.groupby("taxon_name")[trait].agg(["median", "std", "size"])
    return {
        "n_observations": int(len(z)),
        "n_species": int(z.taxon_name.nunique()),
        "total_sd": float(z[trait].std(ddof=1)),
        "median_within_species_sd": float(species_stats["std"].median()),
        "sd_species_medians": float(species_stats["median"].std(ddof=1)),
        "within_species_variance_fraction": within_fraction,
        "between_species_variance_fraction": between_fraction,
    }


def bootstrap_within_fraction(
    df: pd.DataFrame,
    trait: str,
    min_species_n: int,
    repeats: int,
    rng: np.random.Generator,
) -> tuple[float, float]:
    z = df[["taxon_name", trait]].copy()
    z[trait] = pd.to_numeric(z[trait], errors="coerce")
    z = z[np.isfinite(z[trait])]
    counts = z.groupby("taxon_name")[trait].size()
    keep = counts[counts >= min_species_n].index.tolist()
    if len(keep) < 3:
        return np.nan, np.nan
    groups = {s: z.loc[z.taxon_name.eq(s), trait].to_numpy(float) for s in keep}
    vals = []
    for _ in range(repeats):
        sampled = rng.choice(keep, size=len(keep), replace=True)
        chunks = []
        labels = []
        for j, s in enumerate(sampled):
            v = groups[s]
            chunks.append(v)
            labels.extend([j] * len(v))
        y = np.concatenate(chunks)
        lab = np.asarray(labels)
        grand = y.mean()
        total_ss = ((y - grand) ** 2).sum()
        if total_ss <= 0:
            continue
        means = pd.Series(y).groupby(lab).transform("mean").to_numpy()
        vals.append(float(((y - means) ** 2).sum() / total_ss))
    if not vals:
        return np.nan, np.nan
    return float(np.quantile(vals, 0.025)), float(np.quantile(vals, 0.975))


def percentile_rank(s: pd.Series, ascending: bool = True) -> pd.Series:
    return s.rank(pct=True, ascending=ascending, method="average")


def main() -> None:
    a = parse_args()
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.random_state)

    df = pd.read_csv(a.environment, low_memory=False)
    fixed = pd.read_csv(a.spde_fixed)
    stability = pd.read_csv(a.spde_stability)
    hyper = pd.read_csv(a.spde_hyper) if a.spde_hyper and Path(a.spde_hyper).exists() else pd.DataFrame()

    rows = []
    for trait, module in TRAIT_MODULES.items():
        if trait not in df.columns:
            continue
        row = {"trait": trait, "module": module}
        row.update(variance_partition(df, trait, a.min_species_n))
        lo, hi = bootstrap_within_fraction(
            df, trait, a.min_species_n, a.bootstrap_repeats, rng
        )
        row["within_fraction_boot_low"] = lo
        row["within_fraction_boot_high"] = hi

        fx = fixed[(fixed["trait"] == trait) & (fixed["term"] != "Intercept")].copy()
        if len(fx):
            row["median_abs_environment_beta"] = float(fx["mean"].abs().median())
            row["max_abs_environment_beta"] = float(fx["mean"].abs().max())
            cred = fx.get("credible_95_excludes_zero", pd.Series(False, index=fx.index)).astype(bool)
            row["credible_effect_fraction"] = float(cred.mean())
            if "q_bh_posterior_tail_global" in fx:
                row["bh_supported_effect_fraction"] = float(
                    (pd.to_numeric(fx["q_bh_posterior_tail_global"], errors="coerce") < 0.05).mean()
                )
            else:
                row["bh_supported_effect_fraction"] = np.nan
        else:
            row.update(
                median_abs_environment_beta=np.nan,
                max_abs_environment_beta=np.nan,
                credible_effect_fraction=np.nan,
                bh_supported_effect_fraction=np.nan,
            )

        st = stability[stability["trait"] == trait].copy()
        if len(st):
            row["sign_consistent_fraction"] = float(st["sign_consistent"].astype(bool).mean())
            row["mean_groups_per_effect"] = float(pd.to_numeric(st["n_groups"], errors="coerce").mean())
        else:
            row["sign_consistent_fraction"] = np.nan
            row["mean_groups_per_effect"] = np.nan

        if len(hyper):
            hp = hyper[hyper["trait"] == trait].copy()
            range_rows = hp[hp["parameter"].astype(str).str.contains("Range", case=False, na=False)]
            row["median_spatial_range_km"] = (
                float(pd.to_numeric(range_rows["mean"], errors="coerce").median())
                if len(range_rows)
                else np.nan
            )
        rows.append(row)

    res = pd.DataFrame(rows)
    components = {
        "within_rank": "within_species_variance_fraction",
        "env_beta_rank": "median_abs_environment_beta",
        "bh_support_rank": "bh_supported_effect_fraction",
        "sign_stability_rank": "sign_consistent_fraction",
    }
    for new, old in components.items():
        res[new] = percentile_rank(pd.to_numeric(res[old], errors="coerce"), ascending=True)

    # Candidate lability is high when within-species variability and repeatable environmental
    # responsiveness are high. This is not an evolutionary-rate estimate because no phylogeny
    # is included.
    res["candidate_lability_score"] = res[list(components)].mean(axis=1, skipna=True)
    res["candidate_conservatism_score"] = 1.0 - res["candidate_lability_score"]
    res["lability_rank"] = res["candidate_lability_score"].rank(ascending=False, method="min").astype(int)
    q1 = res["candidate_lability_score"].quantile(1 / 3)
    q2 = res["candidate_lability_score"].quantile(2 / 3)
    res["screening_class"] = np.where(
        res["candidate_lability_score"] >= q2,
        "candidate_labile",
        np.where(res["candidate_lability_score"] <= q1, "candidate_conservative", "intermediate"),
    )
    res = res.sort_values("candidate_lability_score", ascending=False)
    res.to_csv(out / "trait_lability_conservatism_screen.csv", index=False)

    module = (
        res.groupby("module", as_index=False)
        .agg(
            n_traits=("trait", "size"),
            mean_candidate_lability=("candidate_lability_score", "mean"),
            median_within_species_variance_fraction=("within_species_variance_fraction", "median"),
            median_abs_environment_beta=("median_abs_environment_beta", "median"),
            mean_bh_supported_effect_fraction=("bh_supported_effect_fraction", "mean"),
            mean_sign_consistent_fraction=("sign_consistent_fraction", "mean"),
        )
        .sort_values("mean_candidate_lability", ascending=False)
    )
    module.to_csv(out / "module_lability_conservatism_summary.csv", index=False)

    metadata = {
        "interpretation": "screening index for candidate phenotypic lability/conservatism, not a phylogenetic evolutionary-rate estimate",
        "components": list(components.values()),
        "bootstrap_repeats": a.bootstrap_repeats,
        "min_species_n": a.min_species_n,
        "random_state": a.random_state,
        "traits": len(res),
    }
    (out / "lability_screen_metadata.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(res[["trait", "module", "candidate_lability_score", "screening_class"]].to_string(index=False))


if __name__ == "__main__":
    main()
