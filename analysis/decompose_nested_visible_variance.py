#!/usr/bin/env python3
"""Nested visible-variance decomposition for the balanced Chapter 1 image atlas.

The analysis partitions observed sums of squares into:
- among assigned species means;
- among photographs/observations within assigned species;
- among detected heads within the same photograph.

The atlas contains one retained photograph per public observation, so photograph and
observation are the same level in this dataset. This is a descriptive decomposition
of visible image variance, not a claim about genetic or developmental variance.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

TRAITS = {
    "orientation_angle_degrees": ("orientation_status", "Orientation"),
    "corolla_lab_lightness": ("colour_status", "Lightness"),
    "corolla_lab_chroma": ("colour_status", "Chroma"),
    "corolla_hue_sin": ("colour_status", "Hue sine"),
    "corolla_hue_cos": ("colour_status", "Hue cosine"),
    "shape_aspect_ratio": ("shape_status", "Aspect ratio"),
    "shape_circularity": ("shape_status", "Circularity"),
    "shape_solidity": ("shape_status", "Solidity"),
    "shape_width_cv": ("shape_status", "Width-profile CV"),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--head-table", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--bootstrap-repeats", type=int, default=2000)
    parser.add_argument("--balanced-repeats", type=int, default=500)
    parser.add_argument("--balanced-photos-per-species", type=int, default=10)
    parser.add_argument("--seed", type=int, default=20260722)
    return parser.parse_args()


def stable_number(value: str, seed: int) -> int:
    payload = f"{seed}|{value}".encode("utf-8")
    return int(hashlib.sha256(payload).hexdigest()[:16], 16)


def nested_ss(work: pd.DataFrame, value: str) -> dict[str, float | int]:
    y = work[value].to_numpy(float)
    grand = float(np.mean(y))
    total_ss = float(np.sum((y - grand) ** 2))

    photo = work.groupby(["taxon_name", "photo_id"], observed=True)[value].agg(["size", "mean"])
    species = work.groupby("taxon_name", observed=True)[value].agg(["size", "mean"])
    species_mean = species["mean"]

    species_ss = float(np.sum(species["size"] * (species["mean"] - grand) ** 2))
    photo_species_means = photo.index.get_level_values("taxon_name").map(species_mean)
    photo_ss = float(
        np.sum(
            photo["size"].to_numpy(float)
            * (photo["mean"].to_numpy(float) - np.asarray(photo_species_means, float)) ** 2
        )
    )

    photo_mean_map = photo["mean"]
    indexed = work.set_index(["taxon_name", "photo_id"])
    residual = indexed[value] - indexed.index.map(photo_mean_map)
    head_ss = float(np.sum(np.square(residual.to_numpy(float))))

    closure = species_ss + photo_ss + head_ss
    if not np.isclose(total_ss, closure, rtol=1e-8, atol=1e-8 * max(1.0, total_ss)):
        raise RuntimeError(f"Nested SS does not close for {value}: total={total_ss}, parts={closure}")

    return {
        "n_heads": int(len(work)),
        "n_photos": int(work["photo_id"].nunique()),
        "n_taxa": int(work["taxon_name"].nunique()),
        "total_ss": total_ss,
        "among_species_ss": species_ss,
        "among_photos_within_species_ss": photo_ss,
        "among_heads_within_photo_ss": head_ss,
        "among_species_fraction": species_ss / total_ss if total_ss > 0 else np.nan,
        "among_photos_within_species_fraction": photo_ss / total_ss if total_ss > 0 else np.nan,
        "among_heads_within_photo_fraction": head_ss / total_ss if total_ss > 0 else np.nan,
        "within_assigned_species_fraction": (photo_ss + head_ss) / total_ss if total_ss > 0 else np.nan,
    }


def species_sufficient(work: pd.DataFrame, value: str) -> pd.DataFrame:
    species = (
        work.groupby("taxon_name", observed=True)[value]
        .agg(["size", "mean"])
        .rename(columns={"size": "n", "mean": "species_mean"})
    )
    photo = work.groupby(["taxon_name", "photo_id"], observed=True)[value].agg(["size", "mean"])
    species_mean_map = species["species_mean"]
    photo_species_means = photo.index.get_level_values("taxon_name").map(species_mean_map)
    photo_part = pd.Series(
        photo["size"].to_numpy(float)
        * (photo["mean"].to_numpy(float) - np.asarray(photo_species_means, float)) ** 2,
        index=photo.index,
    ).groupby(level=0).sum()
    indexed = work.set_index(["taxon_name", "photo_id"])
    residual = indexed[value] - indexed.index.map(photo["mean"])
    head_part = pd.Series(np.square(residual.to_numpy(float)), index=indexed.index).groupby(level=0).sum()
    species["photo_ss"] = photo_part.reindex(species.index, fill_value=0.0)
    species["head_ss"] = head_part.reindex(species.index, fill_value=0.0)
    return species.reset_index()


def bootstrap_fractions(work: pd.DataFrame, value: str, repeats: int, seed: int) -> dict[str, float]:
    sufficient = species_sufficient(work, value)
    n = sufficient["n"].to_numpy(float)
    means = sufficient["species_mean"].to_numpy(float)
    photo_ss = sufficient["photo_ss"].to_numpy(float)
    head_ss = sufficient["head_ss"].to_numpy(float)
    rng = np.random.default_rng(seed)
    rows = np.empty((repeats, 3), dtype=float)
    n_species = len(sufficient)
    for repeat in range(repeats):
        indices = rng.integers(0, n_species, n_species)
        sampled_n = n[indices]
        sampled_means = means[indices]
        grand = float(np.sum(sampled_n * sampled_means) / np.sum(sampled_n))
        species_part = float(np.sum(sampled_n * np.square(sampled_means - grand)))
        photo_part = float(np.sum(photo_ss[indices]))
        head_part = float(np.sum(head_ss[indices]))
        total = species_part + photo_part + head_part
        rows[repeat] = [species_part / total, photo_part / total, head_part / total]
    output: dict[str, float] = {}
    for index, name in enumerate(
        ["among_species", "among_photos_within_species", "among_heads_within_photo"]
    ):
        low, high = np.quantile(rows[:, index], [0.025, 0.975])
        output[f"{name}_ci95_low"] = float(low)
        output[f"{name}_ci95_high"] = float(high)
    within = rows[:, 1] + rows[:, 2]
    low, high = np.quantile(within, [0.025, 0.975])
    output["within_assigned_species_ci95_low"] = float(low)
    output["within_assigned_species_ci95_high"] = float(high)
    return output


def one_head_per_photo(work: pd.DataFrame, seed: int) -> pd.DataFrame:
    ranked = work.copy()
    ranked["stable_rank"] = ranked["annotation_unit_id"].astype(str).map(
        lambda value: stable_number(value, seed)
    )
    return ranked.sort_values(["taxon_name", "photo_id", "stable_rank"]).drop_duplicates(
        ["taxon_name", "photo_id"]
    )


def balanced_photo_sensitivity(
    one_head: pd.DataFrame,
    value: str,
    photos_per_species: int,
    repeats: int,
    seed: int,
) -> dict[str, float | int]:
    groups = [
        part[value].to_numpy(float)
        for _, part in one_head.groupby("taxon_name", observed=True)
        if len(part) >= photos_per_species
    ]
    if len(groups) < 3:
        return {"balanced_n_taxa": len(groups)}
    rng = np.random.default_rng(seed)
    fractions = np.empty(repeats, dtype=float)
    for repeat in range(repeats):
        selected = np.vstack(
            [rng.choice(values, size=photos_per_species, replace=False) for values in groups]
        )
        species_means = selected.mean(axis=1)
        grand = float(species_means.mean())
        species_part = float(
            photos_per_species * np.sum(np.square(species_means - grand))
        )
        within_part = float(np.sum(np.square(selected - species_means[:, None])))
        total = species_part + within_part
        fractions[repeat] = within_part / total if total > 0 else np.nan
    fractions = fractions[np.isfinite(fractions)]
    low, median, high = np.quantile(fractions, [0.025, 0.5, 0.975])
    return {
        "balanced_n_taxa": int(len(groups)),
        "balanced_photos_per_species": int(photos_per_species),
        "balanced_repeats": int(repeats),
        "balanced_within_assigned_species_median": float(median),
        "balanced_within_assigned_species_ci95_low": float(low),
        "balanced_within_assigned_species_ci95_high": float(high),
    }


def make_figure(summary: pd.DataFrame, out_dir: Path) -> None:
    order = list(TRAITS)
    plot = summary.set_index("trait").loc[order]
    y_positions = np.arange(len(plot))
    fig, axis = plt.subplots(figsize=(10, 6.8))
    left = np.zeros(len(plot))
    for column, label in [
        ("among_species_fraction", "Among assigned species means"),
        ("among_photos_within_species_fraction", "Among photographs within species"),
        ("among_heads_within_photo_fraction", "Among heads within photograph"),
    ]:
        values = plot[column].to_numpy(float)
        axis.barh(y_positions, values, left=left, label=label)
        left += values
    axis.set_yticks(y_positions, [TRAITS[trait][1] for trait in order])
    axis.invert_yaxis()
    axis.set_xlim(0, 1)
    axis.set_xlabel("Fraction of total visible image sums of squares")
    axis.set_title("Nested decomposition of visible capitulum-trait variance")
    axis.legend(loc="lower center", bbox_to_anchor=(0.5, -0.24), ncol=3)
    fig.tight_layout()
    fig.savefig(out_dir / "figure_nested_visible_variance.png", dpi=300, bbox_inches="tight")
    fig.savefig(out_dir / "figure_nested_visible_variance.svg", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(args.head_table, low_memory=False)
    required = {"annotation_unit_id", "obs_id", "photo_id", "taxon_name", *TRAITS}
    required.update(status for status, _ in TRAITS.values())
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"Head table missing columns: {missing}")
    photos_per_observation = frame.groupby("obs_id")["photo_id"].nunique()
    if int(photos_per_observation.max()) != 1:
        raise ValueError("Balanced atlas no longer has one retained photo per observation")

    rows = []
    for trait_index, (trait, (status, label)) in enumerate(TRAITS.items()):
        work = frame.loc[
            frame[status].eq("usable"),
            ["annotation_unit_id", "taxon_name", "photo_id", trait],
        ].copy()
        work[trait] = pd.to_numeric(work[trait], errors="coerce")
        work = work.dropna(subset=[trait])
        result = {"trait": trait, "label": label, **nested_ss(work, trait)}
        result.update(
            bootstrap_fractions(
                work,
                trait,
                args.bootstrap_repeats,
                args.seed + trait_index,
            )
        )
        one_head = one_head_per_photo(work, args.seed + 1000 + trait_index)
        one_result = nested_ss(one_head, trait)
        result.update(
            {
                "one_head_per_photo_n_heads": one_result["n_heads"],
                "one_head_per_photo_n_photos": one_result["n_photos"],
                "one_head_per_photo_n_taxa": one_result["n_taxa"],
                "one_head_per_photo_within_assigned_species_fraction": one_result[
                    "within_assigned_species_fraction"
                ],
                "one_head_per_photo_among_species_fraction": one_result[
                    "among_species_fraction"
                ],
            }
        )
        result.update(
            balanced_photo_sensitivity(
                one_head,
                trait,
                args.balanced_photos_per_species,
                args.balanced_repeats,
                args.seed + 2000 + trait_index,
            )
        )
        rows.append(result)

    summary = pd.DataFrame(rows)
    summary.to_csv(out_dir / "nested_visible_variance_summary.csv", index=False)
    make_figure(summary, out_dir)
    metadata = {
        "analysis": "nested descriptive sums-of-squares decomposition",
        "levels": [
            "among assigned species means",
            "among photographs/observations within assigned species",
            "among detected heads within photograph",
        ],
        "one_photo_per_observation": True,
        "bootstrap_unit": "assigned species",
        "bootstrap_repeats": args.bootstrap_repeats,
        "balanced_sensitivity": {
            "one stable-hash head per photograph": True,
            "photos_per_species": args.balanced_photos_per_species,
            "repeats": args.balanced_repeats,
        },
        "interpretation_boundary": (
            "Visible image variance includes biological heterogeneity, viewpoint, "
            "illumination, segmentation and other measurement components; it is not "
            "a genetic variance decomposition."
        ),
    }
    (out_dir / "nested_visible_variance_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    print(
        summary[
            [
                "trait",
                "n_heads",
                "n_photos",
                "n_taxa",
                "among_species_fraction",
                "among_photos_within_species_fraction",
                "among_heads_within_photo_fraction",
                "one_head_per_photo_within_assigned_species_fraction",
                "balanced_n_taxa",
                "balanced_within_assigned_species_median",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
