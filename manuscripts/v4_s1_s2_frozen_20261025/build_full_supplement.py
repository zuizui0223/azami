#!/usr/bin/env python3
"""Rebuild V4 S1-S2 SI tables and 28-pair heatmap from the frozen Actions artifact.

Usage:
 python manuscripts/v4_s1_s2_frozen_20261025/build_full_supplement.py \
   --artifact /path/to/v4-s1-s2-extreme-lego-40062126006.zip \
   --out /path/to/v4-s1-s2-si

Requirements: Python 3, numpy, scipy, statsmodels, matplotlib.
Source: Actions 40062126006 / artifact 12547876704; NO original IDs/coordinates/photos.
All resulting values are observational; do not interpret RV as genetic covariance.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import zipfile

import numpy as np
from statsmodels.stats.multitest import multipletests

EXPECTED_SHA = "6e081dc081bb8ee2311a05ac5b7fcddd3f6309be2fa909e87d432789f4e0e205"
CONSTRUCTS = ("angle", "lightness", "chroma", "hue", "aspect", "circularity", "solidity", "width_profile")
PRIMARY_FAMILIES = ("S1/four_bio18", "S2/four_bio18")


def csv_write(path: Path, rows: list[dict], headers: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def get_archive(archive: Path) -> dict:
    digest = hashlib.sha256(archive.read_bytes()).hexdigest()
    if digest != EXPECTED_SHA:
        raise ValueError(f"Frozen original artifact checksum mismatch: {digest}")
    with zipfile.ZipFile(archive) as z:
        matches = [x for x in z.namelist() if x.endswith("V4_S1_S2_four_vs_vif5_vif10_alltraits.json")]
        if len(matches) != 1:
            raise ValueError("Expected exactly one original full four-model results JSON")
        return json.loads(z.read(matches[0]))


def full_tests(data: dict, out: Path) -> None:
    records = []
    for family, expected_n, expected_holm in zip(PRIMARY_FAMILIES, (80, 128), (1, 3)):
        original = data["multiplicity"][family]["rows"]
        if len(original) != expected_n:
            raise ValueError("Number of biological family tests changed")
        ps = np.asarray([float(x["p_original"]) for x in original])
        holm = multipletests(ps, method="holm")[1]
        bh = multipletests(ps, method="fdr_bh")[1]
        if sum(holm < 0.05) != expected_holm:
            raise ValueError(f"Frozen count failed for {family}")
        if max(abs(x["p_Holm"] - y) for x, y in zip(original, holm)) > 1e-10:
            raise ValueError("Original Holm results do not match independent recomputation")
        if max(abs(x["p_BH"] - y) for x, y in zip(original, bh)) > 1e-10:
            raise ValueError("Original BH results do not match independent recomputation")
        for i, x in enumerate(original, 1):
            low, high = (x.get("ci95_wald") or [None, None])
            records.append(
                {
                    "family": family,
                    "test_number": i,
                    "group": x.get("group"),
                    "coordinate": x.get("coordinate"),
                    "environment": x.get("environment"),
                    "scale": x.get("scale"),
                    "test_kind": x.get("test_kind"),
                    "n": x.get("n"),
                    "beta_std": x.get("beta_std"),
                    "ci95_low": low,
                    "ci95_high": high,
                    "p_original": x["p_original"],
                    "p_BH": x["p_BH"],
                    "p_Holm": x["p_Holm"],
                    "Holm_pass": bool(x["p_Holm"] < 0.05),
                    "trait_cohort_vif_pass": x.get("trait_cohort_vif_pass"),
                }
            )
    assert len(records) == 208
    csv_write(out / "Table_S1_full_208_BIO18_tests.csv", records, list(records[0]))


def full_pairs(data: dict, out: Path) -> tuple[list[dict], list[dict]]:
    cohort = data["S1"]["joint_8_constructs_28_pairs"]
    if (cohort["n"], cohort["taxa"]) != (156683, 174):
        raise ValueError("The 8-construct full-photo cohort does not match the freeze")
    raw = cohort["raw_pairs"]
    adj = cohort["models"]["four_bio18"]["pairs"]
    if len(raw) != 28 or len(adj) != 28:
        raise ValueError("Expected all 28 unordered construct pairs")
    if sum(x["among_rv"] > x["within_rv"] for x in raw) != 21:
        raise ValueError("Raw between-over-within count differs")
    if sum(x["among_rv"] > x["within_rv"] for x in adj) != 24:
        raise ValueError("Adjusted between-over-within count differs")
    records = []
    for a, b in zip(raw, adj):
        if (a["left"], a["right"]) != (b["left"], b["right"]):
            raise ValueError("Pair ordering changed")
        record = {
            "left": a["left"], "right": a["right"],
            "left_module": a["left_module"], "right_module": a["right_module"],
            "n_photos": 156683, "n_taxa": 174,
        }
        for prefix, obj in (("raw", a), ("adjusted", b)):
            for key in ("within_rv", "among_rv", "within_signed_r",
                        "among_signed_r", "among_minus_within"):
                record[f"{prefix}_{key}"] = obj.get(key)
        records.append(record)
    csv_write(out / "Table_S2_all_28_RV_pairs.csv", records, list(records[0]))
    for level, same_expect, cross_expect in (
        ("within", 0.08524521946283355, 0.001419360219707819),
        ("among", 0.19745500497587984, 0.015453102845313172),
    ):
        same = [p[f"{level}_rv"] for p in adj if p["left_module"] == p["right_module"]]
        cross = [p[f"{level}_rv"] for p in adj if p["left_module"] != p["right_module"]]
        if len(same) != 9 or len(cross) != 19:
            raise ValueError("Module pair numbers changed")
        if not np.isclose(np.median(same), same_expect, atol=1e-10):
            raise ValueError("Same-domain median does not match frozen result")
        if not np.isclose(np.median(cross), cross_expect, atol=1e-10):
            raise ValueError("Cross-domain median does not match frozen result")
    return raw, adj


def heatmap(raw: list[dict], adj: list[dict], out: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    matrices = []
    for pairs in (raw, adj):
        for level in ("within", "among"):
            matrix = np.full((8, 8), np.nan)
            for item in pairs:
                i = CONSTRUCTS.index(item["left"])
                j = CONSTRUCTS.index(item["right"])
                matrix[i, j] = matrix[j, i] = float(item[f"{level}_rv"])
            matrices.append(matrix)
    upper = max(0.35, np.quantile(np.concatenate(
        [m[np.isfinite(m)] for m in matrices]), 0.99))
    fig, axes = plt.subplots(2, 2, figsize=(13, 11), constrained_layout=True)
    fig.suptitle("V4 S1: RV geometry of eight photographic head traits", fontsize=14)
    for ax, matrix, heading in zip(
        axes.flat, matrices,
        ("Raw · within taxa", "Raw · among taxa",
         "Four-climate-adjusted · within taxa", "Four-climate-adjusted · among taxa"),
    ):
        img = ax.imshow(np.ma.masked_invalid(matrix), cmap="viridis", vmin=0, vmax=upper)
        ax.set_xticks(range(8), CONSTRUCTS, rotation=45, ha="right", fontsize=8)
        ax.set_yticks(range(8), CONSTRUCTS, fontsize=8)
        for k in (0.5, 3.5):
            ax.axvline(k, lw=1, color="white")
            ax.axhline(k, lw=1, color="white")
        ax.set_title(heading)
    cb = fig.colorbar(img, ax=axes.ravel().tolist(), shrink=0.7)
    cb.set_label("RV magnitude (common scale; diagonal undefined)")
    fig.text(.015, .005, "156,683 photos / 174 taxa; observational photographic RV, not genetic covariance.", fontsize=8)
    fig.savefig(out / "Figure_S1_S1_all_28_RV.png", dpi=200, bbox_inches="tight")
    fig.savefig(out / "Figure_S1_S1_all_28_RV.svg", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    data = get_archive(args.artifact)
    full_tests(data, args.out)
    raw, adj = full_pairs(data, args.out)
    heatmap(raw, adj, args.out)
    (args.out / "SOURCE_SHA256.txt").write_text(
        f"Source artifact 12547876704 sha256: {EXPECTED_SHA}\n"
        "S1 80 + S2 128 tests; Holm/BH independently recalculated and matched.\n"
        "S1 28 raw/adjusted RV pairs independently checked.\n"
        "No private source IDs or locations exported.\n", encoding="utf-8"
    )
    print("Verified, generated two complete aggregate tables and 28-pair heatmap:", args.out)


if __name__ == "__main__":
    main()
