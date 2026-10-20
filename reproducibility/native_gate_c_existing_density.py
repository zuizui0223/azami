"""Read-only independent PR92 density-distance follow-up verification.

This validates ACTUAL historical non-native (46276/259 source) CSVs at their
original Git blob identities. It is not the new native 313350 source, nor a
statistical test of independent species pairs. No images, raw IDs, refits,
raster acquisition, KDE execution or new W2 estimates are performed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HISTORICAL_COMMIT = "6061befaf56ee545bb3b2b4a80ed220d2def1190"
GIT_BLOBS = {
    "pairwise": "5d996a95049840cdd97c18f4ca3760ad8cc950ba",
    "same_taxon": "959ce8ffb8314ef7461b9bc6486561240fd66802",
    "figures_summary": "17315a9f4aa8cf68632e555e016cf06a99b6162c",
    "overlap_summary": "a9f64ccde7817e9b57e807146d13428ca4b91c83",
    "bandwidth": "af68d2affd5fa7917420bcef4c74779a6ef72551",
}
CORE_GROUPS = {"orientation": 33, "colour": 33, "outline": 34}
ALL_PAIR_COUNTS = {"orientation": 1431, "colour": 1431, "outline": 1431,
                   "involucre": 15, "whole": 10}
METRICS = ("density_overlap", "wasserstein_w2", "centred_share_w2_squared")


def validate_blob(path: Path, blob_sha: str) -> None:
    data = path.read_bytes()
    envelope = b"blob " + str(len(data)).encode("ascii") + b"\\0" + data
    if hashlib.sha1(envelope).hexdigest() != blob_sha:
        raise ValueError("Historical PR92 Git blob hash changed or wrong input")


def median(seq) -> float:
    return float(np.median(np.asarray(seq, dtype=float)))


def analyze_historical_density(
    pairwise: pd.DataFrame, same_taxon: pd.DataFrame,
    figures_summary: list, overlap_summary: list, bandwidth: pd.DataFrame
) -> dict:
    """Recompute documented summaries from original taxon-level pair tables."""
    for cols, frame in [
        (("group", "taxon_a", "taxon_b", *METRICS,
          "centroid_distance", "centred_w2"), pairwise),
        (("group", "taxon_name", "same_taxon_w2_p95",
          "centred_same_taxon_w2_p95"), same_taxon),
        (("group", "taxon_a", "taxon_b", "factor", "density_overlap"), bandwidth),
    ]:
        missing = sorted(set(cols) - set(frame))
        if missing:
            raise ValueError("Missing historic source fields: " + ", ".join(missing))
    if len(pairwise) != sum(ALL_PAIR_COUNTS.values()):
        raise ValueError("Historical pairwise row count changed")
    if pairwise.groupby("group").size().to_dict() != ALL_PAIR_COUNTS:
        raise ValueError("Original historical group census changed")
    if set(same_taxon.group) != set(CORE_GROUPS):
        raise ValueError("Unexpected within-source groups")
    if same_taxon.groupby("group").size().to_dict() != CORE_GROUPS:
        raise ValueError("Wrong within-taxon baseline census")
    if same_taxon.duplicated(["group", "taxon_name"]).any():
        raise ValueError("Duplicated historical same-taxon baseline")
    if pairwise.duplicated(["group", "taxon_a", "taxon_b"]).any():
        raise ValueError("Duplicated historical pairwise edge")
    if len(bandwidth) != 900 or bandwidth.duplicated(
        ["group", "taxon_a", "taxon_b", "factor"]
    ).any():
        raise ValueError("Bandwidth test pair identity/replication changed")
    for col in ("same_taxon_w2_p95", "centred_same_taxon_w2_p95"):
        if not np.isfinite(pd.to_numeric(same_taxon[col], errors="coerce")).all():
            raise ValueError("Nonfinite within-taxon reference")
    if set(round(float(k), 8) for k in bandwidth.factor) != {0.75,1.0,1.25}:
        raise ValueError("Wrong bandwidth set")
    expected_figures = {x["group"]:x for x in figures_summary}
    expected_overlap = {x["group"]:x for x in overlap_summary}
    if set(expected_figures) != set(CORE_GROUPS) or set(expected_overlap) != set(ALL_PAIR_COUNTS):
        raise ValueError("Source reports changed their registered groups")
    if not np.isfinite(pairwise[[*METRICS, "centroid_distance", "centred_w2"]].to_numpy(float)).all():
        raise ValueError("Nonfinite PR92 distance/overlap")
    if not pairwise.density_overlap.between(0,1).all():
        raise ValueError("Invalid KDE probability-density overlap")

    b = same_taxon.set_index(["group","taxon_name"])
    report = {}
    for group, taxa_n in CORE_GROUPS.items():
        p = pairwise.loc[pairwise.group.eq(group)].copy()
        left = b.loc[list(zip(p.group, p.taxon_a))]
        right = b.loc[list(zip(p.group, p.taxon_b))]
        # Only original taxa whose 30 disjoint 50+50 same-species W2 draws exist.
        mask = np.isfinite(left.same_taxon_w2_p95.to_numpy(float)) & np.isfinite(
            right.same_taxon_w2_p95.to_numpy(float))
        if not mask.all():
            raise ValueError("Same-taxon reference joined unexpectedly")
        # The all-pair CSV has 54 eligible taxa, so restrict to 33/33/34 baseline taxa.
        # The pairs without baseline taxa are NOT included in the reference comparison.
        # Build the appropriate all-pairs subset *before* baseline indexing.
        report[group] = {
            "historical_between_species_pairs":int(len(p)),
            "full_median_overlap":median(p.density_overlap),
            "full_median_w2":median(p.wasserstein_w2),
            "full_median_centred_w2_squared_share":median(p.centred_share_w2_squared),
            "same_taxon_baseline_taxa":taxa_n,
        }
        prior=expected_overlap[group]
        if (len(p),prior["pairs"]) != (ALL_PAIR_COUNTS[group],ALL_PAIR_COUNTS[group]):
            raise ValueError("Original pair count summary changed")
        if abs(report[group]["full_median_overlap"] - prior["median_pair_density_overlap"]) > 1e-12:
            raise ValueError("Original KDE overlap median not reproduced")
        if abs(report[group]["full_median_w2"] - prior["median_pair_w2"]) > 1e-12:
            raise ValueError("Original W2 median not reproduced")
        if abs(report[group]["full_median_centred_w2_squared_share"] -
               prior["median_pair_centred_share"]) > 1e-12:
            raise ValueError("Original centred W2 summary not reproduced")

        eligible = set(b.loc[group].index)
        u = p.loc[p.taxon_a.isin(eligible) & p.taxon_b.isin(eligible)].copy()
        if len(u) != taxa_n * (taxa_n - 1) // 2:
            raise ValueError("Historical same-taxon comparison does not close all pairs")
        lb=b.loc[list(zip(u.group, u.taxon_a))]
        rb=b.loc[list(zip(u.group, u.taxon_b))]
        w = np.maximum(lb.same_taxon_w2_p95.to_numpy(float),
                       rb.same_taxon_w2_p95.to_numpy(float))
        c = np.maximum(lb.centred_same_taxon_w2_p95.to_numpy(float),
                       rb.centred_same_taxon_w2_p95.to_numpy(float))
        above_full=int((u.wasserstein_w2.to_numpy(float)>w).sum())
        above_centred=int((u.centred_w2.to_numpy(float)>c).sum())
        prior_fig=expected_figures[group]
        if (len(u),above_full,above_centred) != (
            prior_fig["pairs"],prior_fig["between_above_within_p95"],
            prior_fig["centred_above_within_p95"]
        ):
            raise ValueError("Historical same-taxon baseline flags/counts changed")
        report[group].update({
            "pairs_with_same_taxon_reference":int(len(u)),
            "between_w2_above_same_taxon_p95":above_full,
            "centred_w2_above_same_taxon_p95":above_centred,
            "full_w2_exceedance_fraction":above_full/len(u),
            "centred_w2_exceedance_fraction":above_centred/len(u),
        })

        bw = bandwidth.loc[bandwidth.group.eq(group)].pivot(
            index=["taxon_a","taxon_b"],columns="factor",values="density_overlap"
        )
        if bw.shape != (100,3) or set(bw.columns) != {0.75,1.0,1.25}:
            raise ValueError("Historical bandwidth sample changed")
        for h, suffix in ((.75,"075"),(1.25,"125")):
            rho=float(spearmanr(bw[1.0].to_numpy(float),bw[h].to_numpy(float)).statistic)
            if abs(rho-prior_fig["bandwidth_rank_rho_"+suffix])>1e-12:
                raise ValueError("Archived rank bandwidth sensitivity changed")
        for h in (.75,1.0,1.25):
            observed=median(bw[h])
            prior_bandwidth=prior_fig["median_overlap_by_bandwidth"][str(h)]
            if abs(observed-prior_bandwidth)>1e-12:
                raise ValueError("Archived bandwidth overlap median changed")
        report[group]["bandwidth_check"]={
            "pairs":100,
            "factors":[.75,1.0,1.25],
            "median_overlap":[median(bw[h]) for h in (.75,1.0,1.25)],
        }
    return {
        "status":"REPRODUCED_PR92_HISTORICAL_NOT_CURRENT_NATIVE",
        "historical_source_commit":HISTORICAL_COMMIT,
        "source_git_blob_sha1":GIT_BLOBS,
        "previous_source_native_filter":False,
        "previous_source_spatially_thinned_observations":46276,
        "same_species_null":"30 disjoint 50+50 photo-observation subsamples; descriptive p95 reference",
        "modules":report,
        "total_original_between_taxon_pairs":len(pairwise),
        "density_and_distance_inference":[
            "W2 centre and shape parts from medians need not obey replicate-level W2 squared identity.",
            "Taxon-pair fractions are descriptive dependent pairs, not p values.",
            "Image measurement and viewpoint error contribute to within-species W2.",
            "Orientation is 1D, colour and outline 3D: raw overlap rankings are not dimension-matched.",
            "Bandwith stability is rank-based only; absolute overlap changes with smoothing.",
            "Not the 313350-observation native frame; no new native KDE/W2 executed.",
        ],
    }


def main() -> None:
    p=argparse.ArgumentParser(description=__doc__)
    for name in GIT_BLOBS:
        p.add_argument("--"+name.replace("_","-"),type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    paths={k:getattr(a,k) for k in GIT_BLOBS}
    for key,path in paths.items():
        validate_blob(path,GIT_BLOBS[key])
    report=analyze_historical_density(
        pd.read_csv(paths["pairwise"]),
        pd.read_csv(paths["same_taxon"]),
        json.loads(paths["figures_summary"].read_text(encoding="utf-8")),
        json.loads(paths["overlap_summary"].read_text(encoding="utf-8")),
        pd.read_csv(paths["bandwidth"]),
    )
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False)+"\\n",
                      encoding="utf-8")
    print(json.dumps({"status":report["status"],"modules":report["modules"]}))


if __name__=="__main__":
    main()
