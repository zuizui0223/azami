#!/usr/bin/env python3
"""Build the frozen GEB Chapter 1 figure, table and manuscript release.

The builder consumes only the two immutable artifact trees frozen by the Azami
continuous-trait and multilevel workflows. It performs presentation-only
transformations plus one explicitly matched descriptive recomputation: the
nine-primary-endpoint PCA is rebuilt with the same complete-case standardised
SVD used by ``build_geb_v2_trait_geography.py``.

No endpoint definitions, cohorts, multiplicity families, evidence grades or
inferential decisions are changed here.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


MODULE_ORDER = [
    "orientation",
    "colour",
    "shape",
    "involucre_architecture",
    "armature",
]
MODULE_LABEL = {
    "orientation": "Orientation",
    "colour": "Visible colour",
    "shape": "Outline shape",
    "involucre_architecture": "Involucre architecture",
    "armature": "Armature",
}
# Okabe–Ito-derived palette. Redundant labels/markers are retained throughout.
MODULE_COLOUR = {
    "orientation": "#0072B2",
    "colour": "#CC79A7",
    "shape": "#009E73",
    "involucre_architecture": "#E69F00",
    "armature": "#D55E00",
}
ENDPOINT_LABEL = {
    "orientation_image_vertical_angle": "Orientation",
    "corolla_lab_lightness": "Lightness",
    "corolla_lab_chroma": "Chroma",
    "corolla_hue_sin": "Hue sine",
    "corolla_hue_cos": "Hue cosine",
    "corolla_hue": "Hue",
    "capitulum_outline_aspect_ratio": "Aspect ratio",
    "capitulum_outline_circularity": "Circularity",
    "capitulum_outline_solidity": "Solidity",
    "capitulum_width_profile_cv": "Width-profile CV",
    "involucre_length_width_ratio": "Involucre L:W",
    "involucre_apical_taper_ratio": "Apical taper",
    "involucre_basal_taper_ratio": "Basal taper",
    "bract_projection_roughness": "Projection roughness",
    "bract_projection_p95": "Projection P95",
    "bract_projection_maximum": "Projection maximum",
    "bract_spread_fraction": "Spread fraction",
    "bract_projection_peak_density": "Peak density",
    "bract_projection_asymmetry": "Asymmetry",
}
ENDPOINT_ABBREV = {
    "orientation_image_vertical_angle": "Orient.",
    "corolla_lab_lightness": "L*",
    "corolla_lab_chroma": "Chroma",
    "corolla_hue_sin": "Hue sin",
    "corolla_hue_cos": "Hue cos",
    "corolla_hue": "Hue",
    "capitulum_outline_aspect_ratio": "Aspect",
    "capitulum_outline_circularity": "Circular.",
    "capitulum_outline_solidity": "Solidity",
    "capitulum_width_profile_cv": "Width CV",
    "involucre_length_width_ratio": "Inv. L:W",
    "involucre_apical_taper_ratio": "Apical",
    "involucre_basal_taper_ratio": "Basal",
    "bract_projection_roughness": "Roughness",
    "bract_projection_p95": "Proj. P95",
    "bract_projection_maximum": "Proj. max",
    "bract_spread_fraction": "Spread",
    "bract_projection_peak_density": "Peak dens.",
    "bract_projection_asymmetry": "Asymmetry",
}
PREDICTOR_ORDER = ["chelsa_bio01", "chelsa_bio04", "chelsa_bio12", "chelsa_bio15"]
PREDICTOR_LABEL = {
    "chelsa_bio01": "BIO1\nMean temp.",
    "chelsa_bio04": "BIO4\nTemp. seasonality",
    "chelsa_bio12": "BIO12\nAnnual precip.",
    "chelsa_bio15": "BIO15\nPrecip. seasonality",
}
BLOCK_ORDER = [
    "core_thermal",
    "core_precipitation",
    "radiative_atmospheric_drying",
    "mechanical_exposure",
    "growing_season_water_input",
    "climatic_productivity",
]
BLOCK_LABEL = {
    "core_thermal": "Core thermal",
    "core_precipitation": "Core precipitation",
    "radiative_atmospheric_drying": "Radiation + VPD",
    "mechanical_exposure": "Wind",
    "growing_season_water_input": "Growing-season precipitation",
    "climatic_productivity": "Potential NPP",
    "all_process_extension": "All process variables",
}
EVIDENCE_GRADE_LABEL = {
    "A_robust": "A — robust",
    "B_supported_range_sensitive": "B — range-sensitive",
    "B_supported_colour_uncalibrated": "B — colour uncalibrated",
    "B_spatial_context": "B — spatial context",
    "B_among_taxon_pattern": "B — among-taxon pattern",
    "C_exploratory_quality_robust": "C — quality-robust candidate",
    "C_exploratory_image_sensitive": "C — image-sensitive candidate",
    "D_historical_context_only": "D — historical context only",
}
EVIDENCE_GRADE_COLOUR = {
    "A_robust": "#0072B2",
    "B_supported_range_sensitive": "#56B4E9",
    "B_supported_colour_uncalibrated": "#CC79A7",
    "B_spatial_context": "#009E73",
    "B_among_taxon_pattern": "#F0E442",
    "C_exploratory_quality_robust": "#E69F00",
    "C_exploratory_image_sensitive": "#D55E00",
    "D_historical_context_only": "#777777",
}


@dataclass(frozen=True)
class PCAResult:
    scores: pd.DataFrame
    loadings: pd.DataFrame
    explained: np.ndarray


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--continuous-root", required=True, type=Path)
    p.add_argument("--cross-root", required=True, type=Path)
    p.add_argument("--manuscript-dir", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--release-date", default="2026-08-27")
    p.add_argument("--continuous-run", type=int, default=32975451732)
    p.add_argument("--continuous-artifact", type=int, default=9612943217)
    p.add_argument(
        "--continuous-digest",
        default="sha256:101e996b638996a0c5ae79d358bf51293c3585f0e84c4a961b91dcbedf96211e",
    )
    p.add_argument("--cross-run", type=int, default=33035785120)
    p.add_argument("--cross-artifact", type=int, default=9632715852)
    p.add_argument(
        "--cross-digest",
        default="sha256:51e7a26b5bd09e030b67b9342586699abaaf46e630f45b6bb4ee7bfc9152ced6",
    )
    return p.parse_args()


def setup_matplotlib() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.labelsize": 9.5,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "axes.linewidth": 0.8,
            "savefig.bbox": "tight",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def require_file(path: Path) -> Path:
    if not path.is_file():
        raise FileNotFoundError(f"Required frozen source is missing: {path}")
    return path


def read_csv(path: Path, **kwargs: object) -> pd.DataFrame:
    return pd.read_csv(require_file(path), **kwargs)


def as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def stable_pca(wide: pd.DataFrame) -> PCAResult:
    """Match the frozen complete-case, population-SD SVD implementation."""
    complete = wide.dropna().copy()
    if complete.shape[0] < 2 or complete.shape[1] < 2:
        raise ValueError("PCA requires at least two complete taxa and two endpoints")
    x = complete.to_numpy(float)
    means = x.mean(axis=0)
    sds = x.std(axis=0, ddof=0)
    keep = np.isfinite(sds) & (sds > 0)
    if keep.sum() < 2:
        raise ValueError("PCA has fewer than two variable endpoints")
    columns = np.asarray(complete.columns)[keep]
    z = (x[:, keep] - means[keep]) / sds[keep]
    u, s, vt = np.linalg.svd(z, full_matrices=False)
    explained = (s**2) / np.sum(s**2)
    n_pc = min(3, vt.shape[0])
    scores = pd.DataFrame(
        u[:, :n_pc] * s[:n_pc],
        index=complete.index,
        columns=[f"PC{i+1}" for i in range(n_pc)],
    )
    scores.index.name = complete.index.name or "taxon_name"
    scores = scores.reset_index()
    loadings = pd.DataFrame({"endpoint_id": columns})
    for i in range(n_pc):
        loadings[f"PC{i+1}_loading"] = vt[i, :]
    return PCAResult(scores=scores, loadings=loadings, explained=explained[:n_pc])


def module_rank(module: str) -> int:
    try:
        return MODULE_ORDER.index(module)
    except ValueError:
        return len(MODULE_ORDER)


def endpoint_label(endpoint_id: str, abbreviated: bool = False) -> str:
    mapping = ENDPOINT_ABBREV if abbreviated else ENDPOINT_LABEL
    return mapping.get(endpoint_id, endpoint_id.replace("_", " "))


def endpoint_module_map(contract: pd.DataFrame) -> dict[str, str]:
    mapping = dict(zip(contract["endpoint_id"], contract["module"]))
    mapping["corolla_hue"] = "colour"
    return mapping


def sorted_units(units: Iterable[str], module_map: dict[str, str]) -> list[str]:
    return sorted(
        set(units),
        key=lambda value: (module_rank(module_map.get(value, "")), endpoint_label(value)),
    )


def symmetric_matrix(
    pairs: pd.DataFrame,
    units: Sequence[str],
    *,
    left: str = "left",
    right: str = "right",
    value: str = "value",
) -> np.ndarray:
    idx = {name: i for i, name in enumerate(units)}
    out = np.full((len(units), len(units)), np.nan, dtype=float)
    for row in pairs[[left, right, value]].itertuples(index=False):
        a, b, val = row
        if a not in idx or b not in idx:
            continue
        i, j = idx[a], idx[b]
        out[i, j] = float(val)
        out[j, i] = float(val)
    return out


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.08,
        1.04,
        f"({label})",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontweight="bold",
        fontsize=11,
    )


def save_figure(fig: plt.Figure, stem: Path) -> list[Path]:
    stem.parent.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []
    for suffix, kwargs in [
        (".png", {"dpi": 300}),
        (".pdf", {}),
        (".svg", {}),
    ]:
        path = stem.with_suffix(suffix)
        fig.savefig(path, facecolor="white", **kwargs)
        outputs.append(path)
    plt.close(fig)
    return outputs
