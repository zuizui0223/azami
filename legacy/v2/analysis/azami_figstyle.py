#!/usr/bin/env python3
"""Single source of truth for Chapter 1 submission figures.

This module fixes final manuscript geometry, typography, colour semantics,
endpoint labels/order and export behaviour.  It intentionally contains no
analysis logic: figure builders may reshape frozen outputs but must not fit or
select models here.

The current review-DOCX usable width is 6.27 in (A4 page, 1 in left/right
margins).  Figures are therefore built at that exact width and inserted at 6.27
in without Word-side scaling.  A production profile can be added later if the
publisher supplies a different final-width specification.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import matplotlib as mpl
import matplotlib.pyplot as plt

__all__ = [
    "WIDTH", "MAX_HEIGHT_IN", "DPI_RASTER", "FONT", "C", "SCALE_CLASS",
    "STABILITY", "MODULE_COLOUR", "Endpoint", "ENDPOINTS",
    "PCA_ENDPOINT_KEYS", "S7_MATRIX_KEYS", "order", "labels", "use",
    "figure", "panel", "panel_letter", "check", "savefig",
]

# ---------------------------------------------------------------------------
# 1. Geometry
# ---------------------------------------------------------------------------

WIDTH = {
    "single": 3.15,   # internal single-column target
    "onehalf": 4.53,  # internal intermediate target
    "double": 6.27,   # exact usable width of the review DOCX
}
MAX_HEIGHT_IN = 7.60
DPI_RASTER = 600
SAVE_VECTOR = True
PDF_DATE = datetime(2026, 9, 2, tzinfo=timezone.utc)

# ---------------------------------------------------------------------------
# 2. Typography
# ---------------------------------------------------------------------------

FONT = {
    "family": ["Arial", "Helvetica", "DejaVu Sans"],
    "panel": 9.0,
    "axis": 8.0,
    "tick": 7.0,
    "legend": 7.0,
    "annot": 7.0,
    "footnote": 6.5,
}

# ---------------------------------------------------------------------------
# 3. Colour semantics (Okabe-Ito)
# ---------------------------------------------------------------------------

C = {
    "black": "#000000",
    "orange": "#E69F00",
    "skyblue": "#56B4E9",
    "green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermil": "#D55E00",
    "purple": "#CC79A7",
    "grey": "#BFBFBF",
    "lightgrey": "#ECECEC",
    "rule": "#808080",
}

SCALE_CLASS = {
    "both_scales": dict(color=C["purple"], hatch=None, label="both scales"),
    "within_only": dict(color=C["blue"], hatch=None, label="within only"),
    "among_only": dict(color=C["vermil"], hatch=None, label="among only"),
    "neither": dict(color=C["lightgrey"], hatch=None, label="neither"),
    "not_comparable": dict(color="white", hatch="///", label="not comparable"),
}

STABILITY = {
    "stable": dict(color=C["blue"], hatch=None, label="direction stable"),
    "unstable": dict(color=C["orange"], hatch="\\\\", label="direction unstable"),
}

MODULE_COLOUR = {
    "orientation": C["blue"],
    "colour": C["purple"],
    "display": C["skyblue"],
    "shape": C["green"],
    "involucre": C["orange"],
    "projection": C["yellow"],
    "surface": C["vermil"],
}

# ---------------------------------------------------------------------------
# 4. Endpoint registry
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Endpoint:
    key: str
    group: str
    long: str
    short: str
    measured: bool
    joint_hue: bool = False


# Exact endpoint IDs from v2_full27_endpoint_inventory.csv.  The ordering is
# presentation order, not a redefinition of the frozen ontology.
ENDPOINTS: tuple[Endpoint, ...] = (
    Endpoint("orientation_image_vertical_angle", "orientation", "orientation angle", "orient.", True),

    Endpoint("corolla_lab_lightness", "colour", "corolla lightness", "L*", True),
    Endpoint("corolla_lab_chroma", "colour", "corolla chroma", "C*", True),
    Endpoint("corolla_hue_sin", "colour", "hue sine", "hue sin", True, joint_hue=True),
    Endpoint("corolla_hue_cos", "colour", "hue cosine", "hue cos", True, joint_hue=True),
    Endpoint("corolla_hue", "colour", "joint corolla hue", "hue", True, joint_hue=True),
    Endpoint("corolla_purple_pixel_fraction", "colour", "purple pixel fraction", "purple", False),
    Endpoint("corolla_redmagenta_pixel_fraction", "colour", "red-magenta pixel fraction", "red-mag.", False),
    Endpoint("corolla_white_pixel_fraction", "colour", "white pixel fraction", "white", False),
    Endpoint("corolla_yellow_pixel_fraction", "colour", "yellow pixel fraction", "yellow", False),

    Endpoint("visible_floret_fraction", "display", "visible floret fraction", "floret", False),

    Endpoint("capitulum_outline_aspect_ratio", "shape", "outline aspect ratio", "aspect", True),
    Endpoint("capitulum_outline_circularity", "shape", "outline circularity", "circular.", True),
    Endpoint("capitulum_outline_solidity", "shape", "outline solidity", "solidity", True),
    Endpoint("capitulum_width_profile_cv", "shape", "width-profile CV", "width CV", True),

    Endpoint("involucre_length_width_ratio", "involucre", "involucre L/W", "inv. L/W", True),
    Endpoint("involucre_apical_taper_ratio", "involucre", "apical taper", "apical", True),
    Endpoint("involucre_basal_taper_ratio", "involucre", "basal taper", "basal", True),

    Endpoint("bract_projection_asymmetry", "projection", "projection asymmetry", "asym.", True),
    Endpoint("bract_projection_maximum", "projection", "projection maximum", "proj. max", True),
    Endpoint("bract_projection_p95", "projection", "projection p95", "proj. p95", True),
    Endpoint("bract_projection_peak_density", "projection", "projection-peak density", "peaks", True),
    Endpoint("bract_projection_roughness", "projection", "projection roughness", "rough.", True),
    Endpoint("bract_spread_fraction", "projection", "spread fraction", "spread", True),

    Endpoint("involucre_surface_edge_density", "surface", "surface edge density", "edge", True),
    Endpoint("involucre_surface_high_frequency_energy", "surface", "surface high-frequency energy", "high-freq.", True),
    Endpoint("involucre_surface_lbp_entropy", "surface", "surface LBP entropy", "LBP", True),
    Endpoint("involucre_surface_specular_fraction", "surface", "surface specular fraction", "specular", True),
)

# PCA uses 18 measured endpoint columns: hue sine/cosine are separate and the
# four validation-only surface signals are excluded.
PCA_ENDPOINT_KEYS: tuple[str, ...] = tuple(
    e.key for e in ENDPOINTS
    if e.measured
    and e.key != "corolla_hue"
    and e.group != "surface"
)

# The repeated-observation S7 matrices use the same biological subset but hue
# is one circular inferential unit, leaving 17 units.
S7_MATRIX_KEYS: tuple[str, ...] = tuple(
    e.key for e in ENDPOINTS
    if e.measured
    and e.key not in {"corolla_hue_sin", "corolla_hue_cos"}
    and e.group != "surface"
)


def order(*, split_hue: bool = False, measured_only: bool = False) -> list[str]:
    """Return canonical endpoint / analysis-unit order.

    split_hue=True returns sine and cosine as separate measured columns and
    omits the synthetic joint-hue unit.  split_hue=False returns joint hue and
    omits sine/cosine, matching the 26-row analysis-unit denominator when
    measured_only=False.
    """
    keys: list[str] = []
    for endpoint in ENDPOINTS:
        if endpoint.joint_hue:
            if split_hue and endpoint.key == "corolla_hue":
                continue
            if not split_hue and endpoint.key in {"corolla_hue_sin", "corolla_hue_cos"}:
                continue
        if measured_only and not endpoint.measured:
            continue
        keys.append(endpoint.key)
    return keys


def labels(keys: Iterable[str], *, short: bool = False) -> list[str]:
    lut = {e.key: (e.short if short else e.long) for e in ENDPOINTS}
    values = list(keys)
    missing = [key for key in values if key not in lut]
    if missing:
        raise KeyError(f"endpoint keys not in figure registry: {missing}")
    return [lut[key] for key in values]

# ---------------------------------------------------------------------------
# 5. Matplotlib defaults
# ---------------------------------------------------------------------------


def use(*, grid: bool = False) -> None:
    """Apply submission figure defaults.

    constrained_layout is deliberately not enabled globally because several
    existing figures use explicit GridSpec / add_axes placement.  Individual
    figures may opt in after visual QA.
    """
    mpl.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": DPI_RASTER,
        "savefig.bbox": "standard",
        "savefig.pad_inches": 0.0,
        "savefig.transparent": False,
        "savefig.facecolor": "white",

        "font.family": "sans-serif",
        "font.sans-serif": FONT["family"],
        "font.size": FONT["annot"],
        "axes.titlesize": FONT["panel"],
        "axes.titleweight": "bold",
        "axes.titlelocation": "left",
        "axes.labelsize": FONT["axis"],
        "xtick.labelsize": FONT["tick"],
        "ytick.labelsize": FONT["tick"],
        "legend.fontsize": FONT["legend"],
        "legend.title_fontsize": FONT["legend"],
        "mathtext.default": "regular",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",

        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.6,
        "axes.edgecolor": "#333333",
        "axes.labelcolor": "black",
        "axes.axisbelow": True,
        "axes.grid": grid,
        "grid.color": "#D9D9D9",
        "grid.linewidth": 0.4,
        "grid.alpha": 1.0,
        "axes.prop_cycle": mpl.cycler(color=[
            C["blue"], C["orange"], C["green"], C["purple"],
            C["vermil"], C["skyblue"], C["yellow"], C["black"],
        ]),

        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.color": "#333333",
        "ytick.color": "#333333",

        "lines.linewidth": 1.0,
        "lines.markersize": 3.5,
        "lines.markeredgewidth": 0.5,
        "patch.linewidth": 0.5,
        "hatch.linewidth": 0.6,
        "errorbar.capsize": 1.5,

        "legend.frameon": False,
        "legend.handlelength": 1.4,
        "legend.handletextpad": 0.5,
        "legend.columnspacing": 1.2,
        "legend.borderaxespad": 0.3,
    })


def figure(*, width: str = "double", height: float | None = None,
           aspect: float | None = None, **kwargs):
    if width not in WIDTH:
        raise ValueError(f"width must be one of {sorted(WIDTH)}")
    w = WIDTH[width]
    if height is None:
        height = w * (0.75 if aspect is None else aspect)
    if height > MAX_HEIGHT_IN:
        warnings.warn(
            f"height {height:.2f} in exceeds {MAX_HEIGHT_IN:.2f} in; redesign panels before submission",
            stacklevel=2,
        )
    return plt.figure(figsize=(w, height), **kwargs)


def panel(ax, letter: str, title: str = "", *, pad: float = 4.0) -> None:
    text = f"({letter})" + (f" {title}" if title else "")
    ax.set_title(text, loc="left", fontsize=FONT["panel"], fontweight="bold", pad=pad)


def panel_letter(fig, ax, letter: str, *, dx: float = -0.02, dy: float = 0.01) -> None:
    bbox = ax.get_position()
    fig.text(
        bbox.x0 + dx, bbox.y1 + dy, f"({letter})",
        fontsize=FONT["panel"], fontweight="bold", ha="left", va="bottom",
    )


def check(fig, *, width: str = "double") -> None:
    w, h = fig.get_size_inches()
    target = WIDTH[width]
    if abs(w - target) > 0.01:
        raise ValueError(
            f"figure width is {w:.3f} in but {width} target is {target:.3f} in; "
            "build at final width and do not use bbox_inches='tight'"
        )
    if h > MAX_HEIGHT_IN:
        warnings.warn(
            f"figure height {h:.2f} in exceeds {MAX_HEIGHT_IN:.2f} in",
            stacklevel=2,
        )


def savefig(fig, stem: str, *, width: str = "double",
            outdir: str | Path = "figures", vector: bool | None = None,
            metadata: dict | None = None) -> list[Path]:
    """Write fixed-canvas 600-dpi PNG plus vector PDF."""
    check(fig, width=width)
    destination = Path(outdir)
    destination.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []

    png = destination / f"{stem}.png"
    fig.savefig(png, dpi=DPI_RASTER, bbox_inches=None, pad_inches=0.0, facecolor="white")
    written.append(png)

    if SAVE_VECTOR if vector is None else vector:
        pdf = destination / f"{stem}.pdf"
        pdf_metadata = {
            "Title": stem,
            "Author": "",
            "Creator": "Azami Chapter 1 canonical figure style",
            "CreationDate": PDF_DATE,
            "ModDate": PDF_DATE,
        }
        if metadata:
            pdf_metadata.update(metadata)
        fig.savefig(pdf, bbox_inches=None, pad_inches=0.0, facecolor="white", metadata=pdf_metadata)
        written.append(pdf)
    return written


if __name__ == "__main__":
    use()
    assert len(order(split_hue=False)) == 26
    assert len(order(split_hue=True, measured_only=True)) == 22
    assert len(PCA_ENDPOINT_KEYS) == 18
    assert len(S7_MATRIX_KEYS) == 17
    assert "corolla_hue" in S7_MATRIX_KEYS
    assert "corolla_hue_sin" not in S7_MATRIX_KEYS
    assert all("surface" not in key for key in S7_MATRIX_KEYS)
    print(
        "analysis_units=26 measured_endpoints=22 "
        f"pca_endpoints={len(PCA_ENDPOINT_KEYS)} s7_units={len(S7_MATRIX_KEYS)}"
    )
