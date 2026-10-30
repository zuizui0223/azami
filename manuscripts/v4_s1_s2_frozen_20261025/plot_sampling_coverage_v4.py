#!/usr/bin/env python3
"""Render Figure 2 from the source-verified, disclosure-limited V4 coverage.

This is an offline presentation renderer. It reads only the public five-degree
cell counts, never source records, photographs, identifiers or point locations.
It changes neither the analytical cohorts nor their original inference. The
primary climate-analogue analysis used one-degree taxon cells; this figure uses
five-degree display cells for all six analytical cohorts.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.ticker import FixedFormatter, FixedLocator, NullLocator
import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DEFAULT_SOURCE = ROOT / "analysis/v4/results/20261030_SAMPLING_COVERAGE_SOURCE_VERIFIED.json"
BASEMAP = ROOT / "legacy/v3/reproducibility/figures/source/natural_earth_110m_land_outline.json"
STEM = "figure2_sampling_coverage"
WIDTH_MM, HEIGHT_MM = 168.0, 196.0
MIN_FONT_PT = 8.0
BLACK = "#111111"
LAND = "#777f86"
GRID = "#c9cdd1"
COLOURS = ("#d8e7f2", "#8eb9d6", "#4e91bd", "#19659a", "#083f69")
EARTH_MEAN_RADIUS_KM = 6371.0088
SCALE_BAR_KM = 10000.0
EXPECTED_BASEMAP_SHA256 = "0f73d385f7d97e02aec1234e0a89fb46109e3fda7b31da7d6ec977094e15f9f4"
ARCHIVED_BASEMAP_CRLF_SHA256 = "c517cff9ef3e3aa6eb89998c6d1e047d741f16be83767a94999574f3b6f3f636"
# Original-source CI 41020812565; an absent or changed verified source must
# never produce a final figure.
EXPECTED_COVERAGE_SHA256: str | None = "ca6d9c787c19758ebed68be53b7370d0fbc8a81b1b10235676c66619df04dc50"
PANEL_ORDER = (
    "s1_angle", "h3c_primary", "arvense_native", "arvense_introduced",
    "vulgare_native", "vulgare_introduced",
)
TITLES = {
    "s1_angle": "Native taxa · image angle",
    "h3c_primary": "Climate-analogue joint traits",
    "arvense_native": r"$\it{C.\ arvense}$ · native",
    "arvense_introduced": r"$\it{C.\ arvense}$ · introduced",
    "vulgare_native": r"$\it{C.\ vulgare}$ · native",
    "vulgare_introduced": r"$\it{C.\ vulgare}$ · introduced",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def provenance_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def read_verified_source(source: Path) -> dict:
    """Accept the final CI bytes and the same public-schema validator as CI."""
    if EXPECTED_COVERAGE_SHA256 is None:
        raise RuntimeError("Awaiting the verified original-source coverage JSON; no figure rendered")
    if sha256(source) != EXPECTED_COVERAGE_SHA256:
        raise ValueError("Figure 2 verified coverage checksum changed")
    if str(ROOT) not in sys.path:
        sys.path.insert(0, str(ROOT))
    from analysis.v4.submission_sampling_coverage import validate_coverage

    data = json.loads(source.read_text(encoding="utf-8"))
    validate_coverage(data)
    verification = data["source_verification"]
    expected = {
        "native_source_records": 313350,
        "focal_source_records": 297622,
        "h3_joint_eligible_records": 157822,
        "h3_joint_eligible_taxa": 242,
        "h3_primary_taxa": 63,
        "h3_primary_taxon_cells": 7117,
        "h3_secondary_taxa": 31,
        "h3_secondary_taxon_cells": 814,
        "original_inference_and_correction_families_unchanged": True,
    }
    if any(verification.get(key) != value for key, value in expected.items()):
        raise ValueError("Figure 2 original-source verification contract changed")
    if not all(p["shown_cells"] > 0 for p in data["panels"]):
        raise ValueError("A required Figure 2 cohort has no disclosable cells")
    return data


def read_basemap(path: Path) -> dict:
    if sha256(path) != EXPECTED_BASEMAP_SHA256:
        raise ValueError("Figure 2 archived land outline checksum changed")
    archived_bytes = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\n", b"\r\n")
    if hashlib.sha256(archived_bytes).hexdigest() != ARCHIVED_BASEMAP_CRLF_SHA256:
        raise ValueError("Figure 2 historical basemap checksum changed")
    land = json.loads(path.read_text(encoding="utf-8"))
    if land.get("source") != "Natural Earth 1:110m land, public domain":
        raise ValueError("Figure 2 basemap provenance changed")
    return land


def cell_polygon(lat_bin: int, lon_bin: int) -> np.ndarray:
    """Sample each cell edge at 1° before the native Mollweide transform.

    These vertices are boundaries of already public coarse cells, not restored
    source coordinates. Subdivision preserves curved meridians at high latitude.
    """
    south, west = -90.0 + 5.0 * lat_bin, -180.0 + 5.0 * lon_bin
    t = np.arange(6, dtype=float)
    points = np.concatenate((
        np.column_stack((west + t, np.full(6, south))),
        np.column_stack((np.full(5, west + 5.0), south + t[1:])),
        np.column_stack((west + 5.0 - t[1:], np.full(5, south + 5.0))),
        np.column_stack((np.full(5, west), south + 5.0 - t[1:])),
    ))
    return np.deg2rad(points)


def coastline_segments(land: dict) -> list[np.ndarray]:
    """Split the antimeridian and sample the outline at intervals of at most 1°."""
    segments = []
    for part in land["parts"]:
        points = np.asarray(part, dtype=float)
        cuts = np.flatnonzero(np.abs(np.diff(points[:, 0])) > 180.0) + 1
        for segment in np.split(points, cuts):
            if len(segment) < 2:
                continue
            sampled = []
            for start, stop in zip(segment[:-1], segment[1:]):
                steps = max(1, int(np.ceil(np.abs(stop - start).max())))
                sampled.extend(np.linspace(start, stop, steps, endpoint=False))
            sampled.append(segment[-1])
            segments.append(np.deg2rad(sampled))
    return segments


def colour_ticks(maximum: int) -> list[int]:
    ticks = []
    value = 5
    while value <= maximum:
        ticks.append(value)
        value *= 10
    # Include the observed maximum without crowding two nearby end labels.
    if maximum > ticks[-1]:
        if len(ticks) > 1 and math.log10(maximum / ticks[-1]) < 0.4:
            ticks[-1] = maximum
        else:
            ticks.append(maximum)
    return ticks


def render(data: dict, land: dict, output: Path) -> dict:
    panels = {panel["panel"]: panel for panel in data["panels"]}
    maximum = max(row["records"] for panel in panels.values() for row in panel["cells"])
    normalizer = LogNorm(vmin=5, vmax=max(5.000001, maximum))
    colourmap = LinearSegmentedColormap.from_list("V4_coverage_blue", COLOURS, N=256)
    outlines = coastline_segments(land)
    metadata = {"panels": [], "text_bounds_checked": False}
    texts = []
    axes = []

    with plt.rc_context({
        "font.family": "DejaVu Sans", "font.size": MIN_FONT_PT,
        "mathtext.fontset": "dejavusans", "text.color": BLACK,
        "axes.labelcolor": BLACK, "axes.edgecolor": LAND,
        "xtick.color": "#50555a", "ytick.color": "#50555a",
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "svg.hashsalt": "azami-figure2-v4-20261030", "savefig.facecolor": "white",
    }):
        fig = plt.figure(figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4), dpi=150)

        def text(x, top, content, *, size=MIN_FONT_PT, **kwargs):
            artist = fig.text(x / WIDTH_MM, 1.0 - top / HEIGHT_MM, content,
                              fontsize=size, va="top", **kwargs)
            texts.append(artist)
            return artist

        for index, name in enumerate(PANEL_ORDER):
            row_index, column = divmod(index, 2)
            panel = panels[name]
            base = 3.0 + 53.0 * row_index
            left = 10.0 + 80.0 * column
            top, width, height = base + 13.0, 72.0, 36.0
            title_x = 4.5 + 82.0 * column
            text(title_x, base, f"({chr(97 + index)}) {TITLES[name]}",
                 size=8.4, fontweight="bold")
            taxa = f"{panel['total_taxa']} taxa · " if index < 2 else ""
            text(title_x, base + 5.0, f"{taxa}n = {panel['total_records']:,}")
            ax = fig.add_axes([left / WIDTH_MM, 1.0 - (top + height) / HEIGHT_MM,
                               width / WIDTH_MM, height / HEIGHT_MM],
                              projection="mollweide")
            axes.append(ax)
            ax.set_facecolor("white")
            ax.set_xticks(np.deg2rad([-120, -60, 0, 60, 120]))
            ax.set_xticklabels(["120°W", "", "0°", "", "120°E"], fontsize=MIN_FONT_PT)
            ax.set_yticks(np.deg2rad([-60, -30, 0, 30, 60]))
            ax.set_yticklabels(["60°S", "", "0°", "", "60°N"], fontsize=MIN_FONT_PT)
            ax.set_longitude_grid_ends(80)
            ax.grid(True, color=GRID, linewidth=0.35, zorder=1)
            ax.spines["geo"].set_linewidth(0.45)
            ax.tick_params(length=0)
            for label in ax.get_xticklabels():
                label.set_bbox({"facecolor": "white", "edgecolor": "none",
                                "alpha": 0.78, "pad": 0.1})
            values = np.asarray([cell["records"] for cell in panel["cells"]])
            # Apply the native geographic transform pointwise to the sampled
            # geometry. Otherwise GeoAxes inserts 75 redundant vertices into
            # every 1° edge, inflating vector exports without improving them.
            projected = [ax.transProjection.transform(cell_polygon(cell["lat_bin"], cell["lon_bin"]))
                         for cell in panel["cells"]]
            display_transform = ax.transAffine + ax.transAxes
            patches = PolyCollection(
                projected,
                facecolors=colourmap(normalizer(values)), edgecolors="none",
                linewidths=0, antialiased=False, zorder=2, transform=display_transform)
            ax.add_collection(patches)
            ax.add_collection(LineCollection(
                [ax.transProjection.transform(segment) for segment in outlines],
                colors=LAND, linewidths=0.35, zorder=3, transform=display_transform))
            # Retain coordinate labels above marks; their equatorial locations
            # are identical in every panel, with a quiet white backing.
            ax.xaxis.set_zorder(4)
            ax.yaxis.set_zorder(4)
            summary = {key: panel[key] for key in (
                "panel", "total_records", "total_taxa", "total_cells",
                "shown_records", "suppressed_records", "shown_cells", "suppressed_cells")}
            summary.update({"letter": chr(97 + index),
                            "map_box_mm": [left, top, width, height],
                            "record_sum_of_rendered_cells": int(values.sum()),
                            "number_of_rendered_cell_polygons": len(values)})
            metadata["panels"].append(summary)

        text(10.0, 166.0, "Records per 5° cell (log scale)")
        colour_ax = fig.add_axes([10.0 / WIDTH_MM, 1.0 - 175.5 / HEIGHT_MM,
                                  84.0 / WIDTH_MM, 2.5 / HEIGHT_MM])
        sm = matplotlib.cm.ScalarMappable(norm=normalizer, cmap=colourmap)
        bar = fig.colorbar(sm, cax=colour_ax, orientation="horizontal")
        bar.solids.set_rasterized(False)
        bar.solids.set_edgecolor("face")  # Avoid PDF viewer seams between colour bins.
        ticks = colour_ticks(maximum)
        bar.locator = FixedLocator(ticks)
        bar.formatter = FixedFormatter([f"{value:,}" for value in ticks])
        bar.update_ticks()
        bar.ax.xaxis.set_minor_locator(NullLocator())
        bar.ax.tick_params(labelsize=MIN_FONT_PT, length=2, width=0.4, pad=2)
        bar.outline.set_linewidth(0.45)

        # Derive the detached common scale from the actual equatorial transform,
        # not from an assumed map width. It is valid at latitude 0° only.
        fig.canvas.draw()
        angular_span = SCALE_BAR_KM / EARTH_MEAN_RADIUS_KM
        endpoints = axes[0].transData.transform(
            np.array([[-angular_span / 2.0, 0.0], [angular_span / 2.0, 0.0]]))
        endpoints = fig.transFigure.inverted().transform(endpoints)
        scale_width = (endpoints[1, 0] - endpoints[0, 0]) * WIDTH_MM
        centre, scale_top = 132.0, 174.3
        start, stop = centre - scale_width / 2.0, centre + scale_width / 2.0
        text(centre, 166.0, "Equatorial scale", ha="center")
        fig.lines.append(plt.Line2D(
            [start / WIDTH_MM, stop / WIDTH_MM],
            [1.0 - scale_top / HEIGHT_MM] * 2,
            transform=fig.transFigure, color=BLACK, linewidth=0.8))
        for x in (start, centre, stop):
            fig.lines.append(plt.Line2D(
                [x / WIDTH_MM] * 2,
                [1.0 - (scale_top - 0.65) / HEIGHT_MM,
                 1.0 - (scale_top + 0.65) / HEIGHT_MM],
                transform=fig.transFigure, color=BLACK, linewidth=0.65))
        text(start, 177.0, "0", ha="center")
        text(stop, 177.0, "10,000 km", ha="center")
        text(WIDTH_MM / 2.0, 186.0,
             "Mollweide equal-area projection · 5° display grid", ha="center")
        text(WIDTH_MM / 2.0, 191.0,
             "Cells with <5 records omitted; panel totals overlap.", ha="center")

        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        all_texts = texts + list(colour_ax.get_xticklabels())
        for ax in axes:
            all_texts.extend(ax.get_xticklabels() + ax.get_yticklabels())
        visible = [artist for artist in all_texts if artist.get_visible() and artist.get_text()]
        for artist in visible:
            box = artist.get_window_extent(renderer)
            if (box.x0 < -0.5 or box.y0 < -0.5 or
                    box.x1 > fig.bbox.width + 0.5 or box.y1 > fig.bbox.height + 0.5):
                raise ValueError(f"Figure 2 text outside canvas: {artist.get_text()}")
        if min(artist.get_fontsize() for artist in visible) < MIN_FONT_PT:
            raise ValueError("Figure 2 typography fell below its minimum size")
        output.mkdir(parents=True, exist_ok=True)
        for extension in ("png", "pdf", "svg"):
            save_metadata = {"Creator": "V4 source-verified sampling-coverage renderer"}
            if extension == "pdf":
                save_metadata.update({"CreationDate": None, "ModDate": None})
            elif extension == "svg":
                save_metadata["Date"] = None
            fig.savefig(output / f"{STEM}.{extension}", dpi=600,
                        metadata=save_metadata, facecolor="white")
        plt.close(fig)

    metadata.update({
        "dimensions_mm": [WIDTH_MM, HEIGHT_MM], "raster_dpi": 600,
        "minimum_text_pt": MIN_FONT_PT, "text_bounds_checked": True,
        "projection": "Mollweide (equal-area), central meridian 0°",
        "graticule_degrees": {"longitude": 60, "latitude": 30},
        "labelled_longitudes_degrees": [-120, 0, 120],
        "labelled_latitudes_degrees": [-60, 0, 60],
        "display_cell_degrees": 5, "minimum_records_per_displayed_cell": 5,
        "colour_scale": {"measure": "pooled records per 5° latitude–longitude cell",
                         "normalization": "logarithmic, shared across all panels",
                         "minimum": 5, "maximum": maximum, "ticks": ticks,
                         "colours": list(COLOURS), "not_area_normalized": True},
        "scale_bar": {"distance_km": SCALE_BAR_KM, "valid_at_latitude_degrees": 0,
                      "spherical_mean_radius_km": EARTH_MEAN_RADIUS_KM,
                      "longitude_span_degrees": math.degrees(angular_span),
                      "length_mm": scale_width,
                      "computed_from_actual_native_projection_transform": True},
        "primary_h3_analysis_cell_degrees": 1,
        "blank_cell_meaning": "Not displayed: absent or fewer than five records; not evidence of biological absence.",
        "cohort_totals_are_not_additive": True,
        "map_is_analytical_sampling_coverage_not_occurrence_range": True,
    })
    return metadata


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--basemap", type=Path, default=BASEMAP)
    parser.add_argument("--output-dir", type=Path, default=HERE / "figures")
    parser.add_argument("--manifest", type=Path, default=HERE / "figure2_render_manifest.json")
    args = parser.parse_args()
    data = read_verified_source(args.source)
    land = read_basemap(args.basemap)
    presentation = render(data, land, args.output_dir)
    manifest = {
        "schema": "V4_FIGURE2_SAMPLING_COVERAGE_RENDER_20261030",
        "status": "SOURCE_VERIFIED_COVERAGE_RENDERED",
        "source": {"path": provenance_path(args.source),
                   "sha256": sha256(args.source), "schema": data["schema"],
                   "public_schema_and_original_cohort_contract_validated": True,
                   "ci_run_id": 41020812565,
                   "ci_commit": "1fa733a463a6e546a69b8386eddbed6b2f1d047d",
                   "ci_artifact_id": 12961447511,
                   "ci_artifact_zip_sha256": "93924e9b36eee790d263ee6f5ca33859a7719fb428a9eb4f22198d3b9df22f8e"},
        "basemap": {"path": provenance_path(args.basemap),
                    "sha256": sha256(args.basemap), "source": land["source"],
                    "archived_crlf_sha256": ARCHIVED_BASEMAP_CRLF_SHA256,
                    "archived_checksum_verified_after_lf_to_crlf": True,
                    "line_ending_note": "Git working-tree LF; historical figure report used CRLF.",
                    "role": "Presentation land outline only; no political boundaries",
                    "rights": "Public domain", "downloaded_by_renderer": False},
        "renderer": {"path": provenance_path(Path(__file__)),
                     "sha256": sha256(Path(__file__)),
                     "matplotlib_version": matplotlib.__version__,
                     "numpy_version": np.__version__, "network_used": False},
        "presentation": presentation,
        "transformations": [
            "Decode public five-degree cell indices to cell boundaries, with no point-level locations.",
            "Project cell vertices pointwise with native Matplotlib Mollweide, then affine-place them; sample boundaries at one-degree intervals.",
            "Colour each displayed cell by its unchanged pooled record count with one shared logarithmic scale.",
            "Keep all disclosed cells; do not smooth, interpolate counts, impute suppressed cells or refit models.",
            "Show full analytical cohort totals in panel titles and retain displayed/suppressed accounting here.",
            "Project the archived Natural Earth land outlines; split antimeridian crossings and subdivide geographic segments to at most one degree.",
        ],
        "outputs": {extension: {"path": provenance_path(args.output_dir / f"{STEM}.{extension}"),
                                "sha256": sha256(args.output_dir / f"{STEM}.{extension}")}
                    for extension in ("png", "pdf", "svg")},
    }
    args.manifest.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
                             encoding="utf-8")
    print(json.dumps({"status": manifest["status"],
                      "source_sha256": manifest["source"]["sha256"],
                      "panels": [{k: p[k] for k in ("panel", "total_records", "shown_records", "suppressed_records")}
                                 for p in presentation["panels"]]}, indent=2))


if __name__ == "__main__":
    main()
