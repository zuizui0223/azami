#!/usr/bin/env python3
"""Compose V4 Figure 1 from checksum-verified archived photographic panels.

This is a presentation-only renderer. It does not download photographs,
recompute measurements, alter an analytical cohort or fit statistical models.
The original image-specific licences and exact source-page verification are
recorded in FIGURE1_PHOTO_CREDITS_20261030.json. All coordinates below refer to
the archived composite, not the original camera images.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle
from PIL import Image


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE = ROOT / "legacy/v3/reproducibility/figures/source"
MANIFEST = HERE / "FIGURE1_PHOTO_CREDITS_20261030.json"
STEM = "figure1_photographic_measurement_pipeline"
WIDTH_MM, HEIGHT_MM = 168.0, 202.0
BLUE = "#19659a"
BLACK = "#111111"
RULE = "#b5bcc2"
EXPECTED = {
    "Figure_1_real_photo_measurement_provenance.png":
        "834e10fb4c32e133988fbcbe6cb98c8ce7d5cedc9941593a31a58adc3b318855",
    "Figure_1_real_photo_measurement_provenance.csv":
        "f69be11eb4496cddfe7720f6371ff19dca85fac4356b43c2b0dcaac26e5e7c22",
    "Figure_1_real_photo_measurement_provenance.json":
        "8935c5842ec3d7e837b5ed005a9587c25bb0f241304b9fc1df2b5fcec409052e",
}
ARCHIVED_CRLF_SHA256 = {
    "Figure_1_real_photo_measurement_provenance.csv":
        "2bd45cbc4eb3de0873808008df12079dd31f3ea3cec0024e840ece65e07952fd",
    "Figure_1_real_photo_measurement_provenance.json":
        "0d95d7f73929440a381b0e6239ffbfe9774f6f275a0900bdd3cf96fb65a611b1",
}
# Left, top, right, bottom; right and bottom are exclusive.
# The top-row crops retain stored detection boxes and their detector scores.
# Detailed panels retain the archived tight head, floral mask, and outline/hull.
CROPS = {
    "a": (100, 580, 642, 942),
    "b": (814, 442, 1492, 894),
    "c": (1574, 486, 2296, 967),
    "d": (129, 1260, 612, 1950),
    "e": (1693, 1260, 2176, 1950),
    "f": (2475, 1260, 2958, 1950),
}


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def read_verified_sources(source: Path, manifest_path: Path):
    checks = {}
    for name, expected in EXPECTED.items():
        data = (source / name).read_bytes()
        observed = sha256(data)
        if observed != expected:
            raise ValueError(f"Figure 1 source checksum changed: {name}")
        check = {"sha256_current_bytes": observed, "byte_checksum_verified": True}
        if name in ARCHIVED_CRLF_SHA256:
            normalized = sha256(data.replace(b"\r\n", b"\n").replace(b"\n", b"\r\n"))
            if normalized != ARCHIVED_CRLF_SHA256[name]:
                raise ValueError(f"Figure 1 archived metadata changed: {name}")
            check.update({
                "sha256_archived_crlf_bytes": normalized,
                "archived_checksum_verified_after_lf_to_crlf": True,
                "line_ending_note": "Git working-tree LF; historical figure report used CRLF.",
            })
        checks[name] = check
    with (source / "Figure_1_real_photo_measurement_provenance.csv").open(
        encoding="utf-8", newline=""
    ) as handle:
        rows = list(csv.DictReader(handle))
    if [r["taxon_name"] for r in rows] != [
        "Cirsium pannonicum", "Cirsium mohavense", "Cirsium subcoriaceum"
    ]:
        raise ValueError("Figure 1 illustrative photo selection changed")
    if any(rows[0][domain + "_status"] != "usable" for domain in
           ("orientation", "colour", "shape")):
        raise ValueError("Detailed example is not usable in all three domains")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest["rights_status"] != "source_photo_pages_and_api_verified":
        raise ValueError("Figure 1 rights manifest has not passed source verification")
    for row, credit in zip(rows, manifest["photographs"], strict=True):
        if (row["photo_id"] != str(credit["photo_id"])
            or row["photo_license_code"] != credit["recorded_license_code"]
            or not credit["source_page_license_verified"]):
            raise ValueError("Photo attribution manifest and archived source disagree")
    image = Image.open(source / "Figure_1_real_photo_measurement_provenance.png").convert("RGB")
    if image.size != (3600, 3154):
        raise ValueError("Figure 1 source composite dimensions changed")
    return image, rows, manifest, checks


def render(image: Image.Image, rows: list[dict], output: Path) -> dict:
    texts = []
    with plt.rc_context({
        "font.family": "DejaVu Sans", "font.size": 8.0,
        "text.color": BLACK, "axes.labelcolor": BLACK,
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "svg.fonttype": "none", "svg.hashsalt": "azami-figure1-v4-20261030",
        "savefig.facecolor": "white",
    }):
        fig = plt.figure(figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4), dpi=150)

        def text(x, y, value, size=8.0, **kwargs):
            artist = fig.text(x / WIDTH_MM, 1 - y / HEIGHT_MM, value,
                              fontsize=size, va="top", **kwargs)
            texts.append(artist)
            return artist

        def line(x1, y1, x2, y2, colour=RULE, width=0.55):
            fig.lines.append(plt.Line2D(
                [x1 / WIDTH_MM, x2 / WIDTH_MM],
                [1 - y1 / HEIGHT_MM, 1 - y2 / HEIGHT_MM],
                transform=fig.transFigure, color=colour, lw=width))

        def arrow(x1, y1, x2, y2, colour=BLUE, width=0.9):
            artist = FancyArrowPatch(
                (x1 / WIDTH_MM, 1 - y1 / HEIGHT_MM),
                (x2 / WIDTH_MM, 1 - y2 / HEIGHT_MM),
                transform=fig.transFigure, arrowstyle="-|>",
                mutation_scale=7, lw=width, color=colour)
            fig.add_artist(artist)

        def photo(panel, x, top, width, height):
            crop = image.crop(CROPS[panel])
            scale = min(width / crop.width, height / crop.height)
            shown_w, shown_h = crop.width * scale, crop.height * scale
            left = x + (width - shown_w) / 2
            y = top + (height - shown_h) / 2
            ax = fig.add_axes([left / WIDTH_MM,
                               1 - (y + shown_h) / HEIGHT_MM,
                               shown_w / WIDTH_MM, shown_h / HEIGHT_MM])
            ax.imshow(crop, interpolation="nearest")
            # Replace legacy raster panel letters with consistent vector labels.
            # The white patch is restricted to the upper-left presentation area,
            # outside the focal head and its measurement overlays.
            patch_w, patch_h = 5.5 / scale, 4.1 / scale
            ax.add_patch(Rectangle((-0.5, -0.5), patch_w, patch_h,
                                   facecolor="white", edgecolor="none", zorder=5))
            label = ax.text(0.7 / scale, 0.45 / scale, f"({panel})", fontsize=8.0,
                            fontweight="bold", color=BLACK, va="top", zorder=6)
            texts.append(label)
            ax.set_xlim(-0.5, crop.width - 0.5)
            ax.set_ylim(crop.height - 0.5, -0.5)
            ax.set_axis_off()
            return {"panel": panel, "source_crop_xyxy": list(CROPS[panel]),
                    "display_box_mm": [left, y, shown_w, shown_h],
                    "legacy_panel_label_replaced": True,
                    "label_patch_mm": [5.5, 4.1],
                    "pixel_interpolation": "nearest"}

        text(5, 4, "Source photographs and capitulum localisation", 10.0, fontweight="bold")
        columns = [5, 59, 113]
        for x, taxon in zip(columns, ("Cirsium pannonicum", "Cirsium mohavense", "Cirsium subcoriaceum")):
            text(x + 25, 12, taxon, 8.5, fontstyle="italic", ha="center")
        crops = [photo(panel, x, 18, 50, 33.3) for panel, x in zip("abc", columns)]
        text(5, 54, "Stored detector boxes identify candidate heads; photographs illustrate the measurement pipeline.", 8.0)
        line(5, 61, 163, 61)

        text(5, 64, "Continuous measurements from the boxed head in (a)", 9.0, fontweight="bold")
        for x, name in zip(columns, ("Orientation", "Colour", "Outline")):
            text(x + 25, 72, name, 9.0, fontweight="bold", ha="center")
        for panel, x in zip("def", columns):
            crops.append(photo(panel, x + 9.5, 79, 31, 44.3))

        # This is a reference direction, not a reconstructed or inferred head axis.
        arrow(12, 111, 12, 84, colour=BLACK, width=0.85)
        text(5.3, 78.6, "Image\ntop", 8.0)
        text(5, 126.5, "Image angle", 8.0)
        text(55, 126.5, f"{float(rows[0]['orientation_angle_degrees']):.1f}°", 8.0, ha="right")
        text(5, 132, "0° top · 90° horizontal", 8.0)
        text(5, 137, "180° bottom", 8.0)

        colour_values = [
            ("Lightness, L*", f"{float(rows[0]['corolla_lab_lightness']):.1f}"),
            ("Chroma, C*", f"{float(rows[0]['corolla_lab_chroma']):.1f}"),
            ("Circular hue", f"{float(rows[0]['corolla_hue_degrees']):.1f}°"),
        ]
        for index, (label, value) in enumerate(colour_values):
            text(59, 126.5 + index * 5.0, label)
            text(109, 126.5 + index * 5.0, value, ha="right")
        outline_values = [
            ("Aspect ratio", f"{float(rows[0]['shape_aspect_ratio']):.2f}"),
            ("Circularity", f"{float(rows[0]['shape_circularity']):.2f}"),
            ("Solidity", f"{float(rows[0]['shape_solidity']):.2f}"),
        ]
        for index, (label, value) in enumerate(outline_values):
            text(113, 126.5 + index * 5.0, label)
            text(163, 126.5 + index * 5.0, value, ha="right")
        text(113, 141.5, "Width-profile CV: Table S2.1", 8.0)
        text(5, 148.5, "Eight constructs; hue is represented by two circular coordinates. Definitions: Table S2.1.", 8.0)
        line(5, 155, 163, 155)

        # Record aggregation is domain-specific: the contributing heads can differ.
        text(84, 159, "Eligible head measurements → photograph-level medians by domain", 8.4,
             fontweight="bold", ha="center")
        text(84, 164.5, "Source links, taxon, location and joint coordinates are retained.", 8.0, ha="center")
        for centre in (30, 84, 138):
            arrow(centre, 171, centre, 175)
        for x, label in zip(columns, ("Among-taxon\nclimate gradients",
                                      "Within-range\ncolour gradients",
                                      "Climate-analogue\njoint resemblance")):
            rect = Rectangle((x / WIDTH_MM, 1 - 186 / HEIGHT_MM),
                             50 / WIDTH_MM, 10 / HEIGHT_MM,
                             transform=fig.transFigure, facecolor="#f3f7fa",
                             edgecolor=RULE, lw=0.55)
            fig.add_artist(rect)
            text(x + 25, 177, label, 8.0, ha="center", linespacing=1.25)

        text(5, 191, "Photo sources: Bobo-X (a, d–f; CC BY 4.0); Jesse Rorabaugh (b; CC0 1.0);", 8.0)
        text(5, 195.5, "Diego Blanco (c; CC0 1.0). Cropped and annotated; source links are in the caption.", 8.0)

        fig.canvas.draw()
        canvas = fig.canvas.get_renderer()
        outside = []
        for artist in texts:
            box = artist.get_window_extent(canvas)
            if (box.x0 < -0.5 or box.y0 < -0.5
                or box.x1 > fig.bbox.width + 0.5 or box.y1 > fig.bbox.height + 0.5):
                outside.append(artist.get_text())
        if outside:
            plt.close(fig)
            raise ValueError(f"Figure 1 text lies outside the canvas: {outside}")
        output.mkdir(parents=True, exist_ok=True)
        for extension in ("png", "pdf", "svg"):
            metadata = {"Creator": "Azami V4 Figure 1 renderer"} if extension != "png" else {"Software": "Azami V4 Figure 1 renderer"}
            if extension == "pdf":
                metadata.update({"CreationDate": None, "ModDate": None})
            elif extension == "svg":
                metadata["Date"] = None
            fig.savefig(output / f"{STEM}.{extension}", dpi=600,
                        metadata=metadata, facecolor="white")
        plt.close(fig)
    return {"width_mm": WIDTH_MM, "height_mm": HEIGHT_MM,
            "minimum_vector_text_size_pt": 8.0, "raster_export_dpi": 600,
            "text_within_canvas": True, "crops": crops}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=SOURCE)
    parser.add_argument("--output-dir", type=Path, default=HERE / "figures")
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    args = parser.parse_args()
    image, rows, manifest, checks = read_verified_sources(args.source_dir, args.manifest)
    render_info = render(image, rows, args.output_dir)
    render_info.update({
        "source_checks": checks,
        "renderer_sha256": sha256(Path(__file__).read_bytes()),
        "output_files": {f"{STEM}.{ext}": sha256((args.output_dir / f"{STEM}.{ext}").read_bytes())
                         for ext in ("png", "pdf", "svg")},
        "numerical_models_rerun": False,
        "original_photographs_downloaded": False,
        "historical_1_1_degree_overlay_used": False,
        "example_width_profile_cv_value_displayed": False,
    })
    manifest["render"] = render_info
    args.manifest.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "rendered", "output": str(args.output_dir),
                      "source_image_sha256_verified": True,
                      "photo_rights_status": manifest["rights_status"],
                      "width_mm": WIDTH_MM, "height_mm": HEIGHT_MM,
                      "minimum_vector_text_size_pt": 8.0}, indent=2))


if __name__ == "__main__":
    main()
