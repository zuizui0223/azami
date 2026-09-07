"""Plot verified all27 cached coverage and eligibility loss, not accuracy."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import numpy as np

from .workflow import ROOT, digest, text_digest

REGISTRY = ROOT / "ch1_global/v2/ontology/ch1_continuous_trait_contract.csv"
CONTRACT = Path(__file__).with_name("technical_coverage_figure_contract.json")
CONDITIONS = Path(__file__).with_name("perturbation_contract.json")
LABELS = [
    "Orientation angle", "Corolla lightness", "Corolla chroma", "Hue sine", "Hue cosine",
    "Visible floral fraction", "White pixel fraction", "Red/magenta pixel fraction",
    "Purple pixel fraction", "Yellow pixel fraction", "Outline aspect ratio",
    "Outline circularity", "Outline solidity", "Width-profile CV", "Involucre length/width",
    "Apical taper", "Basal taper", "Projection roughness", "Projection 95th percentile",
    "Projection maximum", "Bract spread fraction", "Projection-peak density",
    "Projection asymmetry", "Surface edge density", "Surface LBP entropy",
    "Surface high-frequency energy", "Surface specular fraction",
]
CONDITION_LABELS = ["Crop left", "Crop right", "Crop up", "Crop down", "Resolution 75%",
                    "Resolution 50%", "Gamma 0.8", "Gamma 1.2", "Intensity 0.8",
                    "Intensity 1.2", "Warm gain", "Cool gain", "Blur sigma 1"]


def load_data(summary, verification):
    receipt = json.loads(verification.read_text(encoding="utf-8"))
    if receipt["status"] != "SCALAR_PERTURBATION_MEANS_AND_ATTRITION_RECOMPUTED_BY_SQL" or receipt["scalar_rows_checked"] != 1134:
        raise ValueError("A completed scalar arithmetic verification is required")
    if digest(summary) != receipt["summary_sha256"]:
        raise ValueError("Summary hash differs from verified input")
    with summary.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 1218:
        raise ValueError("Full 1218-row summary required")
    with REGISTRY.open(encoding="utf-8", newline="") as handle:
        endpoints = [r["endpoint_id"] for r in csv.DictReader(handle)]
    conditions = [r["id"] for r in json.loads(CONDITIONS.read_text(encoding="utf-8"))["conditions"]]
    selected = [r for r in rows if r["exposure_stratum"] == "all_cached" and r["metric_kind"] == "registered_endpoint"]
    lookup = {(r["metric_id"], r["condition"]): r for r in selected}
    expected = {(e, c) for e in endpoints for c in conditions}
    if len(endpoints) != 27 or len(selected) != len(expected) or set(lookup) != expected:
        raise ValueError("Complete, unique all27 by 14 condition grid required")
    baseline = [lookup[e, "baseline"] for e in endpoints]
    scheduled = {int(r["scheduled_heads"]) for r in selected}
    if len(scheduled) != 1 or next(iter(scheduled)) <= 0:
        raise ValueError("Scheduled-head denominator changed across grid")
    n_heads = scheduled.pop()
    counts = np.array([int(r["baseline_usable_heads"]) for r in baseline])
    if np.any(counts < 0) or np.any(counts > n_heads):
        raise ValueError("Invalid baseline eligibility count")
    loss = np.full((27, 13), np.nan)
    for i, e in enumerate(endpoints):
        for j, c in enumerate(conditions[1:]):
            row = lookup[e, c]
            if int(row["baseline_usable_heads"]) != counts[i]:
                raise ValueError("Original baseline count changed under perturbation")
            raw = row["component_weighted_loss_fraction"]
            if raw:
                value = float(raw)
                if not math.isfinite(value) or not 0 <= value <= 1:
                    raise ValueError("Invalid eligibility loss fraction")
                loss[i, j] = value * 100
    return endpoints, n_heads, counts, loss


def render(summary, verification, out):
    endpoints, n_heads, counts, loss = load_data(summary, verification)
    out = out.resolve()
    if out == ROOT or (ROOT in out.parents and not any(out.is_relative_to(ROOT / p) for p in ("local_data", "outputs"))):
        raise ValueError("Use an external or ignored local output directory")
    out.mkdir(parents=True, exist_ok=False)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "text.color": "#292929", "axes.labelcolor": "#292929",
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    fig = plt.figure(figsize=(10.8, 10.3))
    grid = fig.add_gridspec(1, 2, width_ratios=[1.35, 4.0], left=.257, right=.985,
                           bottom=.237, top=.850, wspace=.11)
    a = fig.add_subplot(grid[0, 0])
    b = fig.add_subplot(grid[0, 1], sharey=a)
    y = np.arange(27)
    a.barh(y, counts / n_heads * 100, height=.7, color="#416C91", edgecolor="#29445C", linewidth=.35)
    for i, n in enumerate(counts):
        a.text(n/n_heads*100+2, i, f"{n:,}", va="center", fontsize=8)
    a.set(xlim=(0, 100), ylim=(26.5, -.5), yticks=y, yticklabels=LABELS,
          xticks=[0, 50, 100], xlabel="Eligible heads (%)")
    a.set_title("(a) Baseline coverage", loc="left", fontsize=10, pad=13)
    a.xaxis.grid(True, color="#DDDDDD", linewidth=.5)
    a.set_axisbelow(True)
    a.tick_params(axis="y", length=0, pad=7)
    for side in ("top", "right"):
        a.spines[side].set_visible(False)
    tones = ["#F3F7FA", "#DBE7F0", "#B8CDDE", "#91AEC6", "#698FAC", "#416C91", "#29445C"]
    cmap = ListedColormap(tones).with_extremes(bad="#DDDDDD")
    bounds = [0, 5, 10, 20, 40, 60, 80, 100]
    im = b.imshow(np.ma.masked_invalid(loss), cmap=cmap, norm=BoundaryNorm(bounds, len(tones)),
                  aspect="auto", interpolation="nearest")
    b.set_title("(b) Loss of baseline eligibility", loc="left", fontsize=10, pad=13)
    b.set_xticks(range(13), CONDITION_LABELS, rotation=55, ha="right", rotation_mode="anchor")
    b.tick_params(axis="y", left=False, labelleft=False)
    b.tick_params(axis="x", length=0, pad=7, labelsize=8)
    for i in range(27):
        for j in range(13):
            v = loss[i, j]
            label = "NA" if np.isnan(v) else ("<1" if 0 < v < 1 else f"{v:.0f}")
            b.text(j, i, label, ha="center", va="center", fontsize=6.8,
                   color="white" if v >= 60 else "#292929")
    for boundary in (9.5, 13.5, 22.5):
        a.axhline(boundary, color="#BBBBBB", linewidth=.55)
        b.axhline(boundary, color="white", linewidth=1.1)
    fig.suptitle("Image-feature coverage under specified perturbations", x=.257, y=.974,
                 ha="left", fontsize=14)
    fig.text(.257, .934, f"Historical local cache: {n_heads:,} detected heads; all 27 registered endpoints.\n"
             "Left: count and share of all heads. Right: component-weighted loss among baseline-eligible heads.",
             ha="left", va="top", fontsize=9, linespacing=1.6)
    cax = fig.add_axes([.495, .110, .455, .015])
    fig.colorbar(im, cax=cax, orientation="horizontal", ticks=bounds)
    cax.set_xlabel("Baseline eligibility lost (%)", fontsize=9)
    fig.text(.025, .037, "Loss includes failure of minimum resolution and other saved quality requirements; raw values remain retained.\n"
             "Specified digital changes are not calibrated camera errors. Coverage and repeatability do not establish accuracy.\n"
             "Source: saved 14-condition pass and independently checked scalar summary; no ecological models in this figure.",
             fontsize=8, va="bottom", linespacing=1.6)
    for ext in ("png", "pdf"):
        fig.savefig(out/f"technical_coverage.{ext}", dpi=300, facecolor="white")
    plt.close(fig)
    report = {"status": "TECHNICAL_COVERAGE_FIGURE_RENDERED_VISUAL_REVIEW_REQUIRED",
              "summary_sha256": digest(summary), "verification_sha256": digest(verification),
              "figure_contract_sha256_text_lf": text_digest(CONTRACT),
              "implementation_sha256_text_lf": text_digest(Path(__file__)),
              "endpoint_order": endpoints, "scheduled_heads": n_heads,
              "outputs": {name: digest(out/name) for name in ("technical_coverage.png", "technical_coverage.pdf")},
              "matplotlib_version": matplotlib.__version__, "ecological_models_executed": False}
    (out/"figure_report.json").write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8", newline="\n")
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--summary", type=Path, required=True)
    p.add_argument("--verification", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    args = p.parse_args()
    print(json.dumps(render(args.summary, args.verification, args.out_dir), indent=2))


if __name__ == "__main__":
    main()
