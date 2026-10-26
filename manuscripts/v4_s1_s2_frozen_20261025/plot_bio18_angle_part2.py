#!/usr/bin/env python3
"""V4 BIO18: main among-taxon result, complete S1 SI, descriptive two-scale SI.

Frozen 24-Oct-2026 source ZIP checksum is mandatory. Plot only originally fitted
rows. The primary figure contains the ONE Holm-supported among-taxon term. SI
provides ALL 80 original native-taxon tests and the raw/BH/Holm distinctions,
including the within-taxon BIO18 angle estimate. It is not a new regression.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
from zipfile import ZipFile
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SOURCE_SHA256 = "6e081dc081bb8ee2311a05ac5b7fcddd3f6309be2fa909e87d432789f4e0e205"
SOURCE_JSON = "V4_S1_S2_four_vs_vif5_vif10_alltraits.json"
EXPECTED = {
    "between_equal_taxa": {
        "beta": .04544684307686103,
        "lo": .022348959442086275,
        "hi": .06854472671163579,
        "n": 168907,
        "raw": .00011507938226118898,
        "bh": .009206350580895119,
        "holm": .009206350580895119,
    },
    "within_taxon": {
        "beta": .05870038416959123,
        "lo": .0027185919353511986,
        "hi": .11468217640383127,
        "n": 168907,
        "raw": .039865107745336185,
        "bh": .39583657274970613,
        "holm": 1.0,
    },
}

def original_frozen_rows(archive: Path) -> tuple[list[dict], dict]:
    if hashlib.sha256(archive.read_bytes()).hexdigest() != SOURCE_SHA256:
        raise ValueError("Frozen V4 source ZIP SHA-256 mismatch")
    with ZipFile(archive) as z:
        all_data = json.loads(z.read(SOURCE_JSON))
    rows = all_data["multiplicity"]["S1/four_bio18"]["rows"]
    if len(rows) != 80 or sum(float(r["p_Holm"]) < .05 for r in rows) != 1:
        raise ValueError("S1 80-test Holm family count changed")
    selected = {r["scale"]: r for r in rows
                if r["group"] == "angle"
                and r["coordinate"] == "orientation_image_vertical_angle"
                and r["environment"] == "chelsa_bio18"}
    if set(selected) != set(EXPECTED):
        raise ValueError("Expected exactly two original S1 scale rows")
    for scale, frozen in EXPECTED.items():
        r = selected[scale]
        actual = {
            "beta": r["beta_std"],
            "lo": r["ci95_wald"][0],
            "hi": r["ci95_wald"][1],
            "n": r["n"],
            "raw": r["p_original"],
            "bh": r["p_BH"],
            "holm": r["p_Holm"],
        }
        for field, expected in frozen.items():
            if abs(float(actual[field]) - expected) > 1e-10:
                raise ValueError(f"Frozen {scale} {field} changed")
    return rows, selected

def export_si(rows: list[dict], selected: dict, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    keys = ("test_number", "group", "coordinate", "environment", "scale",
            "test_kind", "n", "beta_std", "ci95_lower", "ci95_upper",
            "p_original", "p_BH", "p_Holm", "Holm_pass")
    def unpack(r: dict, index: int) -> dict:
        return {
            "test_number": index, "group": r.get("group"),
            "coordinate": r.get("coordinate"), "environment": r.get("environment"),
            "scale": r.get("scale"), "test_kind": r.get("test_kind"),
            "n": r.get("n"), "beta_std": r.get("beta_std"),
            "ci95_lower": (r.get("ci95_wald") or ["", ""])[0],
            "ci95_upper": (r.get("ci95_wald") or ["", ""])[1],
            "p_original": r["p_original"], "p_BH": r["p_BH"],
            "p_Holm": r["p_Holm"], "Holm_pass": r["p_Holm"] < .05,
        }
    with (out / "SI_S1_all_80_original_raw_BH_Holm.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(unpack(r, i) for i, r in enumerate(rows, 1))
    with (out / "SI_S1_BIO18_image_angle_both_scales.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for scale in ("between_equal_taxa", "within_taxon"):
            index = rows.index(selected[scale]) + 1
            w.writerow(unpack(selected[scale], index))

def main_figure(r: dict, dest: Path) -> None:
    b = r["beta_std"]
    lo, hi = r["ci95_wald"]
    fig, ax = plt.subplots(figsize=(9.6, 4.8))
    fig.subplots_adjust(left=.19, right=.95, top=.72, bottom=.28)
    ax.axvline(0, color=".45", linestyle="--", linewidth=1)
    ax.errorbar(b, 0, xerr=[[b-lo], [hi-b]], fmt="o", markersize=9,
                capsize=5, linewidth=1.8, color="#17639c")
    ax.set_yticks([0], ["Among native taxa"])
    ax.set_xlim(-.015, .13)
    ax.set_ylim(-.55, .65)
    ax.grid(axis="x", alpha=.15)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlabel("Standardized image-angle association with BIO18 (pointwise 95% CI)")
    fig.text(.055,.95, "Capitulum presentation and warm-quarter precipitation",
             ha="left", va="top", size=13, weight="bold")
    fig.text(.055,.86,"Original V4 Holm-supported among-taxon association",
             ha="left", va="top", size=10)
    ax.text(.095,.21,f"β={b:+.4f}\n95% CI [{lo:+.4f},{hi:+.4f}]\n"
            f"Holm P={r['p_Holm']:.4f}",ha="left",va="center",fontsize=10)
    fig.text(.055,.11, "168,907 angle-eligible records / 178 taxa; S1 80-test Holm family.\n"
             "BIO18: warmest-quarter mean monthly precipitation, not actual event rainfall.\n"
             "Photograph-relative angle, not measured gravitational head orientation.",
             ha="left", va="top", size=8.6, color=".30")
    for ext in ("png","svg"):
        fig.savefig(dest/f"figure3_BIO18_among_taxa_image_angle.{ext}",
                    dpi=220,bbox_inches="tight")
    plt.close(fig)

def supplementary_figure(rows: dict, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(10.2, 5.4))
    fig.subplots_adjust(top=.77,left=.24,right=.76,bottom=.27)
    ax.axvline(0,color=".4",linestyle="--",linewidth=1)
    for y,(scale,label,marker) in enumerate((
        ("between_equal_taxa","Among native taxa","o"),
        ("within_taxon","Within native taxa","s"))):
        row=rows[scale];b=row["beta_std"];lo,hi=row["ci95_wald"]
        ax.errorbar(b,1-y,xerr=[[b-lo],[hi-b]],fmt=marker,
                    markersize=8,capsize=5,linewidth=1.8,
                    color="#17639c" if y==0 else ".40",
                    markerfacecolor="#17639c" if y==0 else "white")
        status="Holm-supported" if y==0 else "Not Holm-supported"
        ax.text(.12,1-y,
                f"raw P={row['p_original']:.5g}\n"
                f"BH P={row['p_BH']:.4g}\n"
                f"Holm P={row['p_Holm']:.4g}\n{status}",
                va="center",ha="left",fontsize=9)
    ax.set_xlim(-.015,.35)
    ax.set_ylim(-.55,1.55)
    ax.set_yticks([1,0],["Among native taxa","Within native taxa"])
    ax.grid(axis="x",alpha=.14)
    ax.spines[["top","right"]].set_visible(False)
    ax.set_xlabel("Standardized photographic angle × BIO18 coefficient (pointwise 95% CI)")
    fig.text(.055,.95,"SI — BIO18 angle associations before and after multiplicity correction",
             ha="left",va="top",size=13,weight="bold")
    fig.text(.055,.88,"Same original 80-test family; pointwise intervals are NOT simultaneous",
             ha="left",va="top",size=10)
    fig.text(.055,.12,
             "The within-taxon estimate excludes zero in a pointwise interval, but does NOT survive BH or Holm.\n"
             "These are different statistical scales; non-significance is not evidence of a zero slope.\n"
             "The remaining 78 tests and all raw/BH/Holm P values are in SI_S1_all_80_original_raw_BH_Holm.csv.",
             ha="left",va="top",fontsize=8.4,color=".30")
    for ext in ("png","svg"):
        fig.savefig(out/f"SI_Figure_BIO18_image_angle_both_scales.{ext}",
                    dpi=220,bbox_inches="tight")
    plt.close(fig)

def main() -> None:
    parser=argparse.ArgumentParser()
    parser.add_argument("--artifact",required=True,type=Path)
    parser.add_argument("--output-dir",required=True,type=Path)
    parser.add_argument("--si-dir",required=True,type=Path)
    args=parser.parse_args()
    all_rows, selected=original_frozen_rows(args.artifact)
    args.output_dir.mkdir(parents=True,exist_ok=True)
    export_si(all_rows,selected,args.si_dir)
    main_figure(selected["between_equal_taxa"],args.output_dir)
    supplementary_figure(selected,args.si_dir)
    print("Verified original S1 80-family: 1 Holm-positive; main/complete SI exports finished")

if __name__=="__main__":
    main()
