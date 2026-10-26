#!/usr/bin/env python3
"""Figure 4: source-verified native/introduced photographed-colour slopes.

Consumes ONLY the six-row CSV produced by verify_s2_frozen_colour_context.py
from SHA-256-verified original October 24 V4 source artifact 12547876704.

These are range-STRATIFIED regressions, not a range-status interaction, and
standardized beta magnitudes are not guaranteed comparable in physical units.
"""
from __future__ import annotations
import argparse
import csv
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

EXPECTED = [
    ("Cirsium arvense", "native", "lightness", "topo_elevation",
     22936, +0.05626673194073065, 1.0),
    ("Cirsium arvense", "introduced", "lightness", "topo_elevation",
     26171, -0.17373187423223405, 0.0013123349175077526),
    ("Cirsium vulgare", "native", "chroma", "chelsa_rsds_mean",
     17335, +0.1368829067685122, 0.10239759622660725),
    ("Cirsium vulgare", "introduced", "chroma", "chelsa_rsds_mean",
     27200, +0.1345037874025455, 0.014568133833876904),
    ("Cirsium vulgare", "native", "lightness", "chelsa_rsds_mean",
     17335, -0.04880796672586844, 1.0),
    ("Cirsium vulgare", "introduced", "lightness", "chelsa_rsds_mean",
     27200, -0.14280073287481068, 0.019254137758649138)
]

def read_source(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))
    if len(rows) != len(EXPECTED):
        raise ValueError("Exactly six original native/introduced colour source rows required")
    for row, expected in zip(rows, EXPECTED):
        species, status, trait, exposure, n, beta, holm = expected
        actual=(row["species"],row["range_status"],row["trait"],row["exposure"])
        if actual != (species,status,trait,exposure):
            raise ValueError(f"Original source row identity/order changed: {actual}")
        if int(row["n"]) != n or not math.isclose(float(row["beta_std"]), beta, abs_tol=1e-10):
            raise ValueError("Frozen coefficient / n modified")
        if not math.isclose(float(row["holm_p"]), holm, abs_tol=1e-9):
            raise ValueError("Frozen Holm value modified")
        lo,hi=float(row["ci_lower"]),float(row["ci_upper"])
        if not(lo < beta < hi):
            raise ValueError("Source CI invalid")
        if row["range_status"]=="introduced" and holm>=.05:
            raise ValueError("Expected all three original introduced Holms to pass")
    if not (math.isclose(float(rows[2]["bh_p"]), .017760016689709657, abs_tol=1e-9)
            and float(rows[2]["holm_p"])>.05):
        raise ValueError("Native vulgare chroma BH/Holm interpretation changed")
    return rows

def plot(rows: list[dict[str,str]], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    fig,axs=plt.subplots(3,1,figsize=(10.8,9.0),sharex=True)
    titles=[
        r"$Cirsium\ arvense$ — image lightness × elevation",
        r"$Cirsium\ vulgare$ — image chroma × solar radiation",
        r"$Cirsium\ vulgare$ — image lightness × solar radiation",
    ]
    # Native: hollow square and dashed CI; introduced: filled circle + solid CI.
    for i,ax in enumerate(axs):
        ax.axvline(0,lw=1.1,linestyle="--",color="0.35",zorder=0)
        ax.set_title(titles[i],loc="left",fontsize=12,pad=9)
        for j,(y,label,marker) in enumerate(((0,"Native","s"),(1,"Introduced","o"))):
            r=rows[2*i+j]
            b=float(r["beta_std"])
            low, high=float(r["ci_lower"]),float(r["ci_upper"])
            hlm=float(r["holm_p"])
            color="0.28" if j==0 else "#1767a5"
            ax.errorbar(b,y,xerr=[[b-low],[high-b]],fmt=marker,
                        color=color,markerfacecolor="white" if j==0 else color,
                        markeredgecolor=color,markersize=7,capsize=4,
                        linewidth=2,elinewidth=1.7,zorder=3)
            stat=("Holm < 0.05" if hlm<.05 else
                  ("BH < 0.05 only" if i==1 and j==0 else "Not Holm-supported"))
            ax.text(.345,y,f"n={int(r['n']):,}  |  {stat}",
                    va="center",ha="left",fontsize=9,color="0.25")
        ax.set_yticks([0,1],["Native","Introduced"])
        ax.invert_yaxis()
        ax.set_ylim(1.65,-.65)
        ax.set_xlim(-.31,.62)
        ax.grid(axis="x",alpha=.18)
        ax.spines[["top","right"]].set_visible(False)
    axs[-1].set_xlabel("Within-stratum standardized trait–environment coefficient β (pointwise 95% CI)")
    fig.suptitle("Frozen V4 · Original colour–environment associations by range",
                 x=.06,y=.995,ha="left",fontsize=14)
    fig.text(.055,.02,
             "24 Oct 2026 source artifact 12547876704 (SHA verified upstream); 3 paired endpoints, 6 original fits.\n"
             "Differences between separate coefficients have NOT been tested; ranges differ in sampled climate and photographic support.\n"
             "These are photographic CIELAB metrics, not chemical pigment measurements or native/introduced rainfall niches.",
             fontsize=8.5,ha="left",va="bottom",color="0.34")
    fig.subplots_adjust(left=.17,right=.73,bottom=.18,top=.91,hspace=.47)
    for ext in ("png","svg"):
        fig.savefig(output_dir/f"figure4_S2_native_introduced_colour_coefficients.{ext}",
                    dpi=220,bbox_inches="tight")
    plt.close(fig)

def main() -> None:
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path,required=True,
                   help="SHA-verified six-row CSV from verify_s2_frozen_colour_context.py")
    p.add_argument("--output-dir",type=Path,required=True)
    args=p.parse_args()
    plot(read_source(args.source),args.output_dir)
    print("Verified source colour-pair figure saved:",args.output_dir)

if __name__=="__main__":
    main()
