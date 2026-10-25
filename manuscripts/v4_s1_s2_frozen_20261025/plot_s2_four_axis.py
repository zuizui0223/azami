"""Render two GEB S2 figures from source-real aggregated, privacy-safe evidence.

Source JSON: analysis/v4/results/20261025_S2_FOUR_AXIS_SOURCE_VERIFIED_SUMMARY.json
CIs are descriptive source-cell bootstrap intervals, not independent invasions.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NAMES = ("Cirsium arvense","Cirsium vulgare")
AXES = ("BIO18","BIO05","RSDS","elevation")
LABELS = ("BIO18 · warm-quarter precipitation",
          "BIO05 · warmest-month temperature",
          "RSDS · radiation","Elevation")

def verified(path: Path) -> dict:
    d=json.loads(path.read_text(encoding="utf-8"))
    if d["status"]!="SOURCE_REAL_EXECUTION_SUCCEEDED":
        raise ValueError("Source real V4 S2 data absent")
    # Protect original input identity; published summaries never export IDs.
    if d["provenance"]["original_20261024_frozen_four_model_member_sha256"]!="3738292e0a64da183dfe219c854c16cda199a7f9769f4928da8ab07657d60b70":
        raise ValueError("Original four-model V4 freeze checksum mismatch")
    if d["provenance"].get("corrected_joint_result_sha256")!="800f2f5dc37d7156324a4756f372306c3f8b28cb00c88d683e86cd01c2c64063":
        raise ValueError("Corrected unit-standardized real-source result not verified")
    if not all(d["species"][sp]["distance_unit_invariant_native_sd_standardized"] for sp in NAMES):
        raise ValueError("Heterogeneous-unit niche distance input rejected")
    if set(d["species"])!=set(NAMES):
        raise ValueError("Unexpected thistle source strata")
    return d

def save(fig, root: Path, basename: str):
    for ext in ("png","svg"):
        fig.savefig(root/(basename+"."+ext),bbox_inches="tight",dpi=230)
    plt.close(fig)

def plot_niche(d: dict, out: Path):
    fig,ax=plt.subplots(figsize=(11.0,6.0))
    fig.subplots_adjust(left=.26,right=.97,top=.86,bottom=.24)
    y=np.arange(4)[::-1]
    for i,sp in enumerate(NAMES):
        v=d["species"][sp]
        mid=np.array([v["equal_cell_mean_intro_minus_native_in_native_sd"][z] for z in AXES])
        lo=np.array([v["equal_cell_bootstrap_95"][z][0] for z in AXES])
        hi=np.array([v["equal_cell_bootstrap_95"][z][1] for z in AXES])
        ax.errorbar(mid,y+(.16 if i==0 else -.16),
                    xerr=np.vstack((mid-lo,hi-mid)),fmt="o",capsize=3,
                    label=sp.replace("Cirsium ","C. "))
    ax.axvline(0,ls="--",lw=1)
    ax.set_yticks(y,LABELS)
    ax.set_xlim(-.8,1.75)
    ax.set_ylabel("")
    ax.set_xlabel("Introduced − native, equal spatial-cell mean (native-cell SD)")
    ax.set_title("Native–introduced climatic occupancy is scale and species dependent")
    ax.legend(loc="upper right",frameon=False)
    fig.text(.03,.115,"399 spatial-cell bootstrap resamples; NOT independent introduction replicates.",fontsize=9)
    fig.text(.03,.075,"BIO18 source-record medians increase +48.2 / +49.7 mm/month, but",fontsize=9)
    fig.text(.03,.04,"the cell-weighted mean changes direction in C. arvense; estimator AND weighting differ.",fontsize=9)
    save(fig,out,"GEB_S2_four_axis_cell_contrast")

def plot_angle(d: dict, out: Path):
    fig,ax=plt.subplots(figsize=(11.3,4.6))
    fig.subplots_adjust(left=.35,right=.97,top=.84,bottom=.22)
    names=[];betas=[];lower=[];upper=[]
    for sp in NAMES:
        v=d["species"][sp]["camera_angle_four_axis_analog"]
        for label,key,ci in (
            ("Range status effect","adjusted_intro_effect_degrees","adjusted_intro_effect_95"),
            ("BIO18 × range interaction","BIO18_intro_interaction_degrees_per_joint_sd","interaction_95")
        ):
            names.append(sp.replace("Cirsium ","C. ")+" · "+label)
            betas.append(float(v[key]))
            lower.append(float(v[ci][0]))
            upper.append(float(v[ci][1]))
    yy=np.arange(4)[::-1];b=np.asarray(betas)
    ax.axvline(0,ls="--",lw=1)
    ax.errorbar(b,yy,xerr=np.vstack((b-np.asarray(lower),np.asarray(upper)-b)),fmt="o",capsize=4)
    ax.set_yticks(yy,names)
    ax.set_xlim(-13,13)
    ax.set_xlabel("Photograph-frame angle estimate (degrees, pointwise 95% CI)")
    ax.set_title("No resolved conditional head-angle contrast in common four-axis envelope")
    fig.text(.035,.09,"Common rectangular 5–95% climate support, NOT matched field experiment.",fontsize=9)
    fig.text(.035,.045,"Intervals include zero; this is NOT evidence of equivalent or genetically unchanged orientation.",fontsize=9)
    save(fig,out,"GEB_S2_angle_four_axis_range")

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--source",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    a=parser.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    d=verified(a.source)
    plot_niche(d,a.out);plot_angle(d,a.out)

if __name__=="__main__":
    main()
