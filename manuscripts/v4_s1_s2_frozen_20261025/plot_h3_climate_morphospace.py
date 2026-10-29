#!/usr/bin/env python3
"""GEB Figure 5: local climate-analogue 3D photographic floral resemblance.

Three aggregate-only panels:
 a) climate-matched floral distance advantage versus within-taxon permutation
 b) conditional Gaussian-KDE L2 overlap difference versus permutation reference
 c) the separate, unresolved climate-breadth / floral-breadth relationship.

No taxon names, observation IDs, individual photographs or exact coordinates.
The original S1 80 and S2 128 Holm families are untouched; the H3c Holm-2
values correct only two exploratory within-taxon permutation measures.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

SOURCE_BREADTH=Path("analysis/v4/results/20261029_H3_CLIMATE_TRAIT_HYPERVOLUME_SOURCE_VERIFIED.json")
SOURCE_ANALOGUE=Path("analysis/v4/results/20261029_H3C_CLIMATE_ANALOGUE_TRAIT_OVERLAP_SOURCE_VERIFIED.json")
EXPECTED_H3A_RUN=40841007115
EXPECTED_H3A_ARTIFACT=12883021095
EXPECTED_H3C_RUN=40847696537
EXPECTED_H3C_ARTIFACT=12885860903

def read_sources(breadth,analogue):
    a=json.loads(Path(breadth).read_text(encoding="utf-8"))
    b=json.loads(Path(analogue).read_text(encoding="utf-8"))
    if (a["status"]!="VERIFIED_SOURCE_DERIVED_TWO_HYPOTHESES_EXPLORATORY"
        or a["verified_execution"]["run_id"]!=EXPECTED_H3A_RUN
        or a["verified_execution"]["artifact_id"]!=EXPECTED_H3A_ARTIFACT
        or a["source"]["quality_eligible_matched_observations"]!=157822
        or b["status"]!="SOURCE_DERIVED_EXPLORATORY_NOT_ORIGINAL_V4_HOLM"
        or b["source"]["run_id"]!=EXPECTED_H3C_RUN
        or b["source"]["artifact_id"]!=EXPECTED_H3C_ARTIFACT
        or b["source"]["quality_eligible_matched_photographs"]!=157822):
        raise ValueError("H3_FIGURE_SOURCE_VERIFICATION_FAILED")
    for key,n in (("1degree",63),("5degree",31)):
        ca=a["analyses"][key];co=b["results"][key]
        if ca["taxa"]!=n or co["n_taxa"]!=n:
            raise ValueError("H3_FIGURE_TAXON_DENOMINATOR_CHANGED")
        if co["climate_analogue_conditional_KDE_L2_overlap"]["status"]!="estimated":
            raise ValueError("H3_CONDITIONAL_DENSITY_NOT_ESTIMATED")
    return a,b


def plot(a,b,output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    fig,axs=plt.subplots(1,3,figsize=(15.7,4.8),
                         gridspec_kw={"width_ratios":[1.08,1.13,1.15]})
    keys=("1degree","5degree")
    labels=("1° cells (63 taxa)","5° cells (31 taxa)")
    ys=[1,0]
    for i,key in enumerate(keys):
        item=b["results"][key]
        observed=item["primary_trait_distance_advantage_of_climate_matched_cells"]
        null=item["null_mean"]
        p=item["holm2_distance_vs_overlap_p"]
        axs[0].plot([null,observed],[ys[i],ys[i]],lw=1.4,color="0.60",zorder=1)
        axs[0].plot(null,ys[i],marker="s",markersize=7,color="0.55",markerfacecolor="white")
        axs[0].plot(observed,ys[i],marker="o",markersize=7,color="#19659a")
        axs[0].text(.99,ys[i]+.20,f"Holm P={p:.3f}",transform=axs[0].get_yaxis_transform(),
                    ha="right",fontsize=9)
        o=item["climate_analogue_conditional_KDE_L2_overlap"]
        observed_overlap=o["conditional_minus_all_cells_overlap"]
        null_overlap=o["null_mean"]
        axs[1].plot([null_overlap,observed_overlap],[ys[i],ys[i]],lw=1.4,color="0.60",zorder=1)
        axs[1].plot(null_overlap,ys[i],marker="s",markersize=7,color="0.55",markerfacecolor="white")
        axs[1].plot(observed_overlap,ys[i],marker="o",markersize=7,color="#19659a")
        axs[1].text(.99,ys[i]+.20,f"Holm P={o['holm2_distance_vs_overlap_p']:.3f}",
                    transform=axs[1].get_yaxis_transform(),ha="right",fontsize=9)
        result=a["analyses"][key]["breadth"]
        estimate=result["trait_logvolume_slope_per_climate_logvolume"]
        lo,hi=result["HC3_95CI"]
        axs[2].errorbar(estimate,ys[i],xerr=[[estimate-lo],[hi-estimate]],
                        fmt="o",color="#19659a",capsize=3,markersize=6,zorder=2)
        axs[2].text(.99,ys[i]+.20,f"P={result['sector_conditioned_residual_permutation_p_exploratory']:.3f}",
                    transform=axs[2].get_yaxis_transform(),ha="right",fontsize=9)
    for ax in axs:
        ax.set_yticks(ys,labels)
        ax.set_ylim(-.5,1.6)
        ax.axvline(0,ls="--",color="0.55",lw=1)
        ax.grid(axis="x",alpha=.13)
    axs[0].set_title("(a) Local 3D trait alignment",loc="left",fontweight="bold")
    axs[0].set_xlabel("Distance advantage in matched climates")
    axs[1].set_title("(b) Conditional trait-density overlap",loc="left",fontweight="bold")
    axs[1].set_xlabel("Matched minus all-cell KDE L2 overlap")
    axs[2].set_title("(c) Climatic versus floral breadth",loc="left",fontweight="bold")
    axs[2].set_xlabel("Climate-breadth slope (pointwise HC3 CI)")
    axs[0].set_xlim(-.025,.16)
    axs[1].set_xlim(-.18,-.04)
    axs[2].set_xlim(-.53,.60)
    fig.suptitle("Climatic matching of local floral phenotypes versus total morphospace breadth",
                 y=1.005,fontsize=13)
    fig.text(.09,-.022,
             "Blue circles: observed; hollow squares: permutation-null mean (a,b). "
             "KDE overlap compares different cell subsets, not true hypervolume intersection.",
             fontsize=9,color="0.3")
    fig.tight_layout(rect=[0,.035,1,.96],w_pad=2.0)
    path=Path(output);path.mkdir(parents=True,exist_ok=True)
    for ext in ("png","svg"):
        fig.savefig(path/f"figure5_H3_climatic_breadth_floral_morphospace.{ext}",
                    dpi=260,bbox_inches="tight")
    plt.close(fig)
    return path


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path,default=SOURCE_BREADTH)
    p.add_argument("--analogue-source",type=Path,default=SOURCE_ANALOGUE)
    p.add_argument("--output-dir",type=Path,required=True)
    args=p.parse_args()
    a,b=read_sources(args.source,args.analogue_source)
    plot(a,b,args.output_dir)


if __name__=="__main__":
    main()
