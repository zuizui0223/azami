#!/usr/bin/env python3
"""GEB Figure 5: local climate-analogue 3D photographic floral resemblance.

The main figure shows climate-matched floral distance advantage versus the
within-taxon permutation mean at 1 and 5 degrees. The original three-panel
distance, Gaussian-KDE overlap and climate/floral-breadth display is retained
as a separate supplementary export. All displays use aggregate sources only.

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
        or b["verified_execution"]["github_actions_run"]!=EXPECTED_H3C_RUN
        or b["verified_execution"]["aggregate_artifact_id"]!=EXPECTED_H3C_ARTIFACT
        or b["source"]["quality_eligible_matched_observations"]!=157822):
        raise ValueError("H3_FIGURE_SOURCE_VERIFICATION_FAILED")
    for key,n in (("1degree",63),("5degree",31)):
        ca=a["analyses"][key];co=b["results"][key]
        if ca["taxa"]!=n or co["n_taxa"]!=n:
            raise ValueError("H3_FIGURE_TAXON_DENOMINATOR_CHANGED")
        if co["climate_analogue_conditional_KDE_L2_overlap"]["status"]!="estimated":
            raise ValueError("H3_CONDITIONAL_DENSITY_NOT_ESTIMATED")
    return a,b


def plot(a,b,output):
    """Render the H3c distance result alone under the legacy Figure 5 name."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    blue="#19659a"
    neutral="#69737c"
    with plt.rc_context({"font.family":"DejaVu Sans","font.size":10,
                         "axes.labelcolor":"#28323c","text.color":"#28323c",
                         "xtick.color":"#47525c","ytick.color":"#28323c"}):
        fig,ax=plt.subplots(figsize=(9.4,5.3))
        fig.subplots_adjust(left=.25,right=.96,bottom=.33,top=.75)
        fig.text(.04,.95,"Similar local climates, closer head phenotypes",
                 fontsize=16,fontweight="bold",va="top")
        fig.text(.04,.885,"Image-referenced head angle · CIELAB chroma · 2D outline aspect ratio",
                 fontsize=10.5,color=neutral)

        keys=("1degree","5degree")
        ys=(1,0)
        labels=("1° cells · 63 taxa\nPrimary analysis",
                "5° cells · 31 taxa\nSpatial sensitivity")
        for key,y in zip(keys,ys):
            item=b["results"][key]
            observed=item["primary_trait_distance_advantage_of_climate_matched_cells"]
            null=item["null_mean"]
            p=item["holm2_distance_vs_overlap_p"]
            ax.plot(null,y,marker="s",markersize=8,markeredgewidth=1.5,
                    color=neutral,markerfacecolor="white",zorder=3)
            ax.plot(observed,y,marker="o",markersize=9,color=blue,zorder=3)
            ax.annotate(f"{observed:+.5f}",(observed,y),xytext=(0,13),
                        textcoords="offset points",ha="center",fontsize=11,
                        fontweight="bold",color=blue)
            ax.annotate(f"Holm-2 P = {p:.3f}",(observed,y),xytext=(0,-21),
                        textcoords="offset points",ha="center",fontsize=10)
            ax.annotate(f"{null:+.5f}",(null,y),xytext=(0,-21),
                        textcoords="offset points",ha="center",fontsize=10,
                        color=neutral,bbox={"facecolor":"white","edgecolor":"none","pad":1})

        ax.set_yticks(ys,labels)
        ax.tick_params(axis="y",length=0,pad=16,labelsize=10.5)
        ax.tick_params(axis="x",length=4,pad=7)
        ax.set_ylim(-.55,1.55)
        ax.set_xlim(-.009,.115)
        ax.set_xticks([0,.025,.05,.075,.1],
                      ["0","0.025","0.050","0.075","0.100"])
        ax.axvline(0,ls="--",color=neutral,lw=1,zorder=1)
        ax.grid(axis="y",color="#e4e8eb",lw=.8,zorder=0)
        for side in ("top","right","left"):
            ax.spines[side].set_visible(False)
        ax.spines["bottom"].set_color("#89929a")
        ax.set_xlabel("Distance advantage from climate matching (standardised 3D distance)",
                      fontsize=10.5,labelpad=13)
        handles=[Line2D([],[],marker="o",linestyle="none",markersize=8,
                        color=blue,label="Observed"),
                 Line2D([],[],marker="s",linestyle="none",markersize=7,
                        markeredgewidth=1.5,color=neutral,markerfacecolor="white",
                        label="Within-taxon permutation mean")]
        fig.legend(handles=handles,loc="upper left",bbox_to_anchor=(.245,.84),
                   ncol=2,frameon=False,fontsize=10,handletextpad=.6,columnspacing=2)
        fig.text(.25,.17,"All cell pairs − climate-matched pairs; positive values indicate more similar heads.",
                 fontsize=9.5,color=neutral)
        fig.text(.04,.115,"Exploratory Holm-2 correction retains the original distance + KDE-overlap family; 199 permutations.",
                 fontsize=9,color=neutral)
        fig.text(.04,.068,"Statistical support depends on spatial grain, traits and regions; sensitivity analyses are reported in SI.",
                 fontsize=9,color=neutral)
        path=Path(output);path.mkdir(parents=True,exist_ok=True)
        for ext in ("png","svg"):
            fig.savefig(path/f"figure5_H3_climatic_breadth_floral_morphospace.{ext}",
                        dpi=260,bbox_inches="tight")
        plt.close(fig)
    return path


def plot_supplementary(a,b,output):
    """Preserve the original distance, overlap and breadth context panels."""
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
             "KDE overlap compares different cell subsets, not true hypervolume intersection. "
             "Original H3c results; 1° Monte Carlo and trait/region sensitivity are reported in SI.",
             fontsize=9,color="0.3")
    fig.tight_layout(rect=[0,.035,1,.96],w_pad=2.0)
    path=Path(output);path.mkdir(parents=True,exist_ok=True)
    for ext in ("png","svg"):
        fig.savefig(path/f"SI_H3_context_distance_overlap_breadth.{ext}",
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
    plot_supplementary(a,b,args.output_dir)


if __name__=="__main__":
    main()

