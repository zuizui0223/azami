#!/usr/bin/env python3
"""Figure 5: source-verified aggregate-only climatic vs floral volume comparison.

Two deliberately separate inference panels: (A) H3a climate occupancy
log-volume vs 3D photographed morphospace log-volume, with HC3 pointwise CI;
(B) H3b cross-sector climate-centroid vs floral-centroid distance rank
correlation, displayed without a fake CI. Species are the permutation unit.
No private photographs, taxon identities, observer IDs or coordinates are used.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

SOURCE=Path("analysis/v4/results/20261029_H3_CLIMATE_TRAIT_HYPERVOLUME_SOURCE_VERIFIED.json")
RUN=40841007115
ARTIFACT=12883021095


def read_source(path):
    data=json.loads(Path(path).read_text(encoding="utf-8"))
    if (data["status"]!="VERIFIED_SOURCE_DERIVED_TWO_HYPOTHESES_EXPLORATORY"
        or data["verified_execution"]["run_id"]!=RUN
        or data["verified_execution"]["artifact_id"]!=ARTIFACT
        or data["source"]["source_observations"]!=313350
        or data["source"]["quality_eligible_matched_observations"]!=157822):
        raise ValueError("UNVERIFIED_H3_ORIGINAL_SOURCE_AGGREGATE")
    for key,expected in (("1degree",63),("5degree",31)):
        item=data["analyses"][key]
        if (item["taxa"]!=expected or item["status"]!="source_matched_exploratory"
            or item["breadth"]["status"]!="estimated"
            or item["similarity"]["status"]!="estimated"):
            raise ValueError("H3_ELIGIBILITY_OR_RESULT_MISMATCH")
    return data


def plot(data, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    keys=["1degree","5degree"]
    labels=["1° cells","5° cells"]
    fig,ax=plt.subplots(1,2,figsize=(10.4,4.6),gridspec_kw={"width_ratios":[1.3,1]})
    y=[1,0]
    first=ax[0]
    for i,key in enumerate(keys):
        info=data["analyses"][key]
        fit=info["breadth"]
        est=fit["trait_logvolume_slope_per_climate_logvolume"]
        lo,hi=fit["HC3_95CI"]
        first.errorbar(est,y[i],xerr=[[est-lo],[hi-est]],fmt="o",capsize=4,
                       color="0.23",markersize=7,lw=1.8)
        first.text(0.98,y[i]+.21,
                   f"n={info['taxa']} taxa; perm. P={fit['sector_conditioned_residual_permutation_p_exploratory']:.3f}",
                   transform=first.get_yaxis_transform(),ha="right",fontsize=9,color="0.38")
    first.axvline(0,color="0.5",ls="--",lw=1)
    first.set_yticks(y,labels)
    first.set_ylim(-.65,1.65)
    first.set_xlabel("Slope: climate breadth → floral breadth")
    first.set_title("a  Multivariate breadth",loc="left",fontweight="bold",fontsize=12)
    first.grid(axis="x",alpha=.14)
    second=ax[1]
    for i,key in enumerate(keys):
        info=data["analyses"][key]
        fit=info["similarity"]
        est=fit["climate_distance_vs_trait_distance_spearman"]
        second.scatter(est,y[i],s=63,color="0.25",zorder=2)
        second.text(.98,y[i]+.21,
                   f"{fit['cross_sector_taxon_pairs']:,} pairs; perm. P={fit['within_sector_taxon_label_permutation_p_exploratory']:.3f}",
                   transform=second.get_yaxis_transform(),ha="right",fontsize=9,color="0.38")
    second.axvline(0,color="0.5",ls="--",lw=1)
    second.set_yticks(y,labels)
    second.set_ylim(-.65,1.65)
    second.set_xlim(-.34,.34)
    second.set_xlabel("Spearman ρ: climate vs floral distance")
    second.set_title("b  Cross-region similarity",loc="left",fontweight="bold",fontsize=12)
    second.grid(axis="x",alpha=.14)
    fig.suptitle("Does climatic generalism correspond to more floral diversity?",y=1.0,fontsize=13)
    fig.text(.08,-.015,"Ellipsoidal volume proxy, not true KDE hypervolume. "
                       "Intervals are pointwise; tests exploratory and not part of S1/S2 Holm families.",
             fontsize=8.5,color="0.3")
    fig.tight_layout(rect=[0,.025,1,.95],w_pad=2)
    out=Path(output);out.mkdir(parents=True,exist_ok=True)
    for ext in ("png","svg"):
        fig.savefig(out/("figure5_H3_climatic_breadth_floral_morphospace."+ext),
                    dpi=250,bbox_inches="tight")
    plt.close(fig)
    return out


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path,default=SOURCE)
    p.add_argument("--output-dir",type=Path,required=True)
    args=p.parse_args()
    plot(read_source(args.source),args.output_dir)


if __name__=="__main__":
    main()
