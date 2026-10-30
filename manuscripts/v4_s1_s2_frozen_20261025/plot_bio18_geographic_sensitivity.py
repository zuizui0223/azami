"""Generate GEB main Figure 3: original BIO18 slope and geographic dependence.

Input is ONLY original-source verified public aggregate JSON from the
successful 2026-10-29 confidential-source Actions run. No original observation
coordinates, photo data or taxon identity are read or plotted.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

import numpy as np


EXPECTED = 0.04544684307686103
SOURCE_RUN = 40837107465
SOURCE_ARTIFACT = 12881655588
DEFAULT = Path("analysis/v4/results/20261029_S1_BIO18_GEOGRAPHY_SOURCE_VERIFIED.json")


def extract(data):
    meta=data["source"]
    if (data["status"]!="EMPIRICAL_BIO18_SOURCE_VERIFIED_NOT_IN_ORIGINAL_HOLM_FAMILY"
        or data["n"] !=168907 or data["taxa"]!=178
        or meta["original_secure_workflow_run"]!=SOURCE_RUN
        or meta["public_artifact_id"]!=SOURCE_ARTIFACT
        or abs(meta["original_beta_verified_in_runner"]-EXPECTED)>1e-10):
        raise ValueError("UNVERIFIED_SOURCE_BIO18_GEO_AGGREGATE")
    final=[]
    for weighting, prefix in (("photo_weighted","Photographs"),
                              ("equal_5degree_cell","Equal 5° cells")):
        group=data["geographic_fits"][weighting]
        models=group["global_geographic_adjustments"]
        for key, name in (("none","No geographic controls"),
                          ("longitude_sector","3 longitude sectors"),
                          ("spherical_degree2","Spherical geographic basis")):
            fit=models[key]
            if fit["status"]!="estimated" or fit["n_taxa"]!=178:
                raise ValueError("GLOBAL_GEO_FIT_MISSING")
            final.append((f"{prefix}: {name}",fit))
        east=group["leave_one_predeclared_sector_out"]["E060_E180"]
        if east["status"]!="estimated" or east["heldout_n_taxa"]<20:
            raise ValueError("EAST_SECTOR_HOLDOUT_NOT_ESTIMABLE")
        final.append((f"{prefix}: omit E60–E180 sector",east))
    for label, fit in final:
        beta=fit["standardised_beta"]
        lo, hi=fit["HC3_pointwise_95CI"]
        if not (np.isfinite(beta) and np.isfinite(lo) and np.isfinite(hi)
                and lo <= beta <= hi):
            raise ValueError("INVALID_GEO_SENSITIVITY_INTERVAL")
    if abs(final[0][1]["standardised_beta"]-EXPECTED)>3e-7:
        raise ValueError("FROZEN_S1_BETA_NOT_RECOVERED")
    return final


def plot(data_path, out_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    data=json.loads(Path(data_path).read_text(encoding="utf-8"))
    result=extract(data)
    names=[label for label, _ in result]
    y=np.arange(len(result))[::-1]
    estimates=np.array([m["standardised_beta"] for _,m in result])
    intervals=np.array([m["HC3_pointwise_95CI"] for _,m in result])
    err=np.array([estimates-intervals[:,0],intervals[:,1]-estimates])
    fig,ax=plt.subplots(figsize=(9.4,5.0))
    ax.errorbar(estimates,y,xerr=err,fmt="o",capsize=3.5,lw=1.6)
    ax.axvline(0,linestyle="--",linewidth=1)
    ax.axhline(3.5,linewidth=.65)
    ax.set_yticks(y,names)
    ax.set_xlabel("BIO18 among-taxon coefficient (source-standardised)")
    ax.set_ylim(-.75,7.75)
    ax.grid(axis="x",alpha=.18)
    ax.set_title("Geographic sensitivity of the global precipitation–capitulum association")
    fig.text(.12,.015,"HC3 pointwise 95% intervals; exploratory geographic adjustments, not additional Holm tests.",
             fontsize=9)
    fig.tight_layout(rect=[0,.06,1,1])
    target=Path(out_dir)
    target.mkdir(parents=True,exist_ok=True)
    for extension in ("png","svg"):
        fig.savefig(target/f"figure3_BIO18_geographic_sensitivity.{extension}",
                    dpi=250, bbox_inches="tight")
    plt.close(fig)
    return target


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path,default=DEFAULT)
    p.add_argument("--output-dir",type=Path,required=True)
    args=p.parse_args()
    plot(args.source,args.output_dir)


if __name__=="__main__":
    main()
