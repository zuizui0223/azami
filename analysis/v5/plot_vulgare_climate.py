"""Plot V5 aggregate counts and medians, reusing the fixed native occupancy mask."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import PatchCollection
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle


def plot(source, grid_path, destination, scope):
    grid = np.load(grid_path)
    all_tiles = pd.read_csv(source / "vulgare_tiles.csv")
    valid_chroma = all_tiles.loc[all_tiles.n >= 5, "corolla_lab_chroma_median"]
    color_scale = Normalize(2 * np.floor(valid_chroma.min()/2), 2 * np.ceil(valid_chroma.max()/2))
    tiles = all_tiles.query("scope == @scope")
    profile = pd.read_csv(source / "vulgare_profiles.csv").query("scope == @scope")
    expected = profile.n.sum()
    assert (tiles.groupby("view").n.sum() == expected).all()
    assert (tiles.groupby("view").rare_n.sum() == profile.loc[profile.rare, "n"].sum()).all()
    plt.rcParams.update({"font.family": "Meiryo", "axes.spines.top": False,
                         "axes.spines.right": False, "pdf.fonttype": 42})
    fig, axs = plt.subplots(2, 2, figsize=(10, 8), layout="constrained")
    for row, view in enumerate(["pca", "annual"]):
        g = tiles.query("view == @view and n >= 5").copy()
        for ax, value, cmap, norm, label in zip(
            axs[row], ["climate_rarity", "corolla_lab_chroma_median"],
            ["YlGnBu", "magma"], [Normalize(0, 1), color_scale],
            ["在来域での気候の珍しさ（順位 0–1）", "花冠の彩度 C*"]):
            rects = []
            for r in g.itertuples():
                if view == "pca":
                    ix, iy = int(r.x)*4, int(r.y)*4
                    x, y = grid["x"][ix], grid["y"][iy]
                    w, h = grid["x"][ix+4]-x, grid["y"][iy+4]-y
                else:
                    x, y, w, h = r.x, r.y, 3, 300
                rects.append(Rectangle((x, y), w, h))
            pc = PatchCollection(rects, cmap=cmap, norm=norm, edgecolor="white", linewidth=.5)
            pc.set_array(g[value].to_numpy())
            ax.add_collection(pc)
            ax.autoscale_view()
            if view == "pca":
                xx, yy = [(grid[k][1:]+grid[k][:-1])/2 for k in ["x","y"]]
                ax.contour(xx, yy, grid["native95"].T, levels=[.5], colors="#137744", linewidths=1.3)
                for (ix, iy), tag in [((12, 5), "A"), ((6, 8), "B")]:
                    if ((g.x == ix) & (g.y == iy)).any():
                        ax.text(grid["x"][ix*4+2], grid["y"][iy*4+2], tag, ha="center",
                                va="center", fontsize=8, fontweight="bold",
                                bbox=dict(facecolor="white", edgecolor="none", pad=.7))
                ax.set(xlabel=f"気候 PC1（{grid['variance'][0]:.1%}）",
                       ylabel=f"気候 PC2（{grid['variance'][1]:.1%}）")
            else:
                ax.set(xlabel="年平均気温（°C）", ylabel="年降水量（mm）")
            fig.colorbar(pc, ax=ax, shrink=.78, label=label)
        axs[row,0].text(.02,.97,f"表示 {g.n.sum():,} / {expected:,} セル",va="top",
                         transform=axs[row,0].transAxes,fontsize=8,bbox=dict(facecolor="white",edgecolor="none",alpha=.9))
    axs[0,0].set_title("在来域ではどの気候が珍しいか")
    axs[0,1].set_title("移入先のどの気候で彩度が高いか")
    fig.suptitle("Cirsium vulgare  ｜  " + ("移入域全体" if scope=="all" else "両域に共通する背景気候"),fontsize=16)
    fig.supxlabel("各区画は移入域の1°セルの中央値。5セル以上を表示。緑線は在来の気候占有の主要95%境界。", fontsize=9)
    for ext in ["png","pdf"]:
        fig.savefig(destination / f"vulgare_{scope}.{ext}", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    for arg in ["source", "grid", "destination"]:
        p.add_argument(arg, type=Path)
    args = p.parse_args()
    args.destination.mkdir(parents=True, exist_ok=True)
    for scope in ["all","common"]:
        plot(args.source, args.grid, args.destination, scope)
