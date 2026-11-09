"""Display native and introduced cell means in the four rainfall combinations."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TRAITS = {"orientation_angle_degrees_median":"下向き角度（°）",
          "corolla_lab_chroma_median":"花冠の彩度 C*",
          "shape_aspect_ratio_median":"頭花の投影縦横比"}


def plot(folder):
    g = pd.read_csv(folder/"vulgare_rain_groups.csv").query("scope=='whole' and q==0.5")
    plt.rcParams.update({"font.family":"Meiryo", "pdf.fonttype":42})
    fig, axes = plt.subplots(3,2, figsize=(7.5,8), layout="constrained")
    for row,(col,label) in enumerate(TRAITS.items()):
        low,high = g[col].min(),g[col].max()
        for inv, ax in enumerate(axes[row]):
            cells = g[g.introduced==inv].set_index("group")
            ids = np.array([[2,0],[3,1]])
            values = cells[col].loc[ids.ravel()].to_numpy().reshape(2,2)
            ax.imshow(values,origin="lower",cmap="Blues",vmin=low-(high-low)*.8,vmax=high,aspect="auto")
            for (y,x),k in np.ndenumerate(ids):
                ax.text(x,y,f"{cells.loc[k,col]:.2f}\n(n={cells.loc[k,'cells']:.0f})",
                        ha="center",va="center",color="white" if values[y,x]>(low+high)/2 else "#131e33")
            ax.set_xticks([0,1],["少雨","多雨"])
            ax.set_yticks([0,1],["季節差 小","季節差 大"])
            ax.set_title(("在来域" if inv==0 else "移入域")+"　"+label,fontsize=11)
    fig.suptitle("Cirsium vulgare  ｜  降水の組合せと頭花形質",fontsize=15)
    fig.supxlabel(f"境界は在来セルの中央値（年降水量 {g.rain.iloc[0]:.1f} mm、BIO15 {g.season.iloc[0]:.1f}）。\n数値は1°セルの平均。濃淡は各行の2図で共通。",fontsize=9)
    for ext in ["png","pdf"]:
        fig.savefig(folder/f"rain_head_traits.{ext}",dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("folder", type=Path)
    plot(p.parse_args().folder)
