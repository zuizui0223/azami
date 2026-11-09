"""Verify expansion and contrast estimates and display tested rarity slopes."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

ROOT=Path(__file__).parent/"results"
OUT=ROOT/"expansion_consistency"


def main():
    shifts=pd.read_csv(OUT/"rarity_gradient.csv")
    models=pd.read_csv(OUT/"range_gradient.csv")
    support=pd.read_csv(OUT/"gradient_support.csv")
    assert len(shifts)==8 and len(models)==216 and len(support)==8
    assert (support.ess<=support.n_cells+1e-9).all()
    comp=models.pivot(index=["species","coarsen","model","trait","axis"],columns="term",values="beta")
    np.testing.assert_allclose(comp.introduced-comp.native,comp.difference,atol=1e-9)
    assert (shifts.holm>=shifts.p-1e-12).all() and (models.holm>=models.p-1e-12).all()
    font="Meiryo" if "Meiryo" in font_manager.get_font_names() else "Noto Sans CJK JP"
    plt.rcParams.update({"font.family":[font,"DejaVu Sans"],"font.size":9})
    d=models.query("coarsen==4 and model=='occupancy_rank' and term!='difference'")
    fig,axes=plt.subplots(1,2,figsize=(8.4,3.3),layout="constrained",sharex=True,sharey=True)
    for ax,(sp,g) in zip(axes,d.groupby("species")):
        for k,col,name,dy in [("native","#277d83","在来域",-.10),("introduced","#bd642d","移入域",.10)]:
            vv=g[g.term==k].set_index("trait").loc[["angle","chroma","outline"]]
            ax.errorbar(vv.beta,np.arange(3)+dy,xerr=[vv.beta-vv.lo,vv.hi-vv.beta],label=name,c=col,fmt="o",capsize=2)
        ax.set(title=sp.replace("Cirsium","C."),yticks=range(3),yticklabels=["向き","彩度","画像輪郭"])
        ax.axvline(0,color=".65",lw=1)
        ax.spines[["top","right"]].set_visible(False)
    axes[0].invert_yaxis()
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc="outside lower center",ncols=2,frameon=False)
    fig.supxlabel("在来で珍しい気候へ向かうときの形質勾配（共通気候・標準化係数）",fontsize=10)
    fig.savefig(OUT/"range_gradients.png",dpi=170)
    plt.close(fig)

if __name__=="__main__":main()
