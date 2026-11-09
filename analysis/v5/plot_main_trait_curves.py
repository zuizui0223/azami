"""Plot observed climate-bin means and marginal range-specific fits from V5."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent / "results"
DEST = ROOT / "expansion_consistency"


def main():
    df = pd.read_csv(DEST / "main_rarity_curves.csv")
    gradients = pd.read_csv(DEST / "range_gradient.csv").query(
        "coarsen==4 and model=='occupancy_rank' and term!='difference'")
    plt.rcParams.update({"font.family": "Meiryo", "pdf.fonttype": 42,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(2, 3, figsize=(10.5, 6.3), sharex=True, sharey="col",
                             layout="constrained")
    for row, species in enumerate(["Cirsium arvense","Cirsium vulgare"]):
        for col,(trait,label) in enumerate([("angle","下向き角度（°）"),("chroma","花冠の彩度 C*"),
                                             ("outline","頭花輪郭の縦横比")]):
            ax = axes[row,col]
            for k, term, name, color in [(0,"native","在来域","#326798"),(1,"introduced","移入域","#bf532c")]:
                d = df.query("species==@species and trait==@trait and introduced==@k").sort_values("rarity")
                assert len(d)>=2 and d.cells.min()>=5 and d.rarity.between(0,1).all()
                g = gradients.query("species==@species and trait==@trait and term==@term").iloc[0]
                np.testing.assert_allclose(np.diff(d.fitted)/np.diff(d.rarity),
                                           g.beta*g.native_trait_sd/g.axis_sd,rtol=1e-7)
                ax.plot(d.rarity,d.fitted,color=color,label=name,linewidth=1.9)
                ax.fill_between(d.rarity,d.lo,d.hi,color=color,alpha=.12,lw=0)
                ax.scatter(d.rarity,d.observed,color=color,s=12,marker="o" if k==0 else "^")
            ax.set_xlim(0,1.01)
            ax.grid(axis="y",color="#e6e6e6",linewidth=.7)
            if row==0:
                ax.set_title(label)
            if col==0:
                ax.set_ylabel(species.replace("Cirsium","C."),fontsize=12,fontstyle="italic")
    axes[0,0].legend(frameon=False,fontsize=9)
    fig.suptitle("在来域で占有がまれな気候に向かうと、頭花はどう変わるか",fontsize=16)
    fig.supxlabel("在来域の気候占有の希少度（0＝よく占める、1＝まれ）\n点：気候階級ごとの加重平均　　線：分布域ごとの推定平均　　帯：点別95%信頼区間",fontsize=10)
    for ext in ["png","pdf"]:
        fig.savefig(DEST/f"main_trait_curves.{ext}",dpi=180)
    plt.close(fig)


if __name__=="__main__":
    main()
