"""Validate both representations and display the full set of niche gradients."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib import font_manager

root=Path(__file__).parent/"results"
out=root/"record_overlap"


def main():
    d=pd.read_csv(out/"record_density.csv")
    g=pd.read_csv(out/"novel_gradient.csv")
    frozen=pd.read_csv(root/"whole_range_overlap/whole_head_overlap.csv")
    keys=["species","dimensions","h","group"]
    median=d.query("limits=='original' and representation=='medians' and dimensions>=2").set_index(keys).D.sort_index()
    frozen=frozen.rename(columns={"bandwidth":"h"})
    frozen["group"]=frozen.group.str.replace("introduced_","")
    np.testing.assert_allclose(median,frozen.set_index(keys).D.sort_index(),atol=1e-9)
    assert len(g)==288 and len(d)==360
    assert (d.D.between(0,1)).all()
    slopes=g.pivot_table(index=["species","regional","climate","axis","trait"],columns="term",values="beta")
    np.testing.assert_allclose(slopes.novel-slopes.analogue,slopes.difference,atol=1e-10)
    show=d.query("limits=='original' and h==.5 and dimensions==3")
    primary=g.query("regional and term=='novel' and climate=='niche'")
    for condition,name in [(show,"primary_overlap"),(primary,"primary_slopes")]:
        condition.to_csv(out/(name+".csv"),index=False)
    font="Meiryo" if "Meiryo" in font_manager.get_font_names() else "Noto Sans CJK JP"
    plt.rcParams.update({"font.family":[font,"DejaVu Sans"],"font.size":9,"pdf.fonttype":42})
    fig,axs=plt.subplots(2,2,figsize=(9.4,5.8),layout="constrained",width_ratios=[1,1.4])
    axes=["bio01","bio04","bio12","bio15","rsds_mean"]
    vmax=max(abs(primary.beta))
    norm=TwoSlopeNorm(vmin=-vmax,vcenter=0,vmax=vmax)
    for row,(sp,data) in enumerate(show.groupby("species",sort=True)):
        left,right=axs[row]
        for k,(name,label,color) in enumerate([("medians","セル中央値","#858585"),("records","セル内全記録","#14716d")]):
            values=data.query("representation==@name").set_index("group").loc[["all","analogue","novel"]]
            left.scatter(np.arange(3)+(k-.5)*.18,values.D,s=65,c=color,label=label)
        left.set(xticks=range(3),xticklabels=["全移入","類似気候","新規気候"],ylim=(0,1),ylabel="在来頭花との重複度 D",
                 title=sp.replace("Cirsium","C."))
        left.spines[["top","right"]].set_visible(False)
        if row==0:left.legend(frameon=False,fontsize=8,loc="lower left")
        heat=primary[primary.species==sp].pivot(index="trait",columns="axis",values="beta").loc[["angle","chroma","outline"],axes]
        right.imshow(heat,cmap="BrBG",norm=norm,aspect="auto")
        for (i,j),value in np.ndenumerate(heat.to_numpy()):
            right.text(j,i,f"{value:+.2f}",ha="center",va="center",color="white" if abs(value)>.6*vmax else "black")
        right.set(xticks=range(5),xticklabels=["平均気温","気温季節性","年降水","降水季節性","日射"],
                  yticks=range(3),yticklabels=["向き","彩度","輪郭"],title="新規気候内の地域差調整後の標準化係数")
    fig.savefig(out/"overlap_and_climate.png",dpi=180)
    fig.savefig(out/"overlap_and_climate.pdf")
    plt.close(fig)

if __name__=="__main__":main()
