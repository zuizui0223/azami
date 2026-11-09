"""Check both directional partitions, then draw the introduced-range comparison."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

out=Path(__file__).parent/"results/invasion_design"


def main():
    data=pd.read_csv(out/"availability_partition.csv")
    heads=pd.read_csv(out/"availability_heads.csv")
    balanced=pd.read_csv(out/"climate_balanced_heads.csv")
    columns=["shared_core","shared_outside_other_core","outside_other_background","outside_own_background"]
    np.testing.assert_allclose(data[columns].sum(axis=1),1,atol=1e-9)
    assert len(data)==36 and len(balanced)==20
    assert heads.D.dropna().between(0,1).all() and balanced.D_balanced.between(0,1).all()
    primary=data.query("background==.99 and occ_core==.95 and side=='introduced'").set_index("species")
    font="Meiryo" if "Meiryo" in font_manager.get_font_names() else "Noto Sans CJK JP"
    plt.rcParams.update({"font.family":[font,"DejaVu Sans"],"font.size":9})
    fig,(ax,bx)=plt.subplots(2,1,figsize=(8.5,5.8),layout="constrained")
    left=np.zeros(2)
    for k,label,color in zip(columns,["在来の占有中心域","共通背景で在来の占有中心域外","在来の主要背景外","移入の主要背景外"],
                              ["#488785","#dfa251","#96657d","#bababa"]):
        v=100*primary[k].to_numpy()
        ax.barh([1,0],v,left=left,color=color,label=label,height=.5)
        for i,val in enumerate(v):
            if val>8:ax.text(left[i]+val/2,1-i,f"{val:.0f}%",va="center",ha="center")
        left+=v
    names=[s.replace("Cirsium","C.") for s in primary.index]
    ax.set(yticks=[1,0],yticklabels=names,xlim=(0,100),xlabel="移入側の背景補正占有密度（%）")
    ax.legend(ncol=2,loc="upper center",bbox_to_anchor=(.5,-.35),frameon=False)
    for j,(key,label,color) in enumerate(zip(["occupied_native_core","expansion_in_common_climate","outside_native_background"],
          ["在来の占有中心域","共通背景で在来の占有中心域外","在来の主要背景外"],["#488785","#ce852c","#96657d"])):
        g=heads.query("trait=='joint' and group==@key").set_index("species").loc[primary.index]
        bx.scatter(g.D,np.array([1,0])+(1-j)*.18,color=color,s=54)
        for i,r in enumerate(g.itertuples()):
            if np.isfinite(r.D):bx.text(r.D+.012,1-i+(1-j)*.18,f"{r.D:.2f}（{r.cells}セル）",va="center",fontsize=8)
    bx.set(yticks=[1,0],yticklabels=names,ylim=(-.4,1.4),xlim=(0,1.1),xlabel="在来域全体との頭花分布の重複度 D")
    for a in (ax,bx):a.spines[["top","right"]].set_visible(False)
    fig.savefig(out/"availability_and_heads.png",dpi=170)
    plt.close(fig)

if __name__=="__main__":main()
