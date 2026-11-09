"""Check common-cohort comparisons and plot temperature slopes for all endpoints."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

root=Path(__file__).parent/"results"
out=root/"adjusted_exposures"
names={"month_mean":"観察月・平均","month_max":"観察月・最高","previous_max":"前月・最高"}
def main():
    d=pd.read_csv(out/"adjusted_slopes.csv")
    assert len(d)==756
    assert d.groupby("species").novel_cells.nunique().eq(1).all()
    assert d.groupby("species").native_like_cells.nunique().eq(1).all()
    old=pd.read_csv(root/"record_overlap/novel_gradient.csv").query("climate=='head_month'")
    keys=["species","trait","axis","term","regional"]
    chk=d.query("climate=='month_mean'").merge(old,on=keys,validate="one_to_one")
    assert len(chk)==108
    for name in ["beta","lo","hi","holm_model"]:
        if name=="holm_model":
            np.testing.assert_allclose(chk.holm_model,chk.holm,atol=1e-6)
        else:
            np.testing.assert_allclose(chk[name+"_x"],chk[name+"_y"],atol=1e-6)
    pivot=d.pivot(index=["species","climate","regional","trait","axis"],columns="term",values="beta")
    np.testing.assert_allclose(pivot.novel-pivot.analogue,pivot.difference,atol=1e-9)
    temp=d[d.climate.isin(names)&d.axis.isin(["tavg","tmax","prior_tmax"])&~d.regional]
    temp.to_csv(out/"thermal_contrasts.csv",index=False)
    font="Meiryo" if "Meiryo" in font_manager.get_font_names() else "Noto Sans CJK JP"
    plt.rcParams.update({"font.family":[font,"DejaVu Sans"],"font.size":9})
    fig,axs=plt.subplots(2,3,figsize=(10,4.8),layout="constrained",sharex=True,sharey=True)
    for row,(sp,g) in enumerate(temp.groupby("species")):
        for col,(trait,label) in enumerate([("angle","向き"),("chroma","彩度 C*"),("outline","画像輪郭")]):
            ax=axs[row,col]
            for t,offset,color,tag in [("analogue",-.1,"#698a98","在来類似気候"),("novel",.1,"#bd632e","新規気候")]:
                v=g[(g.trait==trait)&(g.term==t)].set_index("climate").loc[list(names)]
                ax.errorbar(v.beta,np.arange(3)+offset,xerr=[v.beta-v.lo,v.hi-v.beta],fmt="o",color=color,capsize=2,label=tag)
            ax.axvline(0,color=".65",lw=.8)
            ax.set(yticks=range(3),yticklabels=list(names.values()),title=sp.replace("Cirsium","C.")+" · "+label)
            ax.spines[["top","right"]].set_visible(False)
    axs[0,0].invert_yaxis()
    handles,labels=axs[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,ncol=2,loc="outside lower center",frameon=False)
    fig.supxlabel("同じ降水・日射条件での温度勾配（標準化係数・点別95%区間）",fontsize=10)
    fig.savefig(out/"thermal_contrasts.png",dpi=170)
    plt.close(fig)
if __name__=="__main__":main()
