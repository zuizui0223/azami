from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

root=Path(__file__).parent/"results"
out=root/"nonoverlap"
def main():
    data=pd.read_csv(out/"excess_coupling.csv")
    d=data.query("h==.5 and side=='introduced'")
    joint=d.query("trait=='joint_PC12'").merge(pd.read_csv(root/"paired_overlap/whole.csv").query("shape and h==.5"),on="species")
    np.testing.assert_allclose(joint.climate_D_x,joint.climate_D_y)
    np.testing.assert_allclose(joint.trait_D,joint.head_D)
    np.testing.assert_allclose(data.head_excess_in_climate_excess-data.head_excess_in_climate_overlap,data.difference,atol=1e-10)
    assert len(data)==60 and data.head_excess_in_climate_excess.between(0,1).all()
    font="Meiryo" if "Meiryo" in font_manager.get_font_names() else "Noto Sans CJK JP"
    plt.rcParams["font.family"]=[font,"DejaVu Sans"]
    fig,axes=plt.subplots(1,2,figsize=(8.5,3.5),layout="constrained",sharey=True)
    for ax,(sp,g) in zip(axes,d.groupby("species")):
        v=g.set_index("trait").loc[["joint_PC12","angle","chroma","outline"]]
        for j,(key,label,color) in enumerate([("head_excess_in_climate_overlap","気候の重なり部分","#467e96"),
                                                ("head_excess_in_climate_excess","移入側の気候超過部分","#c57a3f")]):
            ax.bar(np.arange(4)+(j-.5)*.35,100*v[key],width=.35,label=label,color=color)
        ax.set(title=sp.replace("Cirsium","C."),xticks=range(4),xticklabels=["頭花全体","向き","彩度","輪郭"])
        ax.spines[["right","top"]].set_visible(False)
    axes[0].set_ylabel("頭花が移入側の超過部分に属する重み（%）")
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc="outside lower center",ncols=2,frameon=False)
    fig.savefig(out/"excess_overlap.png",dpi=175)
    plt.close(fig)
if __name__=="__main__":main()
