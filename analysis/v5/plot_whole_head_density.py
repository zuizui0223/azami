"""Full range overlap, with the same native head reference for introduced novelty."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

root=Path(__file__).parent/"results"
out=root/"whole_range_overlap"
species=["Cirsium arvense","Cirsium vulgare"]
def mass_mask(z,coverage):
    sorted_z=np.sort(z.ravel())[::-1]
    return z>=sorted_z[np.searchsorted(sorted_z.cumsum(),coverage*z.sum())]


def panel(ax,x,y,a,b):
    masks=mass_mask(a,.95).astype(int)+2*mass_mask(b,.95).astype(int)
    ax.pcolormesh(x,y,masks.T,cmap=ListedColormap(["white","#8dc7cf","#f1b786","#827c9c"]),
                  vmin=-.5,vmax=3.5,shading="flat")
    for z,color in [(a,"#256675"),(b,"#a95024")]:
        ax.contour((x[1:]+x[:-1])/2,(y[1:]+y[:-1])/2,mass_mask(z,.5).T.astype(float),
                   levels=[.5],colors=[color],linewidths=1.2)
    i,j=np.where(masks>0)
    ax.set_xlim(x[max(i.min()-2,0)],x[min(i.max()+3,len(x)-1)])
    ax.set_ylim(y[max(j.min()-2,0)],y[min(j.max()+3,len(y)-1)])
    ax.spines[["top","right"]].set_visible(False)


def main():
    font="Meiryo" if "Meiryo" in font_manager.get_font_names() else "Noto Sans CJK JP"
    plt.rcParams.update({"font.family":[font,"DejaVu Sans"],"font.size":9,"pdf.fonttype":42})
    table=pd.read_csv(out/"whole_head_overlap.csv").query("dimensions==3 and bandwidth==.5")
    fig,axes=plt.subplots(2,3,figsize=(11,6.7),layout="constrained")
    for row,sp in enumerate(species):
        key=sp.replace(" ","_")
        climate=np.load(root/"niche_overlap"/(key+".npz"))
        ax=axes[row,0]
        panel(ax,climate["x"],climate["y"],climate["native"],climate["introduced"])
        D=np.minimum(climate["native"],climate["introduced"]).sum()
        ax.set(title=f"{sp.replace('Cirsium','C.')}   気候 D = {D:.3f}",xlabel="気候 PC1",ylabel="気候 PC2")
        head=np.load(out/(key+"_head_densities.npz"))
        x,y=head["edges_0"],head["edges_1"]
        n=head["density_0"]
        for col,idx,group in [(1,1,"introduced_all"),(2,3,"introduced_novel")]:
            b=head["density_"+str(idx)]
            value=table.query("species==@sp and group==@group").D.item()
            np.testing.assert_allclose(value,np.minimum(n,b).sum())
            ax=axes[row,col]
            panel(ax,x,y,n.sum(2),b.sum(2))
            ax.set(title=("全移入域" if col==1 else "移入域の新規気候")+f"   頭花 D = {value:.3f}",
                   xlabel="頭花の角度（°）",ylabel="花冠の彩度 C*")
        xlim=axes[row,1].get_xlim();ylim=axes[row,1].get_ylim()
        xlim2=axes[row,2].get_xlim();ylim2=axes[row,2].get_ylim()
        for ax in axes[row,1:]:
            ax.set(xlim=(min(xlim[0],xlim2[0]),max(xlim[1],xlim2[1])),ylim=(min(ylim[0],ylim2[0]),max(ylim[1],ylim2[1])))
    fig.legend(handles=[Patch(color=c,label=l) for c,l in [("#8dc7cf","在来域のみ"),("#f1b786","移入側のみ"),("#827c9c","重なり")]],
               loc="outside lower center",ncols=3,frameon=False)
    for ext in ["png","pdf"]:
        fig.savefig(out/("range_overlap."+ext),dpi=175)
    plt.close(fig)
if __name__=="__main__":main()
