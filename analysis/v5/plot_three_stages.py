"""Visual summary of the prespecified contrasts and consistency checks."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT=Path(__file__).parent/"results/three_stages"


def main():
    signed=pd.read_csv(ROOT/"v5_signed_pairs.csv")
    cross=pd.read_csv(ROOT/"v5_trait_convergence.csv").query("metric=='joint'")
    same=pd.read_csv(ROOT/"v5_within_analogue.csv")
    transfer=pd.read_csv(ROOT/"native_transfer_test.csv")
    shifts=pd.read_csv(ROOT/"native_transfer_components.csv")
    alignment=pd.read_csv(ROOT/"native_transfer_alignment.csv")
    observed=pd.read_csv(ROOT.parent/"configuration/matched_centroids.csv")
    assert len(signed)==72 and len(transfer)==18 and len(alignment)==6
    keys=["species","model","trait"]
    summed=shifts.groupby(keys).predicted_contribution.sum().sort_index()
    np.testing.assert_allclose(summed,transfer.set_index(keys).predicted.sort_index(),atol=1e-7)
    check=transfer.query("model=='annual'").merge(observed,on=["species","trait"])
    np.testing.assert_allclose(check.observed,check.novel_minus_analogous,atol=1e-8)
    assert (alignment.nonzero_pairs<=alignment.pairs).all()
    fig=plt.figure(figsize=(12,4.2),layout="constrained")
    axes=fig.subplots(1,3,width_ratios=[.9,1.15,1.2])
    for table,name,color,offset in [(cross,"Across taxa","#276884",-.07),(same,"Within taxa","#a7502e",.07)]:
        axes[0].scatter(np.arange(2)+offset,table.excess,label=name,color=color,s=55)
    axes[0].axhline(0,color=".6",lw=1)
    axes[0].set(xticks=[0,1],xticklabels=["1°","5°"],xlim=(-.45,1.45),xlabel="Geographic cell size",ylabel="Extra reduction in head distance",
                title="Resemblance in similar climates")
    axes[0].legend(frameon=False,fontsize=9)
    g=signed.query("degrees==1 and sample=='all_cells' and trait!='outline'")
    heat=g.pivot(index="axis",columns=["scale","trait"],values="beta")
    heat=heat.loc[["annual_rain","mean_temperature","radiation"],
                  [(s,t) for s in ["pair_mean","within_pair"] for t in ["angle","chroma"]]]
    ax=axes[1];ax.imshow(heat,vmin=-.65,vmax=.65,cmap="BrBG",aspect="auto")
    for (i,j),v in np.ndenumerate(heat.values):
        ax.text(j,i,f"{v:+.2f}",ha="center",va="center",color="white" if abs(v)>.4 else "black")
    ax.set(yticks=range(3),yticklabels=["Rain","Temp.","Sun"],xticks=range(4),xticklabels=["Angle","C*","Angle","C*"],
           xlabel="Taxon means      Within pair",title="Signed slopes (1° cells)")
    ax.axvline(1.5,color="white",lw=3)
    colors=["#276884","#a7502e"]
    for (sp,t),c in zip(transfer.query("model=='head_month' and trait!='outline'").groupby("species"),colors):
        t=t.set_index("trait").loc[["angle","chroma"]]
        s=observed.query("species==@sp").set_index("trait").loc[["angle","chroma"]]
        sd=s.novel_minus_analogous/s.difference_native_sd
        for field,style in [("observed","-"),("predicted","--")]:
            x,y=t[field]/sd
            axes[2].plot([0,x],[0,y],ls=style,color=c,marker="o",markevery=[1],
                         label=sp.replace("Cirsium","C.")+" "+field)
    axes[2].axvline(0,color=".75",lw=.6);axes[2].axhline(0,color=".75",lw=.6)
    axes[2].set(xlabel="Head angle shift (native SD)",ylabel="Chroma shift (native SD)",
                title="Novel minus native-like climate")
    axes[2].legend(frameon=False,fontsize=7)
    fig.savefig(ROOT/"three_stages.png",dpi=175)
    plt.close(fig)

if __name__=="__main__":
    main()

