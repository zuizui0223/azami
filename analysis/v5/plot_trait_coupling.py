"""Check aligned region pairs and plot overlap against climate rarity."""
from pathlib import Path
import pandas as pd
import numpy as np
from itertools import permutations
from scipy.stats import rankdata
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
root=Path(__file__).parent/"results"
out=root/"continuous_overlap"

def main():
    pairs=pd.read_csv(out/"trait_regional_pairs.csv")
    results=pd.read_csv(out/"trait_coupling.csv")
    primary=pd.read_csv(root/"paired_overlap/pairs.csv").query("shape")
    keys=["species","native_block","intro_block"]
    checked=pairs.query("design=='regional'").rename(columns={"native_region":"native_block","intro_region":"intro_block"}).merge(primary,on=keys,validate="one_to_one")
    assert len(checked)==213
    np.testing.assert_allclose(checked.climate_D_x,checked.climate_D_y,atol=1e-12)
    assert len(results)==108*1 and len(pairs)==228
    assert pairs.native_climate_rarity.between(0,1).all()
    font="Meiryo" if "Meiryo" in font_manager.get_font_names() else "Noto Sans CJK JP"
    plt.rcParams["font.family"]=[font,"DejaVu Sans"]
    fig,axes=plt.subplots(2,2,figsize=(9.6,6.4),layout="constrained")
    exact=[]
    for sp,group in pairs.query("design=='all_native_reference'").groupby("species"):
        indices=np.array(list(permutations(range(len(group)))))
        for test in ["overlap","novelty_direction","rarity_overlap","rarity_direction"]:
            x=rankdata(group.native_climate_rarity if test.startswith("rarity") else 1-group.climate_D)
            x-=x.mean()
            for trait in ["angle","chroma","outline"]:
                y=rankdata(1-group[trait+"_D"] if test.endswith("overlap") else group[trait+"_shift"])
                y-=y.mean()
                rho=(x@y)/(np.linalg.norm(x)*np.linalg.norm(y))
                prob=np.mean(np.abs(y[indices]@x)>=abs(y@x)-1e-9)
                exact.append(dict(species=sp,test=test,trait=trait,rho=rho,p_exact=prob))
    summary=pd.DataFrame(exact)
    summary["holm_exact"]=summary.groupby("test").p_exact.transform(lambda v:multipletests(v,method="holm")[1])
    summary.to_csv(out/"all_native_exact.csv",index=False)
    for row,(sp,data) in enumerate(checked.groupby("species",sort=True)):
        ax=axes[row,0]
        ax.scatter(1-data.climate_D_x,1-data.head_D,s=14,color="#287b78",alpha=.5)
        rho=data.climate_D_x.corr(data.head_D,method="spearman")
        ax.set(xlim=(0,1),ylim=(0,1),xlabel="気候分布のずれ（1−D）",ylabel="頭花分布のずれ（1−D）",title=sp.replace("Cirsium","C.")+f"   ρ={rho:.2f}")
        ref=pairs.query("species==@sp and design=='all_native_reference'")
        ax=axes[row,1]
        for key,label,color,marker in [("angle","向き","#bb5d2d","o"),("chroma","彩度","#26666e","s"),("outline","輪郭","#806787","^")]:
            ax.scatter(ref.native_climate_rarity,1-ref[key+"_D"],label=label,color=color,marker=marker,s=25)
        ax.axhline(0,color=".65",lw=.9)
        ax.set(xlim=(0,1),ylim=(0,1),xlabel="在来域を基準とした移入地域の気候希少度",ylabel="個別形質分布のずれ（1−D）",title=f"在来全域と移入{len(ref)}地域")
        if row==0:ax.legend(loc="best",frameon=False,ncol=3,fontsize=8)
    fig.savefig(out/"overlap_and_rarity.png",dpi=180)
    plt.close(fig)

if __name__=="__main__":main()
