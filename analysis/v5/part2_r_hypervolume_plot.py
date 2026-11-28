"""Draw only model-generated set-operation points, not observations.
All overlap metrics are from the full fitted dimensions, not the projections.
The four priority fits have exactly the original R primary settings.
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT=Path("work/submission")
OUT=ROOT/"part2_r_hypervolume"
SPECIES=("Cirsium arvense","Cirsium vulgare")
COMPONENTS=("Unique_1","Intersection","Unique_2")
COLORS={"Unique_1":"#568B71","Intersection":"#89909C","Unique_2":"#934C9D"}
VOLCOL={"Unique_1":"reference_unique","Intersection":"intersection","Unique_2":"comparison_unique"}

def inputs():
    stats=pd.read_csv(OUT/"comparisons.csv")
    pts=pd.read_csv(OUT/"modeled_set_points.csv")
    prim=stats.loc[stats.scenario=="primary"].copy()
    assert len(prim)==4 and not prim.duplicated(["species","space","design"]).any()
    assert set(prim.species)==set(SPECIES)
    assert (prim.bandwidth==.5).all() and (prim.coverage==.95).all()
    for a,b in [("reference_volume","reference_unique"),("comparison_volume","comparison_unique")]:
        np.testing.assert_allclose(prim[a],prim[b]+prim.intersection,rtol=1e-7,atol=1e-7)
    np.testing.assert_allclose(prim.union,prim.reference_volume+prim.comparison_volume-prim.intersection,rtol=1e-7,atol=1e-7)
    np.testing.assert_allclose(prim.jaccard,prim.intersection/prim.union,rtol=1e-8)
    np.testing.assert_allclose(prim.sorensen,2*prim.intersection/(prim.reference_volume+prim.comparison_volume),rtol=1e-8)
    assert set(pts.component).issubset(COMPONENTS)
    return prim,pts

def select(stats,pts,sp,space,design):
    s=stats[(stats.species==sp)&(stats.space==space)&(stats.design==design)]
    assert len(s)==1
    r=s.iloc[0]
    p=pts[(pts.species==sp)&(pts.space==space)&(pts.design==design)]
    assert len(p)>0
    cols=[f"x{i+1}" for i in range(int(r.dimensions))]
    assert np.isfinite(p[cols].to_numpy()).all()
    # The R export caps each set independently. Rebalance to approximately
    # equal sampling density per unit volume, so small unique sets are not
    # visually over-represented. The plot is not an observation-density plot.
    rates=[]
    for k in COMPONENTS:
        n=int((p.component==k).sum());v=float(r[VOLCOL[k]])
        if v>0:
            assert n>0
            rates.append(n/v)
        else:
            assert n==0
    rate=min(rates)
    rng=np.random.default_rng(20261128)
    result={}
    for k in COMPONENTS:
        z=p.loc[p.component==k,cols].to_numpy()
        n=len(z)
        if not n:continue
        take=min(n,max(1,int(np.floor(rate*float(r[VOLCOL[k]])))))
        result[k]=z[rng.choice(n,take,replace=False)]
    return r,result

def legend(fig,labels):
    h=[Line2D([0],[0],marker="o",linestyle="",color=COLORS[k],markersize=6,label=label)
       for k,label in zip(COMPONENTS,labels)]
    fig.legend(handles=h,loc="lower center",bbox_to_anchor=(.5,.045),ncol=3,frameon=False,fontsize=9)

def save(fig,name):
    fig.savefig(OUT/(name+".png"),dpi=320,bbox_inches="tight",facecolor="white")
    fig.savefig(OUT/(name+".svg"),bbox_inches="tight",facecolor="white")
    plt.close(fig)

def main():
    stats,pts=inputs()
    # Original cell-level file is read only in the private runner to fit the
    # climate projection. Neither the file nor its rows are plotted/exported.
    private=pd.read_csv(ROOT/"private/r_hypervolume_inputs.csv")
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"axes.spines.top":False,
                         "axes.spines.right":False,"svg.fonttype":"none"})
    fig=plt.figure(figsize=(10.6,6.1))
    fig.subplots_adjust(left=.04,right=.97,top=.83,bottom=.20,wspace=.22)
    rows=[]
    trait_sets=[]
    for sp in SPECIES:
        trait_sets.append(select(stats,pts,sp,"trait","native_novel"))
    full=np.vstack([a for _,m in trait_sets for a in m.values()])
    # Common 3D native-SD axes, no truncation or ellipsoid approximation.
    lo=full.min(0);hi=full.max(0);pad=(hi-lo)*.04
    for k,(sp,(r,m)) in enumerate(zip(SPECIES,trait_sets)):
        ax=fig.add_subplot(1,2,k+1,projection="3d")
        for c in COMPONENTS:
            if c in m:
                a=m[c]
                ax.scatter(a[:,0],a[:,1],a[:,2],s=2.5,alpha=.36,
                           color=COLORS[c],linewidths=0,depthshade=False,rasterized=True)
                rows.append(dict(species=sp,space="trait",component=c,plotted_points=len(a)))
        ax.set(xlim=(lo[0]-pad[0],hi[0]+pad[0]),
               ylim=(lo[1]-pad[1],hi[1]+pad[1]),
               zlim=(lo[2]-pad[2],hi[2]+pad[2]),
               xlabel="Head angle",ylabel="Floral chroma",zlabel="Outline")
        ax.set_box_aspect(hi-lo+2*pad)
        ax.view_init(elev=22,azim=-57)
        ax.set_title(sp,fontstyle="italic",fontsize=13,pad=6)
        shared=r.intersection/r.comparison_volume
        ax.text2D(.5,-.10,f"Jaccard = {r.jaccard:.3f}   |   Shared / novel volume = {shared:.1%}\n"
                  f"Native n = {int(r.reference_n):,}; novel-climate n = {int(r.comparison_n):,}",
                  transform=ax.transAxes,ha="center",va="top",fontsize=9)
    fig.suptitle("Three-dimensional capitulum hypervolumes",fontsize=17,fontweight="bold",y=.98)
    fig.text(.5,.914,"Gaussian hypervolume, 95% probability region; native vs introduced novel climates",ha="center",fontsize=10)
    fig.text(.5,.105,"Axes: species-specific native standard deviations. Points sampled from estimated set volumes, not individual heads.",
             ha="center",fontsize=8)
    legend(fig,("Native only","Intersection","Novel-climate only"))
    save(fig,"Fig_trait_hypervolume_native_novel")
    fig,axs=plt.subplots(1,2,figsize=(10.3,5.0))
    fig.subplots_adjust(left=.085,right=.98,top=.79,bottom=.23,wspace=.30)
    for ax,sp in zip(axs,SPECIES):
        r,m=select(stats,pts,sp,"climate_raw","native_all")
        x=private.loc[(private.species==sp)&(private.space=="climate_raw")&
                      (private.design=="native_all")&(private.group=="native"),
                      [f"x{i}" for i in range(1,6)]].to_numpy()
        assert len(x)==int(r.reference_n)
        center=x.mean(0)
        _,s,vt=np.linalg.svd(x-center,full_matrices=False)
        basis=vt[:2].copy()
        for v in basis:
            if v[np.argmax(np.abs(v))]<0:v*=-1
        share=s[:2]**2/(s*s).sum()
        for c in COMPONENTS:
            if c in m:
                a=(m[c]-center)@basis.T
                ax.scatter(a[:,0],a[:,1],s=3,alpha=.32,color=COLORS[c],linewidths=0,rasterized=True)
                rows.append(dict(species=sp,space="climate_raw",component=c,plotted_points=len(a)))
        ax.set(xlabel=f"Native climate PC1 ({share[0]:.0%})",ylabel=f"Native climate PC2 ({share[1]:.0%})")
        ax.set_aspect("equal",adjustable="datalim")
        ax.set_title(sp,fontstyle="italic",fontsize=13)
        ax.text(.02,.98,f"5-D volume ratio = {r.volume_ratio:.2f}\n5-D Jaccard = {r.jaccard:.3f}",
                transform=ax.transAxes,va="top",fontsize=9,
                bbox=dict(facecolor="white",edgecolor="none",alpha=.9,pad=3))
    fig.suptitle("Five-dimensional occupied-climate hypervolumes",fontsize=16,fontweight="bold",y=.98)
    fig.text(.5,.90,"2-D projection of 5-D modeled sets; raw occupied climates, not background-corrected",ha="center",fontsize=10)
    fig.text(.5,.115,"Colors describe membership in the original 5-D sets; projected overlap is not a 5-D overlap estimate.",
             ha="center",fontsize=8)
    legend(fig,("Native only","Intersection","Introduced only"))
    save(fig,"Fig_climate_hypervolume_raw_projection")
    pd.DataFrame(rows).to_csv(OUT/"plot_point_counts.csv",index=False)
    print("Audited set-volume identities and rendered actual package-generated point clouds.",flush=True)

if __name__=="__main__":
    main()
