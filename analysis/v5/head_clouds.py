"""Draw the native reference and exact introduced cell pairs used in S2."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from analysis.v5.novel_head_overlap import T

GROUPS = ["Native", "Introduced / native-like climate", "Introduced / novel climate"]
COLORS = ["#85939d", "#14768a", "#ce562e"]
AXES = ["Head angle (degrees)", "Corolla chroma (C*)", "Outline index"]


def draw(clouds, folder, observations=False):
    species = sorted(clouds)
    points = np.vstack([group for sp in species for group in clouds[sp]])
    assert np.isfinite(points).all() and len(species) == 2
    lo, hi = points.min(0), points.max(0)
    pad = .04*(hi-lo)
    limits = list(zip(lo-pad, hi+pad))
    interactive = make_subplots(rows=1, cols=2, specs=[[{"type":"scene"}]*2],
                                subplot_titles=[s.replace("Cirsium","C.") for s in species])
    fig = plt.figure(figsize=(11.5,8.1),layout="constrained")
    layout = fig.add_gridspec(3,2,height_ratios=[.16,1,.68])
    legend = fig.add_subplot(layout[0,:])
    legend.axis("off")
    legend.legend(handles=[Line2D([],[],marker="o",ls="",color=c,label=g) for g,c in zip(GROUPS,COLORS)],
                  loc="center",ncols=3,frameon=False)
    for i,sp in enumerate(species):
        ax = fig.add_subplot(layout[1,i],projection="3d")
        proj = fig.add_subplot(layout[2,i])
        groups = clouds[sp]
        for j, (values, color, label) in enumerate(zip(groups,COLORS,GROUPS)):
            x,y,z = values.T
            ax.scatter(x,y,z,s=(1.5 if observations else 7) if j==0 else (4 if observations else 17),color=color,alpha=(.035 if observations else .15) if j==0 else (.14 if observations else .73),
                       edgecolors="none",depthshade=False)
            proj.scatter(x,y,s=(1.5 if observations else 9) if j==0 else (4 if observations else 21),color=color,alpha=(.035 if observations else .17) if j==0 else (.14 if observations else .7),edgecolors="none")
            if not observations:
                interactive.add_trace(go.Scatter3d(x=x,y=y,z=z,mode="markers",name=label,legendgroup=label,
                showlegend=i==0,marker=dict(size=2.2 if j==0 else 4,color=color,opacity=.18 if j==0 else .8),
                hovertemplate="Angle %{x:.1f}°<br>C* %{y:.2f}<br>Outline %{z:.3f}<extra>%{fullData.name}</extra>"),
                row=1,col=i+1)
        ax.set(xlabel=AXES[0],ylabel=AXES[1],zlabel=AXES[2],xlim=limits[0],ylim=limits[1],zlim=limits[2])
        ax.view_init(elev=22,azim=-57)
        ax.set_title(sp.replace("Cirsium","C.")+"\n"+" / ".join(str(len(g)) for g in groups)+(" observations" if observations else " cells"),style="italic")
        proj.set(xlim=limits[0],ylim=limits[1],xlabel=AXES[0],ylabel=AXES[1])
        proj.spines[["top","right"]].set_visible(False)
        scene = dict(aspectmode="cube")
        for key, label, lim in zip(["xaxis","yaxis","zaxis"],AXES,limits):
            scene[key]=dict(title=label,range=list(lim))
        interactive.update_scenes(**scene,row=1,col=i+1)
    fig.savefig(folder/("observed_heads.png" if observations else "head_clouds.png"),dpi=175,bbox_inches="tight")
    plt.close(fig)
    if observations:
        return
    interactive.update_layout(height=640,margin=dict(l=10,r=10,t=80,b=20),
        title_text="Head traits across climates · one point per 1° cell",
        legend=dict(orientation="h",x=.5,xanchor="center",y=-.02))
    interactive.write_html(folder/"head_clouds.html",include_plotlyjs=True,config={"displaylogo":False,"responsive":True})

def observations_in_cells(frame, cells):
    valid = frame.loc[frame[T[0]].between(0,180)].dropna(subset=T).copy()
    valid["i"] = np.minimum(179,np.floor(valid.latitude+90)).astype(int)
    valid["j"] = np.floor((valid.longitude+180)%360).astype(int)
    selected = valid.merge(cells[["i","j"]],on=["i","j"],validate="many_to_one")
    assert len(selected) == cells.photos.sum()
    return selected[T].to_numpy()
