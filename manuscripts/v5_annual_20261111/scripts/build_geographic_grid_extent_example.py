#!/usr/bin/env python3
"""Publication/poster schematic of an explicit 1° vs 5° angular geographic extent.

Illustrative Japanese grid bounds, *not* a map of sampled V5 observations.
The underlying climate raster is CHELSA 30 arc-seconds; V5 aggregates photo-
linked values to taxon × 1°/5° cells. East–west span varies with latitude.

Run: python manuscripts/v5_annual_20261111/scripts/build_geographic_grid_extent_example.py
Requires: matplotlib, basemap, pyproj.
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from mpl_toolkits.basemap import Basemap
from pyproj import Geod

plt.rcParams.update({
    "font.family": "Noto Sans CJK JP", "axes.unicode_minus": False,
    "pdf.fonttype": 42, "svg.fonttype": "none"
})
NAVY, TEAL, PURPLE = "#173754", "#008e83", "#714aa0"
OUT = Path("manuscripts/v5_annual_20261111/figures/generated")
g = Geod(ellps="WGS84")
def distance(a,b,c,d): return g.inv(a,b,c,d)[2]/1000

def main():
    distances = {
        "one_north_south_km": distance(135,35,135,36),
        "one_east_west_km": distance(135,35.5,136,35.5),
        "five_north_south_km": distance(135,35,135,40),
        "five_east_west_km": distance(135,37.5,140,37.5)
    }
    assert round(distances["one_north_south_km"],1)==110.9
    assert round(distances["one_east_west_km"],1)==90.7
    assert round(distances["five_north_south_km"],1)==554.9
    assert round(distances["five_east_west_km"],1)==442.1

    fig=plt.figure(figsize=(15.3,5.4),facecolor="white")
    a=fig.add_axes([.035,.15,.29,.72])
    b=fig.add_axes([.37,.15,.37,.72])
    labels=fig.add_axes([.76,.16,.235,.72])
    labels.axis("off")
    for ax,bounds in ((a,(128,30,145,45)),(b,(134.9,34.9,140.1,40.1))):
        x1,y1,x2,y2=bounds
        m=Basemap(projection="cyl",llcrnrlon=x1,llcrnrlat=y1,
                  urcrnrlon=x2,urcrnrlat=y2,resolution="i",ax=ax)
        m.fillcontinents(color="#eaf1ef",lake_color="white",zorder=0)
        m.drawcoastlines(color="#98a9b2",linewidth=.6,zorder=1)
        ax.add_patch(Rectangle((135,35),5,5,fill=False,
                               edgecolor=PURPLE,linewidth=3,zorder=5))
    a.set_title("日本列島での位置（例示）",loc="left",fontsize=14,fontweight="bold",color=NAVY)
    a.set_xticks([130,135,140,145],[f"{x}°E" for x in [130,135,140,145]])
    a.set_yticks([30,35,40,45],[f"{x}°N" for x in [30,35,40,45]])
    for n in range(135,141): b.axvline(n,color=TEAL,linewidth=1,alpha=.7)
    for n in range(35,41): b.axhline(n,color=TEAL,linewidth=1,alpha=.7)
    b.add_patch(Rectangle((135,35),1,1,facecolor="#9bd7cf",alpha=.45,
                          edgecolor=TEAL,linewidth=2.5,zorder=6))
    b.set_title("5°枠の内部：1°格子と座標境界",loc="left",fontsize=14,fontweight="bold",color=NAVY)
    b.set_xticks(list(range(135,141)),[f"{x}°E" for x in range(135,141)])
    b.set_yticks(list(range(35,41)),[f"{x}°N" for x in range(35,41)])
    b.annotate("1°例：35–36°N × 135–136°E",xy=(135.5,35.5),
               xytext=(137.3,39.2),ha="center",fontsize=10,color=TEAL,
               arrowprops=dict(arrowstyle="->",color=TEAL,linewidth=2),
               bbox=dict(boxstyle="round,pad=.3",fc="white",ec="#acddd5"))
    for pos,label,col,size,strong in [
        (.97,"緯度経度で示す具体範囲",NAVY,14,True),
        (.80,"1°セル（青緑）",TEAL,13,True),
        (.69,"35–36°N / 135–136°E",NAVY,10.5,False),
        (.59,"南北 110.9 km",NAVY,12,False),
        (.50,"東西 90.7 km",NAVY,12,False),
        (.39,"5°セル（紫）",PURPLE,13,True),
        (.28,"35–40°N / 135–140°E",NAVY,10.5,False),
        (.18,"南北 554.9 km",NAVY,12,False),
        (.09,"東西 442.1 km",NAVY,12,False),
    ]:
        labels.text(.01,pos,label,va="top",transform=labels.transAxes,
                    color=col,fontsize=size,fontweight="bold" if strong else "normal")
    fig.text(.02,.03,"地理範囲の例示であり、観測点の表示ではない。東西幅は緯度依存。"
             "V5は撮影地点に対応する原CHELSA 30秒角（約1 km）気候値を分類群×1°／5°の地理セルに中央値集約。",
             fontsize=9,color="#6b7784")
    OUT.mkdir(parents=True,exist_ok=True)
    for ext in ("png","svg","pdf"):
        fig.savefig(OUT/f"figure_grid_geographic_extent_example.{ext}",
                    dpi=300,bbox_inches="tight",pad_inches=.1)
    plt.close(fig)
    print(distances)

if __name__=="__main__":
    main()
