from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager
import pandas as pd

def render(path):
    plt.rcParams.update({'font.family':['Meiryo' if 'Meiryo' in font_manager.get_font_names() else 'Noto Sans CJK JP','DejaVu Sans'],'pdf.fonttype':42,'font.size':10})
    points=pd.read_csv(path/'pairs.csv')
    tests=pd.read_csv(path/'regional.csv');whole=pd.read_csv(path/'whole.csv')
    fig,axs=plt.subplots(1,2,figsize=(9.1,5.25),sharex=True,sharey=True)
    fig.subplots_adjust(top=.79,bottom=.24,left=.1,right=.98,wspace=.08)
    for ax,(sp,g) in zip(axs,points[points['shape']].groupby('species')):
        t=tests[(tests.species==sp)&tests['shape']&(tests.h==.5)].iloc[0]
        ax.scatter(g.climate_D,g.head_D,s=22,alpha=.52,c='#26617a',edgecolors='none')
        ax.set_title(sp,fontsize=13,fontfamily='DejaVu Sans',fontstyle='italic')
        ax.text(.04,.96,'ρ = '+format(t.rho,'+.2f')+'    補正 P = '+format(t.holm_primary,'.3f'),
                transform=ax.transAxes,va='top',fontsize=10)
        ax.set(xlim=(0,1),ylim=(0,1),xlabel='気候ニッチの重複度 D',aspect='equal')
        ax.grid(color='#e5e7e8',lw=.6);ax.set_axisbelow(True)
    axs[0].set_ylabel('頭花形質空間の重複度 D')
    fig.suptitle('在来地域と移入地域の気候と頭花は、どの程度一緒に重なるか',fontsize=16,y=.94)
    fig.text(.5,.10,'各点：在来域の1地域 × 移入域の1地域。両空間に同じ写真・同じセルを使用。',ha='center',fontsize=9)
    fig.text(.5,.047,'D は 0 で重なりなし、1 で密度分布が一致。P は地域単位の置換に基づく。',ha='center',fontsize=9)
    for ext in ['pdf','png']:fig.savefig(path/('paired_overlap.'+ext),dpi=185,bbox_inches='tight')
    plt.close(fig)
if __name__=='__main__':render(Path(sys.argv[1]))
