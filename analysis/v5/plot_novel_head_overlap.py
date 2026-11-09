from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import font_manager
from statsmodels.stats.multitest import multipletests

def render(root):
    tab=pd.read_csv(root/'novel_head_overlap.csv')
    d=tab[(tab.native95_mass==.95)&(tab.n_traits==3)&(tab.bandwidth==.5)]
    assert len(d)==4 and np.allclose(tab.D_native_analogous-tab.D_native_novel,tab.drop_in_overlap)
    assert ((tab[['D_native_analogous','D_native_novel','D_analogous_novel']]>=0)&(tab[['D_native_analogous','D_native_novel','D_analogous_novel']]<=1)).all().all()
    tested=d[d.design=='geographic_pairs'];assert np.allclose(tested.holm_primary,multipletests(tested.p_swap,method='holm')[1])
    plt.rcParams.update({'font.family':['Meiryo' if 'Meiryo' in font_manager.get_font_names() else 'Noto Sans CJK JP','DejaVu Sans'],'pdf.fonttype':42,'font.size':10})
    fig,axs=plt.subplots(1,2,figsize=(8.5,4.6),sharex=True)
    fig.subplots_adjust(top=.76,bottom=.23,left=.14,right=.98,wspace=.45)
    for ax,(sp,g) in zip(axs,d.groupby('species')):
        for j,model in enumerate(['all_cells','geographic_pairs']):
            r=g[g.design==model].iloc[0]
            pts=[r.D_native_analogous,r.D_native_novel]
            ax.plot(pts,[j,j],color='.6',lw=1.5,zorder=1)
            for k,(v,col) in enumerate(zip(pts,['#30749b','#b95227'])):
                ax.scatter(v,j,s=65,c=col,zorder=2,label=(['移入先：在来類似気候','移入先：新規気候'][k] if j==0 else None))
                ax.text(v,j+(.16 if k else -.17),format(v,'.3f'),ha='center',va='center',color=col,fontsize=10)
        ax.set(xlim=(0,1),ylim=(1.45,-.47),xlabel='在来域の頭花との重複度 D')
        ax.set_yticks([0,1],['全セル','地理を揃える'])
        ax.set_title(sp,fontfamily='DejaVu Sans',fontstyle='italic',fontsize=13)
        ax.spines[['top','right']].set_visible(False);ax.grid(axis='x',alpha=.2)
    fig.suptitle('新規気候の頭花は、在来域の頭花から離れるか',fontsize=16,y=.96)
    fig.legend(*axs[0].get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.5,.06),ncol=2,frameon=False)
    fig.text(.5,.02,'頭花空間：角度・写真の彩度・投影輪郭。右ほど、在来域の分布に近い。',ha='center',fontsize=9)
    for ext in ['png','pdf']:fig.savefig(root/('native_head_overlap.'+ext),dpi=185,bbox_inches='tight')
    plt.close(fig)
if __name__=='__main__':render(Path(sys.argv[1]))
