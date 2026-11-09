"""Plot overlap of the 95% native and introduced background-corrected densities."""
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from matplotlib import font_manager

def main(path):
    japanese='Meiryo' if 'Meiryo' in font_manager.get_font_names() else 'Noto Sans CJK JP'
    plt.rcParams.update({'font.family':[japanese,'DejaVu Sans'],'pdf.fonttype':42,'font.size':10})
    colors=['#ffffff','#609dc5','#ecac6c','#6b537c']
    fig,axs=plt.subplots(1,2,figsize=(9.4,5.1))
    fig.subplots_adjust(top=.80,bottom=.25,left=.09,right=.97,wspace=.25)
    for ax,sp in zip(axs,['Cirsium arvense','Cirsium vulgare']):
        d=np.load(path/(sp.replace(' ','_')+'.npz'))
        aa,bb=d['native95'],d['introduced95'];a,b=d['native'],d['introduced']
        assert a.shape==b.shape==(80,80) and np.isclose(a.sum(),1) and np.isclose(b.sum(),1)
        group=aa.astype(int)+2*bb.astype(int)
        x,y=d['x'],d['y'];v=d['variance']
        ax.pcolormesh(x,y,group.T,cmap=ListedColormap(colors),vmin=-.5,vmax=3.5,shading='flat',rasterized=True)
        # Limit blank margins while retaining every occupied 95% cell.
        ix,iy=np.where(group>0)
        ax.set_xlim(x[max(0,ix.min()-2)],x[min(len(x)-1,ix.max()+3)])
        ax.set_ylim(y[max(0,iy.min()-2)],y[min(len(y)-1,iy.max()+3)])
        ax.set_aspect('equal',adjustable='box')
        ax.set_title(sp+'   D = '+format(np.minimum(a,b).sum(),'.3f'),fontfamily='DejaVu Sans',fontsize=12)
        ax.set_xlabel('気候 PC1（'+format(100*v[0],'.0f')+'%）')
        ax.set_ylabel('気候 PC2（'+format(100*v[1],'.0f')+'%）')
        ax.spines[['top','right']].set_visible(False)
        ax.tick_params(labelsize=8)
    fig.suptitle('在来域と移入域の気候ニッチの重なり',fontsize=17,y=.97)
    fig.legend([Patch(facecolor=c) for c in colors[1:]],['在来域のみ','移入域のみ','両域で重複'],
               ncol=3,loc='lower center',bbox_to_anchor=(.5,.085),frameon=False,fontsize=11)
    fig.text(.5,.042,'色面は各分布域の背景補正済み出現密度の95%領域。D は密度全体の重複度。',ha='center',fontsize=9)
    fig.text(.5,.005,'気候軸は各種の分布域内で共通。種間で座標を直接比較しない。',ha='center',fontsize=8)
    for ext in ['png','pdf']:
        fig.savefig(path/('V5_S2_niche_overlap.'+ext),dpi=190,bbox_inches='tight')
    plt.close(fig)

if __name__=='__main__':main(Path(sys.argv[1]))
