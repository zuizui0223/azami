"""Draw the species-specific predictive and direction-of-change checks."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def plot(root):
    t=pd.read_csv(root/'transfer.csv')
    r=pd.read_csv(root/'rain_displacement.csv')
    fig,axs=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    groups=['native_spatial_cv','introduced_analogous','introduced_novel']
    colors=['#365b70','#66854a','#a35131']
    for idx,(sp,g) in enumerate(t.groupby('species')):
        g=g.set_index('group')
        for j,name in enumerate(groups):
            ax=axs[0];z=g.loc[name]
            ax.scatter(z.q2,idx+(j-1)*.2,color=colors[j],label=name.replace('_',' ') if idx==0 else None)
    axs[0].set_yticks(range(2),['C. arvense','C. vulgare'])
    axs[0].set_title('Native model transfer')
    axs[0].set_xlabel('Q² relative to native mean')
    axs[0].legend(fontsize=7,loc='best')
    select=r[(r.group=='introduced_novel')&r.comparison.isin(['wetter','drier'])]
    for idx,z in enumerate(select.itertuples()):
        axs[1].plot([z.lower,z.upper],[idx,idx],color='#666')
        axs[1].scatter(z.angle_estimate,idx,color='#365b70' if z.comparison=='wetter' else '#aa6944')
    axs[1].set_yticks(range(len(select)),[sp.split()[-1]+' / '+d for sp,d in zip(select.species,select.comparison)])
    axs[1].set_title('Novel sites versus native neighbours')
    axs[1].set_xlabel('Difference in photographed angle (degrees)')
    for ax in axs:ax.axvline(0,color='.65',lw=.9);ax.invert_yaxis()
    fig.savefig(root/'transfer_summary.png',dpi=190)
    plt.close(fig)
if __name__=='__main__':plot(Path(sys.argv[1]))
