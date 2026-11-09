"""Show the tested signs of the two most directly motivated ecological predictions."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def main(path):
    t=pd.read_csv(path/'climate_trait_direction.csv')
    fig,axes=plt.subplots(1,2,figsize=(9.6,4.5),layout='constrained')
    for ax,(trait,axis,title) in zip(axes,[('angle','bio18','Rain and photographed angle'),('chroma','rsds_mean','Radiation and photographic chroma')]):
        g=t[(t.panel=='S1')&(t.trait==trait)&(t.axis==axis)]
        for j,s in enumerate(['between','within']):
            for k,geo in enumerate(['False','True']):
                v=g[(g.scope==s)&(g.geo==geo)].iloc[0]
                ax.errorbar(v.beta,j+.14*(2*k-1),xerr=[[v.beta-v.lo],[v.hi-v.beta]],color=['#245976','#ad582f'][k],fmt='o',capsize=3,label=['Climate only','Climate + geography'][k] if j==0 else None)
        ax.axvline(0,color='.55',lw=1)
        ax.set_title(title);ax.set_yticks([0,1],['Among taxa','Within taxa'])
        ax.set_ylim(1.55,-.5);ax.set_xlabel('Standardized slope (prediction: positive)')
    axes[1].legend(loc='lower left',fontsize=8)
    fig.savefig(path/'directional_tests.png',dpi=180);plt.close(fig)
if __name__=='__main__':main(Path(sys.argv[1]))
