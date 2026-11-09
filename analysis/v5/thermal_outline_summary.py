"""Export source-linked focal comparisons and plot the archived outline index."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def summarise(root):
    s1=pd.read_csv(root/'v5_focal_no_elevation_gate.csv')
    s2=pd.read_csv(root/'trait_gradients.csv')
    slopes=s2[s2.coverage.eq(.95)].copy()
    assert len(slopes)==10
    assert np.allclose(s2.slope_new-s2.slope_sim,s2.slope_difference)
    assert (s2.sim_lower<=s2.slope_sim).all() and (s2.sim_upper>=s2.slope_sim).all()
    facet=s1[s1.spec.eq('conditional')]
    assert len(facet)==20
    facet.to_csv(root/'S1_hypothesis_slopes.csv',index=False)
    slopes.to_csv(root/'S2_hypothesis_slopes.csv',index=False)
    fig,axs=plt.subplots(1,2,figsize=(11,3.6),layout='constrained')
    for ax in axs:ax.axvline(0,color='.65',linewidth=.8)
    ci=s1[s1.trait.eq('aspect')&~s1.geo]
    for i,s in enumerate(ci.itertuples()):
        axs[0].errorbar(s.beta,i,xerr=[[s.beta-s.lower],[s.upper-s.beta]],fmt='o',color='#285872',capsize=3)
    axs[0].set_yticks(range(len(ci)),[('Temp.' if a=='bio05' else 'Radiation')+' / '+b for a,b in zip(ci.axis,ci.scale)])
    ci=slopes[slopes.trait.eq('aspect')]
    for i,s in enumerate(ci.itertuples()):
        for offset,kind,color in [(-.12,'sim','#285872'),(.12,'new','#b24d34')]:
            b=getattr(s,'slope_'+kind);lo=getattr(s,kind+'_lower');hi=getattr(s,kind+'_upper')
            axs[1].errorbar(b,i+offset,xerr=[[b-lo],[hi-b]],fmt='o',color=color,capsize=3,
                            label=('Analogous' if kind=='sim' else 'Novel') if i==0 else None)
    axs[1].set_yticks(range(len(ci)),[s.split()[-1]+' / '+('temp.' if a=='bio05' else 'radiation') for s,a in zip(ci.species,ci.axis)])
    axs[0].set_title('S1 / native taxa');axs[1].set_title('S2 / introduced sites')
    axs[0].set_xlabel('Standardized index slope');axs[1].set_xlabel('Index change per climate SD')
    axs[1].legend(fontsize=8)
    for ax in axs:ax.invert_yaxis()
    fig.suptitle('Archived 2D outline index (sensitive to image rotation)')
    fig.savefig(root/'outline_climate.png',dpi=190)
    plt.close(fig)

if __name__=='__main__':summarise(Path(sys.argv[1]))
