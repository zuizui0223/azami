"""Render the three V5 questions from the source-recovered aggregate tables."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def render(path):
    focal=pd.read_csv(path/'v5_focal_no_elevation_gate.csv')
    s2=pd.read_csv(path/'s2_interactions.csv')
    s3=pd.read_csv(path/'s3_analogues.csv')
    fig,axes=plt.subplots(1,3,figsize=(13.8,4.25),layout='constrained')
    f=focal[(focal.spec=='conditional')&focal.trait.isin(['angle','chroma'])]
    order=[('angle','between'),('angle','within'),('chroma','between'),('chroma','within')]
    for i,(t,s) in enumerate(order):
        for geo,color,shift in [(False,'#135d79',-.11),(True,'#ae5931',.11)]:
            v=f[(f.trait==t)&(f.scale==s)&(f.geo==geo)].iloc[0]
            axes[0].errorbar(v.beta,i+shift,xerr=[[v.beta-v.lower],[v.upper-v.beta]],fmt='o',ms=4,
                capsize=2,color=color,label=('Climate' if not geo else 'Climate + geography') if i==0 else None)
    axes[0].set_yticks(range(4),['Angle / between','Angle / within','Chroma / between','Chroma / within'])
    axes[0].set_title('1. Taxonomic scales');axes[0].set_xlabel('Standardized slope')
    axes[0].legend(fontsize=8)
    f=s2[(s2.spec=='warm') & (((s2.trait=='angle')&(s2.axis=='bio18'))|((s2.trait=='chroma')&(s2.axis=='rsds_mean')))]
    for i,v in enumerate(f.itertuples()):
        axes[1].errorbar(v.beta,i,xerr=[[v.beta-v.lower],[v.upper-v.beta]],fmt='o',capsize=2,color='#135d79')
    axes[1].set_yticks(range(len(f)),[a.split()[-1]+' / '+b for a,b in zip(f.species,f.trait)])
    axes[1].set_title('2. Range comparison');axes[1].set_xlabel('Introduced minus native slope')
    f=s3[s3.threshold==.9]
    for i,v in enumerate(f.itertuples()):
        axes[2].scatter(v.excess,i,color='#135d79')
    axes[2].set_yticks(range(len(f)),[str(a)+'° / '+('heads matched' if b else 'climate matched') for a,b in zip(f.deg,f.reverse)])
    axes[2].set_title('3. Cross-taxon resemblance');axes[2].set_xlabel('Distance reduction above null')
    for ax in axes:
        ax.axvline(0,color='.5',linewidth=.7);ax.invert_yaxis()
    fig.savefig(path/'V5_results.png',dpi=200)
    plt.close(fig)

if __name__=='__main__':
    import sys
    render(Path(sys.argv[1]))
