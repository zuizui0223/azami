from pathlib import Path
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def render(root):
    contrasts=pd.read_csv(root/'trait_contrasts.csv')
    cv=pd.read_csv(root/'trait_covariation.csv')
    fig,axes=plt.subplots(2,3,figsize=(11,5.8),layout='constrained')
    palette=['#275e85','#a55135']
    for j,(trait,label) in enumerate([('angle','Photo angle (degrees)'),('chroma','Photographic chroma (C*)'),('aspect','Outline height / width')]):
        ax=axes[0,j]
        sub=contrasts[(contrasts.coverage==.95)&(contrasts.trait==trait)].sort_values('species')
        for i,s in enumerate(sub.itertuples()):
            if s.status=='estimated':
                ax.errorbar(s.beta,i,xerr=[[s.beta-s.lower],[s.upper-s.beta]],fmt='o',color=palette[i],capsize=3)
        ax.set_yticks(range(len(sub)),[s.split()[-1] for s in sub.species]);ax.axvline(0,lw=.8,color='.5')
        ax.set_title(label);ax.set_xlabel('Novel minus analogous climate')
    for j,(x,y) in enumerate([('angle','aspect'),('angle','chroma'),('aspect','chroma')]):
        ax=axes[1,j];sub=cv[(cv.predictor==x)&(cv.outcome==y)&cv.climate_date_adjusted].sort_values(['species','novel'])
        for i,s in enumerate(sub.itertuples()):
            ax.errorbar(s.beta,i,xerr=[[s.beta-s.lower],[s.upper-s.beta]],fmt='o',color=palette[int(s.novel)],capsize=3)
        ax.set_yticks(range(len(sub)),[a.split()[-1]+' / '+('novel' if b else 'analogous') for a,b in zip(sub.species,sub.novel)])
        ax.set_title(x.title()+' and '+y);ax.set_xlabel('Climate-adjusted standardized slope')
        ax.axvline(0,lw=.8,color='.5')
    for ax in axes.ravel():ax.invert_yaxis()
    fig.savefig(root/'S2_novelty_traits.png',dpi=190);plt.close(fig)
if __name__=='__main__':render(Path(sys.argv[1]))
