"""Brownian-GLS sensitivity for the signed trait hypotheses on the nuclear subset."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import statsmodels.api as sm
from Bio import Phylo
from terminal_approach import rooted_ingroup,maps

def run(root):
    source=pd.read_csv(root/'nuclear_input.csv')
    rows=[]
    for file in sorted((root.parent/'two_axes/inputs').glob('*.nwk')):
        for low in [True,False]:
            d=source if low else source[source.S1_source_taxon!='Cirsium esculentum']
            tree=rooted_ingroup(Phylo.read(file,'newick'),source.run,list(d.run))
            tree.root.branch_length=0
            k,*_=maps(tree,list(d.run))
            x=d[['chelsa_bio18','chelsa_bio05','chelsa_rsds_mean']]
            x=sm.add_constant((x-x.mean())/x.std(ddof=0))
            for trait,terms in [('orientation_image_vertical_angle',['chelsa_bio18']),('corolla_lab_chroma',['chelsa_bio05','chelsa_rsds_mean'])]:
                y=d[trait].to_numpy();y=(y-y.mean())/y.std()
                fit=sm.GLS(y,x,sigma=k).fit()
                for term in terms:
                    lo,hi=fit.conf_int().loc[term]
                    rows.append(dict(tree='ASTRAL' if 'ASTRAL' in file.name else 'ML',taxa=len(d),trait=trait,term=term,beta=fit.params[term],lo=lo,hi=hi,p_pointwise=fit.pvalues[term]))
    result=pd.DataFrame(rows);result.to_csv(root/'signed_nuclear_sensitivity.csv',index=False)
    print(result[result.taxa.eq(24)].round(3).to_string(index=False))
if __name__=='__main__':run(Path(sys.argv[1]))
