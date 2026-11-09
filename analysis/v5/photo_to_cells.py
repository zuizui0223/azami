"""Account for observations, taxon-cells and unique map cells in the V5 paired tests."""
import numpy as np
import pandas as pd
from analysis.v5.trait_convergence import root
from analysis.v4.s1_climate_trait_hypervolume import load_source,geo_aggregate,standardise_cells
from analysis.v4 import s1_h3_climate_analogue_alignment as h


def main():
    source=root/"private/decoded"
    frame,meta=load_source(source/"observations.csv.gz",source/"measurement_receipt.json")
    result=[]
    for degree in [1,5]:
        cells,_=geo_aggregate(frame,degree)
        cells,_=standardise_cells(cells)
        cells=cells.reset_index(drop=True)
        taxa=h.eligible_taxa(cells)
        indices=[g.index.to_numpy() for _,g in cells.groupby("accepted_key",sort=True)]
        rng=np.random.default_rng(20261109)
        union=set()
        analogous=set()
        all_combos=[];near_combos=[]
        for _ in range(24):
            picks=[rng.choice(len(ix),8,replace=False) for ix in indices]
            x,y=[np.stack([t[k][s] for t,s in zip(taxa,picks)]) for k in ["x","y"]]
            selected=np.stack([ix[s] for ix,s in zip(indices,picks)])
            p,stat=h.matched_pairs(x,y,np.array([t["sector"] for t in taxa]))
            assert stat["status"]=="qualified"
            endpoints=np.unique(np.r_[p["i"],p["j"]])
            union.update(selected[endpoints].ravel())
            edge,left,right=np.nonzero(p["mask"])
            analogous.update(selected[p["i"][edge],left])
            analogous.update(selected[p["j"][edge],right])
            near_combos.append(len(edge))
            all_combos.append(len(p["i"])*64)
        row=dict(degrees=degree,source_obs=meta["source_observations"],quality_obs=len(frame),
                 quality_taxa=frame.accepted_key.nunique(),eligible_taxa=len(taxa),draws=24,
                 sampled_taxon_cells_per_draw=8*len(taxa),
                 mean_all_combinations=np.mean(all_combos),mean_analogue_combinations=np.mean(near_combos))
        for key,g in [("eligible",cells),("compared",cells.loc[sorted(union)]),("analogue",cells.loc[sorted(analogous)])]:
            row[key+"_taxon_cells"]=len(g)
            row[key+"_map_cells"]=len(g.drop_duplicates(["_latcell","_loncell"]))
            row[key+"_obs"]=int(g.n_photos.sum())
        result.append(row)
    pd.DataFrame(result).to_csv(root/"public/v5_photo_to_cells.csv",index=False)

if __name__=="__main__":main()
