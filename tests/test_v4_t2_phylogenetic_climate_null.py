"""T2 climate-source and conditional phylogenetic null validation tests."""
from io import StringIO
from Bio import Phylo
import numpy as np
import pytest
from analysis.v4 import t2_phylogenetic_climate_null as m

def example():
    names=[f"Cirsium_{i:02d}" for i in range(20)]
    y={name:int(i in (2,4,6,8,10)) for i,name in enumerate(names)}
    climate={name:{"BIO12":300+11*i+3*(i%5)**2,
        "BIO1":15+i*2.2+4*((i+2)%3),
        "Elevation":100+i*27+30*((i+1)%4)}
        for i,name in enumerate(names)}
    t=Phylo.read(StringIO("("+",".join(f"{z}:0.5" for z in names)+");"),"newick")
    return names,y,climate,t

def test_residualization_and_scale_invariance():
    names,y,climate,_=example()
    states,resid,info=m.climate_contrast(names,y,climate)
    assert info["n_tips"]==20 and info["n_nodding"]==5
    assert np.allclose(resid.mean(axis=0),0,atol=1e-10)
    for row in climate.values():
        row["BIO12"]*=2
        row["BIO1"]-=3000
    _,again,_=m.climate_contrast(names,y,climate)
    assert np.allclose(resid,again,atol=1e-10)

def test_binary_tree_null_root_and_two_correction():
    names,y,climate,tree=example()
    states,resid,_=m.climate_contrast(names,y,climate)
    observed=m.signed_contrast(states[None,:],resid)[0]
    shuffle=m.simulate_fixed_prevalence_shuffle(20,5,resid,1024,11)
    infer=m.inference_from_null(observed,shuffle,{"status":"calibrated"},200)
    assert all(a>=b for a,b in zip(infer["two_axis_maxT_adjusted_p"],
                                    infer["raw_empirical_two_sided_p"]))
    for root in ("erect","stationary"):
        null,d=m.conditional_mk_null(tree,y,resid,root_mode=root,seed=8,
            target_draws=100,max_proposals=20000,min_accepted=30)
        assert null.shape[1]==2 and d["observed_tip_prevalence"]==5
        assert d["accepted_occupancy_matched"]>=len(null)
        assert d["total_simulated_tip_histories"]<=20000

def test_matched_original_taxa_required():
    names,y,climate,_=example()
    del climate[names[0]]
    with pytest.raises(ValueError,match="Invalid common"):
        m.climate_contrast(names,y,climate)
