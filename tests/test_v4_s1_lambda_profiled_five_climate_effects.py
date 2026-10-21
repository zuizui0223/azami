"""Same S1 measured angle under independently fitted lambda, forced lambda0/1."""
from io import StringIO

from Bio import Phylo
import numpy as np
import pytest

from analysis.v4 import s1_lambda_profiled_five_climate_effects as eff
from analysis.v4.joint_s1_t2_phylogeny import EXPOSURES
from analysis.v4.s1_continuous_phylogenetic_signal import LambdaProfile, covariance_for


def fixture(n=28):
    rng=np.random.default_rng(337)
    tips=[f"Cirsium_t{i:02d}" for i in range(n)]
    raw=np.array([rng.normal(size=5) for _ in range(n)])
    raw[:,0]*=10.;raw[:,1]*=100.;raw[:,2]=raw[:,2]*200+800
    raw[:,3]=raw[:,3]*15+50;raw[:,4]=raw[:,4]*2+9
    y=70.+0.01*raw[:,2]+rng.normal(0,.35,n)
    rows=[{
       "tip":name,"accepted_key":str(i),
       "image_mean_deg":float(y[i]),"image_median_deg":float(y[i]+.5),
       "S1_five_climate_means":dict(zip(EXPOSURES,raw[i].tolist())),
       "source_concept_collision":False,
       "n_eligible_images":20,
    } for i,name in enumerate(tips)]
    a=",".join(f"{name}:0.5" for name in tips[:n//2])
    b=",".join(f"{name}:0.5" for name in tips[n//2:])
    tree=Phylo.read(StringIO("(("+a+"):1,("+b+"):1);"),"newick")
    return rows,tree


def test_profiled_BIO12_S1_image_degree_units_and_CIs():
    rows,tree=fixture()
    fits=eff.model_one_tree(rows,tree,nuggets=(.05,))
    assert len(fits)==6, "2 continuous traits x 3 lambda models"
    assert {r["fit"] for r in fits}==set(eff.FIXED_MODEL_CLASS)
    assert all(r["source_five_climate_variables"]==list(EXPOSURES) for r in fits)
    lambda0=next(r for r in fits if r["fit"]=="independent_lambda0"
                 and r["response"]=="species_mean")
    beta=lambda0["coefficients"]["chelsa_bio12"]
    assert beta["increment"]==100.
    assert beta["effect_image_deg_per_original_increment"]>0
    assert beta["conditional_t_95pct"][0]>0
    assert lambda0["lambda"]==0.
    assert all(r["not_gravity_head_angle"] for r in fits)


def test_effect_equivariant_to_climate_units_after_original_increment_fixed():
    rows,tree=fixture()
    x=eff.model_one_tree(rows,tree,nuggets=(.05,))
    for row in rows:
        row["S1_five_climate_means"]["chelsa_bio12"]*=2
    y=eff.model_one_tree(rows,tree,nuggets=(.05,))
    a=next(z for z in x if z["fit"]=="independent_lambda0" and z["response"]=="species_mean")
    b=next(z for z in y if z["fit"]=="independent_lambda0" and z["response"]=="species_mean")
    # Double the original physical BIO12 values => half per 100-unit effect.
    assert b["coefficients"]["chelsa_bio12"]["effect_image_deg_per_original_increment"]==pytest.approx(
       a["coefficients"]["chelsa_bio12"]["effect_image_deg_per_original_increment"]/2,abs=1e-8)


def test_lambda_bounds_conditional_nugget_and_climate_crosswalk():
    rows,tree=fixture()
    from analysis.v4.s1_continuous_phylogenetic_signal import model_design
    _,X=model_design(rows,True)
    covariance=covariance_for(tree,rows)
    p=LambdaProfile(covariance,X,.05)
    raw=np.array([r["image_mean_deg"] for r in rows])
    sd=np.std([[r["S1_five_climate_means"][v] for v in EXPOSURES]
               for r in rows],axis=0)
    with pytest.raises(ValueError,match="Lambda outside"):
        eff.estimate_coefficients(p,raw,1.2,sd)
    with pytest.raises(ValueError,match="standardization"):
        eff.estimate_coefficients(p,raw,0.,np.zeros(5))
    assert eff.estimate_coefficients(p,raw,0.,sd)["chelsa_bio15"]["increment"]==10.
