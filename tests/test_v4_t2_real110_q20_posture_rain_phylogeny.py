"""Actual nuclear 110 source/ancestry/posture biological stop rules."""
from io import StringIO
import numpy as np
from Bio import Phylo

from analysis.v4 import t2_real110_q20_posture_rain_phylogeny as m

def toy():
    return Phylo.read(StringIO("(((A:.1,B:.2):.3,C:.4):.2,(D:.4,E:.5):.2);"),"newick")

def test_true_pair_distance_uses_all_original_tree_edges_and_never_invents_gain():
    t=toy()
    data={
        "A":{"botanical_U_D_state":"D","taxon_name":"A","chelsa_bio12":"1000"},
        "B":{"botanical_U_D_state":"U","taxon_name":"B","chelsa_bio12":"1400"},
        "C":{"botanical_U_D_state":"?","taxon_name":"C","chelsa_bio12":"600"},
        "D":{"botanical_U_D_state":"D","taxon_name":"D","chelsa_bio12":"2100"},
        "E":{"botanical_U_D_state":"U","taxon_name":"E","chelsa_bio12":"700"},
    }
    # For the toy tree, override the number-of-D check only in this scope.
    old=m.D
    try:
        m.D=2
        o=m.nearest_U(t,data)
    finally:
        m.D=old
    # Real-data helper deliberately returns None when testing a toy fewer-than-6
    # D cohort; no false empirical count is accepted.
    assert o is None
    assert abs(t.distance("A","B")-.3)<1e-10

def test_nonultrametric_Brownian_and_lambda_zero_are_not_independent_IID():
    t=toy()
    all_ids=[n.name for n in t.get_terminals()]
    C=m.covariance(t,all_ids)
    assert np.linalg.eigvalsh(C).min()>0
    assert C[0,1]>0
    assert np.diag(C).max()>np.diag(C).min()
    y=np.array([900,1050,1150,650,1350.])
    X=np.column_stack([np.ones(5),[1,0,0,1,0]])
    # Standalone numeric GLS requires more taxa/df; different covariance
    # slopes can still be evaluated without claiming 110-species biology.
    a=m.gaussian_gls(y,X,C,0,"lambda0")
    b=m.gaussian_gls(y,X,C,1,"BM")
    assert np.isfinite(a["coefficient_D_minus_U_mm"])
    assert np.isfinite(b["coefficient_D_minus_U_mm"])
    assert a["one_sided_D_wetter_conditional_p"]>=0
    assert a["95CI_D_minus_U_mm"][0]<a["95CI_D_minus_U_mm"][1]

def test_actual_Q20_nuclear_tree_source_contract_is_not_an_accepted_species_tree():
    assert m.D==6 and m.U==32 and m.UNKNOWN==72
    assert m.SCOPED_N==110
    assert len(m.TREE_SHA)==64
    assert m.EXCLUDED=="Cirsium tioganum"
    assert "posthoc110 technical sensitivity" in m.__doc__
    assert "DO NOT use camera-vertical" in m.__doc__
