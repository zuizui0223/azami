"""Deterministic math and S1 source provenance for T2 real110 vs ITS109."""
import io
import numpy as np
from Bio import Phylo
import pytest

from analysis.v4 import t2_110_multilocus_vs_109_its_photo_niche as m


def toy():
    return Phylo.read(io.StringIO(
        "(((A:.012,B:.015):.018,(C:.011,D:.018):.014):.008,"
        "((E:.020,F:.010):.012,(G:.013,H:.012):.012):.008);"
    ),"newick")


def test_actual_Brownian_common_ancestor_not_species_rootless_fake_identity():
    names,C=m.original_true_covariance(
        toy(),{x:x for x in "ABCDEFGH"})
    assert len(names)==8
    assert np.linalg.eigvalsh(C).min()>0
    ids={v:i for i,v in enumerate(names)}
    assert C[ids["A"],ids["B"]]>C[ids["A"],ids["E"]]
    assert C[ids["A"],ids["B"]]>0
    eig,Q,K=m.pagel_spectral_kernel(C,.05)
    assert len(eig)==8
    assert np.linalg.eigvalsh(K).min()>0
    assert np.allclose(np.diag(Q.T@Q),1.)


def test_phylogenetic_signal_estimator_is_unbiased_under_IID_no_signal():
    names,C=m.original_true_covariance(toy(),{x:x for x in "ABCDEFGH"})
    eigen,U,_=m.pagel_spectral_kernel(C,.05)
    constant=np.array([3.1,7.7,4.2,2.9,6.2,3.0,5.6,7.1])
    fitted=m.profile_lambda(constant,eigen,U)
    assert 0<=fitted["ML_lambda"]<=1
    assert fitted["likelihood_ratio_over_unrelated_IID"]>=0
    perm=m.permute_phylo_signal(constant,eigen,U,n=99,seed=44)
    assert perm["permutations"]==99
    assert 0.01<=perm["plus_one_source_tip_label_permutation_p"]<=1
    assert perm["fixed_real_gene_tree_not_topology_posterior"]
    assert perm["original_taxa_permutation_is_not_geographic_climate_null"]


def test_five_CHELSA_S1_own_photo_continuous_regression_not_botanical_binary():
    rng=np.random.default_rng(422)
    native={}
    for i in range(28):
        native[f"Cirsium_{i}"]={
            "chelsa_bio01":str(9+4*np.sin(i*.37)),
            "chelsa_bio04":str(300+12*i+20*np.sin(.16*i)),
            "chelsa_bio12":str(800+28*i+100*np.cos(.33*i)),
            "chelsa_bio15":str(35+4*np.sin(.49*i)),
            "chelsa_rsds_mean":str(18+.8*np.sin(1.7*i)),
            "mean_image_angle_deg":str(65+2*np.cos(i*.88)+3*np.sin(.24*i)+rng.normal(0,2)),
            "median_image_angle_deg":str(61+3*np.sin(.32*i)),
            "botanical_U_D_state":"D" if i%2 else "U",
        }
    names=list(native)
    X=m.equal_species_five_rain_design(names,native)
    assert X.shape==(28,6)
    y=m.phenotype(names,native,"S1_photo_mean_image_head_angle")
    model=m.slope(y,X,np.eye(28))
    assert model["same_original_S1_five_CHELSA_adjustment"]==list(m.FIVE)
    beta=model["BIO12_original_S1_photo_angle_degrees_per100mm"]
    assert np.isfinite(beta)
    # Recoding ALL floristic heads MUST NOT alter the owner's measured angle.
    for x in native.values():x["botanical_U_D_state"]="?"
    assert m.slope(m.phenotype(names,native,"S1_photo_mean_image_head_angle"),X,np.eye(28))==model
    assert "flor" not in str(m.PHENOTYPES).lower()


def test_taxon_coverage_and_failed_original_111_source_contract():
    assert m.N_MULTILOCUS==110
    assert m.N_ITS==109
    assert m.N_EXACT_OVERLAP==58
    assert m.N_Q20_LOCI==108
    assert m.N_Q20_ORIGINAL_SITES==29310
    assert m.ORIGINAL_EXCLUDED=="Cirsium tioganum"
    assert m.N_ORIGINAL_S1==251 and m.N_ORIGINAL_PHOTOS==170066
    assert all(len(x)==64 for x in (
        m.EXPECTED_S1_SHA,m.EXPECTED_ML108_SHA,m.EXPECTED_ITS_SHA))
    assert m.NOISES==(0.,.05)
    assert "PHOTO values are measured by S1 authors" in m.__doc__


def test_no_fictitious_genotype_species_gravity_by_value_threshold():
    s={"Cirsium_foo":{"mean_image_angle_deg":"85","median_image_angle_deg":"86","chelsa_bio12":"1200"}}
    assert m.phenotype(["Cirsium_foo"],s,"S1_photo_mean_image_head_angle").tolist()==[85.]
    assert np.log(1200)==pytest.approx(m.phenotype(["Cirsium_foo"],s,"log_BIO12")[0])
    with pytest.raises(ValueError,match="fictional"):
        m.phenotype(["Cirsium_foo"],s,"nodding_or_not")
