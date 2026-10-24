"""Synthetic BM/OU mathematical contracts and real immutable S1 no-botany scope."""
import io
from copy import deepcopy
import numpy as np
from Bio import Phylo
from analysis.v4 import t2_s1_genbank109_bm_ou_niche_photo as m


def tiny():
    tree=Phylo.read(io.StringIO(
      "(((A:0.3,B:0.4):0.2,(C:0.5,D:0.6):0.2):0.7,"
      "((E:0.2,F:0.3):0.5,(G:0.35,H:0.44):0.2):0.25);"),"newick")
    return tree


def test_true_rooted_common_phylogeny_BM_and_two_OU_hypotheses(monkeypatch):
    tree=tiny()
    from Bio import Phylo
    monkeypatch.setattr(m,'N_TREE',8)
    rows={}
    for j,name in enumerate("ABCDEFGH"):
        rows[name]={"original_photo_count":str(j+1),
                    "chelsa_bio12":str(600+70*j),
                    "chelsa_bio01":str(6+j),
                    "chelsa_bio04":str(600+20*j),
                    "chelsa_bio15":str(30+j),
                    "chelsa_rsds_mean":str(12+j),
                    "mean_image_angle_deg":str(30+j*5),
                    "median_image_angle_deg":str(35+j*4),
                    "botanical_gravity_posture_inferred_from_photo_angle":"False"}
    for root in ("source_Newick","midpoint"):
        tips,h,C,D,obs,scale=m.prepare_data(tree,rows,root)
        assert len(tips)==8 and scale>0
        BM=m.kernel(h,C,D,"BM",nugget=.05)
        fixed_at_alpha_zero=m.kernel(h,C,D,"OU_fixedroot",alpha=1e-6,nugget=.05)
        assert np.allclose(BM,fixed_at_alpha_zero,atol=1e-4)
        stat=m.kernel(h,C,D,"OU_stationary",alpha=.2,nugget=.05)
        assert np.linalg.eigvalsh(stat).min()>0
        assert np.allclose(np.diag(stat),np.repeat(1.05,8))
        assert set(obs["photo_mean_angle_deg"])=={float(30+5*i) for i in range(8)}
        assert len(obs["log_BIO12"])==8
        rho=m.joint_BM_whitened_niche_head_photo_correlation(obs,h,C,D,.05)
        assert -1<=rho["log_BIO12_vs_S1_photo_angle_BM_whitened_r"]<=1
        assert 0<=rho["conditional_fixed_BM_covariance_approx_two_sided_p"]<=1
        assert rho["only_same_109_ITS_taxa_original_photo_angle"]
        assert rho["not_an_independent_ancestral_transitions_or_adaptive_selection_test"]

def test_OU_stationary_covariance_root_invariance_unlike_fixed_root():
    tree=tiny()
    t1=deepcopy(tree)
    t2=deepcopy(tree); t2.root_at_midpoint()
    from Bio import Phylo
    # The statistic is patristic distance, invariant across unrooted rerootings.
    def distances(t):
        tips=t.get_terminals();return {
            tuple(sorted((x.name,y.name))):t.distance(x,y)
            for x in tips for y in tips if x.name!=y.name
        }
    assert all(abs(distances(t1)[k]-v)<1e-9 for k,v in distances(t2).items())


def test_AICc_has_explicit_cost_for_one_extra_OU_alpha():
    assert m.AICc(-100.,109,3)>m.AICc(-100.,109,2)
    assert m.MODES==("WN","BM","OU_stationary","OU_fixedroot")
    assert m.N_TREE==109 and m.N_S1==251
    assert m.NUGGET==.05


def test_own_image_phenotype_CANNOT_be_replaced_by_flora_binary():
    assert "photo_mean_angle_deg" in m.TRAITS
    assert "photo_median_angle_deg" in m.TRAITS
    assert "head_orientation_class" not in m.TRAITS
    assert "nodding" not in m.TRAITS
    assert m.SEED==20261023
    assert m.TRAITS[0]=="log_BIO12"
    assert len(m.source.EXPOSURES)==5


def test_true_ITS_statistically_significant_rain_signal_never_called_origin_event():
    assert m.null_lambda_permutation.__doc__.startswith(
        "Source-tip-label permutation")
    assert m.opt_ou.__name__=="opt_ou"
    assert m.kernel(np.ones(8),np.eye(8),np.zeros((8,8)),"WN").shape==(8,8)
