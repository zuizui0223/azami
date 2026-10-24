"""No fabricated 161-tip species tree or categorical flora U/D in S1 T2."""
from io import StringIO
import numpy as np
from Bio import Phylo
import pytest

from analysis.v4 import t2_s1_real_nuclear_phylogenetic_clustering as m


def toy():
    return Phylo.read(StringIO(
      "(((A:1,B:1):1,(C:1,D:1):1):1,((E:1,F:1):1,(G:1,H:1):1):1);"),"newick")


def test_closest_tip_distance_is_rooted_tree_metric_not_traits():
    t=toy()
    near,d=m.nearest_phylo_neighbor(t)
    assert list(near)==[1,0,3,2,5,4,7,6]
    assert float(d[0,1])==2.0
    assert np.isinf(d[0,0])


def test_real_continuous_scalar_detects_permuted_clade_similarity():
    arr=[10.,10.,30.,30.,50.,50.,70.,70.]
    near,_=m.nearest_phylo_neighbor(toy())
    out=m.trait_clustering(arr,near,"test-mean_image_angle_deg",reps=999)
    assert out["observed_closest_tip_mean_absolute_difference"]==0
    assert 0<out["single_sided_clustering_p"]<.05
    assert out["observed_to_null_ratio_below1_indicates_clustering"]==0
    assert out["n_permutations"]==999


def test_permutation_replay_is_deterministic():
    vals=np.array([20,30,50,60,80,100,120,140.])
    near,_=m.nearest_phylo_neighbor(toy())
    a=m.trait_clustering(vals,near,"same-cohort-mean_image_angle_deg",reps=999)
    b=m.trait_clustering(vals,near,"same-cohort-mean_image_angle_deg",reps=999)
    assert a==b


def test_fixed_scientific_limits_from_original_real_sources():
    assert m.N_S1==251 and m.PHOTO_S1==170066
    assert m.TREE_COHS["NA42_real_authored_nuclear"]==42
    assert m.TREE_COHS["ITS109_GenBank_nuclear_rDNA"]==109
    assert m.TREE_COHS["Q20_posthoc110_candidate_loci"]==110
    assert m.TREE_COHS["same58_ITS_gene_tree"]==58
    assert m.TREE_COHS["same58_Q20_candidate_tree"]==58
    assert m.MIN_PHOTOS==(1,10,30)
    assert m.ANGLE_ENDPOINTS==("mean_image_angle_deg","median_image_angle_deg")
    assert m.NA63_GIT_BLOB_SHA=="6417a38793427bc7c8562dc0e2cab59b74d8c316"
