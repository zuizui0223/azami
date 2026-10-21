"""Photo-angle versus independent botanical posture (not gravity calibration)."""
import numpy as np
import pytest
from analysis.v4 import s1_angle_external_morphology_check as a
from analysis.v4.joint_s1_t2_phylogeny import EXPOSURES

def test_image_ranking_auc_with_ties_and_direction():
    assert a.rank_auc([20.,40.,70.,90.],[0,0,1,1])==1.
    assert a.rank_auc([20.,40.,70.,90.],[1,1,0,0])==0.
    assert a.rank_auc([50.,50.,50.,50.],[0,0,1,1])==.5
    with pytest.raises(ValueError,match="No botanical-state"):
        a.rank_auc([20.,40.],[0,0])

def fixture():
    rng=np.random.default_rng(7)
    rows=[]
    for i in range(46):
        NA=i<42
        is_D=(i in (0,11,20,42,43,44))
        climate=rng.normal(size=5)
        rows.append(dict(
            accepted_key=str(3000+i),
            source_taxon=f"Cirsium simulated_{i:02d}",
            source=a.NA if NA else a.EA,
            state="nodding" if is_D else "erect",
            s1_five_climate_centroids=dict(zip(EXPOSURES,climate.tolist())),
            s1_image_orientation_median=125. if is_D else 65.,
        ))
    return rows

def test_finite_prevalence_exact_46_source_AUC_and_null():
    result=a.source_blocked_exact_auc(fixture())
    assert result["n_S1_angle_taxa"]==46
    assert result["n_botanical_nodding"]==6
    assert result["source_prevalence_exact_assignments"]==45920
    assert result["pooled_mixture_AUC_descriptive_not_source_adjusted"]==1
    assert result["within_source_pair_weighted_AUC"]==1
    assert result["within_source_equal_weight_AUC"]==1
    assert 0<result["within_source_pair_weighted_upper_p"]<.05
    assert result["source_strata"]["North_American_published_author"]["source_specific_total_possible_assignments"]==11480
    assert result["source_strata"]["East_Asian_EAzami_curated"]["source_specific_total_possible_assignments"]==4
    assert result["external_binary_states_not_derived_from_image_angle"] is True
    assert result["cannot_convert_194_image_angle_to_true_nodding_independent_gain_counts"] is True

def test_no_picture_angle_as_gravity_and_range_guard():
    rows=fixture()
    rows[0]["s1_image_orientation_median"]=199
    with pytest.raises(ValueError,match="range"):
        a.source_blocked_exact_auc(rows)
