"""Registered computational robustness checks for native Gates B/C/D.

Test files use SYNTHETIC images and observer labels only. Real numerical
estimates must come from the source-hashed private Actions executor.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from analysis.v4.native_density_robustness import (
    source_scores, pair_metrics, taxon_paired_difference,
    bandwidth_sensitivity, deletion_sensitivity, analyse,
)
from analysis.v4.native_gates_ae import CORE9


def mock_candidate(taxa=8,per_taxon=110):
    rng=np.random.default_rng(20261021)
    rows=[]
    for i in range(taxa):
        for j in range(per_taxon):
            angle=j/5+i*3+rng.normal(scale=.7)
            hue=(i*.3+j*.07) % (2*np.pi)
            row={
                "accepted_key":str(100+i),
                "obs_id":f"synthetic-{i:03d}-{j:04d}",
                "user_id":f"observer-{i}-{j%24}",
                "_cell":int((j%12)+i*20),
                "latitude":float(-35 + i*.4 + j*.025),
                "longitude":float(-165+((j*3+i*31)%330)),
                "orientation_image_vertical_angle":angle,
                "corolla_lab_lightness":float(40+i*.2+rng.normal(scale=2)),
                "corolla_lab_chroma":float(25+i*.1+rng.normal(scale=3)),
                "corolla_hue_sin":float(np.sin(hue)),
                "corolla_hue_cos":float(np.cos(hue)),
                "capitulum_outline_aspect_ratio":float(1.3+i*.01+rng.normal(scale=.08)),
                "capitulum_outline_circularity":float(.7+i*.001+rng.normal(scale=.04)),
                "capitulum_outline_solidity":float(.85+rng.normal(scale=.02)),
                "capitulum_width_profile_cv":float(.2+i*.003+rng.normal(scale=.02)),
            }
            rows.append(row)
    assert set(CORE9)<=set(rows[0])
    return pd.DataFrame(rows)


def test_30_method_components_share_identical_pair_reference_universe():
    f=mock_candidate(8,110)
    taxa=sorted(f.accepted_key.unique())
    scores,retained=source_scores(f,taxa)
    assert set(scores)=={"presentation","colour","outline"}
    assert scores["colour"][taxa[0]].shape==(110,3)
    assert 0<retained["colour"]<=1
    summary,private=pair_metrics(scores,taxa,reps=3,null_reps=4,n=10)
    assert len(summary)==3
    assert all(x["pairs_with_both_taxa_null"]==28 for x in summary.values())
    assert all(x["reference_uses_pair_median_not_first_rep_only"] for x in summary.values())
    assert all(summary[x]["rarefaction_replicates"]==3 for x in summary)
    assert len(private["colour"]["pair_index"])==28
    assert np.array_equal(private["colour"]["matched"],private["outline"]["matched"])
    contrast=taxon_paired_difference(private,taxa,f,bootstrap=40)
    assert contrast["shared_species_pairs"]==28
    assert contrast["reference_supported_taxa"]==8
    assert contrast["taxon_cluster_bootstrap_valid"]>=38
    assert sum(contrast["pair_discordance"].values())==28
    assert contrast["exclude_two_most_reported_taxa"]["remaining_pairs"]==15
    assert len(contrast["taxon_cluster_bootstrap_q025_q50_q975"])==3
    assert 0<=contrast["taxon_cluster_bootstrap_positive_fraction"]<=1
    assert not any("observer-" in str(v) for v in contrast.values())


def test_paired_kde_bandwidth_changes_absolute_overlap_not_cohort_names():
    f=mock_candidate(8,110)
    taxa=sorted(f.accepted_key.unique())
    scores,_=source_scores(f,taxa)
    result=bandwidth_sensitivity(scores,taxa,pairs_sample=5,reps=3,draws=32)
    assert result["status"]=="completed_matched_bandwidth_screen"
    assert result["same_taxon_pair_indices_for_all_three_modules"]
    for k,m in result["modules"].items():
        assert m["shared_sampled_taxon_pairs"]==5
        assert m["paired_rarefaction_replicates"]==3
        assert m["KDE_bandwidth_multipliers"]==[.75,1.,1.25]
        assert m["same_pairs_and_same_30_based_subsampling_seed_across_bandwidths"]
        assert all(0<=x<=1 for x in m["median_pairwise_KDE_overlap"])
        assert m["rank_spearman_vs_base"]["0.75"] is None or (
            -1<=m["rank_spearman_vs_base"]["0.75"]<=1
        )
    with pytest.raises(ValueError,match="Bandwidth numerical"):
        bandwidth_sensitivity(scores,taxa,pairs_sample=1,reps=2,draws=30)


def test_predeclared_observer_and_spatial_deletion_fail_closed_with_8_taxa():
    f=mock_candidate(8,110)
    taxa=sorted(f.accepted_key.unique())
    scores,_=source_scores(f,taxa)
    result=deletion_sensitivity(f,scores,taxa,reps=2,null_reps=2)
    assert result["status"]=="five_a_priori_source_deletion_analyses"
    assert result["original_observers_removed_count"]>0
    assert len(result["scenarios"])==5
    assert all(x["status"].startswith("not_comparable")
               for x in result["scenarios"].values())
    assert all("remaining_candidate_rows" in x
               for x in result["scenarios"].values())


def test_primary_robustness_never_silently_changes_40_taxon_contract():
    f=mock_candidate(8,110)
    with pytest.raises(ValueError,match="registered native taxa changed"):
        analyse(f,reps=6,null_reps=8,bootstrap=40)
    forty=mock_candidate(40,49)
    with pytest.raises(ValueError,match="fixed 50 rarefaction"):
        analyse(forty,reps=6,null_reps=8,bootstrap=40)
    forty=mock_candidate(40,50)
    with pytest.raises(ValueError,match="Insufficient declared robustness"):
        analyse(forty,reps=5,null_reps=8,bootstrap=40)
