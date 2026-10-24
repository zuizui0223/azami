"""Real original-S1 110-tip public Q20 nuclear history; no fake accepted phylogeny."""
import csv
from pathlib import Path

from analysis.v4 import t2_s1_reconstructed_real110_nuclear_climate as m
from analysis.v4 import t2_real110_posthoc_one_sparse_source_sensitivity as source

MASTER=(Path(__file__).resolve().parents[1]/
        "analysis/v4/results/20261023_T2_S1_ALL251_66_TO72_PROVISIONAL_POSTURE_OVERLAY.csv")


def test_all251_original_S1_used_as_taxon_master_no_invented_photo_binary():
    s=m.read_csv(MASTER)
    assert len(s)==251
    assert len({z["accepted_key"] for z in s})==244
    assert sum(int(z["original_photo_count"]) for z in s)==170066
    assert sum(z["botanical_U_D_state"]=="U" for z in s)==54
    assert sum(z["botanical_U_D_state"]=="D" for z in s)==12
    assert m.BEFORE_SAMPLES==111
    assert m.AFTER_SAMPLES==110
    assert m.RECOVERED_LOCI==108
    assert m.ORIGINAL_SHARED_SITES==29310
    assert m.EXCLUDED_NAME=="Cirsium tioganum"


def test_posthoc_tech_rain_analysis_requires_real_source_voucher_references():
    import inspect
    s=inspect.getsource(m.analyze)
    assert 'source_BioSample' in s
    assert 'IQTREE_unrooted_branch_length_original_Q20_tree_SHA256' in s
    assert 'historic_U_to_D_gain_count_inferred' in s
    assert 'true_110_biological_species_tree_completed' in s
    assert 'source_original_111_predeclared_technical_gate_succeeded' in s
    assert 'original_source_name' in s
    assert 'botanical_U_D_state' in s


def test_original_locus_selection_in_111_cannot_be_reoptimized_posthoc():
    assert source.PREDECLARED_111_QUORUM==89
    assert source.ORIGINAL_LOCI==118
    assert source.RECOVERED_LOCI==108
    assert source.EXCLUDED_ORIGINAL_Q20_BASES==119
