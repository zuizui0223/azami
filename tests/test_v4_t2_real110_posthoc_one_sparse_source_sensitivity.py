"""Synthetic regression: preserve original111 failure and one-tip posthoc110 source mask.

No fake genotypes are substituted into the actual source workflow. These tests
verify posthoc choice, 89/111 site quorum, original 108 loci and full source scope.
"""
from collections import Counter

import pytest

from analysis.v4 import t2_real110_posthoc_one_sparse_source_sensitivity as m


def original111_case():
    sources={}
    seqs={}
    broad={}
    sample_rows=[]
    named=[
        "Cirsium conspicuum","Cirsium scopulorum",
        "Cirsium flodmanii","Cirsium tracyi","Cirsium vulgare",
    ]
    for i in range(110):
        tip=f"WCVP{206710000+i}"
        taxon=named[i] if i<len(named) else f"Cirsium synthetic{i:03d}"
        suffix="C" if i%2 else "A"
        bases=suffix*200+"A"*(m.ORIGINAL_SHARED_SITES-200)
        seqs[tip]=bases
        sources[tip]={
            "tip":tip,"taxon_name":taxon,
            "source_original_run":f"SRR77000{i:04d}",
            "source_BioSample":f"SAMN88000{i:04d}",
            "n_real_read_pairs":120000
        }
        broad[tip]={"n_Q20_bases_in_80pct_site_alignment":len(bases)}
        sample_rows.append({"phylogenetic_tip":tip})
    badtip=m.EXCLUDED_TIP
    badseq="A"*m.EXCLUDED_ORIGINAL_Q20_BASES + "N"*(m.ORIGINAL_SHARED_SITES-m.EXCLUDED_ORIGINAL_Q20_BASES)
    seqs[badtip]=badseq
    sources[badtip]={
        "tip":badtip,"taxon_name":m.EXCLUDED_NAME,
        "source_original_run":"SRR99999991",
        "source_BioSample":"SAMN99999991",
        "n_real_read_pairs":120000
    }
    broad[badtip]={"n_Q20_bases_in_80pct_site_alignment":m.EXCLUDED_ORIGINAL_Q20_BASES}
    sample_rows.append({"phylogenetic_tip":badtip})
    all_genes=[]
    for i in range(118):
        selected=i<108
        sites=271+(1 if i<42 else 0) if selected else 0
        all_genes.append({
            "locus":f"CP{i:04d}","original_reference_copy":f"lett-CP{i:04d}",
            "retain_min100_Q20_positions":selected,
            "n_sites_same_reference_called_at_least80pct_taxa":sites,
        })
    assert sum(g["n_sites_same_reference_called_at_least80pct_taxa"] for g in all_genes if g["retain_min100_Q20_positions"])==m.ORIGINAL_SHARED_SITES
    data=list(sources.values())
    return sample_rows,data,all_genes,seqs,broad,{"technical_ML_alignment_source_gate_passed":False}


def test_one_known_low_depth_specimen_only_retained_as_missing():
    a=original111_case()
    seq,rows,bounds,receipt=m.sensitivity(*a)
    assert len(seq)==110 and m.EXCLUDED_TIP not in seq
    assert len(rows)==110 and len(bounds)==108
    assert sum(int(z["end_1based"])-int(z["start_1based"])+1 for z in bounds)==29310
    assert receipt["original_111_reference_position_selection_recomputed_on110"] is False
    assert receipt["original111_technical_gate_passed"] is False
    assert receipt["original_111_predeclared_89_of111_quorum_unmodified"]==89
    assert receipt["posthoc110_actual_parsimony_informative_sites"]>=50
    assert set(m.REQUIRED_NODDING)<=set(x["taxon_name"] for x in rows)
    assert not receipt["copy_number_paralogy_polyploidy_review_completed"]
    assert not receipt["global_EAzami294_true_species_tree_accepted"]


def test_two_original_low_coverage_samples_must_fail_not_be_dropped():
    a=list(original111_case())
    a[4]["WCVP206710000"]["n_Q20_bases_in_80pct_site_alignment"]=900
    with pytest.raises(ValueError,match="Exactly one genuine"):
        m.sensitivity(*a)


def test_original_observed_sites_cannot_be_reselected_posthoc():
    a=list(original111_case())
    # Original 89/111 mask is not recomputed after selecting 110.
    a[2][0]["n_sites_same_reference_called_at_least80pct_taxa"]-=1
    with pytest.raises(ValueError,match="source locus/position mask changed"):
        m.sensitivity(*a)
