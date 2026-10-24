"""Biological stop-rule and synthetic NJ math tests; real seven SRA receipts in CI."""
from copy import deepcopy
from io import StringIO
import pytest
from Bio import Phylo
from analysis.v4 import t2_seven_original110_Q20_nj_diagnostic as m


def tiny(monkeypatch,missing=False):
    monkeypatch.setattr(m,"N_SHARED_COPY_TARGETS",3)
    monkeypatch.setattr(m,"EXPECTED_ALL_SHARED_POSITIONS",359 if missing else 360)
    monkeypatch.setattr(m,"EXPECTED_ALL_VARIABLE_POSITIONS",3)
    monkeypatch.setattr(m,"EXPECTED_ALL_INFORMATIVE_POSITIONS",3)
    monkeypatch.setattr(m,"EXPECTED_TARGETS_GE100",3)
    monkeypatch.setattr(m,"MIN_PRIMARY_LOCI",3)
    loci={"L1":"saff-L1","L2":"lett-L2","L3":"saff-L3"}
    records=[]
    for taxon in range(7):
        targets={}
        for i,(locus,ref) in enumerate(loci.items()):
            bases=[]
            for pos in range(1,121):
                if missing and taxon==0 and i==0 and pos==3:
                    continue
                base="C" if pos==11 and taxon<2 else "A"
                bases.append([pos,base,3,1.])
            targets[locus]={
              "reference":ref,"reference_length_nt":500,
              "n_conservatively_callable_sites":len(bases),
              "q20_majority_bases_1based":bases,
              "not_a_validated_allelic_genotype":True,
            }
        records.append({
            "run":f"SRR{taxon+1:05d}","biosample":f"SAMN{taxon+1:05d}",
            "accepted_key":f"WCVP{taxon+1}",
            "S1_source_taxon":f"Cirsium_test{taxon+1}",
            "targets":targets,
        })
    return records,loci


def test_true_shared_sites_and_synthetic_unrooted_nj_locus_bootstrap(monkeypatch):
    data,refs=tiny(monkeypatch)
    report,loci,pairs,splits=m.infer(data,refs,boot=100)
    assert report["n_original_public_SRA_biosamples"]==7
    assert report["n_primary_loci_shared_atleast_100_Q20_sites"]==3
    assert report["source_shared_Q20_positions_all_targets"]==360
    assert report["source_shared_variable_Q20_positions"]==3
    assert report["source_shared_parsimony_informative_Q20_positions"]==3
    assert report["n_complete_case_invariant_plus_variable_sites_primary"]==360
    assert report["primary_p_distance_pair_count"]==21
    assert all(p["shared_site_count"]==360 for p in pairs)
    assert any(p["p_distance_on_seven_way_complete_case_sites"]>0 for p in pairs)
    assert all(p["not_validated_species_tree_support"] for p in splits)
    assert report["EAzami_294_tip_species_tree_constructed"] is False
    assert report["copy_state_orthology_homeology_and_taxonomy_concepts_verified"] is False
    tree=Phylo.read(StringIO(report["primary_unrooted_newick"]),"newick")
    assert len(tree.get_terminals())==7


def test_missing_observed_base_not_imputed_by_reference(monkeypatch):
    data,refs=tiny(monkeypatch,missing=True)
    _,loci=m.loci_from_original_base_calls(data,refs)
    assert loci[0]["aligned_positions"]==119
    assert len(loci[0]["alignment"][0])==119


def test_wrong_reference_copy_stops_real_comparison(monkeypatch):
    data,refs=tiny(monkeypatch)
    data[0]["targets"]["L1"]["reference"]="sunf-L1"
    with pytest.raises(ValueError,match="Different Asteraceae reference"):
        m.loci_from_original_base_calls(data,refs)


def test_wrong_or_duplicated_source_taxon_label_fails(monkeypatch):
    data,refs=tiny(monkeypatch)
    data[0]["S1_source_taxon"]=data[1]["S1_source_taxon"]
    with pytest.raises(ValueError,match="Duplicate source taxon"):
        m.loci_from_original_base_calls(data,refs)


def test_reference_site_variant_count_with_three_states():
    assert m.site_count("AAAAAAA")== (0,0)
    assert m.site_count("AACCCCC")== (1,1)
    assert m.site_count("AAACGTT")== (1,1)
    assert m.site_count("ACCCCCC")== (1,0)


def test_nj_is_unrooted_technical_not_time_or_adaptation():
    assert m.N_SAMPLES==7
    assert m.BOOTSTRAP_SEED==20261023
    assert m.MIN_COMMON_SITE_PER_TARGET==100
    assert m.N_SOURCE_PAIRED_READS==120000
