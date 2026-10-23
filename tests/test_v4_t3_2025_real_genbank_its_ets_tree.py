"""Original GenBank ITS / ETS are LINKED rDNA, not two independent loci."""
import io
from pathlib import Path
import numpy as np
import pytest
from Bio import Phylo
from analysis.v4 import t3_2025_real_genbank_its_ets_tree as m

def fake_alignment(n=12,seq_len=320):
    out={}
    for j in range(n):
        x=list("ACGT"*(seq_len//4))
        for p in range(45):
            x[p*2]="C" if j%2 else "A"
        out[f"Parasenecio fake{j:02}"]="".join(x)
    return out


def test_synthetic_actual_locus_alignment_keeps_original_observed_site_occupancy():
    actual,source=m.trim_marker(
        fake_alignment(),occupancy=.55,min_tips=8,
        min_information=12,kind="ITS")
    assert source["status"]=="GENBANK_SOURCE_RIBOSOMAL_REGION_MSA_QC_PASS"
    assert source["source_tips"]==12 and source["retained_sites"]==320
    assert source["parsimony_informative_sites"]>=40
    assert len(actual)==12
    assert source["ambiguous_real_genbank_bases_are_N_NOT_reference_fills"]


def test_rDNA_site_occupancy_excludes_unsupported_fake_reference_bases():
    a=fake_alignment()
    for n in list(a)[:8]:
        a[n]=a[n][:180]+"-"*140
    filtered,source=m.trim_marker(a,occupancy=.55,min_tips=8,
                                   min_information=8,kind="ETS")
    assert source["status"]=="GENBANK_SOURCE_RIBOSOMAL_REGION_MSA_QC_PASS"
    assert source["retained_sites"]==180
    assert all(len(seq)==180 for seq in filtered.values())
    assert all("N" not in seq for seq in filtered.values())


def test_site_informative_requires_distinct_observed_nucleotides():
    a={name:"ACGT"*80 for name in fake_alignment()}
    result,receipt=m.trim_marker(a,occupancy=.55,min_tips=8,
                                  min_information=10,kind="ITS")
    assert not result
    assert receipt["status"]=="INSUFFICIENT_HOMOLOGOUS_REAL_RDNA_SITES"


def test_botanical_direct_source_does_not_force_erect_states_onto_unknowns():
    states=m.original_states()
    assert states["Cremanthodium campanulatum"]["posture"]=="D"
    for species in ("Cremanthodium lineare","Parasenecio latipes",
                    "Parasenecio cyclotus","Parasenecio roborowskii",
                    "Parasenecio lidjiangensis"):
        assert states[species]["posture"]=="D"
    assert states["Emilia sonchifolia"]["posture"]=="V"
    assert states["Cremanthodium humile"]["posture"]=="?"
    assert states["Xenophyllum esquilachense"]["posture"]=="U"
    assert "Ligularia duciformis" not in states


def test_frozen_published_163_rows_not_163_biological_species():
    from analysis.v4.t3_public_2025_its_ets_genbank_reconstruct import validate_original
    original=m.rows(m.SOURCE_TABLE_PATH)
    validate_original(original)
    assert len(original)==163
    assert len({r["original_published_scientific_name"] for r in original})==154
    assert m.MIN_ITS_TIPS==100
    assert m.MIN_ETS_TIPS==50
    assert 0<m.OCCUPANCY_ITS<1 and 0<m.OCCUPANCY_ETS<1
