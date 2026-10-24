"""No borrowed tree: only accession-verified S1 original source labels can enter."""
import csv
import json
from pathlib import Path
import pytest
from Bio import Phylo

from analysis.v4 import t2_s1_genbank_its_tree_rain as m


def sample(n=16):
    r={}
    for i in range(n):
        chars=list(("ACGT"*80))
        for j in range(25):
            chars[3*j]=(lambda z:"ACGT"[(("ACGT".index(z)+i%3)%4)])(chars[3*j])
        r[f"Cirsium_species{i:02}"]="".join(chars)
    return r


def test_identical_columns_and_source_tip_count_have_real_variable_sites():
    s=sample()
    result, meta=m.trim_alignment(s)
    assert meta["status"]=="PASS_TO_GENBANK_ITS_IQTREE_EXPLORATORY_ML"
    assert len(result)==16
    assert meta["selected_real_ITS_alignment_columns"]==320
    assert meta["selected_variable_nucleotide_sites"]>=10
    assert meta["selected_parsimony_informative_sites"]>=10
    assert all(len(seq)==320 for seq in result.values())
    assert not meta["single_copy_nuclear_orthology_validated"]


def test_missing_ITS_comparability_does_not_gain_taxa_or_fake_gene():
    result,meta=m.trim_alignment(sample(8))
    assert not result and meta["status"]=="BLOCKED_TOO_FEW_EXACT_S1_GENBANK_ITS_TAXA"
    a=sample()
    for k in a:
        a[k]="ACGT"*80
    result,meta=m.trim_alignment(a)
    assert not result and meta["status"]=="BLOCKED_ALIGNMENT_WITH_TOO_LITTLE_REAL_ITS_INFORMATION"


def test_partial_nonoverlap_is_not_force_aligned():
    a=sample()
    a["Cirsium_species00"]="N"*280+"ACGT"*10
    result,meta=m.trim_alignment(a)
    assert not result
    assert meta["status"]=="BLOCKED_ONE_OR_MORE_SELECTED_ITS_SPECIES_HAS_TOO_LITTLE_OVERLAP"


def test_only_full251_all_sources_success_can_be_used_for_ml(tmp_path):
    n=12
    accessions=[]
    for i in range(n):
        taxon=f"Cirsium species{i:02}"
        accessions.append({"original_source_name":taxon,
                           "source_original_species":taxon,
                           "S1_species_to_GB_record_exact_name":"True",
                           "chelsa_bio12":"1000","chelsa_bio01":"18",
                           "botanical_U_D_state":"D" if i<4 else "U"})
    msa=tmp_path/"mafft.fa"
    msa.write_text("".join(f">{x['original_source_name'].replace(' ','_')}\nACGT\n"
                           for x in accessions))
    receipt=tmp_path/"receipt.json"
    receipt.write_text(json.dumps({
        "source_original_S1_master_names":251,
        "source_original_WCVP_accepted_keys":244,
        "scope":"S1_PREFIX_PILOT",
        "all251_queries_completed_without_network_failures":False,
    }))
    with pytest.raises(ValueError,match="Cannot construct full-S1-source"):
        m.audit_accessions(accessions,receipt,msa)
    r=json.loads(receipt.read_text())
    r["scope"]="S1_ALL251"
    r["all251_queries_completed_without_network_failures"]=True
    receipt.write_text(json.dumps(r))
    with pytest.raises(ValueError,match="Malformed or missing"):
        m.audit_accessions(accessions,receipt,msa)


def test_root_direction_and_full_multilocus_species_tree_not_claimed():
    assert m.MIN_TAXA_FOR_ITS_ML>=10
    assert m.MIN_PHYLO_POSTURE_STATES==(3,5)
    assert m.MAX_MISSING_ITS_TIP_FRAC<0.6
