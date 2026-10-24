"""Offline exact-name/grafting-boundary tests for S1×Moreyra nuclear reads."""
import pandas as pd
import pytest

from analysis.v4.t2_all244_eazami_moreyra_read_crosswalk import (
    analyse,name_representation,
)


def synthetic_roster():
    rows=[]
    for i in range(244):
        name=("Cirsium_arvense" if i==0 else
              "Cirsium_vulgare" if i==1 else f"Cirsium_test{i:03d}")
        coll=191<=i<194
        if coll:
            name+="|Cirsium_extra"+str(i)
        rows.append({
            "accepted_key":str(1000+i),
            "source_taxon_names":name,
            "n_source_orientation_photos":169823 if i==0 else 1,
            "tree194_covered":i<194,
            "strict_phylogenetic_eligible":i<191,
            "source_name_collision":coll,
            "s2_focal":i in (0,1),
        })
    return pd.DataFrame(rows)


def test_taxon_coverage_and_s2_original_read_name_match():
    reads=[
        {"tree_code":"Cirsium arvense","run":"SRR1","biosample":"SAMN1"},
        {"tree_code":"Cirsium vulgare 1","run":"SRR2","biosample":"SAMN2"},
        {"tree_code":"Cirsium test193","run":"SRR3","biosample":"SAMN3"},
    ]
    frame,summary=analyse(synthetic_roster(),reads)
    assert len(frame)==244
    assert summary["n_S1_taxa_with_public_raw_read_name_match"]==2
    assert summary["n_S1_taxa_without_grafted_tree_but_with_Moreyra_public_reads"]==0
    assert summary["source_S1_strict_clean_tree"]==191
    a=frame.loc[frame.original_source_taxon_names=="Cirsium arvense"].iloc[0]
    v=frame.loc[frame.original_source_taxon_names=="Cirsium vulgare"].iloc[0]
    assert a["moreyra_read_name_match_status"]=="public_original_tree_code_exact_name"
    assert v["moreyra_read_name_match_status"]=="public_original_tree_code_numbered_replicate_only"
    assert not a["sample_name_match_is_validated_tree_tip"]
    assert all(not row["source_taxonomic_concept_equivalence_validated"] for row in frame.to_dict("records"))


def test_infraspecific_names_not_collapsed():
    assert name_representation("Cirsium vulgare 3")==("Cirsium vulgare 3","Cirsium vulgare")
    assert name_representation("Cirsium japonicum var. australe")==(
        "Cirsium japonicum var. australe","Cirsium japonicum var. australe")


def test_missing_source_species_or_total_photos_blocks():
    v=synthetic_roster()
    v.loc[0,"n_source_orientation_photos"]=1
    with pytest.raises(ValueError,match="S1 original source cohort"):
        analyse(v,[])
