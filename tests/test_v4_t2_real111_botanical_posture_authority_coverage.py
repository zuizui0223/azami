"""Exact botanical-source U/D coverage, original full251, and real111 sample tests."""
from collections import defaultdict
import pytest

from analysis.v4 import t2_real111_botanical_posture_authority_coverage as m


def test_authority_agreement_unknown_and_conflict_preserved():
    sources=defaultdict(list)
    m.add(sources,"Cirsium test0","U","original_authority","https://example.org/U","exact description","U")
    m.add(sources,"Cirsium test0 var. other","D","different_variety","https://example.org/var","variety","D")
    assert m.infer(sources,"Cirsium_test0")["botanical_U_D_state"]=="U"
    assert m.infer(sources,"Cirsium test0 var. other")["botanical_U_D_state"]=="D"
    m.add(sources,"Cirsium test0","D","independent_conflict","https://example.org/D","same exact taxon name","D")
    item=m.infer(sources,"Cirsium_test0")
    assert item["botanical_U_D_state"]=="?"
    assert item["botanical_state_status"]=="conflicting_independent_botanical_sources"
    assert not item["botanical_U_D_same_public_voucher_verified"]

    m.add(sources,"Cirsium test1","1","published_HeadPosition","https://example.org/NA","NA63","1")
    m.add(sources,"Cirsium test1","?","uncertain","https://example.org/und","ambiguity","?")
    assert m.infer(sources,"Cirsium_test1")["botanical_U_D_state"]=="D"
    assert m.infer(sources,"Cirsium_test1")["any_explicit_unresolved_source"]


def test_original251_and_real111_taxon_name_exact_join_preserves_all_photos():
    original=[]
    for i in range(244):
        key=str(200000000+i)
        for k in range(2 if i<7 else 1):
            original.append({
                "accepted_key":key,
                "S1_original_source_name":f"Cirsium_test{i:03d}_{k}" if i<7 else f"Cirsium_test{i:03d}",
                "original_S1_photo_count":"680",
            })
    assert len(original)==251
    original[0]["original_S1_photo_count"]="66"
    # To make 170066: 250 * 680 + 66
    nuclear=[]
    seen=set()
    for row in original:
        key=row["accepted_key"]
        if key in seen or len(nuclear)>=111:
            continue
        seen.add(key)
        nuclear.append({
            "accepted_key":key,
            "taxon_name":row["S1_original_source_name"].replace("_"," "),
            "run":"SRR"+str(310000000+len(nuclear)),
            "biosample":"SAMN"+str(410000000+len(nuclear)),
            "phylogenetic_tip":"WCVP"+key,
        })
    source=defaultdict(list)
    m.add(source,"Cirsium test000 0","U","NMNS","https://example.org/a","source","U")
    m.add(source,"Cirsium test001 0","D","NA","https://example.org/b","source","1")
    result,report=m.join(original,nuclear,source)
    assert len(result)==251
    assert report["source_original_S1_accepted_keys"]==244
    assert report["source_accepted_name_collisions"]==7
    assert report["source_real_120k_nuclear_sample_tips"]==111
    assert report["nuclear111_with_nodding"]==1
    assert report["nuclear111_with_erect"]==1
    assert report["nuclear111_without_verified_taxon_concept_state"]==109
    assert report["molecular_BioSample_photophenotype_same_voucher_assessed"] is False
    assert result[0]["original_source_name"]=="Cirsium test000 0"


def test_unknown_original_named_taxon_is_not_reassigned_by_binomial():
    index=defaultdict(list)
    m.add(index,"Cirsium japonicum","U","broad","https://x","broad concept","U")
    assert m.infer(index,"Cirsium japonicum var. horridum")["botanical_U_D_state"]=="?"
    assert m.infer(index,"Cirsium japonicum")["botanical_U_D_state"]=="U"


def test_all_original_external_eazami_sources_have_frozen_git_blob_contracts():
    assert len(m.SOURCES)==4
    assert [n for *_,n in m.SOURCES.values()]==[22,2,20,13]
    assert m.NA_BOTANICAL==63 and m.NA_NODDING==5
    assert m.ORIGINAL_UNITS==251 and m.ORIGINAL_IMAGES==170066
    assert m.MEASURED_NUCLEAR==111


def test_real_111_sra_sample_is_not_same_voucher_botanical_measurement():
    o=[{"accepted_key":"10","S1_original_source_name":"Cirsium_test",
        "original_S1_photo_count":"170066"}]
    n=[{"accepted_key":"10","taxon_name":"Cirsium test","run":"SRR10",
        "biosample":"SAMN10","phylogenetic_tip":"WCVP10"}]
    with pytest.raises(ValueError,match="S1 full251"):
        m.join(o,n,defaultdict(list))
