"""T3 Chinese species Flora head-down D admission and genuine 2025 nuclear ITS source."""
from pathlib import Path
from analysis.v4 import t3_131_chinese_flora_new_Ligularia_nodding as m

ROOT=Path(__file__).resolve().parents[1]/"analysis"/"v4"/"results"


def test_chinese_original_flora_says_HEAD_nodding_not_stem():
    r=m.csv_rows(ROOT/"20261024_T3_LIGULARIA_PURDOMII_ORIGINAL_CHINESE_FLORA_HEAD.csv")
    assert len(r)==1
    x=r[0]
    assert x["taxon"]=="Ligularia purdomii"
    assert x["chinese_species_head_phrase"]==m.EXPECTED_FOC_HEAD
    assert x["head_posture"]=="D"
    assert x["english_FOC_describes_erect_stem_not_head"]=="yes"
    assert x["does_english_FOC_assert_nodding"]=="no"
    assert x["double_blind_head_review_completed"]=="no"
    assert x["DNA_and_flora_same_flowering_specimen_verified"]=="no"
    assert x["independent_U_to_D_origin_or_rain_utility_demonstrated"]=="no"


def test_published_2025_genbank_has_actual_nuclear_ITS_no_ETS_for_this_nodding_species():
    a=m.csv_rows(ROOT/"20261023_T3_TAIWANIA2025_ORIGINAL163_ITS_ETS_GENBANK_ACCESSIONS.csv")
    assert len(a)==163
    matches=[z for z in a if z["original_published_scientific_name"]==m.FLORA]
    assert len(matches)==1
    assert matches[0]["ITS_GenBank_accession"]==m.ITS_ACC
    assert matches[0]["ETS_GenBank_accession"]==""
    assert matches[0]["ITS_and_ETS_independent_unlinked_gene_loci"]=="no"


def test_genuine_original_fourtip_LCA_is_not_a_declared_two_tip_monophyletic_species_pair():
    assert m.REQUIRED_LCA_TIPS=={
        "Cremanthodium_campanulatum",
        "Ligularia_yunnanensis",
        "Ligularia_curvisquama",
        "Ligularia_purdomii",
    }
    assert len(m.REQUIRED_LCA_TIPS)==4
    assert m.SHELTER=="Cremanthodium campanulatum"
    assert m.N_GENUINE_GB_ITS_ETS==131
    assert m.EXPECTED_FOC_HEAD.endswith("下垂")
    assert len(m.TREE_SHA)==64
    assert m.FRPS_CHINESE_URL.startswith("https://www.iflora.cn/")
    assert m.FOC_ENGLISH_URL.startswith("https://www.efloras.org/")
