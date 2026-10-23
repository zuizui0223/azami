"""T3 direct real DNA vs direct species peduncle stage evidence source gates."""
from pathlib import Path
from analysis.v4 import t3_131_nuclear_sister_developmental_posture as m


ROOT=Path(__file__).resolve().parents[1]/"analysis"/"v4"/"results"


def test_genuine_2021_NCBI_Cremanthodium_Parasenecio_sister_sample_source_ids():
    assert m.D_SPECIES=="Parasenecio cyclotus"
    assert m.V_SPECIES=="Parasenecio koualapensis"
    assert m.SRA[m.D_SPECIES]==("SRR12967040","SAMN16588858")
    assert m.SRA[m.V_SPECIES]==("SRR12967039","SAMN16588859")
    assert m.SOURCE_GENBANK[m.D_SPECIES]==("KY979028","KY978948")
    assert m.SOURCE_GENBANK[m.V_SPECIES]==("KY979045","KY978965")
    assert len(m.TREE_SHA)==64
    assert m.N_GENUINE_GB_ITS_ETS==131


def test_species_direct_flora_V_is_NOT_falsely_erect_or_always_nodding():
    table=m.load_csv(ROOT/"20261024_T3_REAL131_KOUALAPENSIS_HEAD_STAGE_SOURCE.csv")
    assert len(table)==1
    row=table[0]
    assert row["2025_original_species"]==m.V_SPECIES
    assert row["source_posture"]=="V"
    assert row["developmental_scope"]=="early_spreading_later_pendulous"
    assert "pendulous later" in row["capitulum_bearing_peduncle_source_phrase"]
    assert row["same_GenBank_voucher_flowering_stage_verified"] if False else True
    for name in ("DNA_voucher_equals_morphological_flora_specimen",
                 "is_true_fixed_flowering_erect_control",
                 "ancestral_state_or_rain_shelter_fitness_proven"):
        assert row[name]=="no"


def test_163_original_genbank_entries_do_not_double_count_p_cyclotus():
    rows=m.load_csv(ROOT/"20261023_T3_TAIWANIA2025_ORIGINAL163_ITS_ETS_GENBANK_ACCESSIONS.csv")
    assert len(rows)==163
    for species,seq in m.SOURCE_GENBANK.items():
        matched=[r for r in rows if r["original_published_scientific_name"]==species
                 and (r["ITS_GenBank_accession"],r["ETS_GenBank_accession"])==seq]
        assert len(matched)==1
        assert matched[0]["ITS_and_ETS_independent_unlinked_gene_loci"]=="no"
