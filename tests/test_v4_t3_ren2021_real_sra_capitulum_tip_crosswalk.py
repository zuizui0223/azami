"""Real 48 SRA source accession and plant head-stage noninvented morphology."""
import csv
import io
from pathlib import Path
import pytest
from analysis.v4 import t3_ren2021_real_sra_capitulum_tip_crosswalk as m

ROOT=Path(__file__).resolve().parents[1]/"analysis"/"v4"/"results"
ORIG=ROOT/"20261023_T3_SENECIONEAE_POSTURE_SOURCE_CANDIDATES.csv"
EXTRA=ROOT/"20261023_T3_REN2021_TRUE_SRA_ORIGINAL_POSTURE_SUPPLEMENT.csv"

def fake_runinfo():
    targets=list(m.FOCAL_SOURCE_STATES)
    samples=[]
    for i,name in enumerate(targets):
        samples.append({"run_accession":m.SOURCE_FOCAL_BY_SRR[name],
             "scientific_name":name,
             "sample_accession":f"SAMN165888{i+10}",
             "experiment_accession":f"SRX94200{i+10}",
             "fastq_ftp":f"ftp.sra.ebi.ac.uk/vol1/fastq/SRR129/fake/{i}_1.fastq.gz",
             "fastq_bytes":"23423"})
    for i in range(43):
        genus=("Ligularia" if i%2 else "Parasenecio")
        samples.append({"run_accession":f"SRR12969{i:03}",
             "scientific_name":f"{genus} fake_source_species{i}",
             "sample_accession":f"SAMN16580{i:03}",
             "experiment_accession":f"SRX9500{i}",
             "fastq_ftp":f"ftp.sra.ebi.ac.uk/vol1/fastq/SRR129/fake/{i}_1.fastq.gz",
             "fastq_bytes":"23423"})
    # Preserve genuine Cremanthodium four, not a synthetic fifth
    assert len(samples)==48
    return samples

def test_source_direct_additions_keep_species_critical_head_stage_unknown():
    original=m.csv_rows(ORIG)
    additional=m.csv_rows(EXTRA)
    assert len(original)==10
    assert len(additional)==5
    assert {r["taxon"]:r["source_posture"] for r in additional}==m.FOCAL_SOURCE_STATES
    humile=next(x for x in additional if x["taxon"]=="Cremanthodium humile")
    assert humile["source_posture"]=="?"
    emilia=next(x for x in additional if x["taxon"]=="Emilia sonchifolia")
    assert emilia["source_posture"]=="V"
    assert emilia["ontogenetic_direction"]=="D_before_to_U_later"
    assert m.SOURCE_FOCAL_BY_SRR["Cremanthodium campanulatum"]=="SRR12967014"


def test_synthetic_exact_48_has_only_three_direct_nodding_and_no_erect():
    rows=fake_runinfo()
    scored, counts=m.validate_ena(rows,m.csv_rows(ORIG),m.csv_rows(EXTRA))
    assert len(scored)==48
    assert counts=={"D":3,"V":1,"?":44}
    c={x["original_species_name_exact_from_ENA"]:x for x in scored}
    assert c["Cremanthodium humile"]["T3_documented_head_state_D_U_V_or_unknown"]=="?"
    assert c["Emilia sonchifolia"]["ontogenetic_direction_not_phylogenetic_event"]=="D_before_to_U_later"
    assert sum(x["DNA_sequence_reads_or_original_phylo_Newick_downloaded"] for x in scored)==0
    assert all(not x["identical_DNA_voucher_and_flowering_phenotype_confirmed"] for x in scored)


def test_refuse_vouchered_nonhomologous_or_fabricated_source_codes():
    rows=fake_runinfo()
    old=m.csv_rows(ORIG)
    extra=m.csv_rows(EXTRA)
    rec=[x.copy() for x in extra]
    rec[3]["source_posture"]="D"  # C. humile has only generic heads nodding
    with pytest.raises(ValueError,match="source codes changed"):
        m.validate_ena(rows,old,rec)
    bad=[x.copy() for x in rows]
    next(x for x in bad if x["scientific_name"]=="Cremanthodium campanulatum")["run_accession"]="SRR00000000"
    with pytest.raises(ValueError,match="original real SRA"):
        m.validate_ena(bad,old,extra)


def test_original48_source_is_independent_of_Cirsium_Cardueae_tree():
    assert m.NCBI_PRJNA=="PRJNA672812"
    assert m.EXPECTED_SOURCE_RUNS==48
    assert m.EXPECTED_CREMANTHODIUM==4
    assert m.OTHER_TRIBE_PROOF_ONLY_IF_CIRS_CIRSIUM_INDEPENDENT_TREE
    assert m.SOURCE_FOCAL_BY_SRR["Cremanthodium thomsonii"]=="SRR12967057"
