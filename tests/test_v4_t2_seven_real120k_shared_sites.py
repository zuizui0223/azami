"""No fake genotype creation; empirical site intersection only after each source run."""
import pytest
from analysis.v4 import t2_seven_real120k_shared_sites as m


def test_q20_majority_is_not_reference_imputation_or_heterozygote_forcing():
    assert m.majority_base_call([]) is None
    assert m.majority_base_call(["", "N", "A", "A", "A"])==("A",3,1.0)
    assert m.majority_base_call(["A","A","G"]) is None
    assert m.majority_base_call(["C","C","C","C","T","C","C","C","C","C"])[0]=="C"
    assert m.majority_base_call(["A","A"]) is None


def test_original_110_reference_copy_gate(monkeypatch):
    monkeypatch.setattr(m,"N_SHARED_TARGETS",2)
    rows=[]
    for i in range(241):
        locus=f"L{i}"
        enabled=i<2
        rows.append({
            "locus":locus,
            "deep7_all_Q20_dp3_ge25pct_same_ref":str(enabled),
            "best_reference_copy_identical_all7":str(enabled),
            "original_ref_names":f"saff-{locus}",
        })
    selected=m.targets(rows)
    assert selected=={"L0":"saff-L0","L1":"saff-L1"}
    rows[1]["original_ref_names"]="lett-WRONG"
    with pytest.raises(ValueError,match="reference"):
        m.targets(rows)


def test_all_seven_real_public_sra_voucher_ids_not_duplicated(monkeypatch):
    monkeypatch.setattr(m,"N_ORIGINAL_RUNS",7)
    source=[];chosen=[]
    for i,run in enumerate(sorted(m.SEVEN)):
        src=dict(run=run,biosample=f"SAMN{i}",accepted_key=f"K{i}",
                 S1_source_taxon=f"Cirsium_test_{i}")
        source.append(src)
        chosen.append(dict(run=run,biosample=f"SAMN{i}",accepted_key=f"K{i}",
                           source_taxon=f"Cirsium_test_{i}",same_run_120k_dp3_ge50pct="100"))
    assert m.source_sample(source,chosen,source[0]["run"])==source[0]
    chosen[0]["biosample"]="SAMNWRONG"
    with pytest.raises(ValueError,match="BioSample"):
        m.source_sample(source,chosen,source[0]["run"])


def demo_seven(monkeypatch):
    monkeypatch.setattr(m,"N_SHARED_TARGETS",2)
    contigs={"L1":"saff-L1","L2":"lett-L2"}
    seven=[]
    for i,run in enumerate(sorted(m.SEVEN)):
        a="G" if i<3 else "A"
        one=[[pos,(a if pos==3 else "A"),4,1.] for pos in range(1,101)]
        two=[[pos,"C",5,1.] for pos in range(1,121 if i<6 else 91)]
        seven.append({
            "run":run,"biosample":f"SAMN{i}","accepted_key":f"K{i}",
            "n_source_pairs":120000,
            "n_original_shared_reference_targets":2,
            "targets":{
                "L1":{"reference":"saff-L1","reference_length_nt":350,
                      "not_a_validated_allelic_genotype":True,
                      "q20_majority_bases_1based":one},
                "L2":{"reference":"lett-L2","reference_length_nt":250,
                      "not_a_validated_allelic_genotype":True,
                      "q20_majority_bases_1based":two}
            }
        })
    return contigs,seven


def test_real_identical_reference_site_overlap_and_parsimony_information(monkeypatch):
    contigs,samples=demo_seven(monkeypatch)
    for s in samples:
        assert m.validate_single(s,s,contigs)
    loci,result=m.aggregate(samples,contigs)
    assert len(loci)==2
    one=next(x for x in loci if x["locus"]=="L1")
    two=next(x for x in loci if x["locus"]=="L2")
    assert one["n_shared_Q20_majority_sites_all7"]==100
    assert one["n_shared_variable_sites_all7"]==1
    assert one["n_shared_parsimony_informative_sites_all7"]==1
    assert two["n_shared_Q20_majority_sites_all7"]==90
    assert result["total_shared_Q20_majority_positions_all7"]==190
    assert result["total_shared_variable_positions_all7"]==1
    assert result["n_targets_all7_shared_Q20_majority_atleast100_sites"]==1
    assert not result["species_tree_inferred"]


def test_site_unvalidated_source_phasing_stops_invalid_base_calls(monkeypatch):
    contigs,samples=demo_seven(monkeypatch)
    samples[0]["targets"]["L1"]["q20_majority_bases_1based"][0][1]="N"
    with pytest.raises(ValueError,match="Malformed"):
        m.validate_single(samples[0],samples[0],contigs)


def test_missing_one_real_run_never_imputes_its_alleles(monkeypatch):
    contigs,samples=demo_seven(monkeypatch)
    with pytest.raises(ValueError,match="All seven"):
        m.aggregate(samples[:-1],contigs)


def test_never_substitute_four_reads_for_three_phred_q20_independent_bases():
    assert m.MIN_SUPPORT==3
    assert m.Q20==20
    assert m.MIN_MAJOR_FRACTION==.9
