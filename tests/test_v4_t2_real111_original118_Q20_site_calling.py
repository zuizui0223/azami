"""Original111 real nuclear 118-copy Q20 site/source safeguards (synthetic data)."""
from collections import defaultdict
import pytest
from analysis.v4 import t2_real111_original118_Q20_nucleotide_calls as s
from analysis.v4 import t2_real111_original118_shared_Q20_aggregate as a


def source_fixture(monkeypatch):
    monkeypatch.setattr(s,"N_SOURCE",13)
    monkeypatch.setattr(s,"N_MEASURED",11)
    monkeypatch.setattr(s,"N_TARGETS",3)
    monkeypatch.setattr(s,"QUORUM",9)
    names=["Cirsium conspicuum","Cirsium scopulorum"]+[f"Cirsium test{i}" for i in range(9)]
    source=[]
    samples=[]
    for i,name in enumerate(names):
        key=str(206000001+i);run=f"SRR{25200000+i}";bio=f"SAMN{34200000+i}"
        source.append({"run":run,"biosample":bio,"accepted_key":key,
                       "S1_source_taxon":name})
        samples.append({"run":run,"biosample":bio,"accepted_key":key,
                        "phylogenetic_tip":"WCVP"+key,"taxon_name":name})
    source.extend([
       {"run":"SRR999000001","biosample":"SAMN999000001","accepted_key":"707010001",
        "S1_source_taxon":"Cirsium absent"},
       {"run":"SRR999000002","biosample":"SAMN999000002","accepted_key":"707010002",
        "S1_source_taxon":"Cirsium absentother"}
    ])
    target=[{"locus":f"L{i}","original_compositae_reference":f"saff-L{i}",
             "empirically_supported_source_taxa":"9",
             "minimum_source_taxa_required":"9"} for i in range(3)]
    return source,samples,target


def test_predeclared_real_runs_are_one_per_biosample_with_two_nodding_anchors(monkeypatch):
    source,samples,targets=source_fixture(monkeypatch)
    entries,refs=s.source_inputs(source,samples,targets)
    assert len(entries)==11 and len(refs)==3
    assert {x[0]["taxon_name"] for x in entries}>=s.BOTANICAL_NODDING
    source[0]["biosample"]="SAMN_wrong"
    with pytest.raises(ValueError,match="BioSample"):
        s.source_inputs(source,samples,targets)


def test_one_copy_quorum_fixed_before_seeing_basecalls(monkeypatch):
    original,samples,targets=source_fixture(monkeypatch)
    targets[0]["empirically_supported_source_taxa"]="8"
    with pytest.raises(ValueError,match="lacks 89"):
        s.source_inputs(original,samples,targets)


def test_full111_shard_preserves_every_original_sample_without_sequence_imputation(monkeypatch):
    original,samples,targets=source_fixture(monkeypatch)
    monkeypatch.setattr(s,"N_SHARDS",1)
    def mock(one,org,contigs,original241,reference,scratch):
        return {
            "source_original_run":one["run"],
            "tip":one["phylogenetic_tip"],
            "n_Q20_majority_called_sites":50,
            "n_targets_at_least100_Q20_called_sites":0,
            "loci":{l:{"reference":ref} for l,ref in contigs.items()},
        }
    result=s.collect(original,samples,targets,[f"L{i}" for i in range(241)],None,None,0,runner=mock)
    assert result["requested"]==11 and result["success"]==11
    assert result["failure"]==0
    assert not result["genus_species_phylogeny_inferred"]


def test_real_original_positions_keep_uncalled_N_and_not_replace_reference(monkeypatch):
    monkeypatch.setattr(a,"MIN_SITES_PER_LOCUS",3)
    samples=[]
    for i in range(5):
        seen=[]
        for j in range(1,6):
            if i==4 and j==5:continue
            val="C" if j==3 and i in (2,3) else "A"
            seen.append([j,val,3,.95])
        samples.append({"tip":f"WCVP{i}","loci":{
            "L0":{"reference":"saff-L0","reference_length_nt":50,
                  "q20_majority_bases_1based":seen}}})
    genes,alignment,breadth,summary=a.align_usable_sites(samples,{"L0":"saff-L0"},quorum=4)
    assert genes[0]["n_sites_same_reference_called_at_least80pct_taxa"]==5
    assert genes[0]["n_variable_sites_in_80pct_called"]==1
    assert genes[0]["n_parsimony_informative_sites_in_80pct_called"]==1
    assert alignment["WCVP4"][-1]=="N"
    assert alignment["WCVP0"]=="AAAAA"
    assert summary["strict_single_copy_orthology_gate_passed"] is False
    assert summary["orientation_evolution_rainfall_test_permitted"] is False


def test_missing_original_same_reference_contig_raises():
    sample={"tip":"WCVP1","loci":{"L0":{
        "reference":"lett-L0","reference_length_nt":150,
        "q20_majority_bases_1based":[[1,"A",3,1.0]]}}}
    with pytest.raises(ValueError,match="reference-copy identity"):
        a.align_usable_sites([sample],{"L0":"saff-L0"},quorum=1)


def test_unmeasured_genotype_cannot_be_filled_from_reference():
    sample={"tip":"WCVP1","loci":{"L0":{
        "reference":"saff-L0","reference_length_nt":150,
        "q20_majority_bases_1based":[[1,"A",2,.9]]}}}
    with pytest.raises(ValueError,match="Invalid original"):
        a.align_usable_sites([sample],{"L0":"saff-L0"},quorum=1)
