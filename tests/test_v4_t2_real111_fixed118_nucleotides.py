"""Real source sample/read-to-118 fixed original reference call contract checks."""
import json
import pytest

from analysis.v4 import t2_real111_fixed118_nucleotides as m


def fixture(monkeypatch):
    monkeypatch.setattr(m,"SOURCE_RUNS",4)
    monkeypatch.setattr(m,"SAMPLE_COUNT",3)
    monkeypatch.setattr(m,"TARGET_COUNT",2)
    allruns=[]
    samples=[]
    for i in range(4):
        run="SRR"+str(200+i)
        orig={"run":run,"biosample":"SAMN"+str(i),"accepted_key":str(205500+i),
              "S1_source_taxon":"Cirsium_test"+str(i),
              "original_tree_code":"Cirsium test"+str(i)}
        allruns.append(orig)
        if i<3:
            samples.append({
                "run":run,"biosample":orig["biosample"],
                "accepted_key":orig["accepted_key"],
                "taxon_name":"Cirsium test"+str(i),
                "phylogenetic_tip":"WCVP"+orig["accepted_key"]})
    refs=[{
       "locus":name,"original_compositae_reference":"saff-"+name,
       "empirically_supported_source_taxa":"89"
    } for name in ("L1","L2")]
    return allruns,samples,refs


def test_exact_source_specimens_and_predeclared_reference_copies(monkeypatch):
    originals,samples,targets=fixture(monkeypatch)
    rows,selected,refs=m.audit_samples(originals,samples,targets)
    assert len(rows)==4 and len(selected)==3 and len(refs)==2
    assert set(refs)=={"L1","L2"}
    assert set(selected)=={"SRR200","SRR201","SRR202"}
    samples[0]["biosample"]="SAMNFAKE"
    with pytest.raises(ValueError,match="BioSample"):
        m.audit_samples(originals,samples,targets)


def test_no_cross_locus_or_duplicate_reference(monkeypatch):
    orig,samples,targets=fixture(monkeypatch)
    targets[0]["original_compositae_reference"]="lett-Another"
    with pytest.raises(ValueError,match="80%-copy"):
        m.audit_samples(orig,samples,targets)


def test_validate_genuine_bases_and_scientific_missingness(monkeypatch):
    orig,samples,targets=fixture(monkeypatch)
    _,selected,refs=m.audit_samples(orig,samples,targets)
    src=selected["SRR200"]
    item={
        "run":src["run"],"biosample":src["biosample"],
        "accepted_key":src["accepted_key"],
        "taxon_name":src["taxon_name"],
        "phylogenetic_tip":src["phylogenetic_tip"],
        "source_public_ENA_paired_read_count":120000,
        "n_predeclared_80percent_reference_loci":2,
        "observed_original_reference_q20_site_calls":{
            l:{"reference":ref,"reference_length_nt":150,
               "not_a_validated_allelic_genotype":True,
               "q20_majority_bases_1based":[[1,"A",3,1.0],[3,"C",4,.95]]}
            for l,ref in refs.items()},
        "n_real_Q20_majority_callable_positions":4,
        "uncalled_bases_not_filled_by_reference":True,
        "original_single_copy_orthology_verified":False,
        "FASTQ_BAM_exported":False,
    }
    assert m.validate_one(item,src,refs)
    item["observed_original_reference_q20_site_calls"]["L1"][
        "q20_majority_bases_1based"][0][2]=2
    with pytest.raises(ValueError,match="Invalid positionwise"):
        m.validate_one(item,src,refs)


def test_all_three_selected_in_batch_even_when_poor_coverage(monkeypatch):
    originals,samples,targets=fixture(monkeypatch)
    monkeypatch.setattr(m,"SHARDS",3)
    def f(_source,selected,refs,reference,scratch):
        return {
            "run":selected["run"],"biosample":selected["biosample"],
            "accepted_key":selected["accepted_key"],
            "taxon_name":selected["taxon_name"],
            "phylogenetic_tip":selected["phylogenetic_tip"],
            "status":"real_source111_120k_Q20_identical_copy_base_calls",
            "source_public_ENA_paired_read_count":120000,
            "n_predeclared_80percent_reference_loci":2,
            "observed_original_reference_q20_site_calls":{
                l:{"reference":ref,"reference_length_nt":100,
                   "not_a_validated_allelic_genotype":True,
                   "q20_majority_bases_1based":[]} for l,ref in refs.items()},
            "n_real_Q20_majority_callable_positions":0,
            "uncalled_bases_not_filled_by_reference":True,
            "original_single_copy_orthology_verified":False,
            "FASTQ_BAM_exported":False}
    for i in range(3):
        d=m.shard(originals,samples,targets,None,None,i,3,runner=f)
        assert d["requested"]==1 and d["passed"]==1 and d["failed"]==0
        assert d["not_phylogenetic_species_tree"]


def test_real_global_scope_constants_stable():
    assert m.SOURCE_RUNS==125 and m.SAMPLE_COUNT==111
    assert m.TARGET_COUNT==118 and m.READ_PAIRS==120000
    assert m.AMBIGUOUS_ORIGINAL_SOURCE_RUNS=={"SRR30887230"}
