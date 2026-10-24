"""Observed-only 111-source Q20 positions and strict-provenance technical ML gates."""
import json
from pathlib import Path
import pytest
from analysis.v4 import t2_real111_fixed118_nucleotide_alignment as p


def synthetic(monkeypatch):
    monkeypatch.setattr(p,"FULL",6)
    monkeypatch.setattr(p,"LOCI",2)
    monkeypatch.setattr(p,"MIN_SITES_PER_LOCUS",2)
    monkeypatch.setattr(p,"MIN_TAXA_FOR_TECHNICAL_ML",3)
    monkeypatch.setattr(p,"MIN_LOCI_FOR_TECHNICAL_ML",2)
    monkeypatch.setattr(p,"MIN_BASES_FOR_TECHNICAL_ML",5)
    monkeypatch.setattr(p,"MIN_VARIABLE_SITES_FOR_TECHNICAL_ML",1)
    monkeypatch.setattr(p,"NUCLEOTIDE_QUORUM",.8)
    monkeypatch.setattr(p,"CONFLICTING_RUN",{"SRR105"})
    samples=[]
    runs={}
    targets={}
    for l in ["L0","L1"]:
        targets[l]="saff-"+l
    for i in range(6):
        run=f"SRR10{i}"
        name=f"Cirsium species{i}"
        samples.append({
           "run":run,"biosample":f"SAMN10{i}",
           "accepted_key":str(20670000+i),
           "phylogenetic_tip":"WCVP"+str(20670000+i),
           "taxon_name":name})
        per_locus={}
        for locus,ref in targets.items():
            vals=[[1,"A",4,.99],
                  [2,("G" if i>=2 else "C"),5,.95],
                  [3,"T",4,1.0]]
            if i==1:
                vals=vals[:-1]
            per_locus[locus]={
                "reference":ref,"reference_length_nt":100,
                "not_a_validated_allelic_genotype":True,
                "q20_majority_bases_1based":vals}
        runs[run]={
           **samples[-1],"status":"real_source111_120k_Q20_identical_copy_base_calls",
           "source_public_ENA_paired_read_count":120000,
           "n_predeclared_80percent_reference_loci":2,
           "source_run_label_conflict_unresolved":i==5,
           "uncalled_bases_not_filled_by_reference":True,
           "original_single_copy_orthology_verified":False,
           "observed_original_reference_q20_site_calls":per_locus}
    target_rows=[{"locus":key,"original_compositae_reference":val}
                 for key,val in targets.items()]
    return samples,target_rows,runs,targets


def test_sites_intersect_only_called_bases_not_reference(monkeypatch):
    samples,targets,runs,refs=synthetic(monkeypatch)
    qc,loci,full,technical,state=p.site_and_taxa_alignment(runs,refs)
    assert state["original_real120k_SRA_source_runs"]==6
    assert state["original_unambiguous_taxon_sample_count"]==5
    assert state["source_number_for_individual_site"]==4
    assert state["primary_technical_alignment_real_nucleotide_sites"]==6
    assert len(full)==5
    assert "WCVP20670005" not in full
    assert full["WCVP20670001"].count("N")==2
    assert len(loci)==2 and all(x["n_sites_retained_in_phylo_alignment"]==3 for x in loci)
    assert not state["biological_accepted_taxon_concepts_or_orthologs_validated"]
    assert state["technical_source_taxa_with25pct_sites"]==5
    assert state["ML_readiness"]


def test_too_few_observed_sites_is_not_filled_by_reference(monkeypatch):
    samples,targets,runs,refs=synthetic(monkeypatch)
    for r in runs.values():
        for loc in refs:
            r["observed_original_reference_q20_site_calls"][loc]["q20_majority_bases_1based"]=[]
    qc,loci,full,technical,state=p.site_and_taxa_alignment(runs,refs)
    assert state["primary_technical_alignment_real_nucleotide_sites"]==0
    assert not state["ML_readiness"]
    assert all(x["n_sites_retained_in_phylo_alignment"]==0 for x in loci)


def test_all_receipts_source_run_identity_and_reference_exact(tmp_path,monkeypatch):
    samples,targets,runs,refs=synthetic(monkeypatch)
    for i in range(10):
        path=tmp_path/f"shard_{i}.json"
        group=[runs[f"SRR10{i}"]] if i<6 else []
        path.write_text(json.dumps({
            "shard":i,"n_shards":10,"no_original_FASTQ_or_BAM_preserved":True,
            "requested":len(group),"passed":len(group),"failed":0,"outcomes":group}))
    received,copy=p.read_receipts(tmp_path,samples,targets)
    assert len(received)==6 and copy==refs
    runs["SRR100"]["observed_original_reference_q20_site_calls"]["L0"]["reference"]="saff-OTHER"
    (tmp_path/"shard_0.json").write_text(json.dumps({
        "shard":0,"n_shards":10,"no_original_FASTQ_or_BAM_preserved":True,
        "requested":1,"passed":1,"failed":0,"outcomes":[runs["SRR100"]]}))
    with pytest.raises(ValueError,match="same-copy"):
        p.read_receipts(tmp_path,samples,targets)


def test_unmeasured_111_specimen_not_dropped_or_duplicated(monkeypatch):
    samples,targets,runs,refs=synthetic(monkeypatch)
    samples[0]["run"]="SRRFAKE"
    with pytest.raises(ValueError,match="sample inventory"):
        p.source_targets(samples,targets)


def test_predeclared_real_global_thresholds():
    assert p.FULL==111 and p.LOCI==118
    assert p.NUCLEOTIDE_QUORUM==.8
    assert p.MIN_SITES_PER_LOCUS==100
    assert p.MIN_TAXA_FOR_TECHNICAL_ML==60
    assert p.CONFLICTING_RUN=={"SRR30887230"}
