"""Unit checks for genuine but tiny public S2 capture-read pilot evidence."""
import json
from pathlib import Path

import pytest

from analysis.v4 import t2_s2_bounded_public_reads as pilot


def write_sample(folder,source):
    root=folder/source["run"]
    root.mkdir(parents=True)
    receipt={**source,"max_spots":30000,
             "forward_fastq_reads":25000,"paired_fastq_reads":25000}
    (root/"extraction_metrics.json").write_text(json.dumps(receipt),encoding="utf-8")
    (root/"mapq20_primary_count.txt").write_text("198\n",encoding="utf-8")
    (root/"idxstats.tsv").write_text(
        "saff-L001\t500\t150\t10\n"
        "lett-L002\t450\t300\t11\n"
        "sunf-L003\t390\t11\t0\n"
        "*\t0\t0\t45000\n",encoding="utf-8")
    return root


def test_two_real_accessions_and_reference_locus_overlap(tmp_path,monkeypatch):
    monkeypatch.setattr(pilot,"EXPECTED_REFERENCE_LOCI",3)
    for source in pilot.PANEL:
        write_sample(tmp_path,source)
    result=pilot.audit_metrics(tmp_path)
    assert result["run_count"]==2
    assert {r["run"] for r in result["S2_public_target_capture_run_results"]}=={
        "SRR25265555","SRR25265662"}
    assert all(r["mapped_reference_loci_at_least_one_alignment"]==3
               for r in result["S2_public_target_capture_run_results"])
    assert all(not r["nuclear_species_phylogeny_estimated"]
               for r in result["S2_public_target_capture_run_results"])
    assert result["n_raw_fastq_or_bam_files_uploaded"]==0


def test_reference_locus_count_drift_fails(tmp_path,monkeypatch):
    monkeypatch.setattr(pilot,"EXPECTED_REFERENCE_LOCI",4)
    for x in pilot.PANEL:
        write_sample(tmp_path,x)
    with pytest.raises(ValueError,match="locus count"):
        pilot.audit_metrics(tmp_path)


def test_locus_ids_cannot_be_inferred_from_unknown_reference(tmp_path,monkeypatch):
    monkeypatch.setattr(pilot,"EXPECTED_REFERENCE_LOCI",3)
    for x in pilot.PANEL:
        folder=write_sample(tmp_path,x)
        (folder/"idxstats.tsv").write_text(
            "badref-L001\t500\t150\t10\n"
            "lett-L002\t450\t300\t11\n"
            "sunf-L003\t390\t11\t0\n",encoding="utf-8")
    with pytest.raises(ValueError,match="Unknown original"):
        pilot.audit_metrics(tmp_path)


def test_invalid_source_biosample_blocks_promotion(tmp_path,monkeypatch):
    monkeypatch.setattr(pilot,"EXPECTED_REFERENCE_LOCI",3)
    for x in pilot.PANEL:
        folder=write_sample(tmp_path,x)
        data=json.loads((folder/"extraction_metrics.json").read_text())
        data["biosample"]="SAMN_STALE"
        (folder/"extraction_metrics.json").write_text(json.dumps(data))
    with pytest.raises(ValueError,match="Wrong accession"):
        pilot.audit_metrics(tmp_path)


def test_original_reference_not_a_reproduction_of_moreyra_350():
    assert pilot.MAX_SPOTS==30000
    assert pilot.EXPECTED_REFERENCE_LOCI==1061
    assert pilot.ORIGINAL_REFERENCE_BLOB=="4f89e234007f367ffa8aa5e2be536bc44f31f445"
