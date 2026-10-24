"""Audit REAL bounded SRA reads from the two S2 focal Cirsium taxa.

This is a sequencing read-access/mapping pilot, not a 241-locus HybPiper
assembly, accepted 2-tip phylogeny, or a genus-wide species tree. Only a
bounded first-30000 SRA spots per run are exported, and source reads are NOT
included in GitHub Actions artifacts.

Source run and BioSample accessions are pinned to the 2026-10-22 full
EAzami294/Moreyra crosswalk with an explicitly selected non-conflicted
C. vulgare voucher. Original Compositae1061 reference is independently
authenticated in EAzami.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

PANEL = (
    {"taxon":"Cirsium arvense","run":"SRR25265555",
     "biosample":"SAMN34240259","source":"Moreyra2025"},
    {"taxon":"Cirsium vulgare","run":"SRR25265662",
     "biosample":"SAMN34240353","source":"Moreyra2025"},
)
MAX_SPOTS=30000
EXPECTED_REFERENCE_LOCI=1061
ORIGINAL_REFERENCE_BLOB="4f89e234007f367ffa8aa5e2be536bc44f31f445"


def read_idxstats(path):
    """One reference-sequence alignment summary; loci may have >1 reference."""
    raw=Path(path).read_text(encoding="utf-8").splitlines()
    if not raw:
        raise ValueError("BWA/samtools idxstats missing")
    contigs=[]
    observed=set()
    loci=set()
    total_mapped=0
    for line in raw:
        parts=line.rstrip().split("\t")
        if len(parts)!=4:
            raise ValueError("Invalid four-column samtools idxstats line")
        name,length,mapped,unmapped=parts
        if name=="*":
            continue
        if name in observed:
            raise ValueError("Duplicated reference contig in index")
        observed.add(name)
        if not name.startswith(("lett-","saff-","sunf-")) or "-" not in name:
            raise ValueError("Unknown original reference contig header")
        gene=name.split("-",1)[1]
        if not gene or int(length)<50 or min(int(mapped),int(unmapped))<0:
            raise ValueError("Invalid reference locus or samtools counts")
        contigs.append(name)
        if int(mapped)>0:
            loci.add(gene)
            total_mapped+=int(mapped)
    if len(contigs)<EXPECTED_REFERENCE_LOCI or len(set(x.split("-",1)[1] for x in contigs))!=EXPECTED_REFERENCE_LOCI:
        raise ValueError("Real Comp1061 reference locus count different from frozen source")
    return {
        "reference_sequence_records":len(contigs),
        "source_unique_loci":EXPECTED_REFERENCE_LOCI,
        "mapped_reference_loci_at_least_one_alignment":len(loci),
        "mapped_alignments_from_idxstats":total_mapped,
    }


def audit_metrics(folder):
    results=[]
    for target in PANEL:
        path=Path(folder)/target["run"]
        files=["extraction_metrics.json","mapq20_primary_count.txt","idxstats.tsv"]
        for n in files:
            if not (path/n).is_file():
                raise ValueError(f"Missing actual public read pilot evidence: {target['run']} / {n}")
        metric=json.loads((path/"extraction_metrics.json").read_text(encoding="utf-8"))
        for k in ("run","biosample","taxon"):
            if metric.get(k)!=target[k]:
                raise ValueError("Wrong accession / original source BioSample in real read result")
        if metric.get("source")!=target["source"] or int(metric.get("max_spots",0))!=MAX_SPOTS:
            raise ValueError("Not original source or unexpected SRA spot limit")
        paired=int(metric.get("paired_fastq_reads",0))
        reads=int(metric.get("forward_fastq_reads",0))
        if paired<1 or reads<1 or paired!=reads or reads>MAX_SPOTS:
            raise ValueError("Pilot does not contain bounded paired-end public reads")
        hq=int((path/"mapq20_primary_count.txt").read_text().strip())
        if not 0<=hq<=2*reads:
            raise ValueError("Impossible uniquely mapped count")
        summary=read_idxstats(path/"idxstats.tsv")
        if hq>summary["mapped_alignments_from_idxstats"]:
            raise ValueError("High MAPQ primary reads exceed total index mappings")
        results.append({
            "taxon":target["taxon"],
            "run":target["run"],
            "biosample":target["biosample"],
            "source":target["source"],
            "max_spots_attempted":MAX_SPOTS,
            "read_pairs_extracted":reads,
            "high_mapq20_primary_alignments":hq,
            **summary,
            "locus_count_is_alignment_overlap_not_full_HybPiper_recovery":True,
            "partial_first_spots_not_random_sample":True,
            "nuclear_species_phylogeny_estimated":False,
            "original_350_locus_phylogeny_reproduced":False,
        })
    if {r["run"] for r in results}!={x["run"] for x in PANEL}:
        raise ValueError("Missing S2 source accession")
    return {
        "status":"ACTUAL_BOUNDED_PUBLIC_TARGET_CAPTURE_READ_MAPPING_PILOT",
        "run_count":len(results),
        "reference_source":"Siniscalchi original Compositae1061, not augmented Moreyra reference",
        "original_reference_git_blob_sha1":ORIGINAL_REFERENCE_BLOB,
        "source_public_repo":"zuizui0223/EAzami original locus reference recovery",
        "n_raw_fastq_or_bam_files_uploaded":0,
        "full_125_S1_run_phylogeny_estimated":False,
        "full_294_sample_EAzami_tree_estimated":False,
        "full_244_source_taxa_retained_in_existing_ledger":True,
        "S2_public_target_capture_run_results":results,
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-workdir",type=Path,required=True)
    p.add_argument("--out-dir",type=Path,required=True)
    a=p.parse_args()
    result=audit_metrics(a.source_workdir)
    a.out_dir.mkdir(parents=True,exist_ok=True)
    (a.out_dir/"S2_REAL_PUBLIC_READ_MAPPING_30000_SPOTS_RECEIPT.json").write_text(
        json.dumps(result,sort_keys=True,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("T2_S2_REAL_NUCLEAR_CAPTURE_READ_PILOT="+json.dumps(result,sort_keys=True,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
