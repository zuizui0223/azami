"""Freeze whole-source 111 taxon / 241 locus modal-reference tree input candidates.

Empirical input is the COMPLETED all-123-run, 111-selected-taxon Q20
per-locus depth study. Keep ALL 113 Moreyra source-key candidates and master
S1 251 operational source names. Select a reference for each candidate locus
using depth-3x breadth >=25% in at least 80% of 111 distinct source taxa.
The 80% cutoff is fixed BEFORE inspecting full base calls.

This is a *data-readiness* gate only; modal identical original reference
copies cannot certify single-copy orthology and cannot yet justify phylogenetic
species-tree, head-position transition, BM/OU adaptation or chronology claims.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

ORIGINAL_SOURCE_NAMES=251
MASTER_ACCEPTED_KEYS=244
AUTHOR_LOCI=241
MOREYRA_CANDIDATES=113
MEASURED=111
MISSING_SHORT=2
EXACT_REF_QUORUM=.80
MIN_EXPECTED_PHYLO_CANDIDATE_LOCI=30
SOURCE_HEADERS={"accepted_key","original_source_taxon","selected_single_real_run",
                "selected_single_biosample","real120k_depth_measured"}
LOCUS_HEADERS={"locus","n_taxa_Q20_dp3_ge25pct_same_original_reference_copy",
               "modal_Q20_dp3_ge25pct_original_ref",
               "n_taxa_Q20_dp3_ge25pct","n_taxa_Q20_dp3_ge50pct"}
NUCLEAR_POSTURE={"Cirsium conspicuum":"nodding","Cirsium scopulorum":"nodding"}


def read_csv(path):
    with Path(path).open(newline="",encoding="utf-8-sig") as f:
        rows=list(csv.DictReader(f))
    if not rows:
        raise ValueError("Empty source verified original nuclear inventory")
    return rows


def original_tip(taxon_key):
    val=str(taxon_key)
    if not val.isdigit():
        raise ValueError("Taxon accepted key not numeric; fail closed")
    return "WCVP"+val


def preflight(samples,loci):
    if len(samples)!=MOREYRA_CANDIDATES or not SOURCE_HEADERS<=set(samples[0]):
        raise ValueError("S1 original 113 source names / sample columns changed")
    if len(loci)!=AUTHOR_LOCI or not LOCUS_HEADERS<=set(loci[0]):
        raise ValueError("Original author241 locus matrix unavailable")
    if len({r["accepted_key"] for r in samples})!=len(samples):
        raise ValueError("One S1 accepted concept occurs twice in source")
    if len({r["locus"] for r in loci})!=len(loci):
        raise ValueError("Source 241 candidate loci duplicated")
    observed=[]
    missing=[]
    for row in samples:
        measured=str(row["real120k_depth_measured"]).lower()=="true"
        if str(row["real120k_depth_measured"]).lower() not in {"true","false"}:
            raise ValueError("Source public 120k read evidence is ambiguous")
        name=" ".join(row["original_source_taxon"].replace("_"," ").split())
        run=row["selected_single_real_run"]
        bio=row["selected_single_biosample"]
        if measured:
            if not run.startswith("SRR") or not bio.startswith("SAMN"):
                raise ValueError("One 120k source SRA sample lacks real BioSample identity")
            observed.append({
                "accepted_key":row["accepted_key"],"taxon_name":name,
                "phylogenetic_tip":original_tip(row["accepted_key"]),
                "run":run,"biosample":bio,
                "independent_original_public_sample":True,
                "source_name_taxon_concept_verified":False,
            })
        else:
            if run or bio:
                raise ValueError("Unavailable 120k sample must not have fabricated run selection")
            missing.append({"accepted_key":row["accepted_key"],"source_taxon":name,
                            "not_measured_120k_not_missing_gene":True})
    if (len(observed)!=MEASURED or len(missing)!=MISSING_SHORT or
            len({x["run"] for x in observed})!=MEASURED or
            len({x["biosample"] for x in observed})!=MEASURED):
        raise ValueError("Source 111 real BioSamples plus two documented short archives drifted")
    n=math.ceil(EXACT_REF_QUORUM*MEASURED)
    for z in loci:
        for f in ("n_taxa_Q20_dp3_ge25pct_same_original_reference_copy",
                  "n_taxa_Q20_dp3_ge25pct","n_taxa_Q20_dp3_ge50pct"):
            value=int(z[f])
            if not 0<=value<=MEASURED:
                raise ValueError("Invalid nuclear Q20 taxon occupancy")
        if int(z["n_taxa_Q20_dp3_ge25pct_same_original_reference_copy"])>int(z["n_taxa_Q20_dp3_ge25pct"]):
            raise ValueError("Modal original nuclear reference copy count exceeds total taxa")
    chosen=[]
    for row in loci:
        same=int(row["n_taxa_Q20_dp3_ge25pct_same_original_reference_copy"])
        name=row["modal_Q20_dp3_ge25pct_original_ref"]
        if same<n:
            continue
        if not name.startswith(("lett-","saff-","sunf-")) or not name.endswith("-"+row["locus"]):
            raise ValueError("No unambiguous original single reference copy for genomic candidate")
        chosen.append({
           "locus":row["locus"],
           "original_compositae_reference":name,
           "empirically_supported_source_taxa":same,
           "fraction_of_source_taxa_with_3x25pct_same_reference":same/MEASURED,
           "minimum_source_taxa_required":n,
           "single_original_reference_copy_is_not_validated_ortholog":True,
        })
    chosen.sort(key=lambda row:row["locus"])
    named={x["taxon_name"] for x in observed}
    if set(NUCLEAR_POSTURE)-named:
        raise ValueError("Both source independently botanical coded nodding anchors absent")
    status=({
       "schema":"azami_T2_original111_common241_reference_preflight_20261023_v1",
       "original_S1_accepted_taxa_master":MASTER_ACCEPTED_KEYS,
       "original_S1_named_taxa_master":ORIGINAL_SOURCE_NAMES,
       "source_Moreyra_name_candidates":MOREYRA_CANDIDATES,
       "source_120k_BioSamples":MEASURED,
       "source_public_reads_short_no120k":MISSING_SHORT,
       "original_author_no_warning_locus_universe":AUTHOR_LOCI,
       "source_site_depth": "samtools baseQ20 mapQ20, depth>=3, breadth>=25%, same original reference copy",
       "source_taxon_occupancy_gate":EXACT_REF_QUORUM,
       "source_taxon_number_for_gate":n,
       "candidate_80pct_original_reference_loci":len(chosen),
       "minimum_candidates_before_technical_ML":MIN_EXPECTED_PHYLO_CANDIDATE_LOCI,
       "technical_ML_source_locus_gate_passed":len(chosen)>=MIN_EXPECTED_PHYLO_CANDIDATE_LOCI,
       "source_botanical_nodding_author_names_present":sorted(NUCLEAR_POSTURE),
       "no_fastq_or_BAM":True,"no_nuclear_species_tree":True,
       "same_original_reference_not_single_copy_orthology":True,
       "no_inference_of_gravity_nodding_from_image_angle":True,
       "not_an_accepted_global_244_or_294_species_tree":True,
    })
    return observed,missing,chosen,status


def write_csv(path,rows,fields=None):
    with Path(path).open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields or list(rows[0]),extrasaction="raise")
        w.writeheader();w.writerows(rows)


def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("--taxa",type=Path,required=True)
    a.add_argument("--loci",type=Path,required=True)
    a.add_argument("--outdir",type=Path,required=True)
    args=a.parse_args()
    taxa=read_csv(args.taxa)
    loci=read_csv(args.loci)
    samples,missing,genes,summary=preflight(taxa,loci)
    summary["source_sha256"]={
       "taxa":hashlib.sha256(args.taxa.read_bytes()).hexdigest(),
       "loci":hashlib.sha256(args.loci.read_bytes()).hexdigest()}
    args.outdir.mkdir(parents=True,exist_ok=True)
    write_csv(args.outdir/"S1_111_REAL120K_SAMPLE_SOURCE_RUNS.csv",samples)
    write_csv(args.outdir/"S1_TWO_SHORT_SOURCE_TAXA_NOT_GENETIC_ABSENCE.csv",missing)
    if genes:
        write_csv(args.outdir/"AUTHOR241_S1_111_SOURCE_COMMON_REFERENCE_ML_TARGETS.csv",genes)
    (args.outdir/"S1_111_REAL_ORIGINAL241_GLOBAL_PHYLO_PREFLIGHT.json").write_text(
        json.dumps(summary,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("T2_REAL111_SOURCE_PHYLO_PREFLIGHT="+json.dumps(summary,ensure_ascii=False,sort_keys=True),flush=True)


if __name__=="__main__":
    main()
