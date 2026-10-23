"""T3 source-authenticated 2021 Senecioneae target-capture 48-sample gene panel.

Recover *actual* PRJNA672812 read-run scientific names from ENA's public API
(NO downloading ~17GB of FASTQ), and join the original 10 documented head
state candidates plus independently source-checked Cremanthodium and Emilia.
Output a 48-real-DNA-tip native molecular collection with EXPLICIT unknown
capitulum posture, not a fabricated tree or 48 morphologically scored tips.

A focus result is that Ren2021 genuinely sequenced *Cremanthodium
campanulatum*, the exact SPECIES with a published rain/UV shelter posture
experiment, and C. lineare/thomsonii with direct flora-described nodding.
*Emilia sonchifolia*, ALSO a genuine Ren2021 SRA sample, has documented
opposite developmental movement (pendulous before anthesis; erect later).
No species-level U or independently rooted U->D evolutionary gains are
established from these few morphology-matched sample names. DNA specimen
voucher is NOT the experimentally manipulated flowering plant.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import io
import json
from pathlib import Path

import requests

NCBI_PRJNA="PRJNA672812"
ENA_URL="https://www.ebi.ac.uk/ena/portal/api/filereport"
FIELDS=("run_accession","scientific_name","sample_accession","experiment_accession",
        "fastq_ftp","fastq_bytes")
EXPECTED_SOURCE_RUNS=48
EXPECTED_CREMANTHODIUM=4
EXTERNAL_OLD_T3=10
FOCAL_SOURCE_STATES={
    "Cremanthodium campanulatum":"D",
    "Cremanthodium lineare":"D",
    "Cremanthodium thomsonii":"D",
    "Cremanthodium humile":"?",
    "Emilia sonchifolia":"V",
    "Parasenecio latipes":"D",
    "Parasenecio cyclotus":"D",
    "Parasenecio roborowskii":"D",
    "Parasenecio lidjiangensis":"D",
}
SOURCE_FOCAL_BY_SRR={
    "Cremanthodium campanulatum":"SRR12967014",
    "Cremanthodium lineare":"SRR12967012",
    "Cremanthodium thomsonii":"SRR12967057",
    "Cremanthodium humile":"SRR12967013",
    "Emilia sonchifolia":"SRR12967055",
    "Parasenecio latipes":"SRR12967038",
    "Parasenecio cyclotus":"SRR12967040",
    "Parasenecio roborowskii":"SRR12967033",
    "Parasenecio lidjiangensis":"SRR12967037",
}
OTHER_TRIBE_PROOF_ONLY_IF_CIRS_CIRSIUM_INDEPENDENT_TREE=True


def csv_rows(path):
    with Path(path).open(encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f))


def validate_ena(rows, original, supplemental):
    if len(rows)!=EXPECTED_SOURCE_RUNS:
        raise ValueError("Ren et al. 2021 original targeted capture BioProject not exactly 48 runs")
    if {c for r in rows for c in FIELDS if c not in r}:
        raise ValueError("Some published original DNA metadata fields missing")
    if len({r["run_accession"] for r in rows})!=EXPECTED_SOURCE_RUNS or len({
        r["sample_accession"] for r in rows})!=EXPECTED_SOURCE_RUNS:
        raise ValueError("Ren2021 duplicated sequencing runs/BioSamples; no pseudoreplicates")
    if any(not r["run_accession"].startswith("SRR") or
           not r["sample_accession"].startswith(("SAM","SRS")) or
           len(r["scientific_name"].split())<2 or
           not r["fastq_ftp"].startswith("ftp.sra.ebi.ac.uk/") or
           not r["fastq_bytes"] for r in rows):
        raise ValueError("NCBI original target capture SRA run source not authenticated")
    if len(original)!=EXTERNAL_OLD_T3 or len({r["taxon"] for r in original})!=EXTERNAL_OLD_T3:
        raise ValueError("Previous 10 independently sourced original Senecioneae candidates lost")
    if len(supplemental)!=len(FOCAL_SOURCE_STATES) or (
        {r["taxon"]:r["source_posture"] for r in supplemental}!=FOCAL_SOURCE_STATES):
        raise ValueError("Direct species flora and developmental-stage source codes changed")
    if any(not r["independent_flora_or_experiment_url"].startswith("https://")
           for r in supplemental):
        raise ValueError("Missing published original flora/experimental head-position evidence")
    published={r["taxon"]:r for r in original}
    for r in supplemental:
        old=published.get(r["taxon"])
        if old and (old["source_posture"]!=r["source_posture"]
                    or old["same_tip_rooted_nuclear_tree_verified"]!="no"):
            raise ValueError("New verified source must never overwrite previous original posture or pretend ready root")
    names={r["scientific_name"] for r in rows}
    if len({r["scientific_name"] for r in rows})!=EXPECTED_SOURCE_RUNS:
        raise ValueError("Repeated 2021 biological scientific name needs specimen-level resolution")
    if sum(n.startswith("Cremanthodium ") for n in names)!=EXPECTED_CREMANTHODIUM:
        raise ValueError("Expected four original Cremanthodium data-source sample names")
    if not set(FOCAL_SOURCE_STATES)<=names:
        raise ValueError("At least one independently described fixed SRA Senecioneae species not in true BioProject")
    if any(next(r for r in rows if r["scientific_name"]==taxon)["run_accession"]!=srr
           for taxon,srr in SOURCE_FOCAL_BY_SRR.items()):
        raise ValueError("Original real SRA-BioSample taxon anchor changed")
    old_matched=[z for z in original if z["taxon"] in names]
    if len(old_matched)!=1 or old_matched[0]["taxon"]!="Cremanthodium campanulatum":
        raise ValueError("Original ten flora examples do not have the expected single REAL nuclear specimen")
    outputs=[]
    for row in rows:
        name=row["scientific_name"]
        evo=next((x for x in supplemental if x["taxon"]==name),None)
        fixed=evo["source_posture"] if evo else "?"
        if fixed=="V" and evo["ontogenetic_direction"]!="D_before_to_U_later":
            raise ValueError("Developmental head movement cannot masquerade as fixed gain")
        if fixed in ("D","U") and (
                not evo or evo["source_evidence_class"]=="genus_only_not_species_orientation"):
            raise ValueError("Original tribe genus nomenclature cannot impute biological head trait")
        outputs.append({
            "original_nuclear_SRA_run":row["run_accession"],
            "original_nuclear_BioSample":row["sample_accession"],
            "original_nuclear_SRA_experiment":row["experiment_accession"],
            "original_species_name_exact_from_ENA":name,
            "T3_documented_head_state_D_U_V_or_unknown":fixed,
            "capitulum_direct_source_phrase":evo["original_head_morphology_phrase"] if evo else "",
            "capitulum_source_URL":evo["independent_flora_or_experiment_url"] if evo else "",
            "capitulum_documented_life_stage":evo["life_stage_scope"] if evo else "unknown",
            "ontogenetic_direction_not_phylogenetic_event":evo["ontogenetic_direction"] if evo else "none",
            "genuine_public_target_capture_bases_available":True,
            "DNA_sequence_reads_or_original_phylo_Newick_downloaded":False,
            "identical_DNA_voucher_and_flowering_phenotype_confirmed":False,
            "biological_head_state_fully_blind_double_scored":False,
            "nuclear_homeolog_reticulation_and_polyploid_orthology_audited":False,
        })
    counts=Counter(x["T3_documented_head_state_D_U_V_or_unknown"] for x in outputs)
    if counts!={"D":7,"V":1,"?":40}:
        raise ValueError("Source-authenticated 48-tip available morphology coverage changed")
    return outputs,counts


def fetch_public_ena(session=None):
    client=session or requests.Session()
    response=client.get(ENA_URL,params={
        "accession":NCBI_PRJNA,"result":"read_run",
        "fields":",".join(FIELDS),"format":"tsv"},
        headers={"Accept":"text/tab-separated-values",
                 "User-Agent":"azami-T3-original-public-phylogeny-source/1.0"},
        timeout=(15,60))
    response.raise_for_status()
    if len(response.content)>1_000_000:
        raise ValueError("Refusing any public source resembling full raw FASTQ download")
    data=response.content
    rows=list(csv.DictReader(io.StringIO(data.decode("utf-8-sig")),delimiter="\t"))
    if not rows:
        raise ValueError("Actual original Ren2021 ENA BioProject public run metadata unavailable")
    return rows,data

def emit(ena_data,original,supplemental,outdir):
    rows=list(csv.DictReader(io.StringIO(ena_data.decode("utf-8-sig")),delimiter="\t"))
    result,counts=validate_ena(rows,csv_rows(original),csv_rows(supplemental))
    target=Path(outdir);target.mkdir(parents=True,exist_ok=True)
    with (target/"T3_REN2021_PRJNA672812_ORIGINAL_NUCLEAR_SRA_48_TIP_HEAD_STATE_CROSSWALK.csv").open(
        "w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(result[0]),extrasaction="raise")
        w.writeheader();w.writerows(result)
    # Public scientific run accessions and published species names are okay;
    # byte-level FASTQ URLs and original sample coordinates need not be exposed.
    receipt={
        "schema":"azami_T3_Ren2021_PRJNA672812_actual_public48_nuclear_SRA_flora_crosswalk_v1",
        "status":"ORIGINAL_48_REAL_NUCLEAR_SRA_SAMPLES_PROVEN_WITH_SEVEN_NODDING_TWO_GENERA_AND_ONE_STAGE_VARIABLE",
        "source_paper_doi":"10.1016/j.ympev.2021.107202",
        "source_bioproject":"PRJNA672812",
        "source_original_48_SRA_runs":len(rows),
        "source_original_48_distinct_BioSamples":len({z["sample_accession"] for z in rows}),
        "source_original_48_distinct_species_names":len({z["scientific_name"] for z in rows}),
        "Cremanthodium_original_nuclear_SRA_species":EXPECTED_CREMANTHODIUM,
        "independent_original_morphology_10_ever_matching_2021_nuclear_specimens":1,
        "source_reconciled_direct_morphology_states":{"D":7,"V":1,"unknown":40,"U":0},
        "new_direct_flora_nodding_species":["Cremanthodium lineare","Cremanthodium thomsonii",
                                             "Parasenecio latipes","Parasenecio cyclotus",
                                             "Parasenecio roborowskii","Parasenecio lidjiangensis"],
        "experimental_other_tribe_posture_species_in_actual_nuclear_dataset":"Cremanthodium campanulatum",
        "documented_during_development_opposite_direction":"Emilia sonchifolia: pendulous before anthesis, erect later",
        "one_documented_genus_only_not_species_Cremanthodium_humile_retained_unknown":True,
        "no_empirically_scored_fixed_erect_sister_species_in_2021_SRA":True,
        "original_T3_10_species_exact_target_capture_overlap":1,
        "published_nuclear_tree_Newick_directly_retrieved":False,
        "DNA_FASTQ_reads_downloaded":False,
        "rooted_independent_U_to_D_evolutionary_gains_established":False,
        "rain_shelter_adaptive_fitness_claim_for_other_species":False,
        "Cremanthodium_campanulatum_NCBI_BioSample_not_same_experimental_individual":True,
        "2021_nuclear_runinfo_ena_SHA256":hashlib.sha256(ena_data).hexdigest(),
        "original_10_flora_treatments_SHA256":hashlib.sha256(Path(original).read_bytes()).hexdigest(),
        "additional_source_species_stage_treatments_SHA256":hashlib.sha256(Path(supplemental).read_bytes()).hexdigest(),
        "2021_NCBI_BioProject_URL":"https://www.ncbi.nlm.nih.gov/bioproject/PRJNA672812",
        "2021_archive_record_URL":"https://data.mendeley.com/datasets/bkm34wgc3n/1",
    }
    (target/"T3_REN2021_48_ORIGINAL_SRA_DIRECT_TIP_BOTANICAL_CROSSWALK.json").write_text(
        json.dumps(receipt,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print("T3_REAL_NUCLEAR_2021_SRA_FLORA_HEAD_CROSSWALK="+json.dumps(receipt,ensure_ascii=False,sort_keys=True),flush=True)
    return receipt


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--original-10",required=True,type=Path)
    p.add_argument("--added-5",required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    a=p.parse_args()
    rows,data=fetch_public_ena()
    emit(data,a.original_10,a.added_5,a.outdir)

if __name__=="__main__":main()
