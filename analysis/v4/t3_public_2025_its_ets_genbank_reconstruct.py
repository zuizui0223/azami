"""Reconstruct a genuine *exploratory* Senecioneae nuclear ITS+ETS tree.

Accession-source: Liu & Jiang (2025) Taiwania 70:403-410 supplemental
Table S1, 163 public GenBank accession ROWS / 154 distinct taxon names
(9 species repeated across source records). Each accession is retrieved
from NCBI EUtilities by exact ID and its metadata matched to the original
published source biological binomial. No synonym coercion or species name
imputation. Missing/inconsistent NCBI record is missing, NOT a reference base.

ITS and ETS are LINKED tandem ribosomal DNA spacers, NOT independent
single-copy genes: this is not the validated Ren2021 target-capture nuclear
species tree and cannot identify the number of U->D transitions, nor rain
adaptive gains. One independently scored flowering plant is not the same
as any of the sequencing vouchers in this supplement.

This file performs the real small GenBank download / original-source
accession auditing and FASTA preparation, not a mocked tree.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from collections import Counter
from pathlib import Path
import re
import time

from Bio import SeqIO
import requests

SOURCE_URL="https://taiwania.ntu.edu.tw/pdf/tai.2025.70.403.a.pdf"
SOURCE_PAPER="Liu & Jiang 2025 Taiwania 70(3) 403-410"
SOURCE_MARKERS=("ITS","ETS")
ORIGINAL_ROWS=163
ORIGINAL_NAMES=154
EXPECTED_PUBLISHED_ETS_ROWS=90
MAX_FASTA_RESPONSE=1_500_000
BATCH=14
THROTTLE=0.55
MIN_DNA_LEN=200
MAX_DNA_LEN=2000
NCBI_EFETCH="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
# Never disclose/remove missing source names simply because no sequence fetched.


def read_csv(p):
    with Path(p).open(newline="",encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


def validate_original(rows):
    if len(rows)!=ORIGINAL_ROWS or len({r["original_supplement_row"] for r in rows})!=ORIGINAL_ROWS:
        raise ValueError("Official 2025 Taiwania Table S1 163 accession source rows missing")
    if len({r["original_published_scientific_name"] for r in rows})!=ORIGINAL_NAMES:
        raise ValueError("Published original 154 species names not conserved")
    if Counter(bool(r["ETS_GenBank_accession"]) for r in rows)[True]!=EXPECTED_PUBLISHED_ETS_ROWS:
        raise ValueError("Original 90 ETS entries not conserved")
    if set(int(r["original_supplement_row"]) for r in rows)!=set(range(1,164)):
        raise ValueError("Original accession source row IDs must cover every 1..163")
    for r in rows:
        name=r["original_published_scientific_name"]
        if (not name.split()[0][0].isupper()
            or len(name.split())!=2
            or len(name)>100
            or r["original_published_source_url"]!=SOURCE_URL
            or r["ITS_and_ETS_independent_unlinked_gene_loci"]!="no"
            or r["species_and_seq_origin_same_individual_verified"]!="no"):
            raise ValueError("GenBank specimen or unlinked rRNA locus provenance distorted")
        for marker in SOURCE_MARKERS:
            acc=r[f"{marker}_GenBank_accession"]
            if acc and not re.fullmatch(r"[A-Z]{1,3}\d{5,9}",acc):
                raise ValueError("Invalid GenBank source accession")
    return sorted({r[f"{x}_GenBank_accession"] for r in rows for x in SOURCE_MARKERS if r[f"{x}_GenBank_accession"]})


class NCBI:
    def __init__(self,session=None,delay=THROTTLE):
        self.client=session or requests.Session()
        self.delay=delay
        self.last=0.
        self.n_calls=0
        self.n_remote_bytes=0

    def fetch_batch(self,ids):
        if not 1<=len(ids)<=BATCH or any(not re.fullmatch(r"[A-Z]{1,3}\d{5,9}",x) for x in ids):
            raise ValueError("Unbounded/unpinned request for original GenBank source")
        params={"db":"nuccore","id":",".join(ids),"rettype":"gb","retmode":"text",
                "tool":"azami_T3_ITS_ETS_published2025"}
        for attempt in range(5):
            wait=self.last+self.delay-time.monotonic()
            if wait>0:time.sleep(wait)
            self.last=time.monotonic()
            self.n_calls+=1
            try:
                r=self.client.get(NCBI_EFETCH,params=params,
                    headers={"User-Agent":"azami-T3-authenticated-rDNA/1"},
                    timeout=(12,55))
                if r.status_code in (429,500,502,503,504):
                    time.sleep(min(2**attempt,10))
                    continue
                r.raise_for_status()
                if len(r.content)>MAX_FASTA_RESPONSE:
                    raise ValueError("Unexpected bulk GenBank download; reject")
                self.n_remote_bytes+=len(r.content)
                records=list(SeqIO.parse(io.StringIO(r.text),"genbank"))
                if not records:
                    raise ValueError("Original published NCBI ID returned no parseable GenBank record")
                if len(records)>len(ids):
                    raise ValueError("NCBI accession duplication in requested batch")
                return records
            except (requests.RequestException, ValueError):
                if attempt==4:raise
                time.sleep(min(2**attempt,10))
        raise RuntimeError("Unreachable NCBI state")


def check_record(rec,published_taxon,marker,accession):
    original=str(rec.id).split(".")[0]
    organism=" ".join(rec.annotations.get("organism","").split())
    if original!=accession:
        return None,"NCBI_accession_record_ID_mismatch"
    if organism.lower()!=published_taxon.lower():
        return None,"organism_name_differs_from_published_source_name"
    title=rec.description.lower()
    if any(x in title for x in ("chloroplast","plastid","mitochondr","complete genome")):
        return None,"nonnuclear_title"
    seq=str(rec.seq).upper().replace("U","T")
    if len(seq)<MIN_DNA_LEN or len(seq)>MAX_DNA_LEN:
        return None,"unusable_rDNA_sequence_length"
    if any(c not in "ACGTNRYWSKMBDHV" for c in seq) or seq.count("N")>len(seq)*.06:
        return None,"invalid_or_highly_ambiguous_genbank_sequence"
    # Published marker provenance is primary; title may simply say rDNA
    # region without the words 'internal'/'external'.
    info={"accession_version":rec.id,"marker":marker,
          "published_species_name":published_taxon,"NCBI_species_organism":organism,
          "genbank_title":rec.description,
          "n_original_bases":len(seq),"raw_nucleotide_sequence":seq,
          "title_explicit_marker":bool(re.search(
               r"ITS|internal transcribed" if marker=="ITS"
               else r"ETS|external transcribed",title,re.I)),
          "source_sample_same_in_ITS_and_ETS_not_verified":True}
    return info,None


def collect_original(rows,fetcher,*,taxa_limit=None):
    originals=sorted(rows,key=lambda r:int(r["original_supplement_row"]))
    if taxa_limit is not None:
        if not 1<=taxa_limit<=ORIGINAL_ROWS:raise ValueError("Bounded pilot scope")
        originals=originals[:taxa_limit]
    requested=sorted({r[f"{m}_GenBank_accession"] for r in originals
                      for m in SOURCE_MARKERS if r[f"{m}_GenBank_accession"]})
    records={}
    err=[]
    for i in range(0,len(requested),BATCH):
        subset=requested[i:i+BATCH]
        try:
            for g in fetcher.fetch_batch(subset):
                key=g.id.split(".")[0]
                if key not in subset:raise ValueError("NCBI response contained unexpected public sequence ID")
                if key in records:raise ValueError("NCBI repeated the same accession in batch")
                records[key]=g
        except Exception as exc:
            err.append({"batch_start":i,"batch_size":len(subset),
                        "reason":str(exc)[:180],
                        "source_access_unavailable_not_absence":True})
    scored=[]
    for r in originals:
        result={"supplement_source_row":int(r["original_supplement_row"]),
                "published_source_name":r["original_published_scientific_name"]}
        for marker in SOURCE_MARKERS:
            accession=r[f"{marker}_GenBank_accession"]
            result[f"{marker}_original_accession"]=accession
            if not accession:
                result[f"{marker}_source_status"]="not_listed_in_supplement"
                continue
            seq=records.get(accession)
            if seq is None:
                result[f"{marker}_source_status"]="NCBI_unavailable_not_biological_locus_absence"
                continue
            clean,why=check_record(seq,r["original_published_scientific_name"],marker,accession)
            if clean is not None:
                result[f"{marker}_source_status"]="original_GenBank_exact_species_confirmed"
                result[f"{marker}_sequence"]=clean
            else:result[f"{marker}_source_status"]=why
        scored.append(result)
    # Prefer 2025 paper's paired ITS/ETS accession from one source row,
    # then ITS-only rows; do not count 9 repeated source name lines as 9
    # independent species/evolutionary replicates.
    bytax={}
    for row in scored:
        if row["ITS_source_status"]!="original_GenBank_exact_species_confirmed":
            continue
        taxon=row["published_source_name"]
        a=bool(row.get("ETS_sequence"))
        if taxon not in bytax or (
            (a, len(row["ITS_sequence"]["raw_nucleotide_sequence"]),
             -row["supplement_source_row"])>
            (bool(bytax[taxon].get("ETS_sequence")),
             len(bytax[taxon]["ITS_sequence"]["raw_nucleotide_sequence"]),
             -bytax[taxon]["supplement_source_row"])):
            bytax[taxon]=row
    outcome={
        "source_status":"T3_2025_AUTHENTIC_GENBANK_ITS_ETS_RETRIEVAL_PILOT" if taxa_limit else "T3_2025_AUTHENTIC_ALL163_ORIGINAL_ACCESSION_RETRIEVAL",
        "original_163_source_rows":len(rows),
        "original_154_unique_published_taxa":len({r["original_published_scientific_name"] for r in rows}),
        "source_rows_queried":len(scored),
        "remote_original_accessions_queried":len(requested),
        "n_original_genbank_sequences_fetched":len(records),
        "source_failure_batches_not_ITS_absence":err,
        "n_unique_taxa_with_exact_published_organism_ITS":len(bytax),
        "n_unique_taxa_with_exact_published_organism_both_ITS_and_ETS":sum(bool(r.get("ETS_sequence")) for r in bytax.values()),
        "n_selected_taxa_without_ETS_source":sum(not bool(r.get("ETS_sequence")) for r in bytax.values()),
        "original_ITS_ETS_rDNA_spacers_physically_linked_not_independent_nuclear_orthologs":True,
        "original_2021_48_target_capture_phylogeny_not_reproduced":True,
        "same_DNA_voucher_for_two_marker_accessions_not_proven":True,
        "T3_ancestral_nodding_transition_count_and_rain_fitness_not_inferred":True,
        "n_source_its_rejection":dict(Counter(r["ITS_source_status"] for r in scored)),
        "n_source_ets_rejection":dict(Counter(r["ETS_source_status"] for r in scored)),
        "ncbi_requests":fetcher.n_calls,
        "ncbi_download_bytes":fetcher.n_remote_bytes,
    }
    return scored,bytax,outcome


def export(rows,taxa,receipt,outdir,source_csv):
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    meta=[]
    for taxon,item in sorted(taxa.items()):
        meta.append({"published_source_species":taxon,
                     "SRA_published_2025_table_row":item["supplement_source_row"],
                     "published_ITS_GenBank_accession":item["ITS_original_accession"],
                     "published_ETS_GenBank_accession":item["ETS_original_accession"],
                     "ETS_source_verified":bool(item.get("ETS_sequence")),
                     "DNA_specimen_vouchers_for_both_markers_confirmed":False,
                     "taxon_tip_label":taxon.replace(" ","_")})
    with (out/"T3_2025_NUCLEAR_rDNA_SOURCE_ACCESSION_TAXON_MAPPING.csv").open(
        "w",newline="",encoding="utf-8") as h:
        columns=["published_source_species","SRA_published_2025_table_row",
                 "published_ITS_GenBank_accession","published_ETS_GenBank_accession",
                 "ETS_source_verified","DNA_specimen_vouchers_for_both_markers_confirmed",
                 "taxon_tip_label"]
        writer=csv.DictWriter(h,fieldnames=columns);writer.writeheader();writer.writerows(meta)
    for marker in SOURCE_MARKERS:
        with (out/f"T3_2025_ORIGINAL_GENBANK_{marker}_SEQUENCE.fasta").open("w",encoding="ascii") as f:
            for species,record in sorted(taxa.items()):
                seq=record.get(f"{marker}_sequence")
                if not seq:continue
                label=species.replace(" ","_")
                f.write(f">{label}\n")
                nuc=seq["raw_nucleotide_sequence"]
                for i in range(0,len(nuc),80):f.write(nuc[i:i+80]+"\n")
    compact=[]
    for z in rows:
        compact.append({k:v for k,v in z.items() if not k.endswith("_sequence")})
    with (out/"T3_2025_ALL163_ORIGINAL_PUBLISHED_SOURCE_ACCESSION_AUDIT.csv").open(
        "w",newline="",encoding="utf-8") as h:
        cols=["supplement_source_row","published_source_name",
              "ITS_original_accession","ITS_source_status","ETS_original_accession","ETS_source_status"]
        w=csv.DictWriter(h,fieldnames=cols,extrasaction="ignore")
        w.writeheader();w.writerows(compact)
    receipt["supplement_original_source_table_sha256"]=hashlib.sha256(Path(source_csv).read_bytes()).hexdigest()
    receipt["published_supplement_url"]=SOURCE_URL
    receipt["ITS_fasta_SHA256"]=hashlib.sha256((out/"T3_2025_ORIGINAL_GENBANK_ITS_SEQUENCE.fasta").read_bytes()).hexdigest()
    receipt["ETS_fasta_SHA256"]=hashlib.sha256((out/"T3_2025_ORIGINAL_GENBANK_ETS_SEQUENCE.fasta").read_bytes()).hexdigest()
    (out/"T3_2025_TRUE_SOURCE_GENBANK_ITS_ETS_RECOVERY_STATUS.json").write_text(
        json.dumps(receipt,indent=2,ensure_ascii=False,sort_keys=True,allow_nan=False)+"\n",
        encoding="utf-8")
    print("T3_NUCLEAR_MARKER_163_SOURCE="+json.dumps(receipt,ensure_ascii=False,sort_keys=True),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--original-2025-table",type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    p.add_argument("--limit",type=int)
    a=p.parse_args()
    original=read_csv(a.original_2025_table)
    validate_original(original)
    client=NCBI()
    scored,selected,receipt=collect_original(original,client,taxa_limit=a.limit)
    export(scored,selected,receipt,a.outdir,a.original_2025_table)
    if receipt["source_failure_batches_not_ITS_absence"]:
        raise RuntimeError("GenBank transport failed, report retained; no claim of missing biological loci")


if __name__=="__main__":main()
