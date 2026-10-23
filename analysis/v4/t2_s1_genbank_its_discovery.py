"""Retrieve actual GenBank nuclear ITS sequences for ALL S1 original Cirsium names.

Scientific identity:
- S1 remains 251 original source-name × WCVP-accepted-key records (244 keys),
  not a randomly selected subset or a reconstructed genus-wide census.
- The NCBI organism of the actual GB accession must match the original
  biological binomial EXACTLY. No blind accepted-name or synonym collapse.
- ITS is multi-copy nuclear ribosomal DNA; this cannot alone determine a
  single-copy species tree, infer introgression, or date adaptive gains.
- No photo EXIF trait enters tree reconstruction. Observed present CHELSA
  climate and independently documented botanical posture are JOINED ONLY
  after the sequence/accession tree is built.
- Source failures are NOT GenBank negative searches; they stay unknown.
NCBI EUtilities <=2.3 requests/second and bounded retry; all 251 taxa, no
user-supplied or guessed email address. See NCBI EUtilities usage guidelines.
"""
from __future__ import annotations
import argparse
from collections import Counter
import csv
from io import StringIO
import hashlib
import json
from pathlib import Path
import re
import time
from urllib.parse import urljoin

from Bio import SeqIO
import requests

S1_N=251
S1_WCVP=244
S1_IMAGE_COUNT=170066
S1_NUCLEAR=111
S1_SCORED=66
S1_SOURCE_BASENAME="20261023_T2_S1_ALL251_66_TO72_PROVISIONAL_POSTURE_OVERLAY.csv"
NCBI_BASE="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
NCBI_TOOL="azamiS1CirsiumPhylogeny"
MAX_HITS_PER_TAXON=8
REQUEST_GAP_SEC=.48
MIN_ITS_BP=350
MAX_ITS_BP=1100
MAX_AMBIGUOUS_FRAC=.035
ITS_PATTERN=re.compile(r"(?:\bITS(?:\s*[-_]?\s*[12])?\b|internal transcribed spacer)",re.I)
CHLORO_PATTERN=re.compile(r"(?:chloroplast|plastid|mitochondr|genome assembly)",re.I)
FIELDS=("accepted_key","original_source_name","original_photo_count",
        "botanical_U_D_state","chelsa_bio12","chelsa_bio01")


def read_s1(path):
    with Path(path).open(newline="",encoding="utf-8-sig") as handle:
        rows=list(csv.DictReader(handle))
    if len(rows)!=S1_N or len({x["accepted_key"] for x in rows})!=S1_WCVP:
        raise ValueError("Original S1 251 names / 244 accepted taxonomy changed")
    if len({(x["accepted_key"],x["original_source_name"]) for x in rows})!=S1_N:
        raise ValueError("Duplicate original S1 taxa")
    if sum(int(x["original_photo_count"]) for x in rows)!=S1_IMAGE_COUNT:
        raise ValueError("Original 170,066 S1 specimens missing")
    if sum(x["botanical_U_D_state"] in ("U","D") for x in rows)!=S1_SCORED:
        raise ValueError("Original independent botanical phenotype census changed")
    if sum(x["nuclear_real120k_BioSample"]=="True" for x in rows)!=S1_NUCLEAR:
        raise ValueError("S1 real nuclear SRA sample count changed")
    for x in rows:
        if (not x["original_source_name"].startswith("Cirsium ")
                or len(x["original_source_name"].split())!=2
                or not x["original_source_name"].split()[1].islower()
                or x["botanical_U_D_state"] not in ("D","U","?")
                or float(x["chelsa_bio12"])<=0):
            raise ValueError("S1 source original taxon/rain/author coding invalid")
    return rows


class Eutils:
    def __init__(self, session=None, delay=REQUEST_GAP_SEC, *,max_calls=1050):
        self.session=session or requests.Session()
        self.delay=delay
        self.calls=0
        self.max_calls=max_calls
        self.next_time=0.

    def get(self,method,params):
        if method not in ("esearch.fcgi","efetch.fcgi"):
            raise ValueError("Unlisted NCBI EUtilities method")
        request={"db":"nuccore","tool":NCBI_TOOL,**params}
        last=None
        for attempt in range(4):
            if self.calls>=self.max_calls:
                raise RuntimeError("NCBI source request budget reached before all species reviewed")
            pause=self.next_time-time.monotonic()
            if pause>0:time.sleep(pause)
            self.next_time=time.monotonic()+self.delay
            self.calls+=1
            try:
                r=self.session.get(NCBI_BASE+method,params=request,
                                   headers={"User-Agent":"azami-S1-public-sequences/1.0"},
                                   timeout=(10,45))
                if r.status_code in (429,500,502,503,504):
                    last=RuntimeError("EUtilities temporary status "+str(r.status_code))
                    time.sleep(min(2**attempt,8))
                    continue
                r.raise_for_status()
                if len(r.content)>2_500_000:
                    raise ValueError("One NCBI response oversized; refuse unverifiable dump")
                return r
            except requests.RequestException as e:
                last=e
                time.sleep(min(2**attempt,8))
        raise RuntimeError("NCBI retrieval not verified after retries: "+str(last)[:140])

    def search(self,name):
        if not name.startswith("Cirsium ") or len(name.split())!=2:
            raise ValueError("Only original S1 two-word biological species names")
        term=(f'"{name}"[Organism] AND '
              '(("internal transcribed spacer"[Title] OR ITS[Title] '
              'OR ITS1[Title] OR ITS2[Title] OR "ribosomal DNA"[Title]) '
              'AND 300:1200[SLEN]) NOT (uncultured[Title] OR metagenom*[Title])')
        r=self.get("esearch.fcgi",{"term":term,"retmode":"json",
                                   "retmax":MAX_HITS_PER_TAXON})
        obj=r.json().get("esearchresult",{})
        if not isinstance(obj.get("idlist"),list) or "count" not in obj:
            raise ValueError("NCBI ESearch missing count/idlist; do not infer absence")
        return obj["idlist"],int(obj["count"]),term

    def fetch(self,uids):
        if not uids:return []
        if len(uids)>MAX_HITS_PER_TAXON or any(not str(x).isdigit() for x in uids):
            raise ValueError("Unknown or oversized source GenBank UID batch")
        r=self.get("efetch.fcgi",{"id":",".join(uids),
                                   "rettype":"gb","retmode":"text"})
        try:
            seqs=list(SeqIO.parse(StringIO(r.text),"genbank"))
        except Exception as e:
            raise ValueError("Retrieved NCBI GenBank INSDC records could not be parsed") from e
        if not seqs:
            raise ValueError("ESearch returned UIDs but EFetch did not return sequence records")
        return seqs


def qualifying_its(record,name):
    organism=" ".join(str(record.annotations.get("organism","")).split())
    if organism.casefold()!=name.casefold():
        return None,"different_NCBI_source_organism"
    title=str(record.description)
    if CHLORO_PATTERN.search(title):
        return None,"organelle_or_genome_title_excluded"
    if not ITS_PATTERN.search(title):
        return None,"ITS_marker_not_explicitly_labeled"
    qualifiers={}
    for feat in record.features:
        if feat.type=="source":
            qualifiers=feat.qualifiers
            break
    if "organelle" in qualifiers:
        return None,"organelle_feature_excluded"
    seq=str(record.seq).upper().replace("U","T")
    if not MIN_ITS_BP<=len(seq)<=MAX_ITS_BP:
        return None,"ITS_length_outside_common_region"
    nonbase=sum(z not in "ACGT" for z in seq)
    if nonbase/len(seq)>MAX_AMBIGUOUS_FRAC:
        return None,"ITS_too_many_ambiguous_bases"
    if any(z not in "ACGTNRYWSKMBDHV" for z in seq):
        return None,"invalid_IUPAC_character"
    accession=str(record.id)
    if not re.fullmatch(r"[A-Z]{1,5}_?\d{4,12}(?:\.\d+)?",accession):
        return None,"unrecognised_public_INSDC_accession"
    # Compare marker annotations and specimen vouchers. Do not assume ITS
    # shares orthology with all copies in a polyploid or hybrid Cirsium.
    definition=title.lower()
    voucher=("|".join(qualifiers.get("specimen_voucher",[])) or
             "|".join(qualifiers.get("isolate",[])) or "")
    score=(3*("internal transcribed spacer" in definition)
           +2*("5.8" in definition or "its1" in definition and "its2" in definition)
           +1*bool(voucher)
           +min(len(seq),800)/800
           -nonbase/len(seq))
    return {
        "source_original_species":name,
        "GenBank_accession_version":accession,
        "GenBank_original_organism":organism,
        "original_GenBank_definition":title,
        "original_GenBank_record_length_bp":len(seq),
        "source_specimen_voucher_or_isolate":voucher[:150],
        "seq_A_C_G_T_fraction":round(1-nonbase/len(seq),6),
        "S1_species_to_GB_record_exact_name":True,
        "ITS_phylogeny_single_copy_orthology_confirmed":False,
        "polyploidy_ITS_homogenization_evaluated":False,
        "score_for_preferring_vouchered_full_ITS":round(score,6),
        "nucleotide_sequence":seq
    },None


def discover(rows,client,max_species=None):
    scanned=(rows if max_species is None else rows[:max_species])
    ledger=[];selected=[];all_excluded=Counter()
    for n,r in enumerate(scanned,1):
        name=r["original_source_name"]
        status="not_checked";hits=0;selected_row=None;fail=None;term=""
        try:
            uids,hits,term=client.search(name)
            if not uids:status="NCBI_search_no_ITS_candidates"
            else:
                best=[]
                for candidate in client.fetch(uids):
                    ok,reason=qualifying_its(candidate,name)
                    if ok:best.append(ok)
                    else:all_excluded[reason]+=1
                if best:
                    best.sort(key=lambda z:(-z["score_for_preferring_vouchered_full_ITS"],
                                              z["GenBank_accession_version"]))
                    selected_row=best[0]
                    status="exact_species_ITS_accession_verified"
                    selected.append({**{k:r[k] for k in FIELDS},**selected_row,
                                     "all_matching_ncbi_GB_candidates_considered":len(best)})
                else:status="NCBI_ITS_candidates_not_exact_sufficient_annotation"
        except Exception as exc:
            status="NCBI_access_unresolved_not_biological_absence"
            fail=f"{type(exc).__name__}: {str(exc)[:165]}"
        ledger.append({
            **{k:r[k] for k in FIELDS},
            "NCBI_ITS_discovery_status":status,
            "NCBI_candidate_search_count":hits,
            "selected_GenBank_accession_version":(
                selected_row["GenBank_accession_version"] if selected_row else ""),
            "source_access_error":fail or "",
            "original_NCBI_Esearch_term":term,
            "n_original_source_photo_observations":r["original_photo_count"],
            "source_only_GenBank_ITS_not_assembled_single_copy_nuclear_species_tree":True,
        })
        if n%20==0:
            print(f"GENBANK_S1_REVIEW_PROGRESS={n}/{len(scanned)} accepted_exact={len(selected)}",flush=True)
    if len(ledger)!=len(scanned) or len({r["original_source_name"] for r in ledger})!=len(scanned):
        raise ValueError("Original S1 species must each have one complete audit row")
    return ledger,selected,all_excluded


def write_csv(path,rows,headers=None):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    cols=headers or list(rows[0].keys())
    with path.open("w",newline="",encoding="utf-8") as h:
        w=csv.DictWriter(h,fieldnames=cols,extrasaction="raise")
        w.writeheader();w.writerows(rows)


def export(source,ledger,selected,excluded,outdir,*,max_species=None,request_count=0):
    target=Path(outdir);target.mkdir(parents=True,exist_ok=True)
    if len(selected)!=sum(r["NCBI_ITS_discovery_status"]=="exact_species_ITS_accession_verified" for r in ledger):
        raise ValueError("Biological sequence success did not match original S1 ledger")
    write_csv(target/"S1_GENBANK_ITS_ALL_ORIGINAL_NAMES_SOURCE_AUDIT.csv",ledger)
    meta=[{k:v for k,v in row.items() if k!="nucleotide_sequence"} for row in selected]
    write_csv(target/"S1_GENBANK_ITS_EXACT_SPECIES_ACCESSIONS.csv",meta,
              headers=(list(meta[0]) if meta else list(FIELDS)+[
                  "GenBank_accession_version","GenBank_original_organism"]))
    path=target/"S1_GENBANK_ITS_EXACT_ORGANISM_SEQUENCE.fasta"
    with path.open("w",encoding="ascii") as f:
        for entry in selected:
            # IQ-TREE tip is exact SOURCE binomial: accession stored separately.
            label=entry["source_original_species"].replace(" ","_")
            f.write(f">{label}\n")
            seq=entry["nucleotide_sequence"]
            for i in range(0,len(seq),80):f.write(seq[i:i+80]+"\n")
    outcome=Counter(x["NCBI_ITS_discovery_status"] for x in ledger)
    unknown=outcome.get("NCBI_access_unresolved_not_biological_absence",0)
    scope="S1_ALL251" if max_species is None else "S1_PREFIX_PILOT"
    summary={
        "schema":"azami_T2_S1_original251_GenBank_exact_binomial_ITS_source_index_v1",
        "status":("S1_GENBANK_FULL_ORIGINAL_NAMES_SOURCE_DISCOVERY_COMPLETED" if
            max_species is None and unknown==0 else "PILOT_OR_SOURCE_ACCESS_NOT_COMPLETED"),
        "scope":scope,
        "source_original_S1_master_names":S1_N,
        "source_original_WCVP_accepted_keys":S1_WCVP,
        "total_source_photographs_retained":S1_IMAGE_COUNT,
        "source_names_actually_queried":len(ledger),
        "all251_queries_completed_without_network_failures":len(ledger)==S1_N and unknown==0,
        "accepted_exact_species_nuclear_ITS_accessions":len(selected),
        "GenBank_record_status_counts":dict(outcome),
        "rejected_candidate_record_reasons":dict(excluded),
        "NCBI_requests_including_retries":request_count,
        "selected_accessions_are_exact_original_organism_names":True,
        "some_GenBank_sequences_may_span_only_ITS1_or_ITS2":True,
        "conspecific_multi_copy_and_polyploid_ITS_homogeneity_not_audited":True,
        "original_S1_bio12_centroids_are_present_distribution_not_evolutionary_histories":True,
        "head_posture_binary_states_independent_botanical_authorities_only":True,
        "S1_photo_EXIF_angles_not_used_for_GenBank_species_selection":True,
        "whole_S1_251_tip_verified_true_nuclear_species_tree_completed":False,
        "mitochondrial_or_chloroplast_sequences_excluded_from_nuclear_ITS":True,
        "source_file_sha256":hashlib.sha256(Path(source).read_bytes()).hexdigest(),
        "output_accession_ledger_SHA256":hashlib.sha256(
            (target/"S1_GENBANK_ITS_EXACT_SPECIES_ACCESSIONS.csv").read_bytes()).hexdigest(),
        "output_original_ITS_FASTA_SHA256":hashlib.sha256(path.read_bytes()).hexdigest(),
    }
    (target/"S1_GENBANK_SOURCE_COVERAGE_RECEIPT.json").write_text(
        json.dumps(summary,sort_keys=True,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("S1_GENBANK_REAL_ITS_SOURCE=" + json.dumps(summary,sort_keys=True),flush=True)
    return summary


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--s1",type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    p.add_argument("--limit",type=int)
    args=p.parse_args()
    original=read_s1(args.s1)
    if args.limit is not None and not 1<=args.limit<=S1_N:
        raise ValueError("Pilot limit outside original 251")
    client=Eutils()
    ledger,seq,excluded=discover(original,client,args.limit)
    result=export(args.s1,ledger,seq,excluded,args.outdir,
                  max_species=args.limit,request_count=client.calls)
    if result["GenBank_record_status_counts"].get("NCBI_access_unresolved_not_biological_absence",0):
        raise RuntimeError("Public NCBI access failed: not legitimate to call these taxa ITS-negative")


if __name__=="__main__":
    main()
