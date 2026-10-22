"""Probe NCBI SRA run metadata sizes for ALL 125 original source-matched Cirsium runs.

Metadata only; do not download raw reads to GitHub Actions. A verified 244-key
source/125-run ledger is the input. NCBI reports archived run size in MB, NOT
uncompressed FASTQ nor storage needed for multi-locus inference.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
import csv
import io
import json
from pathlib import Path
import time
import urllib.error
import urllib.parse
import urllib.request

BIOPROJECT="PRJNA957074"
RUN_N=125
URL="https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo"
FIELDS=("accepted_key","S1_source_taxon","original_tree_code","biosample","run",
        "experiment","voucher_and_herbarium","scope_class",
        "name_reconciliation_priority","tree_code_vs_sra_name",
        "accepted_concept_verified")


def source_runs(path):
    with Path(path).open(newline="",encoding="utf-8-sig") as f:
        rows=list(csv.DictReader(f))
    if len(rows)!=RUN_N or set(rows[0])!=set(FIELDS):
        raise ValueError("Original Moreyra source-voucher run inventory changed")
    runs=[r["run"] for r in rows]
    if len(set(runs))!=RUN_N or not all(r.startswith("SRR") for r in runs):
        raise ValueError("Expected 125 unique source-verified public SRR runs")
    return rows


def parse_runinfo(text):
    """Official SRA backend runinfo is plain text; CSV or TSV headers are valid."""
    text=text.lstrip("\ufeff").strip()
    if not text or text.startswith("<") or text.startswith("{"):
        raise ValueError("NCBI returned no delimited RunInfo metadata")
    first=text.splitlines()[0]
    delimiter="\t" if first.count("\t")>first.count(",") else ","
    parsed=list(csv.DictReader(io.StringIO(text),delimiter=delimiter))
    if not parsed:
        raise ValueError("NCBI RunInfo returned a header-only response")
    keyset=set(parsed[0])
    runfield="Run" if "Run" in keyset else "run" if "run" in keyset else "acc"
    biofield="BioSample" if "BioSample" in keyset else "biosample"
    sizefield="size_MB" if "size_MB" in keyset else "MBytes"
    if {runfield,biofield,sizefield}-keyset:
        raise ValueError("NCBI backend runinfo lacks run/BioSample/size columns: "+
                         ",".join(sorted(keyset)))
    result={}
    for row in parsed:
        run=(row.get(runfield) or "").strip()
        if not run or not run.startswith(("SRR","ERR","DRR")):
            continue
        if run in result:
            raise ValueError("Duplicate public SRA RunInfo run accession")
        try:
            size=float(row[sizefield])
        except (ValueError,TypeError):
            continue
        if not 0<size<1e8:
            continue
        result[run]={
            "Run":run,"BioSample":(row.get(biofield) or "").strip(),
            "size_MB":size,
            "spots":row.get("spots",""),
            "bases":row.get("bases",""),
            "LibraryStrategy":row.get("LibraryStrategy",""),
            "LibraryLayout":row.get("LibraryLayout",""),
            "Platform":row.get("Platform",""),
            "BioProject":row.get("BioProject",""),
        }
    return result


def download_project_runinfo(requested,fetch=None):
    """Fetch explicit accessions in bounded batches; never full public reads."""
    if fetch is None:
        def fetch(url):
            req=urllib.request.Request(url,headers={
                "User-Agent":"azami-source-verified-public-SRA-run-size-audit/1.0",
                "Accept":"text/plain,text/csv,*/*"})
            with urllib.request.urlopen(req,timeout=90) as response:
                return response.read().decode("utf-8-sig")
    accessions=sorted(set(requested))
    if len(accessions)!=RUN_N or not all(x.startswith("SRR") for x in accessions):
        raise ValueError("Only 125 original public runs can enter the size screen")
    allruns={}
    attempts=[]
    for start in range(0,len(accessions),20):
        batch=accessions[start:start+20]
        url=URL+"?"+urllib.parse.urlencode({"acc":",".join(batch)})
        errors=[]
        rows=None
        for attempt in range(3):
            try:
                rows=parse_runinfo(fetch(url))
                break
            except urllib.error.HTTPError as err:
                errors.append(f"HTTP {err.code}")
                if err.code in (400,401,403,404):
                    break
            except (urllib.error.URLError,TimeoutError,ValueError) as err:
                errors.append(str(err))
            if attempt<2:
                time.sleep(attempt+1)
        if rows is None:
            raise RuntimeError(
                f"NCBI official backend runinfo failed on batch {start//20}: "+
                " / ".join(errors))
        unexpected=set(rows)-set(batch)
        if unexpected:
            raise ValueError("NCBI backend returned unexpected SRA run accessions")
        allruns.update(rows)
        attempts.append({"batch":start//20,"requested":len(batch),
                         "returned":len(rows)})
    if set(allruns)!=set(accessions):
        raise ValueError("NCBI backend lacks requested runs: "+
                         str(sorted(set(accessions)-set(allruns))[:6]))
    return allruns,URL,attempts


def audit(source,records):
    requested={r["run"] for r in source}
    matched=requested & records.keys()
    missing=sorted(requested-records.keys())
    if missing:
        raise ValueError(f"NCBI RunInfo lacks {len(missing)} pinned public runs: {missing[:5]}")
    out=[]
    for row in source:
        src=records[row["run"]]
        if src["BioSample"]!=row["biosample"]:
            raise ValueError("NCBI run/BioSample identifier differs from EAzami published reconciliation")
        if src["BioProject"] and src["BioProject"]!=BIOPROJECT:
            raise ValueError("Unexpected BioProject on source-pinned run")
        out.append({
            "S1_accepted_key":row["accepted_key"],
            "source_taxon":row["S1_source_taxon"],
            "SRA_run":row["run"],
            "BioSample":row["biosample"],
            "size_MB_compressed_SRA_metadata":src["size_MB"],
            "library_strategy":src["LibraryStrategy"],
            "library_layout":src["LibraryLayout"],
            "platform":src["Platform"],
            "source_name_priority":row["name_reconciliation_priority"],
            "source_taxon_name_disagreement":row["tree_code_vs_sra_name"]!="exact",
            "voucher_accepted_concept_verified":False,
        })
    sizes=sorted((r["size_MB_compressed_SRA_metadata"] for r in out))
    total=float(sum(sizes))
    summaries={
        "status":"PUBLIC_NCBI_SRA_SIZE_METADATA_SOURCE_VERIFIED_NO_FASTQ_DOWNLOADED",
        "source_BioProject":BIOPROJECT,
        "source_matched_run_count":len(out),
        "source_matched_accepted_key_count":len({r["S1_accepted_key"] for r in out}),
        "compressed_SRA_size_total_MB":total,
        "compressed_SRA_size_total_GB_decimal":total/1000,
        "compressed_SRA_size_median_MB":float(sizes[len(sizes)//2]),
        "compressed_SRA_size_max_MB":max(sizes),
        "n_runs_with_original_name_disagreement":sum(r["source_taxon_name_disagreement"] for r in out),
        "name_voucher_identity_verified":False,
        "SRA_runinfo_compressed_size_not_FASTQ_or_peak_disk":True,
        "all_244_source_species_have_nuclear_data":False,
        "finished_original_350_locus_tree":False,
        "all_125_fastq_downloaded":False,
    }
    return out,summaries


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--crosswalk",type=Path,required=True)
    p.add_argument("--output-dir",type=Path,required=True)
    a=p.parse_args()
    src=source_runs(a.crosswalk)
    allrun,url,batches=download_project_runinfo([r["run"] for r in src])
    rows,summary=audit(src,allrun)
    summary["NCBI_documented_SRA_backend_RunInfo"]=url
    summary["bounded_size_batches"]=batches
    a.output_dir.mkdir(parents=True,exist_ok=True)
    with (a.output_dir/"S1_MOREYRA_125_PUBLIC_RUN_SIZE_METADATA.csv").open(
            "w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),extrasaction="raise")
        w.writeheader()
        w.writerows(rows)
    (a.output_dir/"S1_MOREYRA_125_RUN_SIZE_SUMMARY.json").write_text(
        json.dumps(summary,ensure_ascii=False,indent=2,sort_keys=True)+"\n",
        encoding="utf-8")
    print("T2_MOREYRA_NCBI_125_RUN_SIZE="+json.dumps(summary,sort_keys=True,
                                               ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
