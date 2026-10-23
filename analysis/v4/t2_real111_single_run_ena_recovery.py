"""Recover ONE genuine Moreyra 120k public accession after bounded ENA gzip prefix failure.

Do not redo 110 already confirmed original Cirsium SRA data. The first 120000
paired reads are unchanged; only the maximum HTTP compressed byte prefix for
SRR30887180 is raised from 64 MiB to 128 MiB. Exact source BioSample,
predeclared 118 original reference copies, baseQ20/mapQ20, >=3 read Q20
majority >=0.90 and biological uncertainty stop gates remain unchanged.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

from analysis.v4 import t2_four_species_30k_120k_nuclear_depth as readops
from analysis.v4 import t2_real111_fixed118_nucleotides as original
from analysis.v4 import t2_moreyra_all125_ena_read_mapping as read_registry

RUN="SRR30887180"
BIO="SAMN44017849"
TAXON="Cirsium flodmanii"
MAX_COMPRESSED_BYTES=128*1024*1024
TRY_MAX=(16*1024*1024,32*1024*1024,64*1024*1024,
         96*1024*1024,128*1024*1024)


def bounded_extended_fetch(row):
    if row["run"]!=RUN or row["biosample"]!=BIO:
        raise ValueError("Targeted rescue must never silently include other biological specimens")
    original_limit=readops.MAX_BYTES
    try:
        readops.MAX_BYTES=TRY_MAX
        return readops.acquire_larger_prefix(row)
    finally:
        readops.MAX_BYTES=original_limit


def recover(source125,selected111,targets,reference,scratch,worker=original.map_and_call):
    byrun,selected,copies=original.audit_samples(source125,selected111,targets)
    if RUN not in selected or RUN not in byrun:
        raise ValueError("Original 111 source SRA count does not include previously successful seven-sample run")
    specimen=selected[RUN]
    if specimen["biosample"]!=BIO or specimen["taxon_name"]!=TAXON:
        raise ValueError("Published 7-source DNA sample no longer matches original taxon label")
    item=worker(byrun[RUN],specimen,copies,reference,scratch,
                acquire=bounded_extended_fetch)
    original.validate_one(item,specimen,copies)
    if item["n_real_Q20_majority_callable_positions"]<=0:
        raise ValueError("Actual recovered source run has no observed Q20 major bases")
    item["rescue_protocol"]="same_real_first120000_pairs_larger_compressed_HTTP_prefix_only"
    item["rescue_max_compressed_bytes_per_mate"]=MAX_COMPRESSED_BYTES
    item["rescue_not_random_reads_or_new_biological_sample"]=True
    item["source_first64MiB_failure_preserved_in_original_10_shards"]=True
    return item


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--original125",type=Path,required=True)
    p.add_argument("--selected111",type=Path,required=True)
    p.add_argument("--targets",type=Path,required=True)
    p.add_argument("--reference",type=Path,required=True)
    p.add_argument("--scratch",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    source=read_registry.load_source(a.original125)
    samples=original.load_csv(a.selected111,original.SOURCE_TAXA_FIELDS)
    copies=original.load_csv(a.targets,original.TARGET_FIELDS)
    result=recover(source,samples,copies,a.reference,a.scratch)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,
                                   allow_nan=False)+"\n",encoding="utf-8")
    print("T2_EXACT_ONE_111_SOURCE_SRA_GZIP_PREFIX_RESCUED="+json.dumps({
        "run":RUN,"biosample":BIO,"source_Q20_bases":
        result["n_real_Q20_majority_callable_positions"]},sort_keys=True),flush=True)


if __name__=="__main__":
    main()
