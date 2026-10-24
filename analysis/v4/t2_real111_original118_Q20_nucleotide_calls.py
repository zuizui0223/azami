"""Run REAL 120k paired public ENA DNA reads on 111 S1-linked taxa and 118 loci.

Source: SHA-verified original EAzami/Moreyra 125-run bridge; actual 111
one-BioSample-per-source-WCVP preflight and 118 original Compositae1061
reference copies supported by >=89/111 specimens at Q20 depth3 25% breadth.
Calls require mappingQ>=20, baseQ>=20, >=3 reads and >=90% majority base,
following the tested prior seven-specimen empirical site-calling contract.

Unobserved positions and possible heterozygotes are missing (no reference
imputation). This is reference-coordinate homolog screening, not copy-number
validated orthologs, a dated phylogeny, or head-angle adaptation.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import shutil
import time

from analysis.v4 import t2_moreyra_all125_ena_read_mapping as original
from analysis.v4 import t2_four_species_30k_120k_nuclear_depth as readops
from analysis.v4 import t2_moreyra_all125_author241_site_depth as depth
from analysis.v4 import t2_seven_real120k_shared_sites as caller

N_SOURCE=125
N_MEASURED=111
N_TARGETS=118
N_SHARDS=10
QUORUM=89
READS=120000
STATUS="ACTUAL_ORIGINAL111_REAL120K_Q20_SITE_CALL_NOT_CERTIFIED_ORTHOLOGY"
BOTANICAL_NODDING={"Cirsium conspicuum","Cirsium scopulorum"}


def read_csv(path):
    with Path(path).open(newline="",encoding="utf-8-sig") as f:
        result=list(csv.DictReader(f))
    if not result:
        raise ValueError("Empty actual nucleotide-source evidence")
    return result


def source_inputs(original_rows,samples,targets):
    if len(original_rows)!=N_SOURCE or len({r["run"] for r in original_rows})!=N_SOURCE:
        raise ValueError("Missing one of all125 original matched Moreyra reads")
    if len(samples)!=N_MEASURED or len({r["run"] for r in samples})!=N_MEASURED:
        raise ValueError("Preflight source 111 BioSamples not distinct")
    if len(targets)!=N_TARGETS:
        raise ValueError("Original 118 nuclear targets not preserved")
    src={r["run"]:r for r in original_rows}
    tip_seen=set()
    biological_seen=set()
    out=[]
    for row in samples:
        run,bio,key=row["run"],row["biosample"],row["accepted_key"]
        if run not in src:
            raise ValueError("Selected S1 source accession absent from Moreyra original125")
        org=src[run]
        if (org["biosample"]!=bio or org["accepted_key"]!=key or
            "WCVP"+key != row["phylogenetic_tip"]):
            raise ValueError("Run, BioSample or original WCVP source concept diverges")
        if row["phylogenetic_tip"] in tip_seen or bio in biological_seen:
            raise ValueError("Duplicated biological specimen/tip across source runs")
        named=lambda x:" ".join(str(x).replace("_"," ").split())
        if named(org["S1_source_taxon"])!=named(row["taxon_name"]):
            raise ValueError("Public Moreyra source name differs from preflight original")
        tip_seen.add(row["phylogenetic_tip"])
        biological_seen.add(bio)
        out.append((row,org))
    names={s["taxon_name"] for s,_ in out}
    if not BOTANICAL_NODDING<=names:
        raise ValueError("Both independent author-scored nodding taxon names lost")
    contigs={}
    for row in targets:
        locus,contig=row["locus"],row["original_compositae_reference"]
        if not contig.startswith(("lett-","saff-","sunf-")) or not contig.endswith("-"+locus):
            raise ValueError("Unverified original single reference sequence for nuclear target")
        if int(row["minimum_source_taxa_required"])!=QUORUM:
            raise ValueError("Source 80% taxon occupancy gate changed")
        if int(row["empirically_supported_source_taxa"])<QUORUM:
            raise ValueError("A source reference copy lacks 89 valid Q20 reference-mapped taxa")
        if locus in contigs:
            raise ValueError("Duplicate original nucleus locus")
        contigs[locus]=contig
    if len(set(contigs.values()))!=N_TARGETS:
        raise ValueError("One reference contig assigned to two original loci")
    return sorted(out,key=lambda r:r[0]["phylogenetic_tip"]),contigs


def source_call(one,original_row,contigs,author241,reference,scratch):
    run=one["run"]
    directory=Path(scratch)/run
    directory.mkdir(parents=True,exist_ok=True)
    try:
        mates,_metadata=readops.acquire_larger_prefix(original_row)
        if len(mates)!=2 or any(len(reads)!=4*READS for reads in mates):
            raise ValueError("This original source run did not provide exactly 120k paired reads")
        qc=readops.map_read_pairs(original_row,reference,author241,mates,directory)
        if qc["input_read_pairs"]!=READS:
            raise ValueError("BWA original run no longer matches claimed 120k reads")
        calls=caller.call_sites_from_pysam(directory/"aligned.sorted.bam",contigs)
        if set(calls)!=set(contigs):
            raise ValueError("One original 118 nuclear ref copy missing")
        for locus,rec in calls.items():
            if rec["reference"]!=contigs[locus] or not rec["not_a_validated_allelic_genotype"]:
                raise ValueError("Site basecalls swapped contig identity or claimed validated genotype")
        summary={
            "source_original_run":run,"source_BioSample":one["biosample"],
            "source_S1_accepted_key":one["accepted_key"],
            "taxon_name":one["taxon_name"],"tip":one["phylogenetic_tip"],
            "n_real_read_pairs":READS,"n_original_common_reference_loci":N_TARGETS,
            "n_Q20_majority_called_sites":sum(z["n_conservatively_callable_sites"] for z in calls.values()),
            "n_targets_at_least100_Q20_called_sites":sum(
                z["n_conservatively_callable_sites"]>=100 for z in calls.values()),
            "original_author_botanical_posture":(
                "nodding" if one["taxon_name"] in BOTANICAL_NODDING else "not_independently_scored_here"),
            "orthology_copy_number_validated":False,
            "gravity_head_orientation_inferred_from_image":False,
            "read_fastq_bam_in_public_artifact":False,
            "loci":calls,
        }
        return summary
    finally:
        shutil.rmtree(directory,ignore_errors=True)


def collect(original_rows,samples,targets,author241,reference,scratch,shard,runner=source_call):
    source,contigs=source_inputs(original_rows,samples,targets)
    if not 0<=shard<N_SHARDS:
        raise ValueError("Shard 0..9 required")
    if len(author241)!=241 or not set(contigs).issubset(author241):
        raise ValueError("Original source filtered 241 gene name universe mismatches 118")
    chosen=[x for i,x in enumerate(source) if i%N_SHARDS==shard]
    if len(chosen) not in (11,12):
        raise ValueError("One real111 source shard has unexpected number of nuclear vouchers")
    outputs=[]
    for i,(sample,org) in enumerate(chosen):
        start=time.monotonic()
        try:
            x=runner(sample,org,contigs,author241,reference,scratch)
            if x["source_original_run"]!=sample["run"] or x["tip"]!=sample["phylogenetic_tip"]:
                raise ValueError("Real source accession changed after sampling")
            if set(x["loci"])!=set(contigs):
                raise ValueError("Missing original Compositae1061 nuclear target")
            x["status"]="real_Q20_majority_reference_site_calls"
            x["elapsed_seconds"]=round(time.monotonic()-start,2)
            outputs.append(x)
            print("T2_REAL111_Q20_SITE_OK="+json.dumps({
                "run":sample["run"],"called":x["n_Q20_majority_called_sites"],
                "loci_100nt":x["n_targets_at_least100_Q20_called_sites"]},
                sort_keys=True),flush=True)
        except Exception as err:
            fail={
                "status":"failed","source_original_run":sample["run"],
                "source_BioSample":sample["biosample"],
                "tip":sample["phylogenetic_tip"],
                "error_type":type(err).__name__,"error":str(err)[:300],
                "elapsed_seconds":round(time.monotonic()-start,2),
            }
            outputs.append(fail)
            print("T2_REAL111_Q20_SITE_FAILED="+json.dumps(fail,sort_keys=True),flush=True)
    if len(outputs)!=len(chosen) or len({x["source_original_run"] for x in outputs})!=len(chosen):
        raise ValueError("One of original selected111 SRA specimen sample outputs lost")
    return {
      "schema":"azami_T2_111_actual120k_118original_Q20_consensus_shard_v1",
      "status":STATUS,"shard":shard,"n_shards":N_SHARDS,
      "n_source_full_S1_WCVP":244,"n_source_S1_named_units":251,
      "n_source_nuclear_BioSamples":N_MEASURED,"n_source_original_loci":N_TARGETS,
      "requested":len(chosen),
      "success":sum(x["status"]=="real_Q20_majority_reference_site_calls" for x in outputs),
      "failure":sum(x["status"]=="failed" for x in outputs),
      "raw_reads_or_BAM_publication":False,
      "single_original_reference_not_confirmed_ortholog":True,
      "genus_species_phylogeny_inferred":False,
      "outcomes":outputs,
    }


def main():
    a=argparse.ArgumentParser(description=__doc__)
    for key in ("--source","--sample-panel","--target-panel","--author241","--reference","--scratch","--out"):
        a.add_argument(key,type=Path,required=True)
    a.add_argument("--shard",type=int,required=True)
    opts=a.parse_args()
    original_rows=original.load_source(opts.source)
    samples=read_csv(opts.sample_panel)
    targets=read_csv(opts.target_panel)
    loci=depth.loci_from_pinned_author_list(opts.author241)
    if not opts.reference.is_file():
        raise ValueError("Source-authenticated original Compositae1061 reference missing")
    output=collect(original_rows,samples,targets,loci,opts.reference,opts.scratch,opts.shard)
    opts.out.parent.mkdir(parents=True,exist_ok=True)
    opts.out.write_text(json.dumps(output,sort_keys=True,ensure_ascii=False,
                                    indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print("T2_REAL111_Q20_CONSENSUS_SHARD="+json.dumps({
        "id":opts.shard,"requested":output["requested"],
        "success":output["success"],"failed":output["failure"]},
        sort_keys=True),flush=True)


if __name__=="__main__":
    main()
