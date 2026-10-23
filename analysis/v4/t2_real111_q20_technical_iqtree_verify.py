"""Validate REAL original111 118-locus masked-Q20 technical ML tree and sample IDs.

The 111 source tips are directly linked to public Moreyra/EAzami verified runs.
Original genomic reference-copy Q20 coordinates were called on 120k paired
public reads, but still require paralogy/homeology and gene tree checks.
Do not convert these data into an accepted 294-tip species tree, clock,
trait transition, or rain adaptive radiation test.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
from Bio import Phylo

N_TIPS=111
N_ORIGINAL=244
N_SOURCE_NAMES=251


def verify(input_manifest,original_samples,tree_file,output):
    summary=json.loads(Path(input_manifest).read_text(encoding="utf-8"))
    with Path(original_samples).open(newline="",encoding="utf-8-sig") as f:
        samples=list(csv.DictReader(f))
    if (len(samples)!=N_TIPS or
        len({r["phylogenetic_tip"] for r in samples})!=N_TIPS or
        summary.get("source_111_independent_public_SRA_BioSamples")!=N_TIPS or
        not summary.get("preliminary_GTRG4_IQTREE_permitted") or
        summary.get("single_copy_and_homeolog_orthology_verified") or
        summary.get("global_EAzami294_species_tree_accepted") or
        summary.get("adaptive_radiation_or_rain_transition_supported")):
        raise ValueError("Source real111 Q20 same-position input is not accepted for technical ML")
    t=Phylo.read(str(tree_file),"newick")
    tips=[z.name for z in t.get_terminals()]
    source={r["phylogenetic_tip"] for r in samples}
    if len(tips)!=N_TIPS or len(set(tips))!=N_TIPS or set(tips)!=source:
        raise ValueError("Actual IQ-TREE ML tree tip IDs not identical to original111 sources")
    if any(z.branch_length is not None and z.branch_length<0 for z in t.find_clades()):
        raise ValueError("Negative ML branch length")
    inferred=sum(z.branch_length is not None for z in t.find_clades())
    if inferred<N_TIPS-1:
        raise ValueError("Too few actual substitution-model branch lengths")
    result={
        "schema":"azami_T2_111source_118original_realQ20_unrooted_GTRG4_ML_technical_v1",
        "status":"EMPIRICAL_REAL_SOURCE111_NUCLEAR_TECHNICAL_ML_NOT_ACCEPTED_SPECIES_TREE",
        "source_original_photo_WCVP_keys":N_ORIGINAL,
        "source_original_S1_named_units":N_SOURCE_NAMES,
        "real_original_SRA_BioSamples":N_TIPS,
        "same_original_reference_nuclear_loci_input":
            summary["all111_selected_original_reference_locus_partitions"],
        "observed_Q20_reference_coordinate_positions_input":
            summary["all111_source_Q20_80pct_reference_positions"],
        "real_source_GTRG4_tree_tips":len(tips),
        "estimated_original_reference_branch_length_edges":inferred,
        "IQTREE_tree_sha256":hashlib.sha256(Path(tree_file).read_bytes()).hexdigest(),
        "branch_lengths_are_substitutions_per_site_not_years":True,
        "actual_public_FASTQ_and_original_reference_Q20_site_calls":True,
        "missing_source_bases_are_N_not_reference_imputed":True,
        "nuclear_copy_identity_and_paralog_homeolog_review_completed":False,
        "author350_locus_reproduction_completed":False,
        "EAzami294_independent_species_tree_completed":False,
        "genus_global_species_tree_accepted":False,
        "true_botanical_nodding_Mk_events_inferred":False,
        "BIO12_BIO15_ancient_rain_adaptation_inferred":False,
        "require_gene_block_topology_support_and_copy_validation":True,
    }
    out=Path(output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False,sort_keys=True)+"\n",encoding="utf-8")
    print("T2_REAL111_PUBLIC_Q20_IQTREE_TECHNICAL_RESULT="+json.dumps(result,sort_keys=True),flush=True)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input-manifest",type=Path,required=True)
    p.add_argument("--samples",type=Path,required=True)
    p.add_argument("--tree",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    verify(a.input_manifest,a.samples,a.tree,a.output)


if __name__=="__main__":
    main()
