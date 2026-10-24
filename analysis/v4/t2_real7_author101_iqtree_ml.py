"""Actual seven publicly sourced SRA Q20 observed sites -> technical IQ-TREE ML.

Independent from existing NJ p-distance method, but NOT independent reads:
both methods use the same 7 original Moreyra public run 120k-pair prefixes.
101 loci with >=100 fully observed, coordinate-orthologous *candidate*
sites is the primary observed-only matrix (28726 sites). All 110 loci
(29270 observed sites) is a sensitivity. No missing positions filled with
ancestral/reference nucleotide; no A/C/G/T ambiguity fabricated. These are
spliced observed sites, NOT intact recovered loci, paralog-validated genes
or an accepted Cirsium species phylogeny.

Two prespecified IQ-TREE 2 ML runs (GTR+G4, 1000 UFBoot/SH-aLRT)
are compared to NJ and the independently published Siniscalchi NA63
nuclear quartet on exact 4 overlap taxon labels.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from hashlib import sha256
from io import StringIO
import json
from pathlib import Path
import re

from Bio import Phylo

from analysis.v4 import t2_seven_original110_Q20_nj_diagnostic as nj
from analysis.v4 import t2_real7_NA63_quartet_crosscheck as independent
from analysis.v4 import t2_seven_real120k_shared_sites as source

EXPECTED_7=7
EXPECTED_101=101
EXPECTED_110=110
EXPECTED_101_SITES=28726
EXPECTED_110_SITES=29270
EXPECTED_110_VARIABLE=542
EXPECTED_110_INFORMATIVE=86
MODEL="GTR+G4"
REPS=1000
SEED=20261023


def read_csv(path):
    with Path(path).open(newline="",encoding="utf-8-sig") as f:
        data=list(csv.DictReader(f))
    if not data:
        raise ValueError("Missing source-authenticated SRA specimen or locus registry")
    return data


def build_matrices(original125,chosen,refs,receipts):
    if len(original125)!=source.N_ORIGINAL_RUNS or len(chosen)!=EXPECTED_7:
        raise ValueError("Original125+7 EAzami SRA biological-sample scope changed")
    if len(refs)!=EXPECTED_110 or len(set(refs.values()))!=EXPECTED_110:
        raise ValueError("Original source 110 unique target reference copies lost")
    records=nj.load_original_receipts(receipts,original125,chosen,refs)
    names,loci=nj.loci_from_original_base_calls(records,refs)
    primary=[x for x in loci if x["eligible_primary"]]
    values={
        "n_public_source_samples":len(names),
        "n_public_original_candidate_targets":len(loci),
        "n_candidate_loci_primary":len(primary),
        "n_original_all_candidate_sites":sum(z["aligned_positions"] for z in loci),
        "n_primary_observed_only_sites":sum(z["aligned_positions"] for z in primary),
        "n_source_variable_all110":sum(z["variable_positions"] for z in loci),
        "n_source_parsimony_informative_all110":
            sum(z["parsimony_informative_positions"] for z in loci),
    }
    wanted=(EXPECTED_7,EXPECTED_110,EXPECTED_101,EXPECTED_110_SITES,
            EXPECTED_101_SITES,EXPECTED_110_VARIABLE,EXPECTED_110_INFORMATIVE)
    got=tuple(values.values())
    if got!=wanted:
        raise ValueError(f"Actual verified seven-sample Q20 site totals drifted: {got} != {wanted}")
    expected_ids={x["run"] for x in records}
    if expected_ids!=source.SEVEN or len({x["biosample"] for x in records})!=EXPECTED_7:
        raise ValueError("Known 7 true BioSample sequence/accession matching changed")
    return names,records,loci,primary,values


def alignment_by_samples(names,loci):
    if not names or len(set(names))!=len(names):
        raise ValueError("Duplicate or empty original species labels")
    full={}
    for i,name in enumerate(names):
        if not re.fullmatch(r"Cirsium_[A-Za-z0-9_]+",name):
            raise ValueError("Unverified author source taxon labels must not enter ML")
        seq="".join(x["alignment"][i] for x in loci)
        if not seq or set(seq)-set("ACGT"):
            raise ValueError("Unknown, unobserved or reference-imputed DNA in complete-case site matrix")
        full[name]=seq
    if len(set(len(x) for x in full.values()))!=1:
        raise ValueError("Not identical observed coordinate sites in all seven biological samples")
    return full


def fasta_payload(sequences):
    if not sequences:
        raise ValueError("No original SRA observed sequences")
    n={len(x) for x in sequences.values()}
    if len(n)!=1:
        raise ValueError("Different alignment lengths")
    return "".join(
        f">{taxon}\n"+"\n".join(seq[i:i+80] for i in range(0,len(seq),80))+"\n"
        for taxon,seq in sorted(sequences.items())
    )


def prepare(source125,chosen,targets,receipts,dest):
    original=read_csv(source125)
    chosen_rows=read_csv(chosen)
    references=source.targets(read_csv(targets))
    names,samples,all_loci,primary,counts=build_matrices(
        original,chosen_rows,references,receipts)
    dest=Path(dest)
    dest.mkdir(parents=True,exist_ok=True)
    artifacts={}
    for label,seq_loci,n_expected in (
        ("primary101",primary,EXPECTED_101_SITES),
        ("all110_sensitivity",all_loci,EXPECTED_110_SITES),
    ):
        alignment=alignment_by_samples(names,seq_loci)
        if len(next(iter(alignment.values())))!=n_expected:
            raise ValueError("No fresh site evidence matching original source totals")
        payload=fasta_payload(alignment)
        name=dest/(label+".fasta")
        name.write_text(payload,encoding="ascii")
        artifacts[label]={
            "input_path":name.name,
            "alignment_fasta_sha256":sha256(payload.encode("ascii")).hexdigest(),
            "n_observed_complete_case_sites":n_expected,
            "n_original_locus_candidates":len(seq_loci),
            "n_source_real_nuclear_samples":len(names),
            "unobserved_or_ambiguous_bases_added":0,
            "sequences_are_discontiguous_observed_Q20_positions":True,
            "loci_verified_orthologous_single_copy":False,
        }
    summary={
      "schema":"azami_v4_seven_real_Q20_101_110_ML_input_contract_v1",
      "source":"original Moreyra source 120k SRA run pairs and EAzami original 241",
      "true_source_sample_taxa":names,
      "true_source_run_ids":[x["run"] for x in samples],
      "true_source_biosamples":[x["biosample"] for x in samples],
      **counts,
      "candidate_alignments":artifacts,
      "predeclared_ML_model":MODEL,
      "predeclared_UFBOOT_and_SHaLRT_replicates":REPS,
      "predeclared_random_seed":SEED,
      "reference_ancestral_site_imputation_applied":False,
      "source_nucleotides_real_not_independently_orthology_validated":True,
      "full_294_EAzami_species_tree_completed":False,
      "seven_sample_output_is_valid_species_phylogeny":False,
    }
    (dest/"REAL7_SOURCE_Q20_ML_INPUT_MANIFEST.json").write_text(
        json.dumps(summary,indent=2,ensure_ascii=False,sort_keys=True)+"\n",encoding="utf-8")
    return summary


def splits_with_support(tree,expected):
    names={z.name for z in tree.get_terminals()}
    if len(names)!=expected or len(tree.get_terminals())!=expected:
        raise ValueError("Maximum-likelihood tree omitted a real SRA taxon")
    result={}
    for clade in tree.get_nonterminals():
        branch={x.name for x in clade.get_terminals()}
        other=names-branch
        if min(len(branch),len(other))<2:
            continue
        group=tuple(sorted(branch))
        alt=tuple(sorted(other))
        key=min(group,alt) if len(group)==len(alt) else group if len(group)<len(alt) else alt
        bs=clade.confidence
        support=float(bs) if bs is not None else None
        if key in result and result[key]!=support and support is not None:
            raise ValueError("Non-consistent support on the same unrooted ML split")
        if key not in result or result[key] is None:
            result[key]=support
    if len(result)!=expected-3:
        raise ValueError("Expected fully binary unrooted ML tree internal splits")
    return result


def quartet_split_four(tree):
    return independent.quartet_split(tree)


def inspect_tree(path,expected_names,author_tree,technical_nj):
    tree=Phylo.read(str(path),"newick")
    names={x.name for x in tree.get_terminals()}
    if names!=set(expected_names):
        raise ValueError("ML sample names diverge from exact original SRA seven names")
    if any((z.branch_length is not None and z.branch_length<0)
            for z in tree.find_clades()):
        raise ValueError("Negative ML inferred substitutions-per-site branch length")
    splits=splits_with_support(tree,EXPECTED_7)
    nj_splits=nj.splits_unrooted(technical_nj)
    if len(nj_splits)!=EXPECTED_7-3:
        raise ValueError("Original source-authentic p-distance NJ split set invalid")
    if len(splits)!=len(nj_splits):
        raise ValueError("Unresolved ML/NJ comparison for same seven SRA samples")
    actual=quartet_split_four(tree)
    published=quartet_split_four(author_tree)
    if published!=["Cirsium_flodmanii","Cirsium_tracyi"]:
        raise ValueError("Source-authenticated NA63 published quartet identity changed")
    return {
       "ML_tree_path":str(path),
       "n_real_sra_tips":len(names),
       "n_unrooted_binary_splits":len(splits),
       "ML_matches_published_independent_NA63_4tip_quartet":actual==published,
       "observed_real_ML_four_tip_quartet_side":actual,
       "published_NA63_four_tip_quartet_side":published,
       "n_ML_splits_same_as_uncorrected_NJ":len(set(splits)&nj_splits),
       "ML_all_internal_unrooted_splits_equal_NJ":set(splits)==nj_splits,
       "splits":[{
            "unrooted_canonical_side":"|".join(side),
            "IQTREE_branch_support_number_if_parsed":support,
       } for side,support in sorted(splits.items())],
       "not_independent_ortholog_family_gene_trees":True,
       "not_species_tree_nor_rooted_nor_dated":True,
    }


def verify(ml_manifest,nj_report,nj_tree,primary_tree,all_tree,output):
    manifest=json.loads(Path(ml_manifest).read_text(encoding="utf-8"))
    old=json.loads(Path(nj_report).read_text(encoding="utf-8"))
    if (manifest.get("n_candidate_loci_primary")!=EXPECTED_101 or
        manifest.get("n_primary_observed_only_sites")!=EXPECTED_101_SITES or
        manifest.get("n_original_all_candidate_sites")!=EXPECTED_110_SITES or
        manifest.get("n_source_parsimony_informative_all110")!=EXPECTED_110_INFORMATIVE
        or old.get("n_original_public_SRA_biosamples")!=EXPECTED_7
        or old.get("source_shared_Q20_positions_all_targets")!=EXPECTED_110_SITES
        or old.get("n_complete_case_invariant_plus_variable_sites_primary")!=EXPECTED_101_SITES
        or not old.get("nuclear_tree_preliminary_only")
        or old.get("EAzami_294_tip_species_tree_constructed")):
        raise ValueError("Original real SRA seven Q20 source provenance absent")
    original=independent.load_na63()
    nj_phy=Phylo.read(str(nj_tree),"newick")
    out={
      "schema":"azami_v4_7source_101locus_observedQ20_iqtree2_maximum_likelihood_technical_v1",
      "status":"REAL_Q20_PUBLIC_SRA_SEVEN_SOURCE_ML_TECHNICAL_NOT_ACCEPTED_PHYLOGENY",
      "source_real_run_ids":manifest["true_source_run_ids"],
      "source_real_biosamples":manifest["true_source_biosamples"],
      "source_real_taxon_labels":manifest["true_source_sample_taxa"],
      "primary_loci":EXPECTED_101,
      "primary_observed_only_Q20_positions":EXPECTED_101_SITES,
      "sensitivity_loci":EXPECTED_110,
      "sensitivity_observed_only_Q20_positions":EXPECTED_110_SITES,
      "all110_original_variable_sites":EXPECTED_110_VARIABLE,
      "all110_original_parsimony_informative_sites":EXPECTED_110_INFORMATIVE,
      "predeclared_model":MODEL,
      "IQTREE_UFBOOT_replicates":REPS,
      "IQTREE_SHaLRT_replicates":REPS,
      "results":{
           "primary_101":inspect_tree(primary_tree,manifest["true_source_sample_taxa"],original,nj_phy),
           "sensitivity_110":inspect_tree(all_tree,manifest["true_source_sample_taxa"],original,nj_phy),
       },
      "intact_orthologous_sequences_or_copy_number_verified":False,
      "full_Cirsium_Moreyra_350_original_phylogeny_recovered":False,
      "EAzami_294sample_global_species_tree_reconstructed":False,
      "evolutionary_ou_bm_or_rain_nodding_selection_supported_by_this_test":False,
      "raw_SRA_fastq_or_BAM_distributed":False,
      "large_gene_region_alignment_in_final_artifact":False,
    }
    output=Path(output);output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(out,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    print("REAL7_GTRG4_MAXIMUM_LIKELIHOOD_NJ_NA63=",json.dumps({
       "primary_quartet_agreement":out["results"]["primary_101"]["ML_matches_published_independent_NA63_4tip_quartet"],
       "sensitivity_quartet_agreement":out["results"]["sensitivity_110"]["ML_matches_published_independent_NA63_4tip_quartet"],
       "primary_splits_match_NJ":out["results"]["primary_101"]["n_ML_splits_same_as_uncorrected_NJ"],
       "sensitivity_splits_match_NJ":out["results"]["sensitivity_110"]["n_ML_splits_same_as_uncorrected_NJ"],
    },sort_keys=True),flush=True)
    return out


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mode",choices=("prepare","verify"),required=True)
    ap.add_argument("--sra-source",type=Path)
    ap.add_argument("--chosen",type=Path)
    ap.add_argument("--targets",type=Path)
    ap.add_argument("--receipts",type=Path)
    ap.add_argument("--out-dir",type=Path)
    ap.add_argument("--ml-manifest",type=Path)
    ap.add_argument("--nj-report",type=Path)
    ap.add_argument("--nj-tree",type=Path)
    ap.add_argument("--ml-primary",type=Path)
    ap.add_argument("--ml-all",type=Path)
    ap.add_argument("--output",type=Path)
    a=ap.parse_args()
    if a.mode=="prepare":
        if not all([a.sra_source,a.chosen,a.targets,a.receipts,a.out_dir]):
            ap.error("prepare requires original SRA sample and seven original Q20 receipt inputs")
        result=prepare(a.sra_source,a.chosen,a.targets,a.receipts,a.out_dir)
        print("REAL7_GTRG4_ML_INPUT_PREPARED=",json.dumps({
           "samples":result["n_public_source_samples"],
           "primary_loci":result["n_candidate_loci_primary"],
           "primary_sites":result["n_primary_observed_only_sites"],
           "sensitivity_sites":result["n_original_all_candidate_sites"],
        },sort_keys=True),flush=True)
    else:
        if not all([a.ml_manifest,a.nj_report,a.nj_tree,a.ml_primary,a.ml_all,a.output]):
            ap.error("verify requires original NJ and two actual IQ-TREE ML Newick inputs")
        verify(a.ml_manifest,a.nj_report,a.nj_tree,a.ml_primary,a.ml_all,a.output)


if __name__=="__main__":
    main()
