"""Actual seven-sample sequence-site p-distance and NJ *technical* nuclear pilot.

Reuses only author-source authenticated seven public SRA 120k read-prefix Q20
majority-base receipts from Actions 39868141961. Complete-case same-reference
nucleotide coordinates across all seven, including INVARIANT sites, become
a technical p-distance screening alignment on the original 110 shared-copy
targets. Distances and locus-bootstrap unrooted neighbor-joining splits are
reported. This is NOT an orthology-validated, rooted or dated Cirsium species
tree, an OU/BM fit, or evidence for repeated rain-driven nodding evolution.

No FASTQ, BAM or unvalidated full-locus consensus is created or reuploaded.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from io import StringIO
import json
from pathlib import Path
import random

from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

from analysis.v4 import t2_seven_real120k_shared_sites as upstream

N_SAMPLES=7
N_SHARED_COPY_TARGETS=110
N_SOURCE_PAIRED_READS=120000
EXPECTED_ALL_SHARED_POSITIONS=29270
EXPECTED_ALL_VARIABLE_POSITIONS=542
EXPECTED_ALL_INFORMATIVE_POSITIONS=86
EXPECTED_TARGETS_GE100=101
MIN_COMMON_SITE_PER_TARGET=100
MIN_PRIMARY_LOCI=30
BOOTSTRAP_REPLICATES=300
BOOTSTRAP_SEED=20261023


def read_csv(path):
    with Path(path).open(newline="",encoding="utf-8-sig") as f:
        rows=list(csv.DictReader(f))
    if not rows:
        raise ValueError("Empty source-verified original target or sample table")
    return rows


def load_original_receipts(root,original125,chosen,targets):
    paths=sorted(Path(root).rglob("SRR*.json"))
    if len(paths)!=N_SAMPLES or {p.stem for p in paths}!=upstream.SEVEN:
        raise ValueError("Missing or duplicated seven real source authenticated nucleotide receipts")
    if len(original125)!=upstream.N_ORIGINAL_RUNS or len(chosen)!=N_SAMPLES:
        raise ValueError("Original 125 SRA source and seven-sample list not preserved")
    rows=[]
    bysource={x["run"]:x for x in original125}
    for file in paths:
        record=json.loads(file.read_text(encoding="utf-8"))
        original=upstream.source_sample(original125,chosen,record["run"])
        upstream.validate_single(record,original,targets)
        if (record["n_source_pairs"]!=N_SOURCE_PAIRED_READS or
                record.get("status")!="real_source120k_Q20_majority_nucleotide_position_pilot" or
                len({r["accepted_key"] for r in rows+[record]})!=len(rows)+1):
            raise ValueError("Source actual run/single accepted key or base-call identity drift")
        rows.append(record)
    if len(rows)!=N_SAMPLES or len({x["biosample"] for x in rows})!=N_SAMPLES:
        raise ValueError("Biological SRA sample mistakenly duplicated")
    return sorted(rows,key=lambda r:r["run"])


def site_count(bases):
    counts=Counter(bases)
    variable=int(len(counts)>1)
    informative=int(sum(v>=2 for v in counts.values())>=2)
    return variable,informative


def loci_from_original_base_calls(records,refs):
    """Do not replace unobserved bases with reference nor pool loci across copies."""
    if len(records)!=N_SAMPLES or len(refs)!=N_SHARED_COPY_TARGETS:
        raise ValueError("Not seven specimens and exactly 110 source copy-matched loci")
    names=[r["S1_source_taxon"].replace(" ","_") for r in records]
    if len(set(names))!=N_SAMPLES:
        raise ValueError("Duplicate source taxon label, cannot infer a seven-taxon contrast")
    results=[]
    for locus in sorted(refs):
        calls=[]
        for record in records:
            entry=record["targets"][locus]
            if entry["reference"]!=refs[locus]:
                raise ValueError("Different Asteraceae reference copy selected across specimens")
            bypos={}
            for pos,base,dp,fraction in entry["q20_majority_bases_1based"]:
                if pos in bypos or base not in "ACGT" or dp<3 or not(.9<=fraction<=1):
                    raise ValueError("Malformed original public read Q20 majority nucleotide call")
                bypos[int(pos)]=base
            calls.append(bypos)
        both=set(calls[0])
        for c in calls[1:]:
            both.intersection_update(c)
        common=sorted(both)
        alignment=["".join(call[pos] for pos in common) for call in calls]
        if len({len(x) for x in alignment})!=1:
            raise ValueError("Variable per-source locus alignment lengths")
        variable=informative=0
        for pos in common:
            a,b=site_count([call[pos] for call in calls])
            variable+=a;informative+=b
        results.append({
            "locus":locus,"reference":refs[locus],
            "aligned_positions":len(common),
            "variable_positions":variable,
            "parsimony_informative_positions":informative,
            "alignment":alignment,
            "eligible_primary":len(common)>=MIN_COMMON_SITE_PER_TARGET,
        })
    if len(results)!=N_SHARED_COPY_TARGETS:
        raise ValueError("Source 110 locus loss")
    return names,results


def mismatch_matrix(data,n):
    """Return 7x7 integer mismatches across fully shared, including invariant, sites."""
    out=[[0]*n for _ in range(n)]
    length=0
    for locus in data:
        aligned=locus["alignment"]
        if len(aligned)!=n or len(set(len(x) for x in aligned))!=1:
            raise ValueError("Incomplete original shared-position alignment")
        L=len(aligned[0]);length+=L
        for i in range(n):
            for j in range(i):
                m=sum(a!=b for a,b in zip(aligned[i],aligned[j]))
                out[i][j]+=m;out[j][i]+=m
    if not length:
        raise ValueError("No observed shared Q20 nucleotide positions")
    return out,length


def distance_nj(names,selected_loci):
    mismatches,n=mismatch_matrix(selected_loci,len(names))
    matrix=[[mismatches[i][j]/n for j in range(i+1)] for i in range(len(names))]
    if any(matrix[i][i]!=0 for i in range(len(names))):
        raise ValueError("Malformed self-distance")
    tree=DistanceTreeConstructor().nj(DistanceMatrix(names,matrix))
    observed={t.name for t in tree.get_terminals()}
    if observed!=set(names) or len(observed)!=len(names):
        raise ValueError("NJ lost a real SRA source taxon")
    return tree,mismatches,n


def splits_unrooted(tree):
    """Undirected splits (two sides, neither leaf/trivial); ignore rooting artifact."""
    names={t.name for t in tree.get_terminals()}
    if len(names)!=N_SAMPLES:
        raise ValueError("Wrong original sample count in pilot tree")
    out=set()
    for clade in tree.get_nonterminals():
        left={t.name for t in clade.get_terminals()}
        right=names-left
        if min(len(left),len(right))<2:
            continue
        side=left if (len(left)<len(right) or
                      len(left)==len(right) and tuple(sorted(left))<tuple(sorted(right))) else right
        out.add(tuple(sorted(side)))
    return out


def infer(records,refs,boot=BOOTSTRAP_REPLICATES,seed=BOOTSTRAP_SEED):
    names,loci=loci_from_original_base_calls(records,refs)
    all_pos=sum(x["aligned_positions"] for x in loci)
    all_var=sum(x["variable_positions"] for x in loci)
    all_inf=sum(x["parsimony_informative_positions"] for x in loci)
    n_ge100=sum(z["eligible_primary"] for z in loci)
    if (all_pos,all_var,all_inf,n_ge100)!=(EXPECTED_ALL_SHARED_POSITIONS,
                                             EXPECTED_ALL_VARIABLE_POSITIONS,
                                             EXPECTED_ALL_INFORMATIVE_POSITIONS,
                                             EXPECTED_TARGETS_GE100):
        raise ValueError(
            "Real same-position seven SRA source Q20 29270/542/86/101 receipt drift")
    primary=[x for x in loci if x["eligible_primary"]]
    if len(primary)<MIN_PRIMARY_LOCI or boot<100:
        raise ValueError("Pilot loci or locus bootstrap count insufficient")
    # A tree on complete-case observed positions is a technical distance
    # sensitivity only, NOT a phylogenomic acceptance.
    tree,mi,sites=distance_nj(names,primary)
    secondary,_,all_sites=distance_nj(names,loci)
    if all_sites!=EXPECTED_ALL_SHARED_POSITIONS:
        raise ValueError("Site count across all original 110 targets changed")
    splits=splits_unrooted(tree)
    rng=random.Random(seed)
    supports=Counter()
    for _ in range(boot):
        draw=[primary[rng.randrange(len(primary))] for _ in range(len(primary))]
        t,_,_=distance_nj(names,draw)
        supports.update(splits_unrooted(t))
    pairs=[]
    for i in range(len(names)):
        for j in range(i):
            pairs.append({
                "source_taxon_A":names[j],
                "source_taxon_B":names[i],
                "shared_site_count":sites,
                "n_mismatching_observed_nucleotides":mi[i][j],
                "p_distance_on_seven_way_complete_case_sites":mi[i][j]/sites,
                "uncorrected_distance_not_a_dated_molecular_clock":True,
            })
    split_rows=[{"unrooted_split_side":"|".join(side),
                 "n_bootstrap_locus_resamples_with_this_split":supports[side],
                 "n_locus_bootstrap_resamples":boot,
                 "technical_NJ_locus_resampling_fraction":supports[side]/boot,
                 "not_validated_species_tree_support":True}
                for side in sorted(splits)]
    stream=StringIO(); Phylo.write(tree,stream,"newick")
    result={
        "schema":"azami_v4_seven_source120k_original110_Q20_complete_case_NJ_diagnostic_v1",
        "status":"EMPIRICAL_REAL_SRA_SEQUENCE_DISTANCE_AND_UNROOTED_NJ_TECHNICAL_NOT_A_VALIDATED_PHYLOGENY",
        "n_original_public_SRA_biosamples":N_SAMPLES,
        "n_original_same_copy_candidate_loci":N_SHARED_COPY_TARGETS,
        "n_primary_loci_shared_atleast_100_Q20_sites":len(primary),
        "n_total_loci_with_shared_calls":len(loci),
        "source_shared_Q20_positions_all_targets":all_pos,
        "source_shared_variable_Q20_positions":all_var,
        "source_shared_parsimony_informative_Q20_positions":all_inf,
        "n_complete_case_invariant_plus_variable_sites_primary":sites,
        "primary_p_distance_pair_count":len(pairs),
        "locus_bootstrap_draws":boot,
        "locus_bootstrap_seed":seed,
        "run_first_prefix_paired_reads_per_source":N_SOURCE_PAIRED_READS,
        "method":"actual 120k FASTQ prefixes, 3x/Q20/90%-majority same-reference bases; seven-way full observed position intersection; >=100 site candidate loci; uncorrected p-distance; Bio.Phylo NJ; 300 locus resamples",
        "nuclear_tree_preliminary_only":True,
        "original_Moreyra_350_locus_tree_recovered":False,
        "EAzami_294_tip_species_tree_constructed":False,
        "copy_state_orthology_homeology_and_taxonomy_concepts_verified":False,
        "no_adaptive_evolution_or_orientation_transition_result":True,
        "no_reused_technical_NJ_as_OU_BM_species_tree":True,
        "fasta_fastq_bam_in_generated_outputs":False,
        "primary_unrooted_newick":stream.getvalue().strip(),
    }
    return result,loci,pairs,split_rows


def write_csv(path,rows,fields=None):
    with path.open("w",newline="",encoding="utf-8") as h:
        writer=csv.DictWriter(h,fieldnames=fields or list(rows[0]),extrasaction="raise")
        writer.writeheader();writer.writerows(rows)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-receipts",type=Path,required=True)
    ap.add_argument("--sra-source",type=Path,required=True)
    ap.add_argument("--chosen",type=Path,required=True)
    ap.add_argument("--targets",type=Path,required=True)
    ap.add_argument("--out-dir",type=Path,required=True)
    args=ap.parse_args()
    target_refs=upstream.targets(read_csv(args.targets))
    samples=load_original_receipts(args.run_receipts,
                 read_csv(args.sra_source),read_csv(args.chosen),target_refs)
    report,loci,pairs,splits=infer(samples,target_refs)
    args.out_dir.mkdir(parents=True,exist_ok=True)
    report["source_SRA_runs"]=[x["run"] for x in samples]
    report["source_SRA_BioSamples"]=[x["biosample"] for x in samples]
    write_csv(args.out_dir/"SEVEN_REAL_Q20_101_LOCUS_PAIRWISE_PDIST.csv",pairs)
    write_csv(args.out_dir/"SEVEN_REAL_Q20_101_LOCUS_NJ_BOOTSTRAP_SPLITS.csv",
              splits,fields=["unrooted_split_side",
                            "n_bootstrap_locus_resamples_with_this_split",
                            "n_locus_bootstrap_resamples",
                            "technical_NJ_locus_resampling_fraction",
                            "not_validated_species_tree_support"])
    write_csv(args.out_dir/"SEVEN_REAL_Q20_110_ORIGINAL_LOCUS_SITE_COUNTS.csv",
              [{k:v for k,v in z.items() if k!="alignment"} for z in loci])
    (args.out_dir/"SEVEN_REAL_Q20_101_LOCUS_NJ_TECHNICAL.nwk").write_text(
          report["primary_unrooted_newick"]+"\n",encoding="utf-8")
    (args.out_dir/"SEVEN_REAL_Q20_101_LOCUS_NJ_DIAGNOSTIC.json").write_text(
          json.dumps(report,indent=2,sort_keys=True,ensure_ascii=False,allow_nan=False)+"\n",
          encoding="utf-8")
    print("T2_SEVEN_REAL_NUCLEAR_SITE_NJ_TECHNICAL="+json.dumps({
        "source_real_SRA_taxa":report["n_original_public_SRA_biosamples"],
        "original_loci":report["n_original_same_copy_candidate_loci"],
        "primary_loci":report["n_primary_loci_shared_atleast_100_Q20_sites"],
        "shared_Q20_sites":report["source_shared_Q20_positions_all_targets"],
        "primary_Q20_positions":report["n_complete_case_invariant_plus_variable_sites_primary"],
        "bootstrap_splits":splits,
    },ensure_ascii=False,sort_keys=True),flush=True)


if __name__=="__main__":
    main()
