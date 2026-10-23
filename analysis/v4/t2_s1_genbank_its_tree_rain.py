"""Build OWN exploratory S1 Cirsium nuclear-ITS tree and test BIO12 by tree tips.

Strictly derived from the 251 actual original S1 source names and GenBank
accession-verified ITS entries, NOT a borrowed NA63/EastAsia20 tree and NOT a
52-tree generic taxonomic graft.

1. Curate MAFFT alignment by measured occupancy and informative positions;
2. IQ-TREE2 is run separately in GitHub Actions on this exact accepted alignment;
3. audit ML tree + all accessions + S1 tip-name WCVP evidence;
4. if >=3 independently botanically nodding and >=5 upright matched tips, compare
   S1 present-day CHELSA BIO12 after temperate covariate under exact same ITS
   tree topology with independent model vs Brownian/Pagel-lambda PGLS.

One nuclear ribosomal ITS copy selected per species is NOT a single-copy
ortholog and may show concerted evolution, paralogy and hybridization. Without
an authenticated OUTGROUP and multiple independent nuclear loci, the exploratory
ML tree's root is arbitrary; any midpoint root in PGLS is explicitly a sensitivity
assumption, NOT an ancestral-state identification or established adaptive gain.
"""
from __future__ import annotations
import argparse
from collections import Counter
from copy import deepcopy
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from Bio import SeqIO, Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from analysis.v4.t2_na63_actual_tree_rain_gls import (
    brownian_shared_ancestry, fit_models, equal_variance_independent_OLS,
)

REQUIRED_SOURCE_UNITS=251
REQUIRED_WCVP=244
ACGT=frozenset("ACGT")
MIN_TAXA_FOR_ITS_ML=12
MIN_ALIGNED_COLUMNS=140
MIN_VARIANT_COLUMNS=10
MIN_PHYLO_POSTURE_STATES=(3,5)
MAX_ALLOWED_GAP_FRAC=.40
MAX_MISSING_ITS_TIP_FRAC=.52


def read_csv(path):
    with Path(path).open(encoding="utf-8-sig",newline="") as h:
        return list(csv.DictReader(h))


def audit_accessions(accessions, original_receipt, msa):
    receipt=json.loads(Path(original_receipt).read_text(encoding="utf-8"))
    if (receipt["source_original_S1_master_names"]!=REQUIRED_SOURCE_UNITS or
        receipt["source_original_WCVP_accepted_keys"]!=REQUIRED_WCVP or
        receipt["scope"]!="S1_ALL251" or
        not receipt["all251_queries_completed_without_network_failures"]):
        raise ValueError("Cannot construct full-S1-source GenBank tree from pilot or failed requests")
    if not accessions:
        return None,[],{"status":"BLOCKED_NO_RELIABLE_PUBLIC_GENBANK_ITS_SEQUENCES"}
    if len(set(x["original_source_name"] for x in accessions))!=len(accessions):
        raise ValueError("Original S1 taxon scored with >1 sequence or duplicate name")
    if any(x["original_source_name"]!=x["source_original_species"] or
           x["S1_species_to_GB_record_exact_name"]!="True" for x in accessions):
        raise ValueError("GenBank accession biological source != exact S1 original binomial")
    byname={x["original_source_name"].replace(" ","_"):x for x in accessions}
    seqs=list(SeqIO.parse(str(msa),"fasta"))
    mapped={}
    reverse=[]
    for record in seqs:
        name=record.id
        if name.startswith("_R_"):
            name=name[3:]
            reverse.append(name)
        if name not in byname or name in mapped:
            raise ValueError("MAFFT reverse orientation or source accession label drift")
        mapped[name]=str(record.seq).upper()
    if len(mapped)!=len(byname) or set(mapped)!=set(byname):
        raise ValueError("MAFFT alignment lost authentic S1 GenBank source taxa")
    lengths={len(x) for x in mapped.values()}
    if len(lengths)!=1:
        raise ValueError("MAFFT output has unequal aligned sequence lengths")
    ncols=lengths.pop()
    if ncols>20000 or ncols<MIN_ALIGNED_COLUMNS:
        raise ValueError("Malformed or missing GenBank ITS sequence alignment")
    return mapped,reverse,{
        "source_S1_unique_GenBank_ITS_taxa":len(byname),
        "S1_original_names_scanned":REQUIRED_SOURCE_UNITS,
        "original_WCVP_taxon_keys":REQUIRED_WCVP,
        "MAFFT_full_alignment_columns":ncols,
        "MAFFT_sequences_reverse_complemented":sorted(reverse),
        "source_GenBank_accessions_SHA256":hashlib.sha256(Path(original_receipt).read_bytes()).hexdigest(),
        "MAFFT_alignment_sha256":hashlib.sha256(Path(msa).read_bytes()).hexdigest(),
    }


def trim_alignment(mapped, min_occupancy=.60):
    if mapped is None or len(mapped)<MIN_TAXA_FOR_ITS_ML:
        return {},{"status":"BLOCKED_TOO_FEW_EXACT_S1_GENBANK_ITS_TAXA"}
    if not .50<=min_occupancy<=.85:
        raise ValueError("Unapproved ITS alignment occupancy")
    names=sorted(mapped)
    raw=list(mapped.values())
    n=len(names);columns=len(raw[0])
    kept=[]
    biallelic=variable=parsimony=0
    for i in range(columns):
        seq=[mapped[name][i] for name in names]
        calls=[b for b in seq if b in ACGT]
        if len(calls)/n<min_occupancy:
            continue
        kept.append(i)
        c=Counter(calls)
        variable+=len(c)>1
        parsimony+=sum(z>=2 for z in c.values())>=2
    if len(kept)<MIN_ALIGNED_COLUMNS or variable<MIN_VARIANT_COLUMNS:
        return {},{
            "status":"BLOCKED_ALIGNMENT_WITH_TOO_LITTLE_REAL_ITS_INFORMATION",
            "retained_columns":len(kept),"variable_site_columns":variable,
            "parsimony_informative_columns":parsimony,
        }
    seqs={}
    for name in names:
        cleaned="".join(mapped[name][i] if mapped[name][i] in ACGT else "N"
                        for i in kept)
        seqs[name]=cleaned
    counts={name:len(seq)-seq.count("N") for name,seq in seqs.items()}
    if min(counts.values())/len(kept) < 1-MAX_MISSING_ITS_TIP_FRAC:
        return {},{
            "status":"BLOCKED_ONE_OR_MORE_SELECTED_ITS_SPECIES_HAS_TOO_LITTLE_OVERLAP",
            "min_sample_real_overlap_fraction":min(counts.values())/len(kept),
            "original_species_with_lowest_bases":min(counts,key=counts.get),
        }
    return seqs,{
        "status":"PASS_TO_GENBANK_ITS_IQTREE_EXPLORATORY_ML",
        "selected_real_ITS_alignment_columns":len(kept),
        "selected_variable_nucleotide_sites":variable,
        "selected_parsimony_informative_sites":parsimony,
        "n_genbank_exact_species_tips":n,
        "min_source_real_NA_T_C_G_sites_per_tip":min(counts.values()),
        "median_real_sites_per_tip":int(np.median(list(counts.values()))),
        "source_missing_bases_become_N_not_GenBank_reference_alleles":True,
        "DNA_alignment_based_on_nuclear_ribosomal_ITS_only":True,
        "single_copy_nuclear_orthology_validated":False,
        "true_species_tree_independently_completed":False,
    }


def prepare(accessions_file,receipt_file,msa_file,outdir):
    accessions=read_csv(accessions_file)
    mapped,reverse,source=audit_accessions(accessions,receipt_file,msa_file)
    trimmed,screen=trim_alignment(mapped)
    d=Path(outdir);d.mkdir(parents=True,exist_ok=True)
    out=d/"S1_ITS_ACTUAL_GENBANK_MAFFT_QC_NUCLEAR_ALIGNMENT.fasta"
    with out.open("w",encoding="ascii") as h:
        for name,seq in sorted(trimmed.items()):
            h.write(f">{name}\n")
            for i in range(0,len(seq),80):h.write(seq[i:i+80]+"\n")
    report={
        "schema":"azami_T2_originalS1_ITS_MAFFT_own_tree_QC_v1",
        "original_registry":source,"sequence_quality":screen,
        "no_Moreyra_or_older_published_backbone_used_for_ITS_tree":True,
        "leaf_tip_original_GenBank_accessions_linked_1_to_1_to_S1_source_names":True,
        "single_NCBI_ITS_not_Cirsium_multilocus_species_tree":True,
        "MIDPOINT_ROOT_IF_LATER_USED_IS_NOT_EVOLUTIONARY_ANCESTRAL_ROOT":True,
        "present_annual_rain_not_historic_rain_at_state_transition":True,
    }
    (d/"S1_GENBANK_ITS_ML_PREPARATION.json").write_text(
        json.dumps(report,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("S1_OWN_GENBANK_ITS_TREE_PREPARATION="+json.dumps(report,sort_keys=True),
          flush=True)
    return screen


def summarize_tree(accessions_file,receipt_file,qc_file,tree_file,outdir):
    rows=read_csv(accessions_file)
    original=json.loads(Path(receipt_file).read_text(encoding="utf-8"))
    qc=json.loads(Path(qc_file).read_text(encoding="utf-8"))
    if qc["sequence_quality"]["status"]!="PASS_TO_GENBANK_ITS_IQTREE_EXPLORATORY_ML":
        raise ValueError("Cannot infer ML tree without authentic ITS alignment QC")
    t=Phylo.read(str(tree_file),"newick")
    n=len(rows)
    tips=[r.name for r in t.get_terminals()]
    names={r["original_source_name"].replace(" ","_") for r in rows}
    if n<MIN_TAXA_FOR_ITS_ML or set(tips)!=names or len(set(tips))!=n:
        raise ValueError("IQ-TREE inferred GenBank tree lacks original authenticated S1 taxa")
    if any(z.branch_length is None or z.branch_length<0 for z in t.find_clades() if z is not t.root):
        raise ValueError("Invalid estimated ML nuclear ITS branch lengths")
    bytaxon={r["original_source_name"].replace(" ","_"):r for r in rows}
    known={tip:bytaxon[tip] for tip in tips if bytaxon[tip]["botanical_U_D_state"] in ("U","D")}
    D=sum(x["botanical_U_D_state"]=="D" for x in known.values())
    U=sum(x["botanical_U_D_state"]=="U" for x in known.values())
    result={
        "schema":"azami_T2_S1_all_original_name_GenBank_ITS_IQTREE_own_phylogeny_rain_v1",
        "status":"EMPIRICAL_GENBANK_ACCESSION_DERIVED_EXPLORATORY_ITS_TREE",
        "source_original_WCVP_accepted_taxon_keys":244,
        "source_S1_original_binomial_name_units":251,
        "GenBank_ITS_ML_tree_source_taxa":n,
        "S1_source_taxa_without_qualifying_GenBank_ITS":251-n,
        "author_scored_botanical_nodding_tips_in_ITS_tree":D,
        "author_scored_botanical_erect_tips_in_ITS_tree":U,
        "author_scored_botanical_unknown_tips_in_ITS_tree":n-D-U,
        "reconstructed_S1_ITS_ML_tree_sha256":hashlib.sha256(Path(tree_file).read_bytes()).hexdigest(),
        "GenBank_source_registry_sha256":hashlib.sha256(Path(accessions_file).read_bytes()).hexdigest(),
        "trimmed_original_nuclear_ITS_queried_loci":1,
        "no_specimen_voucher_ITS_copy_number_validation":True,
        "Cirsium_ITS_gene_history_not_single_copy_gene_orthology":True,
        "tree_unrooted_ML_midpoint_root_only_exploratory":True,
        "flowering_study_climate_BIO12_from_S1_species_photo_source_centroids":True,
        "pheno_state_independent_botanical_source_not_photo_EXIF_angle":True,
        "rain_selection_convergence_or_ancestral_U_to_D_gains_confirmed":False,
        "real_phylogenetic_comparison_performed":False,
    }
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    if D>=MIN_PHYLO_POSTURE_STATES[0] and U>=MIN_PHYLO_POSTURE_STATES[1]:
        try:
            t.root_at_midpoint()
            trimmed=deepcopy(t)
            for tip in list(trimmed.get_terminals()):
                if tip.name not in known: trimmed.prune(tip)
            if trimmed.root.branch_length is not None:
                trimmed.root.branch_length=0.
            tipnames=[x.name for x in trimmed.get_terminals()]
            orig=[known[x] for x in tipnames]
            C=brownian_shared_ancestry(trimmed,tipnames)
            y=np.array([float(x["chelsa_bio12"]) for x in orig])
            trait=np.array([x["botanical_U_D_state"]=="D" for x in orig],float)
            temp=np.array([float(x["chelsa_bio01"]) for x in orig])
            X=np.column_stack((np.ones(len(orig)),trait))
            Xt=np.column_stack((X,(temp-temp.mean())/temp.std()))
            models={
                "no_covariate_BIO12_mm":fit_models(y,X,C),
                "BIO1_adjusted_BIO12_mm":fit_models(y,Xt,C),
            }
            independent={
                "no_covariate_BIO12_mm":equal_variance_independent_OLS(y,X),
                "BIO1_adjusted_BIO12_mm":equal_variance_independent_OLS(y,Xt),
            }
            result["real_phylogenetic_comparison_performed"]=True
            result["ITS_midpoint_root_sensitivity_BIO12_GLS"]=models
            result["IID_species_equal_variance_OLS_same_exact_tips"]=independent
            result["BOTANICAL_tip_climate_D_minus_U_mean_mm"]=float(y[trait==1].mean()-y[trait==0].mean())
            result["phylogenetic_covariance_molecular_midpoint_root_assumption"]=True
            result["root_ancestral_posture_inferred"]=False
        except (ValueError,np.linalg.LinAlgError) as e:
            result["posture_BIO12_status"]="BLOCKED_ITS_PGLS_COVARIANCE_OR_IDENTIFIABILITY"
            result["posture_BIO12_reason"]=str(e)[:250]
    else:
        result["posture_BIO12_status"]="BLOCKED_TOO_FEW_INDEPENDENT_BOTANICAL_U_D_ON_GENBANK_ITS_TREE"
    report=out/"S1_GENBANK_ITS_REAL_TREE_BIO12_COMPARISON.json"
    report.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2,
                                 allow_nan=False)+"\n",encoding="utf-8")
    print("S1_OWN_GENBANK_ITS_TREE_VS_RAIN="+json.dumps({
        "status":result["status"],"tree_tips":n,"n_D":D,"n_U":U,
        "climate_phylogenetic_comparison":result["real_phylogenetic_comparison_performed"],
        "estimated_effect":result.get("BOTANICAL_tip_climate_D_minus_U_mean_mm"),
        "block_reason":result.get("posture_BIO12_reason"),
    },sort_keys=True),flush=True)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--stage",choices=("prepare","report"),required=True)
    p.add_argument("--accessions",type=Path,required=True)
    p.add_argument("--source-receipt",type=Path,required=True)
    p.add_argument("--msa",type=Path)
    p.add_argument("--qc",type=Path)
    p.add_argument("--tree",type=Path)
    p.add_argument("--outdir",type=Path,required=True)
    a=p.parse_args()
    if a.stage=="prepare":
        if a.msa is None:raise ValueError("--msa required for prepare")
        prepare(a.accessions,a.source_receipt,a.msa,a.outdir)
    else:
        if a.qc is None or a.tree is None:raise ValueError("--qc and --tree required")
        summarize_tree(a.accessions,a.source_receipt,a.qc,a.tree,a.outdir)


if __name__=="__main__":
    main()
