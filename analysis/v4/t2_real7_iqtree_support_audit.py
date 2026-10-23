"""Audit actual IQ-TREE 2.0.7 SH-aLRT/UFBoot branch labels on two real 7-SRA Q20 trees.

IQ-TREE outputs SH-aLRT/UFBoot numeric pairs such as )92.6/99:0.001.
Supports describe site-resampled, technical ML inferred on same reference-biased
seven-way complete-case observed Q20 sites, NOT orthology validated or dated
species tree; cross-referenced with earlier whole-locus NJ bootstrap frequencies.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from io import StringIO
import json
from pathlib import Path
import re

from Bio import Phylo

from analysis.v4 import t2_real7_author101_iqtree_ml as ml
from analysis.v4 import t2_real7_NA63_quartet_crosscheck as author

SUPPORT_REGEX=re.compile(r"^\s*(\d+(?:\.\d+)?)/(\d+(?:\.\d+)?)\s*$")
EXPECTED_N_SAMPLES=7
EXPECTED_N_INNER_SPLITS=4
PREFIXES=("primary_101","sensitivity_110")

def support(clade):
    values=[str(x).strip() for x in (clade.name,clade.comment) if x is not None]
    for item in values:
        m=SUPPORT_REGEX.fullmatch(item)
        if m:
            sh=float(m.group(1));uf=float(m.group(2))
            if not 0<=sh<=100 or not 0<=uf<=100:
                raise ValueError("Invalid IQ-TREE SH-aLRT/UFBoot numeric branch pair")
            return sh,uf
    if clade.confidence is not None:
        # Unlabelled numeric confidence cannot be assigned to SH-aLRT or UFBoot.
        return None,None
    return None,None


def read_splits(path):
    raw=Path(path).read_text(encoding="utf-8").strip()
    if not raw.endswith(";"):
        raise ValueError("ML Newick missing final terminator")
    tree=Phylo.read(StringIO(raw),"newick")
    names={x.name for x in tree.get_terminals()}
    if len(names)!=EXPECTED_N_SAMPLES:
        raise ValueError("Unexpected source exact seven sample count")
    result={}
    for clade in tree.get_nonterminals():
        part={t.name for t in clade.get_terminals()}
        comp=names-part
        if min(len(part),len(comp))<2:
            continue
        key=min(tuple(sorted(part)),tuple(sorted(comp)),
                key=lambda x:(len(x),x))
        pair=support(clade)
        if key in result:
            # A Newick 'root' can duplicate the same unrooted split on both
            # child sides; only one of their support values might be labelled.
            if result[key]!=(None,None) and pair!=(None,None) and result[key]!=pair:
                raise ValueError("Conflicting IQ-TREE support labels on one unrooted split")
            if result[key]==(None,None):
                result[key]=pair
        else:
            result[key]=pair
    if len(result)!=EXPECTED_N_INNER_SPLITS:
        raise ValueError("Binary seven-tip unrooted ML split universe changed")
    return tree,result


def read_nj_bootstrap(path):
    with Path(path).open(newline="",encoding="utf-8-sig") as f:
        rows=list(csv.DictReader(f))
    if len(rows)!=EXPECTED_N_INNER_SPLITS:
        raise ValueError("Original 7-tip locus-bootstrap NJ split rows missing")
    out={}
    for row in rows:
        side=tuple(sorted(row["unrooted_split_side"].split("|")))
        if side in out or int(row["n_locus_bootstrap_resamples"])!=300:
            raise ValueError("Original source 300-locus-resampled NJ distinct splits missing")
        x=float(row["technical_NJ_locus_resampling_fraction"])
        if not(0<=x<=1):
            raise ValueError("Malformed original real-sequence NJ locus bootstrap proportion")
        out[side]=x
    return out


def evaluate(primary,sensitivity,nj_bootstrap,verified_ml):
    if verified_ml.get("primary_loci")!=101 or verified_ml.get("sensitivity_loci")!=110 or \
       verified_ml.get("primary_observed_only_Q20_positions")!=28726 or \
       verified_ml.get("sensitivity_observed_only_Q20_positions")!=29270 or \
       verified_ml.get("EAzami_294sample_global_species_tree_reconstructed") or \
       not set(PREFIXES)==set(verified_ml["results"]):
        raise ValueError("Actual 7 SRA ML source verification record mismatch")
    boot=read_nj_bootstrap(nj_bootstrap)
    results={}
    base_splits=None
    for label,path in (("primary_101",primary),("sensitivity_110",sensitivity)):
        tree,splits=read_splits(path)
        if base_splits is None:
            base_splits=set(splits)
        elif set(splits)!=base_splits:
            raise ValueError("Primary 101 and sensitivity 110 ML biological split topology changed")
        names=set(verified_ml["source_real_taxon_labels"])
        if {x.name for x in tree.get_terminals()}!=names:
            raise ValueError("Original-seven-ENA-samples ML terminal labels changed")
        if not verified_ml["results"][label]["ML_matches_published_independent_NA63_4tip_quartet"]:
            raise ValueError("Original source-pinned ML versus independent 4-tip comparison changed")
        if verified_ml["results"][label]["n_ML_splits_same_as_uncorrected_NJ"]!=4:
            raise ValueError("ML/NJ topology was not 4/4 under the original case")
        result=[]
        for split,(sh,uf) in sorted(splits.items()):
            result.append({
                "side":"|".join(split),
                "sh_aLRT_percent":sh,
                "ufboot_site_resample_percent":uf,
                "NJ_300_whole_locus_resample_fraction":boot.get(split),
                "same_biological_specimens_as_NJ":True,
                "difference_site_vs_locus_bootstrap_biological_independence":True,
            })
        results[label]=result
    s=results["primary_101"]
    return {
       "schema":"azami_v4_real7_101v110_IQTREE_support_vs_NJ_whole_locus_bootstrap_v1",
       "source_real_SRA_specimens":EXPECTED_N_SAMPLES,
       "primary_observed_only_loci":101,
       "sensitivity_original_copy_loci":110,
       "original_public_sra_source_independent_NA63_quartet_match":True,
       "primary_and_sensitivity_topology_equal":True,
       "primary_internal_binary_unrooted_splits":len(s),
       "primary_sh_aLRT_and_UFBOOT_both_parseable":
          all(row["sh_aLRT_percent"] is not None and
              row["ufboot_site_resample_percent"] is not None for row in s),
       "n_primary_UFBOOT_ge95":sum(
           x["ufboot_site_resample_percent"] is not None and
           x["ufboot_site_resample_percent"]>=95 for x in s),
       "n_primary_NJ_whole_locus_boot_ge95":sum(
           x["NJ_300_whole_locus_resample_fraction"] is not None and
           x["NJ_300_whole_locus_resample_fraction"]>=.95 for x in s),
       "branches_by_alignment":results,
       "source_public_read_nucleotide_copy_homeology_validated":False,
       "same_empirical_7_run_ML_vs_NJ_not_independent":True,
       "site_bootstrap_not_equal_independent_locus_bootstrap":True,
       "cannot_claim_moreyra_350_reconstruction":True,
       "cannot_claim_head_angle_phylogenetic_adaptation":True,
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--primary-tree",type=Path,required=True)
    p.add_argument("--sensitivity-tree",type=Path,required=True)
    p.add_argument("--nj-locus-support",type=Path,required=True)
    p.add_argument("--ml-verified",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    verified=json.loads(a.ml_verified.read_text(encoding="utf-8"))
    result=evaluate(a.primary_tree,a.sensitivity_tree,
                    a.nj_locus_support,verified)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print("T2_REAL7_IQTREE_UFBOOT_SHaLRT_WHOLE_LOCUS_SUPPORT="+json.dumps({
        "primary_UFBOOT_ge95":result["n_primary_UFBOOT_ge95"],
        "NJ_whole_locus_ge95":result["n_primary_NJ_whole_locus_boot_ge95"],
        "support_parseable":result["primary_sh_aLRT_and_UFBOOT_both_parseable"],
        "branches":result["branches_by_alignment"]["primary_101"],
    },sort_keys=True,ensure_ascii=False),flush=True)
if __name__=="__main__":
    main()
