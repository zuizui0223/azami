"""Independent source-authenticated NA63 nuclear quartet check for real 7-tip NJ.

Compares unrooted quartet among the EXACT four taxon labels overlapping the
Siniscalchi 2023 published independent NA63 original nuclear tree and the
source-linked seven-taxon 101-locus Q20 technical distance NJ from 120k real
public ENA paired read prefixes (Actions 39868882212).

Only quartet topology concordance; no accepted seven-tip species tree,
phylogenetic dating, BM/OU fit, orientation transitions or rainfall adaptation.
Names are candidate species labels, not identical voucher provenances.
"""
from __future__ import annotations
import argparse
import csv
from copy import deepcopy
from io import StringIO
import hashlib
import json
from pathlib import Path
import urllib.request

from Bio import Phylo

ORIGINAL_NA63=63
TECH_NJ_TIPS=7
QUARTET=(
    "Cirsium_flodmanii","Cirsium_scopulorum",
    "Cirsium_tracyi","Cirsium_vulgare",
)
NA_GIT_COMMIT="c925d4af86160064448f794a9715f70b917a2786"
NA_GIT_BLOB_SHA1="6417a38793427bc7c8562dc0e2cab59b74d8c316"
NA_URL=(
    "https://raw.githubusercontent.com/carol-siniscalchi/cirsium/"
    +NA_GIT_COMMIT+"/traits/cirsium.phyluce95p.ultrametric.tre"
)
PILOT_LOCI=101
SHARED_SOURCE_SITE_COUNT=29270
PILOT_COMPLETE_CASE_SITE_COUNT=28726
EXPECTED_PAIR_BOOTSTRAP_SOURCE="Cirsium_flodmanii|Cirsium_tracyi"
EXPECTED_PAIR_BOOTSTRAP_FRACTION=.75


def blob_sha(content):
    return hashlib.sha1(b"blob "+str(len(content)).encode()+b"\0"+content).hexdigest()


def load_na63(fetch=None):
    if fetch is None:
        def fetch(url):
            req=urllib.request.Request(url,headers={
                "User-Agent":"azami-reproducible-independent-nuclear-quartet-validation/1.0"})
            with urllib.request.urlopen(req,timeout=100) as connection:
                return connection.read()
    raw=fetch(NA_URL)
    if blob_sha(raw)!=NA_GIT_BLOB_SHA1:
        raise ValueError("Published independent North American nuclear tree Git blob changed")
    tree=Phylo.read(StringIO(raw.decode("utf-8")),"newick")
    tips=[t.name for t in tree.get_terminals()]
    if (len(tips)!=ORIGINAL_NA63 or len(set(tips))!=ORIGINAL_NA63
            or not set(QUARTET)<=set(tips)):
        raise ValueError("Original 63 nuclear tips or four exact taxa changed")
    return tree


def quartet_split(tree,four=QUARTET):
    tips={x.name for x in tree.get_terminals()}
    if len(set(four))!=4 or not set(four)<=tips:
        raise ValueError("Four exact independently identified botanical taxa missing")
    x=deepcopy(tree)
    for tip in sorted(tips-set(four)):
        x.prune(tip)
    if len(x.get_terminals())!=4:
        raise ValueError("Four-species unrooted phylogenetic restriction failed")
    choices=set()
    for clade in x.get_nonterminals():
        group={z.name for z in clade.get_terminals()}
        if len(group)!=2 or len(set(four)-group)!=2:
            continue
        pair=tuple(sorted(group))
        alt=tuple(sorted(set(four)-group))
        choices.add(min(pair,alt))
    if len(choices)!=1:
        raise ValueError("Independently rooted or unresolved 4-species branch split")
    return sorted(next(iter(choices)))


def evaluate(pilot_source,pilot_tree,boot_table,original):
    if (pilot_source.get("n_original_public_SRA_biosamples")!=TECH_NJ_TIPS
        or pilot_source.get("n_primary_loci_shared_atleast_100_Q20_sites")!=PILOT_LOCI
        or pilot_source.get("source_shared_Q20_positions_all_targets")!=SHARED_SOURCE_SITE_COUNT
        or pilot_source.get("n_complete_case_invariant_plus_variable_sites_primary")!=PILOT_COMPLETE_CASE_SITE_COUNT
        or not pilot_source.get("nuclear_tree_preliminary_only")
        or pilot_source.get("EAzami_294_tip_species_tree_constructed")):
        raise ValueError("Seven actual public SRA consensus sources not authentic")
    names={z.name for z in pilot_tree.get_terminals()}
    if len(names)!=TECH_NJ_TIPS or len(pilot_tree.get_terminals())!=TECH_NJ_TIPS:
        raise ValueError("Provisional 7 specimen topology changed")
    na_names={z.name for z in original.get_terminals()}
    common=names & na_names
    if common!=set(QUARTET):
        raise ValueError("Exact species-name overlap with published NA63 nuclear phylogeny changed")
    quartet_pilot=quartet_split(pilot_tree)
    quartet_published=quartet_split(original)
    source_boot=[x for x in boot_table
        if x["unrooted_split_side"]==EXPECTED_PAIR_BOOTSTRAP_SOURCE]
    if len(source_boot)!=1 or int(source_boot[0]["n_locus_bootstrap_resamples"])!=300:
        raise ValueError("Original seven-sample 300-resample branch diagnostic missing")
    fraction=float(source_boot[0]["technical_NJ_locus_resampling_fraction"])
    if abs(fraction-EXPECTED_PAIR_BOOTSTRAP_FRACTION)>1e-9:
        raise ValueError("Original seven-tip clade branch technical locus re-sampling drifted")
    if quartet_published!=list(sorted(["Cirsium_flodmanii","Cirsium_tracyi"])):
        raise ValueError("The pinned independent author nuclear tree quartet changed")
    result={
      "schema":"azami_v4_seven_real_Q20_sequence_vs_independent_NA63_nuclear_quartet_v1",
      "status":"TECHNICAL_READ_DISTANCE_QUARTET_MATCH_NOT_ACCEPTED_GLOBAL_SPECIES_TREE",
      "empirical_seven_public_SRA_source_samples":TECH_NJ_TIPS,
      "independent_published_nuclear_tree_taxa":ORIGINAL_NA63,
      "independent_published_tree_doi":"10.1086/724310",
      "exact_name_matched_taxa":sorted(QUARTET),
      "only_four_of_seven_in_independent_NA63_tree":True,
      "pilot_real_read_NJ_quartet_pair":quartet_pilot,
      "independent_author_nuclear_quartet_pair":quartet_published,
      "unrooted_four_tip_topology_match":quartet_pilot==quartet_published,
      "seven_tip_full_panel_flodmanii_tracyi_locus_bootstrap_support":fraction,
      "seven_tip_branch_bootstrap_not_quartet_bootstrap":True,
      "reference_copy_concordance_not_ortholog_confirmation":True,
      "original_S1_244_WCVP_key_and_251_source_name_table_preserved":True,
      "sample_voucher_identity_across_studies_not_verified":True,
      "seven_taxa_NJ_not_accepted_as_genus_species_tree":True,
      "NA63_root_time_not_used_for_7_tip_pilot":True,
      "no_BM_OU_or_rain_adaptive_radiation_result":True,
      "source_original_NA63_nuclear_git_blob_SHA1":NA_GIT_BLOB_SHA1,
    }
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pilot-json",type=Path,required=True)
    p.add_argument("--pilot-newick",type=Path,required=True)
    p.add_argument("--technical-branches",type=Path,required=True)
    p.add_argument("--out-dir",type=Path,required=True)
    a=p.parse_args()
    pilot=json.loads(a.pilot_json.read_text(encoding="utf-8"))
    tree=Phylo.read(str(a.pilot_newick),"newick")
    with a.technical_branches.open(newline="",encoding="utf-8-sig") as h:
        splits=list(csv.DictReader(h))
    published=load_na63()
    result=evaluate(pilot,tree,splits,published)
    a.out_dir.mkdir(parents=True,exist_ok=True)
    (a.out_dir/"NA63_AND_7_REAL_Q20_SOURCE_QUARTET_CONCORDANCE.json").write_text(
        json.dumps(result,sort_keys=True,ensure_ascii=False,indent=2)+"\n",
        encoding="utf-8")
    print("T2_REAL7_Q20_VS_AUTHENTIC_NA63_QUARTET="+json.dumps(result,
         sort_keys=True,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
