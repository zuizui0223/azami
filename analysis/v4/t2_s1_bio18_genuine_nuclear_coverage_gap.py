"""Real ORIGINAL S1 251 4-climate BIO18 source coverage on two true nuclear trees.

This is a scientific sampling-bias diagnosis, NOT an additional association
fit. Frozen all-S1 BIO18 top quartile must be set before inspecting nuclear
tree tip labels. Compare original photo means, n photos, and true independent
NA42/JP12 nuclear species tree, never artificial global 54/161-tip tree.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from pathlib import Path
import numpy as np
from analysis.v4 import t2_s1_continuous_genuine_nuclear_two_tree_clustering as dna
from analysis.v4 import t2_s1_fourclimate_real_nuclear_pgls_pic as panel

SOURCE_251_SHA=dna.S1_SHA
NEW_S1_FOUR_SHA=panel.NEW_AGGREGATE_SHA
FOUR=panel.COVARIATES
EXPECTED_CUTOFF=588.2739645
EXPECTED_WET=63

def summarize(original,new,trees,raw_sources):
    original_by={r["original_source_name"]:r for r in original.values()}
    by={r["original_source_name"]:r for r in new.values()}
    if len(original_by)!=251 or len(by)!=251 or set(original_by)!=set(by):
        raise ValueError("Source251 original scientific names mismatch")
    vals=np.array([float(by[n]["chelsa_bio18"]) for n in sorted(by)])
    threshold=float(np.quantile(vals,.75))
    wet={n for n in by if float(by[n]["chelsa_bio18"])>=threshold}
    if len(wet)!=EXPECTED_WET or not np.isclose(threshold,EXPECTED_CUTOFF,rtol=0,atol=1e-4):
        raise ValueError("Original BIO18 highest quartile source definition drifted")
    result={}
    all_matches=set()
    labels_wet=[]
    for label,(tree,matching) in trees.items():
        candidates=set(matching.values())
        if len(candidates)!=(42 if label.startswith("NA") else 12):
            raise ValueError("Genuine strict nuclear tree source tip mismatch")
        if all_matches&candidates:
            raise ValueError("Two source nuclear panels share duplicate S1 species")
        all_matches|=candidates
        region=[by[n] for n in candidates]
        wet_in=candidates&wet
        n10_wet={n for n in wet_in if int(by[n]["n_photo_complete_four_axis"])>=10}
        present=np.array([float(r["chelsa_bio18"]) for r in region])
        rec={"n_source_names":len(candidates),
             "n_original_photo_angles_four_axis_complete":sum(int(r["n_photo_complete_four_axis"]) for r in region),
             "n_global_S1_BIO18_wet_quartile_species":len(wet_in),
             "n_BIO18_wet_global_species_with_at_least_10_original_photos":len(n10_wet),
             "n_BIO18_wet_global_species_with_only_1_original_photo":sum(int(by[n]["n_photo_complete_four_axis"])==1 for n in wet_in),
             "n_photo_observations_among_global_top_quartile_species":sum(int(by[n]["n_photo_complete_four_axis"]) for n in wet_in),
             "BIO18_min_mm":float(np.min(present)),
             "BIO18_mean_mm":float(np.mean(present)),
             "BIO18_max_mm":float(np.max(present)),
             "BIO18_global_quartile_threshold_mm":threshold,
             "has_any_highest_global_quartile_BIO18_species":bool(wet_in),
             "source_names_wet_quartile":sorted(wet_in),
             "taxonomy_and_original_photo_count_matched_without_U_D_or_outcome_selection":True}
        result[label]=rec
        labels_wet+=sorted(wet_in)
    if len(all_matches)!=54 or len(wet&all_matches)!=8:
        raise ValueError("True 54 source multilocus match / 8 wet original taxa drift")
    absent=wet-all_matches
    if len(absent)!=55 or result["NA42_genuine_published_nuclear"]["n_global_S1_BIO18_wet_quartile_species"]!=0:
        raise ValueError("Original 55 BIO18-wet species missing from both true nuclear panels")
    return {
        "schema":"V4_T2_S1_SOURCE251_BIO18_GLOBAL_WETTAIL_TRUE_NA42_JAPAN12_NUCLEAR_GAP_20261024_v1",
        "status":"CONFIRMED_TRUE_NUCLEAR_COVERAGE_GAP_FOR_ORIGINAL_S1_BIO18",
        "S1_original_4axis_photo_observations":169020,
        "S1_original_source_names":251,
        "S1_BIO18_global_wettest_quartile_threshold_mm":threshold,
        "S1_original_global_wettest_quartile_species":len(wet),
        "genuine_nuclear_species_tree_panels":result,
        "distinct_original_species_in_true_nuclear_tree_panels":len(all_matches),
        "n_wet_global_species_present_in_two_true_trees":len(wet&all_matches),
        "n_wet_global_species_MISSING_from_two_true_trees":len(absent),
        "wet_missing_species_original_S1_source_names":sorted(absent),
        "wet_covered_species_original_S1_source_names":sorted(wet&all_matches),
        "same_trait_S1_continuous_image_axis_0_to_180":True,
        "no_flora_U_D_used":True,
        "no_plastid_or_ITS_only_or_failed_Q20_genetree_used":True,
        "source_provenance":raw_sources,
        "sampling_scope_warning":"NA42 contains zero of 63 wettest full S1 BIO18 species, JP12 contains eight (4 only 1 original photo); no genus-wide covariate phylogenetic replication can be rejected or confirmed on these two small disjoint trees.",
        "not_two_independent_replications_of_S1":True,
        "no_ancestral_U_to_D_transition_or_climate_causality_claimed":True,
    }

def evaluate(original251,newfour,output):
    if hashlib.sha256(Path(original251).read_bytes()).hexdigest()!=SOURCE_251_SHA:
        raise ValueError("Unverified original source photograph cohort")
    if hashlib.sha256(Path(newfour).read_bytes()).hexdigest()!=NEW_S1_FOUR_SHA:
        raise ValueError("New real BIO18 original photo species source changed")
    nuclear_sources=[dna.upstream(dna.NA_REPO,dna.NA_REF,dna.NA_TREE),
        dna.upstream(dna.J_REPO,dna.J_REF,dna.J_TREE),
        dna.upstream(dna.J_REPO,dna.J_REF,dna.J_MAP),
        dna.upstream(dna.J_REPO,dna.J_REF,dna.J_MANIFEST)]
    original,trees,sources=dna.read_sources(original251,*nuclear_sources)
    four=panel.new_taxon_source(newfour,original)
    result=summarize(original,four,trees,sources)
    out=Path(output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("S1_TO_T2_REAL_NUCLEAR_BIO18_COVERAGE="+json.dumps({
        "global_S1_taxa":251,"global_top25wet":63,
        "true_nuclear_taxa":54,"NA42_wet":result["genuine_nuclear_species_tree_panels"]["NA42_genuine_published_nuclear"]["n_global_S1_BIO18_wet_quartile_species"],
        "JP12_wet":result["genuine_nuclear_species_tree_panels"]["Japan12_genuine_176_loci_nuclear"]["n_global_S1_BIO18_wet_quartile_species"],
        "wet_missing_from_true_nuclear":55,
        "t2_null_cannot_disprove_biogeographic_S1_positive":True
    },sort_keys=True))
    return result

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--original251",type=Path,required=True)
    a.add_argument("--newfour",type=Path,required=True)
    a.add_argument("--out",type=Path,required=True)
    r=a.parse_args()
    evaluate(r.original251,r.newfour,r.out)
if __name__=="__main__":main()
