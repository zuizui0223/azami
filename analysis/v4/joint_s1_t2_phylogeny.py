"""One source-backed S1 spatial / T2 historical evidence ledger, without re-fitting.

S1: 170,066 observation-level image-vertical angles / 244 native accepted keys
and the *same* five jointly adjusted CHELSA exposures. The 52 GBOTB grafted
trees are S1 placement sensitivity on 169,912 images / 194 accepted keys,
not a sample of independently evolved head postures.

T2: only independently author/authority-scored biological U/D states, never
angle-thresholded S1 photos. Joins are exact source concept to accepted key;
collided WCVP keys are unusable for ancestry until adjudication.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import urllib.request

import numpy as np
import pandas as pd

from analysis.v4.native_integration_audit import eligible

ENDPOINT="orientation_image_vertical_angle"
EXPOSURES=("chelsa_bio01","chelsa_bio04","chelsa_bio12",
           "chelsa_bio15","chelsa_rsds_mean")
T2_SOURCE_EA="zuizui0223/EAzami"
T2_COMMIT_EA="c9870e91604ed1514344ddc8f765eae5237dbf27"
T2_PATH_EA="data/evidence/fdt4_existing_material_join_audit_v1.csv"
T2_BLOB_EA="f12ce3f1665c2a4d3de01cd34dea5741812e2af1"
T2_SOURCE_NA="carol-siniscalchi/cirsium"
T2_COMMIT_NA="c925d4af86160064448f794a9715f70b917a2786"
T2_PATH_NA="traits/traits.csv"
T2_BLOB_NA="fa127f0587b676e6e004ea121424d26069586bcb"
BASE=Path("analysis/v4/results/20261020_native_integration")
TREE=Path("analysis/v4/results/20261020_native_phylogeny")
MODEL_PATH=BASE/"native_trait_models.json"
S1_PATH=BASE/"native_trait_effects.csv"
PGLS_PATH=TREE/"phylogenetic_summary.csv"
MAPPING_PATH=TREE/"taxon_mapping.json"
COVERAGE_PATH=TREE/"coverage.json"


def blob_sha(raw):
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()


def fetch_pinned(repo, ref, path, expected):
    url=f"https://raw.githubusercontent.com/{repo}/{ref}/{path}"
    with urllib.request.urlopen(url,timeout=100) as response:
        raw=response.read()
    if blob_sha(raw)!=expected:
        raise ValueError(f"External evidence drift {path}")
    return list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig")))), hashlib.sha256(raw).hexdigest()


def orientation_state_rows(east_asia, north_america):
    """Accept existing reviewed taxon states only; never relabel S1 photograph angles."""
    if len(east_asia)!=20 or len(north_america)!=63:
        raise ValueError("Fixed published biological state inventories changed")
    rows=[]
    seen=set()
    for r in east_asia:
        name=r["accepted_taxon"].strip()
        state=r["orientation_state"]
        if name in seen:
            raise ValueError("Duplicate East Asian concept")
        seen.add(name)
        if r["orientation_resolved"]=="true" and state in ("D","U"):
            rows.append(dict(source="EAzami_FDT4_existing_authority",
                source_taxon=name, state="nodding" if state=="D" else "erect",
                evidence_level="EAzami_curated_not_new_double_blind",
                source_concept=name))
    seen=set()
    for r in north_america:
        name=r["Species"].replace("_"," ").strip()
        if name in seen:
            raise ValueError("Duplicate North American source name")
        seen.add(name)
        if r["HeadPosition"] not in ("0","1"):
            raise ValueError("Invalid original author morphological HeadPosition")
        rows.append(dict(source="Siniscalchi2023_author_scored",
             source_taxon=name,state="nodding" if r["HeadPosition"]=="1" else "erect",
             evidence_level="published_author_score_not_new_double_blind",
             source_concept=name))
    if len(rows)!=80 or sum(r["source"]=="Siniscalchi2023_author_scored" and
                           r["state"]=="nodding" for r in rows)!=5:
        raise ValueError("Biological state source count changed")
    return rows


def s1_table(s1,phylo,models,coverage):
    expected=tuple(models["selection"]["selected"])
    if expected != EXPOSURES:
        raise ValueError("Frozen five-exposure model differs from preregistered S1")
    if models["input_sha256"]!="342d78e9a8a960e8d9bdcf0bc2cdacd13bee48214a8b962419d5a71d5aef46ed":
        raise ValueError("Frozen native corrected source SHA changed")
    s1=s1.loc[(s1.module=="orientation") & (s1.endpoint==ENDPOINT)]
    phylo=phylo.loc[phylo.endpoint==ENDPOINT]
    terms=[f"{e}__{level}" for e in EXPOSURES for level in ("between","within")]
    if (len(s1)!=10 or len(phylo)!=10 or set(s1.term)!=set(terms) or
            set(phylo.term)!=set(terms) or
            not s1.n.eq(170066).all() or not s1.species.eq(244).all() or
            not phylo.n.eq(169912).all() or not phylo.species.eq(194).all()):
        raise ValueError("Frozen S1/PGLS orientation scope mismatch")
    if (coverage["full_n"],coverage["covered_n"],coverage["full_species"],
            coverage["covered_species"])!=(170066,169912,244,194):
        raise ValueError("Tree-covered source and S1 eligibility mismatch")
    table=[]
    for term in terms:
        primary=s1.loc[s1.term==term].iloc[0]
        spatial=phylo.loc[phylo.term==term].iloc[0]
        increment=float(primary.environmental_increment)
        if abs(float(spatial.increment)-increment)>1e-10:
            raise ValueError("PGLS increment differs from S1 per fixed exposure")
        table.append({
            "term":term, "exposure":term.split("__")[0],
            "scale":term.split("__")[1],
            "increment":increment,
            "S1_full_n":170066, "S1_full_accepted_keys":244,
            "S1_effect_image_degrees":float(primary.trait_change),
            "S1_interval": [float(primary.lower95),float(primary.upper95)],
            "S1_tree_covered_n":169912, "S1_tree_covered_keys":194,
            "placement_scenario_count":52,
            "placement_effect_min":float(spatial.phy_min),
            "placement_effect_max":float(spatial.phy_max),
            "placement_95pct_model_interval_envelope":[float(spatial.lower_envelope),
                                                     float(spatial.upper_envelope)],
            "placement_pointwise_positive_interval_count":int(spatial.positive_intervals),
            "placement_pointwise_negative_interval_count":int(spatial.negative_intervals),
            "supported_across_52_pointwise":(
               "all_positive" if int(spatial.positive_intervals)==52 else
               "all_negative" if int(spatial.negative_intervals)==52 else "no_uniform_pointwise_support"),
            "selection_boundary":"same five co-adjusted predictors, no post-hoc EAzami axis selection",
            "phenotype":"EXIF_image_vertical_angle_NOT_gravitational_head_tilt",
        })
    return table


def source_concept_index(mapping):
    by_key={}
    for row in mapping:
        key=str(row["accepted_key"])
        if key in by_key:
            raise ValueError("Duplicated accepted key in historical matching")
        by_key[key]={"source_names":set(row["source_names"]),"tree_matches":set(row["tree_matches"]),
                     "status":row["status"]}
    return by_key


def actual_native_roster(frame,models,coverage):
    """Read already-decoded private observation layer, output only public taxon summaries."""
    if not frame.obs_id.is_unique:
        raise ValueError("Native observation IDs not unique")
    if not all(col in frame for col in ["taxon_name","accepted_key",ENDPOINT,ENDPOINT+"__eligible",*EXPOSURES]):
        raise ValueError("Decoded source missing exact S1 fields")
    mask=eligible(frame,[ENDPOINT]) & np.isfinite(frame[list(EXPOSURES)]).all(axis=1)
    sample=frame.loc[mask,["accepted_key","taxon_name",ENDPOINT]]
    if len(sample)!=170066 or sample.accepted_key.nunique()!=244:
        raise ValueError("Native 244/170066 S1 population has changed")
    by_key={}
    for key,group in sample.groupby("accepted_key",sort=True):
        names=sorted(set(group.taxon_name.dropna().astype(str).str.strip()))
        data=group[ENDPOINT].to_numpy(float)
        by_key[str(key)]={"source_names":names,"n":int(len(group)),
                          "median_image_angle_deg":float(np.median(data)),
                          "mean_image_angle_deg":float(np.mean(data)),
                          "concept_collision_in_source":len(names)!=1}
    return by_key


def match_sources(states,source,tree_index,*,scope_complete=True):
    tree_tip_names=set()
    for row in tree_index.values():
        tree_tip_names.update(row["tree_matches"])
    result=[]
    for r in states:
        canonical=r["source_taxon"].replace(" ","_")
        matching=[(key,v) for key,v in source.items()
                  if canonical in {x.replace(" ","_") for x in v["source_names"]}]
        item={**r,"s1_validated_taxon_match":False,
              "s1_phylogeny_194_key_match":False,
              "join_status":"not_in_S1_orientation_cohort",
              "accepted_key":None}
        if len(matching)==1:
            key,group=matching[0]
            if len(group["source_names"])!=1:
                item["join_status"]="blocked_WCVP_source_name_collision"
            elif group.get("n",0)<=0:
                item["join_status"]="missing_S1_eligibility"
            else:
                item.update(accepted_key=key,
                    s1_validated_taxon_match=True,
                    s1_observations=group["n"],
                    s1_image_orientation_median=group["median_image_angle_deg"],
                    s1_phylogeny_194_key_match=bool(key in tree_index and
                        canonical in tree_index[key]["tree_matches"]),
                    join_status="exact_S1_and_tree" if
                        key in tree_index and canonical in tree_index[key]["tree_matches"]
                        else "exact_S1_only_tree_uncovered")
        elif len(matching)>1:
            item["join_status"]="ambiguous_multiple_WCVP_accepted_keys"
        result.append(item)
    return result


def build(s1_file,phylo_file,models_file,coverage_file,map_file,ea_rows,na_rows,
          native_frame=None,native_source_sha=None):
    models=json.loads(Path(models_file).read_text())
    coverage=json.loads(Path(coverage_file).read_text())[ENDPOINT]
    s1=pd.read_csv(s1_file)
    phylo=pd.read_csv(phylo_file)
    effects=s1_table(s1,phylo,models,coverage)
    mapping=json.loads(Path(map_file).read_text())
    index=source_concept_index(mapping)
    tree_index={key:val for key,val in index.items() if val["status"]=="matched"}
    if len(tree_index)!=194:
        raise ValueError("Expected exactly 194 accepted-key phylogenetic matches")
    states=orientation_state_rows(ea_rows,na_rows)
    if native_frame is not None:
        if hashlib.sha256(Path(native_source_sha).read_bytes()).hexdigest()!=models["input_sha256"]:
            raise ValueError("Decoded private native input hash does not match frozen S1")
        source=actual_native_roster(native_frame,models,coverage)
        status="empirical_full_244_native_join_source_sha_verified"
    else:
        # This mode has only the historical mapping (306 keys) and does not
        # know which of 306 were among the 244 orientation-eligible keys.
        # Therefore it NEVER reports exact 244-state overlap.
        source={}
        status="source_244_roster_unavailable_do_not_infer_taxon_state_coverage"
    joined=match_sources(states,source,tree_index)
    def count(layer,metric):
        return sum(x["source"]==layer and x.get(metric) is True for x in joined)
    return {
      "status":status,
      "current_primary_question":"rain_and_Cirsium_capitulum_orientation",
      "model_scope":{"native_accepted_keys":306,"orientation_eligible_accepted_keys":244,
                     "orientation_eligible_observations":170066,
                     "placement_tree_covered_keys":194,
                     "placement_tree_covered_observations":169912,
                     "grafted_placement_scenarios":52},
      "fixed_covariates":list(EXPOSURES),
      "frozen_S1_five_exposure_effects":effects,
      "T2_independent_state_inventory":{
         "NA_original_author_states":63,"NA_original_nodding_states":5,
         "EAzami_source_concept_state_rows":20,"EAzami_resolved_states":17,
         "not_all_244_have_validated_botanical_U_D":True,
         "same_photo_derived_state_or_inferred_threshold_used":False,
         "NA_exact_S1_244_matches":count("Siniscalchi2023_author_scored","s1_validated_taxon_match"),
         "NA_exact_S1_194_tree_matches":count("Siniscalchi2023_author_scored","s1_phylogeny_194_key_match"),
         "EastAsia_exact_S1_244_matches":count("EAzami_FDT4_existing_authority","s1_validated_taxon_match"),
         "EastAsia_exact_S1_194_tree_matches":count("EAzami_FDT4_existing_authority","s1_phylogeny_194_key_match"),
         "taxon_source_audit":joined},
      "falsification":{
        "EAzami_BIO15":"positive_EastAsian_U_D_direction_but_no_conditional_single_axis_specificity",
        "S1_BIO12":"higher_annual_precipitation_positive_among_image_angle_244",
        "S1_BIO15":"full_five_covar_interval_crosses_zero_and_all_52_scenario_intervals_cross_zero",
        "not_a_prospective_replication":True,
        "not_an_evolutionary_rate_or_transition_dating_analysis":True,
        "S1_photo_vertical_is_not_botanical_posture":True,
        "historical_52_are_scenarios_not_52_independent_origins":True},
      "publication_routing":"S1 full 244 as observational analysis; one phylogenetic history section with 52 placement robustness plus independent-source historical T2, not separate redundant paper claims",
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out",required=True)
    p.add_argument("--native-table",type=Path)
    p.add_argument("--no-native-diagnostic",action="store_true",
                   help="Allow ledger-only mode: exact 244-source join explicitly unavailable")
    args=p.parse_args()
    if args.native_table is None and not args.no_native_diagnostic:
        raise ValueError("Verified native decoded input required unless explicitly ledger-only")
    ea,ea_sha=fetch_pinned(T2_SOURCE_EA,T2_COMMIT_EA,T2_PATH_EA,T2_BLOB_EA)
    na,na_sha=fetch_pinned(T2_SOURCE_NA,T2_COMMIT_NA,T2_PATH_NA,T2_BLOB_NA)
    source=None
    if args.native_table is not None:
        source=pd.read_csv(args.native_table,low_memory=False)
    result=build(S1_PATH,PGLS_PATH,MODEL_PATH,COVERAGE_PATH,MAPPING_PATH,ea,na,
                 native_frame=source,native_source_sha=args.native_table)
    result["provenance"]={"EAzami_git_commit":T2_COMMIT_EA,
           "EAzami_source_git_blob":T2_BLOB_EA,"EAzami_source_sha256":ea_sha,
           "NA_source_git_commit":T2_COMMIT_NA,"NA_source_git_blob":T2_BLOB_NA,
           "NA_source_sha256":na_sha,
           "S1_model_source_sha256":json.loads(MODEL_PATH.read_text())["input_sha256"],
           "native_decoded_source_authenticated":bool(args.native_table)}
    output=Path(args.out)
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("V4_RAIN_S1_T2_JOINT_SUMMARY="+json.dumps({
        "status":result["status"],
        "model_scope":result["model_scope"],
        "T2_inventory":{k:v for k,v in result["T2_independent_state_inventory"].items()
                        if k!="taxon_source_audit"},
        "BIO12_between":next(x for x in result["frozen_S1_five_exposure_effects"]
                             if x["term"]=="chelsa_bio12__between"),
        "BIO15_between":next(x for x in result["frozen_S1_five_exposure_effects"]
                             if x["term"]=="chelsa_bio15__between"),
    },ensure_ascii=False,sort_keys=True),flush=True)


if __name__=="__main__":
    main()
