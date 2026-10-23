"""Empirical BIO12/BIO15/BIO1 by INDEPENDENT botanical U/D on ALL original 251.

Input: original image-eligible 170066 photographs grouped by ORIGINAL source
taxon name nested within WCVP key, plus author-scored botanical capitulum
position from independent NA63/NMNS/EAzami/Flora of China. NEVER derive
upright/nodding from EXIF camera-relative head angle.

All 251 operational units remain in a complete ledger. Inference is limited
to U/D-authority-scored taxa (66), with 185 genuinely unknown; SRA-source
nuclear subset (38 U/D in 111 source real BioSamples) is separate.
Permutation conditions on *authority source lane*, not ancestral geographic
lineage or phylogeny. Results are associations, NOT independent origins,
rainfall-caused adaptation, or comparative phylogenetic OU/BM inference.
"""
from __future__ import annotations

import argparse
from collections import Counter,defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import random
import statistics

TOTAL=251
ACCEPTED=244
IMAGES=170066
NUCLEAR=111
SCORING={"U":54,"D":12,"?":185}
NUCLEAR_SCORING={"U":32,"D":6,"?":73}
SPLIT_KEYS=7
AXES=("chelsa_bio12","chelsa_bio15","chelsa_bio01")
LABELS={"chelsa_bio12":"annual_precipitation_BIO12_mm",
        "chelsa_bio15":"precipitation_seasonality_BIO15",
        "chelsa_bio01":"temperature_BIO1_original_CHELSA_units"}
PERMUTATIONS=4999
SEED=20261023

def normalized(value):
    return " ".join(str(value or "").replace("_"," ").split())

def csv_rows(path):
    with Path(path).open(encoding="utf-8-sig",newline="") as h:
        data=list(csv.DictReader(h))
    if not data:
        raise ValueError("Source original 251 trait/climate table is empty")
    return data

def merge(original,authority):
    if len(original)!=TOTAL or len(authority)!=TOTAL:
        raise ValueError("All original 251 source-name units required")
    by_auth={(r["accepted_key"],normalized(r["original_source_name"])):r
             for r in authority}
    if len(by_auth)!=TOTAL:
        raise ValueError("Ambiguous independent botanical U/D authority source taxon keys")
    if (len({r["accepted_key"] for r in original})!=ACCEPTED or
        sum(int(r["original_photo_count"]) for r in original)!=IMAGES):
        raise ValueError("Original 244 accepted keys/170066 photographed heads changed")
    occurrences=Counter(r["accepted_key"] for r in original)
    if sum(n==2 for n in occurrences.values())!=SPLIT_KEYS or max(occurrences.values())>2:
        raise ValueError("Seven original source-name conflicts not preserved")
    result=[]
    matched=set()
    for r in original:
        name=normalized(r["original_source_name"])
        key=r["accepted_key"];lookup=(key,name)
        auth=by_auth.get(lookup)
        if auth is None:
            raise ValueError("Botanical U/D authority cannot be joined on EXACT source name")
        matched.add(lookup)
        if int(auth["original_S1_photo_count"])!=int(r["original_photo_count"]):
            raise ValueError("Photo count changed when joining botanical U/D source")
        state=auth["botanical_U_D_state"]
        if state not in ("U","D","?"):
            raise ValueError("Unknown botanical posture code")
        src=auth["botanical_source_names"]
        lane=("NA63_author_scored" if "Siniscalchi2023_NA63_HeadPosition" in src
              else "East_Asia_botanical_authority" if state!="?"
              else "no_authority_or_ambiguous")
        nuclear=auth["nuclear_real120k_sample_in_111"]
        if nuclear not in ("True","False"):
            raise ValueError("False positive real nuclear specimen flag")
        climate={}
        for axis in AXES:
            v=float(r[axis])
            if not math.isfinite(v):
                raise ValueError("Source CHELSA climatic taxon centroid is not finite")
            climate[axis]=v
        if climate["chelsa_bio12"]<=0:
            raise ValueError("Source original BIO12 should be positive")
        result.append({
           "accepted_key":key,
           "original_source_name":name,
           "original_photo_count":int(r["original_photo_count"]),
           "WCVP_key_contains_multiple_original_source_names":occurrences[key]>1,
           "botanical_U_D_state":state,
           "independent_botanical_authority_lane":lane,
           "source_botanical_attribution":src,
           "nuclear_real120k_BioSample":nuclear=="True",
           "camera_relative_image_head_angle_mean_deg":float(r["mean_image_angle_deg"]),
           "botanical_state_not_exif_angle_threshold":True,
           "independent_botanical_source_same_SRA_voucher_verified":False,
           "species_tree_independent_origins_inferred":False,
           **climate,
        })
    if len(matched)!=TOTAL:
        raise ValueError("Original 251 WCVP+source names lost")
    states=Counter(r["botanical_U_D_state"] for r in result)
    nuc=Counter(r["botanical_U_D_state"] for r in result if r["nuclear_real120k_BioSample"])
    if states!=SCORING or nuc!=NUCLEAR_SCORING:
        raise ValueError(f"Botanical scored 54/12/185 or real111 32/6/73 source scope changed: {states} {nuc}")
    if sum(r["nuclear_real120k_BioSample"] for r in result)!=NUCLEAR:
        raise ValueError("Exactly 111 real original public nuclear samples required")
    return result

def statistic(rows,states,axis,lanes):
    """Within-authority standardized D-U difference, weighted by effective n."""
    total=0.;weight=0.
    for inds in lanes.values():
        down=[i for i in inds if states[i]=="D"]
        erect=[i for i in inds if states[i]=="U"]
        if not down or not erect:
            continue
        vals=[rows[i][axis] for i in inds]
        std=statistics.stdev(vals) if len(vals)>1 else 0.
        if std<=0:
            continue
        w=2*len(down)*len(erect)/(len(down)+len(erect))
        diff=(statistics.mean(rows[i][axis] for i in down)-
              statistics.mean(rows[i][axis] for i in erect))/std
        total+=w*diff;weight+=w
    return (total/weight if weight>0 else None),weight

def compare(rows,iterations=PERMUTATIONS,seed=SEED):
    selected=[x for x in rows if x["botanical_U_D_state"] in ("U","D")]
    counts=Counter(r["botanical_U_D_state"] for r in selected)
    all_nuclear=all(x["nuclear_real120k_BioSample"] for x in rows)
    desc={
       "n_original_taxa_in_screen":len(rows),
       "n_authority_scored":len(selected),
       "n_D_nodding":counts["D"],
       "n_U_erect":counts["U"],
       "n_unscored":len(rows)-len(selected),
       "n_in_exact_real111_public_nuclear_sample_panel":sum(
            r["nuclear_real120k_BioSample"] for r in rows),
       "n_with_WCVP_source_alias_collision":sum(
            r["WCVP_key_contains_multiple_original_source_names"] for r in selected),
       "nucleotide_tree_necessary_for_evolutionary_origins":True,
       "nuclear_sample_voucher_same_authority_not_verified":True,
    }
    if counts["D"]<3 or counts["U"]<5:
        return {**desc,"status":"NO_GROUNDED_U_D_COMPARISON_TOO_FEW_STATES"}
    by_lane=defaultdict(list)
    for i,r in enumerate(selected):
        by_lane[r["independent_botanical_authority_lane"]].append(i)
    states=[r["botanical_U_D_state"] for r in selected]
    lanes={name:indices for name,indices in by_lane.items()
           if 0<sum(states[i]=="D" for i in indices)<len(indices)}
    lane_diagnostics={}
    for name,indices in sorted(by_lane.items()):
        d=[selected[i] for i in indices if states[i]=="D"]
        u=[selected[i] for i in indices if states[i]=="U"]
        lane_diagnostics[name]={
           "D":len(d),"U":len(u),
           "mean_BIO12_D_mm":statistics.mean(x["chelsa_bio12"] for x in d) if d else None,
           "mean_BIO12_U_mm":statistics.mean(x["chelsa_bio12"] for x in u) if u else None,
           "mean_BIO12_D_minus_U_mm":(
             statistics.mean(x["chelsa_bio12"] for x in d)-
             statistics.mean(x["chelsa_bio12"] for x in u)) if d and u else None,
           "n_source_photos_D":sum(x["original_photo_count"] for x in d),
           "n_source_photos_U":sum(x["original_photo_count"] for x in u),
        }
    primary={}
    for axis in AXES:
        d=[x[axis] for x in selected if x["botanical_U_D_state"]=="D"]
        u=[x[axis] for x in selected if x["botanical_U_D_state"]=="U"]
        score,w=statistic(selected,states,axis,lanes)
        if score is None:
            return {**desc,"status":"NO_WITHIN_AUTHORITY_OVERLAP",
                    "authority_lane_diagnostics":lane_diagnostics}
        primary[axis]={
           "mean_D":statistics.mean(d),"mean_U":statistics.mean(u),
           "pooled_descriptive_difference_D_minus_U":statistics.mean(d)-statistics.mean(u),
           "within_authority_standardized_difference":score,
           "effective_within_lane_comparison_weight":w,
        }
    if not isinstance(iterations,int) or iterations<99:
        raise ValueError("Permutation reproducibility requires at least 99 label draws")
    rng=random.Random(seed)
    observed={a:abs(primary[a]["within_authority_standardized_difference"]) for a in AXES}
    exceed={a:0 for a in AXES}
    family_exceed={a:0 for a in AXES}
    for _ in range(iterations):
        assign=states[:]
        for ix in lanes.values():
            shuffled=[states[i] for i in ix]
            rng.shuffle(shuffled)
            for i,s in zip(ix,shuffled):
                assign[i]=s
        shuffled_effects={}
        for a in AXES:
            effect,_=statistic(selected,assign,a,lanes)
            shuffled_effects[a]=abs(effect)
        highest=max(shuffled_effects.values())
        for a in AXES:
            exceed[a]+=int(shuffled_effects[a]>=observed[a]-1e-12)
            family_exceed[a]+=int(highest>=observed[a]-1e-12)
    for a in AXES:
        primary[a].update({
           "permutation_two_sided_p":(exceed[a]+1)/(iterations+1),
           "maxT_three_climate_axes_adjusted_p":(family_exceed[a]+1)/(iterations+1),
           "null_conditioning":"shuffle authentic botanical U/D only within authority_source_lane",
        })
    return {
        **desc,"status":"EMPIRICAL_OBSERVATIONAL_BOTANICAL_U_D_VS_CLIMATE_NOT_PHYLOGENETIC_CAUSAL",
        "permutation_draws":iterations,
        "random_seed":seed,
        "axes":{LABELS[a]:primary[a] for a in AXES},
        "authority_lane_diagnostics":lane_diagnostics,
        "authority_lanes_with_both_U_D":len(lanes),
        "geographic_autocorrelation_and_lineage_phylogeny_not_adjusted":True,
        "other_185_original_source_name_phenotypes_not_fabricated":True,
        "OU_BM_and_adaptive_radiation_not_inferred":True,
    }

def analyze(photo,authority,draws=PERMUTATIONS):
    joined=merge(photo,authority)
    groups={
       "all_251_authority_scored":joined,
       "original_at_least5_photos":[r for r in joined if r["original_photo_count"]>=5],
       "original_at_least10_photos":[r for r in joined if r["original_photo_count"]>=10],
       "same_real111_public_SRA_taxa":[r for r in joined if r["nuclear_real120k_BioSample"]],
       "single_source_WCVP_key_only":[r for r in joined if not r["WCVP_key_contains_multiple_original_source_names"]],
    }
    results={name:compare(rows,iterations=draws,seed=SEED)
             for name,rows in groups.items()}
    return joined,{
        "schema":"azami_v4_full251_independent_botanical_UD_present_niche_stratified_permutation_20261023",
        "source":"original 170066 observed S1 image taxa + exact independent NA63/NMNS/EAzami/China botanical literature states",
        "all_original_accepted_WCVP_keys":ACCEPTED,
        "all_original_source_name_units":TOTAL,
        "original_source_photos":IMAGES,
        "known_independent_botanical_U_D_source_names":SCORING["U"]+SCORING["D"],
        "unscored_or_ambiguous_original_taxa":SCORING["?"],
        "n_source_real120k_nuclear_taxa":NUCLEAR,
        "phylogenetic_regimes_and_orientation_gain_count_fitted":False,
        "cases":results,
    }

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--original251",required=True,type=Path)
    p.add_argument("--botanical251",required=True,type=Path)
    p.add_argument("--out-dir",required=True,type=Path)
    p.add_argument("--permutations",type=int,default=PERMUTATIONS)
    a=p.parse_args()
    source=csv_rows(a.original251);botanical=csv_rows(a.botanical251)
    combined,report=analyze(source,botanical,draws=a.permutations)
    report["sources_sha256"]={
       a.original251.name:hashlib.sha256(a.original251.read_bytes()).hexdigest(),
       a.botanical251.name:hashlib.sha256(a.botanical251.read_bytes()).hexdigest(),
    }
    a.out_dir.mkdir(parents=True,exist_ok=True)
    with (a.out_dir/"S1_ALL251_ORIGINAL_PHOTO_CLIMATE_BOTANICAL_UD_JOIN.csv").open(
        "w",encoding="utf-8",newline="") as h:
        writer=csv.DictWriter(h,fieldnames=list(combined[0]))
        writer.writeheader();writer.writerows(combined)
    (a.out_dir/"FULL251_BOTANICAL_U_D_CLIMATE_COMPARISONS.json").write_text(
       json.dumps(report,sort_keys=True,indent=2,ensure_ascii=False,allow_nan=False)+"\n",
       encoding="utf-8")
    concise=report["cases"]["all_251_authority_scored"]
    print("T2_FULL251_REAL_BOTANICAL_U_D_RAIN_BIO12_BIO15="+json.dumps({
         "all":report["all_original_source_name_units"],
         "scored":report["known_independent_botanical_U_D_source_names"],
         "cases":{name:{"scored":case["n_authority_scored"],
                       "down":case["n_D_nodding"],
                       "axes":case.get("axes",{})}
                  for name,case in report["cases"].items()},
        },sort_keys=True,ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()
