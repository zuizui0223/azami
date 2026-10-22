"""All 251 original S1 photo source names x verified 294 nuclear samples x real 241 depth.

The original 170,066 S1 angle observations have 244 WCVP accepted keys but
251 taxon-name-within-accepted-key observational units. Seven keys represent
two source names each. Never merge C. brevicaule and C. irumtiense merely
because WCVP returned the same accepted key. Never label 251 as accepted
biological species, or 294 source BioSamples as a completed nuclear tree.

Real 241 Q20 site-depth evidence only exists for 113 name-matched S1 keys;
coverage is NOT zero on missing samples. EAzami294 original source label
overlap is descriptive voucher candidates, not accepted-name proof.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path
import re

N_SOURCE=251
N_ACCEPTED=244
N_PHOTOS=170066
N_MERGED_KEYS=7
N_Q20_TAXA=113
N_EAZAMI_BIOSAMPLES=294
N_EAZAMI_RUNS=295
FOCAL={"Cirsium_brevicaule","Cirsium_irumtiense"}
SHARED_KEY="206745753"

def read_csv(path):
    with Path(path).open(newline="",encoding="utf-8-sig") as f:
        rows=list(csv.DictReader(f))
    if not rows:
        raise ValueError("Missing source S1, depth or EAzami nuclear input")
    return rows

def name(x):
    return "_".join(str(x or "").strip().split()).replace(" ","_")

def specimen_code_stripped(source):
    """Terminal digits on Moreyra 'tree_code' label, not taxonomic rank."""
    return re.sub(r"_\d+$","",name(source))


def build(original,depth,manifests):
    if (len(original)!=N_SOURCE or
        len({(x["accepted_key"],x["original_source_name"]) for x in original})!=N_SOURCE or
        len({x["accepted_key"] for x in original})!=N_ACCEPTED or
        sum(int(x["original_photo_count"]) for x in original)!=N_PHOTOS):
        raise ValueError("251 source-name /244 accepted key /170066 original image scope drift")
    keys=defaultdict(list)
    for x in original:keys[x["accepted_key"]].append(x)
    if sum(len(v)==2 for v in keys.values())!=N_MERGED_KEYS or any(len(v)>2 for v in keys.values()):
        raise ValueError("Seven original WCVP name-mixed pairs lost")
    if (len(depth)!=N_Q20_TAXA or
        len({x["accepted_key"] for x in depth})!=N_Q20_TAXA):
        raise ValueError("Not all 113 original Moreyra real depth name candidate taxa")
    real={x["accepted_key"]:x for x in depth}
    if len(manifests)!=N_EAZAMI_BIOSAMPLES or \
       len({x["biosample"] for x in manifests})!=N_EAZAMI_BIOSAMPLES:
        raise ValueError("EAzami294 deduplicated biological sample panel missing")
    runs={r for x in manifests for r in x["run_accessions"].split("|") if r}
    if len(runs)!=N_EAZAMI_RUNS:
        raise ValueError("Original 295 independent public accession run set changed")
    samples=defaultdict(list)
    for x in manifests:
        original_name=name(x["analysis_taxon_label"])
        norm=(specimen_code_stripped(original_name)
              if x["source_study"]=="Moreyra2025" else original_name)
        if norm:samples[norm].append(x)
    result=[]
    for x in original:
        key=x["accepted_key"]
        source_name=name(x["original_source_name"])
        collided=len(keys[key])>1
        s=sorted(samples.get(source_name,[]),key=lambda a:a["tip_id"])
        row=real.get(key)
        if collided and row is not None:
            raise ValueError("Source WCVP key multi-name collision assigned one SRA depth record")
        if row is not None:
            sample_name=name(row["source_taxon"])
            if sample_name!=source_name:
                raise ValueError("Wrong real 241 depth source taxon on accepted key")
        result.append({
            "accepted_key":key,
            "S1_original_source_name":source_name,
            "original_S1_photo_count":int(x["original_photo_count"]),
            "S1_image_vertical_angle_deg":float(x["mean_image_angle_deg"]),
            "S1_CHELSA_BIO12_mm":float(x["chelsa_bio12"]),
            "original_WCVP_collapsed_pair":collided,
            "n_EAzami294_name_candidate_samples":len(s),
            "n_EAzami294_name_candidate_Moreyra_samples":
                 sum(y["source_study"]=="Moreyra2025" for y in s),
            "n_EAzami294_name_candidate_Chang_samples":
                 sum(y["source_study"]!="Moreyra2025" for y in s),
            "EAzami294_sample_ids":"|".join(y["tip_id"] for y in s),
            "real_241_site_depth_available":row is not None,
            "Q20_dp3_ge50pct_loci_best_single_run":
                int(row["maximum_per_run_loci_dp3_ge50pct"]) if row else None,
            "Q20_dp3_ge25pct_loci_best_single_run":
                int(row["maximum_per_run_loci_dp3_ge25pct"]) if row else None,
            "missing_site_depth_status":
                "" if row else "not_in_original_125_public_Moreyra_SRA_run_panel",
            "one_to_one_accepted_species_voucher_verified":False,
            "original_294_species_tree_completed":False,
        })
    pairs={z["S1_original_source_name"]:z for z in result
           if z["accepted_key"]==SHARED_KEY}
    if set(pairs)!=FOCAL or len(pairs)!=2:
        raise ValueError("Cirsium brevicaule and irumtiense wrongly collapsed")
    if sum(z["real_241_site_depth_available"] for z in result)!=N_Q20_TAXA:
        raise ValueError("113 original source-identified real 241 depth taxa not conserved")
    summary={
       "schema":"azami_v4_S1_all251_original_source_names_with_EAzami294_real241_coverage_v1",
       "total_original_photo_observations":N_PHOTOS,
       "total_original_S1_WCVP_accepted_keys":N_ACCEPTED,
       "original_S1_source_name_operational_units":N_SOURCE,
       "original_source_mixed_taxon_accepted_keys":N_MERGED_KEYS,
       "original_S1_units_with_real_241_depth":N_Q20_TAXA,
       "original_S1_units_without_real_241_depth":N_SOURCE-N_Q20_TAXA,
       "original_S1_units_with_EAzami294_sample_name_candidate":
          sum(z["n_EAzami294_name_candidate_samples"]>0 for z in result),
       "S1_source_name_matches_are_not_accepted_species_or_voucher_proof":True,
       "C_brevicaule_and_C_irumtiense_separate_source_name_units":{
          k:{"n_photo":z["original_S1_photo_count"],
             "n_eazami294_name_matched_samples":z["n_EAzami294_name_candidate_samples"],
             "n_eazami294_Chang_samples":z["n_EAzami294_name_candidate_Chang_samples"]}
          for k,z in pairs.items()},
       "EAzami_global_nuclear_tree_reconstructed":False,
       "nuclear_missingness_does_not_mean_genetic_absence":True,
       "author241_Q20_read_depth_not_orthology":True,
       "real_species_adaptive_radiation_conclusion":False,
    }
    return result,summary


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--s1-source251",type=Path,required=True)
    p.add_argument("--real241-species",type=Path,required=True)
    p.add_argument("--eazami294-samples",type=Path,required=True)
    p.add_argument("--out-dir",type=Path,required=True)
    a=p.parse_args()
    rows,summary=build(read_csv(a.s1_source251),read_csv(a.real241_species),
                       read_csv(a.eazami294_samples))
    a.out_dir.mkdir(parents=True,exist_ok=True)
    with (a.out_dir/"S1_ALL251_SOURCE_NAME_294_EAZAMI_REAL241_CROSSWALK.csv").open(
        "w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),extrasaction="raise")
        w.writeheader();w.writerows(rows)
    (a.out_dir/"S1_ALL251_REAL241_MATCHING_SUMMARY.json").write_text(
        json.dumps(summary,ensure_ascii=False,sort_keys=True,indent=2)+"\n",
        encoding="utf-8")
    print("S1_ALL251_SOURCE_NAME_REAL241_EAzami294_INTEGRATION="+json.dumps(
          summary,ensure_ascii=False,sort_keys=True),flush=True)


if __name__=="__main__":
    main()
