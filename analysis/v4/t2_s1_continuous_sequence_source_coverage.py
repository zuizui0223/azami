"""S1 continuous PHOTO-angle T2 original nuclear sequence coverage, not binary U/D.

Original source S1 251 camera-relative image head angles and five CHELSA variables.
Actual GenBank ITS109 and technically posthoc (original111 gate FAILED)
Q20-110 sources are independently SHA-validated with the existing V4 parser.
Union161 is a *taxon-name coverage union*, never a fabricated nuclear tree.
"""
import argparse,csv,json
from pathlib import Path
import numpy as np
from analysis.v4.t2_110_multilocus_vs_109_its_photo_niche import (
  EXPECTED_S1_SHA,EXPECTED_ML108_SHA,EXPECTED_ITS_SHA,verify_sources)
from analysis.v4.t2_s1_109_vs251_phylogenetic_sampling_explanation import (
  FIVE,RESPONSES,fit_equal_species_5axis_HC3)

COUNTS={
 "all251":(251,170066,63),
 "its109":(109,162146,16),
 "q20_posthoc110":(110,132762,19),
 "source_union161":(161,165296,27),
 "shared58":(58,129612,8),
 "its_only51":(51,32534,8),
 "q20_only52":(52,3150,11),
 "missing90":(90,4770,36),
}
def prepare(original,technical,receipt,nuclear_tree,its_tree):
    by,t1,t2,mapping,same,_ = verify_sources(original,technical,receipt,nuclear_tree,its_tree)
    its={t.name for t in t2.get_terminals()}
    q=set(mapping.values())
    both=its&q
    universe=set(by)
    groups={
      "all251":universe,"its109":its,"q20_posthoc110":q,
      "source_union161":its|q,"shared58":both,
      "its_only51":its-q,"q20_only52":q-its,
      "missing90":universe-(its|q)}
    if len(both)!=58 or both!=same or (its|q)!=groups["source_union161"]:
        raise ValueError("Original S1 source name intersection changed")
    rainfall=np.array([float(by[t]["chelsa_bio12"]) for t in sorted(universe)])
    threshold=float(np.quantile(rainfall,.75))
    wet={t for t in universe if float(by[t]["chelsa_bio12"])>=threshold}
    if len(wet)!=63 or not np.isclose(threshold,1834.900234745,atol=.005):
        raise ValueError("The full S1 wet quartile changed")
    return by,groups,wet,threshold

def slope_raw(raw,y):
    if len(y)<15:return float("nan")
    scale=raw.std(axis=0)
    if (scale<=1e-9).any():return float("nan")
    z=(raw-raw.mean(axis=0))/scale
    z[:,2]=(raw[:,2]-raw[:,2].mean())/100
    X=np.column_stack([np.ones(len(y)),z])
    if np.linalg.matrix_rank(X)!=6:return float("nan")
    return float(np.linalg.lstsq(X,y,rcond=None)[0][3])

def paired_bootstrap(by,its,union,*,draws=2000,seed=20261024):
    keys=sorted(by)
    raw=np.array([[float(by[k][v]) for v in FIVE] for k in keys])
    y=np.array([float(by[k]["mean_image_angle_deg"]) for k in keys])
    i_mask=np.array([k in its for k in keys])
    u_mask=np.array([k in union for k in keys])
    rng=np.random.default_rng(seed)
    differences=[]
    for _ in range(draws):
        sample=rng.integers(0,len(keys),size=len(keys))
        a=sample[i_mask[sample]]
        b=sample[u_mask[sample]]
        d=slope_raw(raw[b],y[b])-slope_raw(raw[a],y[a])
        if not np.isfinite(d):raise ValueError("Invalid paired taxon-name bootstrap")
        differences.append(d)
    return {"draws":draws,"seed":seed,
      "union_minus_nested_its_slope_q025_median_q975":
        [float(x) for x in np.quantile(differences,[.025,.5,.975])],
      "source_taxon_resampling_not_phylogenetic_CI":True}

def analyze(original,technical,receipt,nuclear_tree,its_tree,*,draws=2000):
    by,groups,wet,threshold=prepare(original,technical,receipt,nuclear_tree,its_tree)
    cohort={}
    for label,names in groups.items():
        rows=[by[t] for t in sorted(names)]
        nphoto=sum(int(r["original_photo_count"]) for r in rows)
        nwet=len(names&wet)
        if (len(rows),nphoto,nwet)!=COUNTS[label]:
            raise ValueError("Strict all-source photo/name/wet records changed at "+label)
        models={}
        for endpoint in RESPONSES:
            result=fit_equal_species_5axis_HC3(rows,endpoint)
            models[endpoint]={"slope_deg_per100mm_BIO12":
               result["BIO12_effect_original_photo_angle_degrees_per100mm"],
                "HC3_95CI":result["HC3_95percent_normal_ci"]}
        cohort[label]={"n_source_names":len(names),"n_original_S1_photo_observations":nphoto,
          "wettest_original_251_quartile_names":nwet,
          "mean_current_BIO12_mm":float(np.mean([float(r["chelsa_bio12"]) for r in rows])),
          "angle_BIO12_equal_taxon_5climate_HC3":models}
    novel=sorted(groups["q20_only52"]&wet)
    if len(novel)!=11 or len(groups["missing90"]&wet)!=36:
        raise ValueError("Source wet addition not original 11/36")
    ledger=[{"original_name":name,"accepted_key":by[name]["accepted_key"],
       "n_S1_original_photos":int(by[name]["original_photo_count"]),
       "S1_native_BIO12_mm":float(by[name]["chelsa_bio12"]),
       "S1_image_angle_mean_deg":float(by[name]["mean_image_angle_deg"]),
       "nuclear_ITS109":name in groups["its109"],
       "nuclear_Q20_110_posthoc_original111_gate_FAILED":name in groups["q20_posthoc110"],
       "top_25percent_BIO12_in_original251":name in wet,
       "not_a_botanical_erect_nodding_label":True}
       for name in sorted(by)]
    return {
      "status":"REAL_S1_PHOTO_CONTINUOUS_SOURCE_COVERAGE_NOT_NUCLEAR_SPECIES_TREE",
      "schema":"T2_S1_photo_real109_110_union251_wet_quartile_20261024_v1",
      "original_full_S1_251_taxon_wet_cutoff_mm":threshold,
      "original_full_S1_251_5_climates":list(FIVE),
      "cohorts":cohort,"new_wet_taxa_on_Q20_candidates_not_in_ITS109":novel,
      "unmatched_high_rainfall_source_names":len(groups["missing90"]&wet),
      "paired_source_selection_bootstrap":paired_bootstrap(
           by,groups["its109"],groups["source_union161"],draws=draws),
      "sources":{"S1_sha256":EXPECTED_S1_SHA,
                 "Q20_tree_sha256":EXPECTED_ML108_SHA,
                 "ITS_tree_sha256":EXPECTED_ITS_SHA},
      "not_an_inferred_161_tip_tree":True,
      "Q20_original111_primary_QC_passed":False,
      "no_botanical_U_D_or_Mk":True,
      "S1_camera_vertical_not_gravity":True,
      "not_original_observation_level_S1_model":True,
      "not_evidence_for_independent_adaptive_evolution":True,
    },ledger

def main():
    p=argparse.ArgumentParser()
    for name in ("original251","source110","receipt","Q20-tree","ITS-tree","out","ledger"):
        p.add_argument("--"+name,required=True,type=Path)
    p.add_argument("--bootstraps",type=int,default=2000)
    a=p.parse_args()
    r,ledger=analyze(a.original251,a.source110,a.receipt,a.Q20_tree,a.ITS_tree,
                     draws=a.bootstraps)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(r,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    with a.ledger.open("w",newline="",encoding="utf-8") as h:
        w=csv.DictWriter(h,fieldnames=list(ledger[0]))
        w.writeheader();w.writerows(ledger)
    print("V4_T2_CONTINUOUS_REAL_SEQUENCE_WET_UNION="+json.dumps({
        k:{"n":v["n_source_names"],"wet":v["wettest_original_251_quartile_names"],
            "BIO12_angle_slope":v["angle_BIO12_equal_taxon_5climate_HC3"]["mean_image_angle_deg"]["slope_deg_per100mm_BIO12"]}
        for k,v in r["cohorts"].items()},sort_keys=True))
if __name__=="__main__":
    main()
