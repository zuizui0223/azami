"""S1/S2-compatible 22-trait inventory -> genuine nuclear branchwise climate tests.

Same original 313350 source-native observations, frozen 251-source S1 identity,
endpoint-specific measurable-photo cohorts and four S1/S2 climate pixels.
Only the genuine NA63 and independent Japan38 176-locus trees are primary.
NEVER read source rows from a public artifact: location/photo/user information
stays inside GitHub Actions' authorised private workspace. Output is aggregate.

S2 introduced-range extension is an *availability gate*, not a fabricated
colour/shape remeasurement from the angle-only S2 four-climate source file.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.v4 import s1s2_multitrait_branch_core as core
from analysis.v4.conditional_ancestry_export import reconstruct
from analysis.v4 import s1_s2_real_bio18_bio05_rsds_elevation as climates
from analysis.v4 import t2_s1_continuous_genuine_nuclear_two_tree_clustering as dna
from analysis.v4 import t2_true_nuclear_S1_continuous_eazami_brownian_fourclimate as true_nuclear
from analysis.v4 import t2_world_s1_fouraxis_eazami_brownian_branch_regime as world

SOURCE_N=313350
SOURCE_ACCEPTED=306
FROZEN_S1_ANGLE_N=170066
FROZEN_S1_NEWFOUR_N=169020
PAIRS=(
    ("presentation", "flower_colour"),
    ("presentation", "whole_head_shape"),
    ("flower_colour", "whole_head_shape"),
    ("whole_head_shape", "involucre_body"),
    ("involucre_body", "phyllary_contour_proxies"),
)


def sha(path):
    h=hashlib.sha256()
    with open(path,"rb") as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b""):
            h.update(chunk)
    return h.hexdigest()


def load_sources(table:Path,reference:Path,four:Path,registry:Path):
    # Every externally held aggregate source is verified by the EXISTING
    # certified 2026-10-24 source reader, not inferred from a filename.
    by,_,old=world.source_panels(reference,four)
    measured=pd.read_csv(table,dtype={"obs_id":str,"accepted_key":str},low_memory=False)
    if (len(measured),measured.accepted_key.nunique())!=(SOURCE_N,SOURCE_ACCEPTED):
        raise ValueError("S1 original accepted native cohort 313350/306 changed")
    if not measured.obs_id.is_unique or measured.obs_id.isna().any():
        raise ValueError("Repeated original observation/source identity")
    req={"obs_id","accepted_key","taxon_name","latitude","longitude",*core.FOUR[:0],
         "chelsa_rsds_mean","chelsa_bio01","chelsa_bio04","chelsa_bio12","chelsa_bio15"}
    if not req <= set(measured):raise ValueError("Native source location/old climate contract missing")
    old_frame=pd.DataFrame(old)
    old_frame["accepted_key"]=old_frame.accepted_key.astype(str)
    original_photo=climates.proper_s1(measured,old_frame)
    if len(original_photo)!=FROZEN_S1_ANGLE_N:
        raise ValueError("S1 original photograph/angle baseline not reproduced")
    selection=core.selection_registry(json.loads(registry.read_text(encoding="utf-8")))
    for name in selection["measured"]:
        if name not in measured or name+"__eligible" not in measured:
            raise ValueError("Original measured endpoint unavailable: "+name)
    measured["source_name"]=measured.taxon_name.astype(str).str.replace(" ","_",regex=False).str.strip()
    names={r["original_source_name"]:str(r["accepted_key"]) for r in old}
    matched=measured.source_name.isin(names)
    if (measured.loc[matched,"accepted_key"] != measured.loc[matched,"source_name"].map(names)).any():
        raise ValueError("Original source species and accepted taxonomic key inconsistent")
    return measured,old_frame,original_photo,pd.DataFrame(list(by.values())),selection


def add_physical_four_environment(frame,source_root):
    # Same pinned-source pixel sampler and physical unit decoders as S1/S2 angle.
    numeric=pd.to_numeric
    lat=numeric(frame.latitude,errors="coerce").to_numpy(float)
    lon=numeric(frame.longitude,errors="coerce").to_numpy(float)
    files={
        "chelsa_bio18":"CHELSA_bio18_1981-2010_V.2.1.tif",
        "chelsa_bio05":"CHELSA_bio05_1981-2010_V.2.1.tif",
        "elevation_m":"wc2.1_30s_elev.tif",
    }
    out=frame.copy()
    raster_receipt={}
    for env,filename in files.items():
        out[env],raster_receipt[env]=climates.raster_sample(
            source_root/filename,lat,lon)
    out["chelsa_rsds_mean"]=numeric(out["chelsa_rsds_mean"],errors="coerce")
    return out, raster_receipt


def verify_original_angle_four(frame,original_photo,four):
    # Detect a silent change of geographic pixels, units or eligible source
    # photos by reproducing the *previously frozen* 169020-angle census.
    df=frame.set_index("obs_id",verify_integrity=True).loc[original_photo.obs_id].reset_index()
    good=np.isfinite(df[list(core.FOUR)].to_numpy(dtype=float)).all(axis=1)
    selected=df.loc[good].copy()
    if len(selected)!=FROZEN_S1_NEWFOUR_N:
        raise ValueError("Previously frozen S1 original four-pixel 169020 gate changed")
    reference=four.set_index("original_source_name")
    if len(reference)!=251 or not reference.index.is_unique:
        raise ValueError("Frozen original source name/physical four-axis roster invalid")
    names=set(reference.index)
    if set(selected.source_name)!=names:
        raise ValueError("The same 251 original S1 source names no longer have four pixels")
    group=selected.groupby("source_name",sort=True)
    counts=group.size()
    for name,n in counts.items():
        if int(n)!=int(reference.at[name,"n_photo_complete_four_axis"]):
            raise ValueError("Original S1 four-climate photo count drift: "+name)
    for k in (*core.FOUR,"orientation_image_vertical_angle"):
        expect="mean_image_angle_deg" if k=="orientation_image_vertical_angle" else k
        observed=group[k].mean().sort_index().to_numpy(float)
        target=pd.to_numeric(reference.loc[counts.index,expect],errors="raise").to_numpy(float)
        if not np.allclose(observed,target,atol=0.001,rtol=0):
            raise ValueError("Original S1 angle/climate pixel mean drift: "+k)
    return {"n_original_S1_angle_photos":len(original_photo),
            "four_pixel_complete":len(selected),"source_names":len(names),
            "result":"PASS_SAME_S1_SITE_PIXEL_SOURCE_NOT_NEW_ENVIRONMENT_PROXY"}


def observed_tip_means(frame:pd.DataFrame, endpoints:tuple[str,...]|list[str],
                       accepted_reference:dict[str,str],*,n_min:int=1)->pd.DataFrame:
    """Every multivariate trait group requires the same source photographs."""
    fields=tuple(endpoints)
    if not fields or n_min<1:raise ValueError("No source measurement endpoints")
    keep=frame.source_name.isin(accepted_reference)
    for endpoint in fields:
        x=pd.to_numeric(frame[endpoint],errors="coerce")
        keep &= frame[endpoint+"__eligible"].eq(True) & np.isfinite(x)
    keep &= np.isfinite(frame[list(core.FOUR)].to_numpy(dtype=float)).all(axis=1)
    sample=frame.loc[keep,["source_name","accepted_key",*fields,*core.FOUR]].copy()
    if sample.empty:
        return pd.DataFrame(columns=["source_name","n_original_observations",*fields,*core.FOUR]).set_index("source_name")
    if (sample.accepted_key.astype(str)!=sample.source_name.map(accepted_reference)).any():
        raise ValueError("Multiple biological source names pooled via an accepted taxon")
    for col in fields:
        sample[col]=pd.to_numeric(sample[col],errors="raise")
    out=sample.groupby("source_name",sort=True).agg(
         n_original_observations=("accepted_key","size"),
         **{c:(c,"mean") for c in (*fields,*core.FOUR)})
    return out.loc[out.n_original_observations>=n_min]


def true_nuclear_panels(reference:Path):
    source=[dna.upstream(dna.NA_REPO,dna.NA_REF,dna.NA_TREE),
            dna.upstream(dna.J_REPO,dna.J_REF,dna.J_TREE),
            dna.upstream(dna.J_REPO,dna.J_REF,dna.J_MAP),
            dna.upstream(dna.J_REPO,dna.J_REF,dna.J_MANIFEST)]
    _,tree_panels,provenance=dna.read_sources(reference,*source)
    out={}
    for label,(raw_tree,mapping) in tree_panels.items():
        names=sorted(mapping.values())
        real=true_nuclear.independent_tree(raw_tree,mapping)
        out[label]={"tree":real,"names":names}
    if {k:len(v["names"]) for k,v in out.items()}!={
        "NA42_genuine_published_nuclear":42,"Japan12_genuine_176_loci_nuclear":12}:
        raise ValueError("Published source-verified authentic nuclear coverage drift")
    return out,provenance


def _test_unit(panel_label,treeinfo,measures,cols,scope,nperm,gate,unit):
    names=sorted(set(treeinfo["names"]) & set(measures.index))
    receipt={"unit":unit,"n_tip":len(names),"min_original_observations":gate,
             "scope":scope,"panel":panel_label,
             "n_original_observations":int(measures.loc[names,"n_original_observations"].sum()) if names else 0}
    if len(names)<core.MIN_TIPS:
        return {**receipt,"status":"STOP_INSUFFICIENT_AUTHENTIC_NUCLEAR_TIPS"}
    tree=world.pruned(treeinfo["tree"],names)
    D,internal=core.branch_operator(tree,names)
    if scope=="internal_edges":D=D[internal]
    if D.shape[0]<3:
        return {**receipt,"status":"STOP_FEW_INDEPENDENT_INTERNAL_EDGES"}
    trait=measures.loc[names,list(cols)].to_numpy(float)
    climate=measures.loc[names,list(core.FOUR)].to_numpy(float)
    try:
        axes=core.alignment_test(D,climate,trait,nperm=nperm,
             seed=core.stable_seed("|".join((panel_label,unit,scope,str(gate)))))
    except ValueError as exc:
        return {**receipt,"status":"STOP_UNIDENTIFIABLE_REAL_TRAIT_OR_CLIMATE",
                "reason":str(exc)}
    output={**receipt,"status":"COMPUTED_REAL_PHOTO_CONTINUOUS_NUCLEAR_BROWNIAN",
            "n_placement_scenarios":1,"n_tip_shuffles":nperm,
            "edge_count_scored":int(D.shape[0]),"axes":axes}
    if scope=="all_edges":
        output["conditional_ancestry"]=reconstruct(
            tree,names,np.column_stack((trait,climate)),[*cols,*core.FOUR])
    return output


def _test_coscore(panel_label,treeinfo,means,group_a,group_b,nperm,scope,gate):
    names=sorted(set(means.index) & set(treeinfo["names"]))
    unit="+".join(group_a)+"__AND__"+"+".join(group_b)
    row={"unit":unit,"panel":panel_label,"scope":scope,"min_original_observations":gate,
         "n_tip":len(names)}
    if len(names)<core.MIN_TIPS:
        return {**row,"status":"STOP_INSUFFICIENT_SAME_PHOTO_MULTITRAIT_TIPS"}
    t=world.pruned(treeinfo["tree"],names)
    D,internal=core.branch_operator(t,names)
    if scope=="internal_edges":D=D[internal]
    if len(D)<3:return {**row,"status":"STOP_FEW_INTERNAL_EDGES"}
    try:
        test=core.co_shift_test(D,means.loc[names,list(group_a)].to_numpy(float),
            means.loc[names,list(group_b)].to_numpy(float),nperm=nperm,
            seed=core.stable_seed("coshift:"+panel_label+unit+scope+str(gate)))
    except ValueError as exc:
        return {**row,"status":"STOP_UNIDENTIFIABLE_PHOTO_MODULE",
                "reason":str(exc)}
    return {**row,"status":"COMPUTED_SAME_PHOTO_MODULE_COSHIFT","n_tip_shuffles":nperm,**test}


def compute(frame,source251,four251,registry,*,reference:Path,nperm:int=9999):
    if nperm<999:raise ValueError("Cannot downgrade frozen permutation resolution")
    reference_keys={r["original_source_name"]:str(r["accepted_key"]) for r in source251.to_dict("records")}
    trees,treesrc=true_nuclear_panels(reference)
    units=core.tests_from_registry(registry)
    results=[]; coverage=[]; pairs=[]
    for unit,trait_fields in units.items():
        available=observed_tip_means(frame,trait_fields,reference_keys,n_min=1)
        coverage.append({"unit":unit,"n_source_taxa_4env":len(available),
            "n_original_source_observations":int(available.n_original_observations.sum()),
            "n_source_taxa_photo10":int((available.n_original_observations>=10).sum())})
        for panel_name,info in trees.items():
            for gate in (1,10):
                sub=available.loc[available.n_original_observations>=gate]
                for scope in ("all_edges","internal_edges"):
                    results.append(_test_unit(panel_name,info,sub,trait_fields,scope,nperm,gate,unit))
    # Holm over all eligible endpoint + block x climate tests in each
    # panel/eligibility/scope family, not one cherry-picked colour or axis.
    for panel_name in trees:
        for gate in (1,10):
            for scope in ("all_edges","internal_edges"):
                rows=[row for row in results if row["panel"]==panel_name and
                      row["min_original_observations"]==gate and row["scope"]==scope and
                      row["status"].startswith("COMPUTED")]
                flat=[(r["unit"]+"|"+axis,r["axes"][axis]["conditional_tip_shuffle_fraction"])
                     for r in rows for axis in core.FOUR]
                adjusted=core.holm_family(flat)
                for row in rows:
                    row["family_resolution"]=core.permutation_resolution(nperm,len(flat))
                    for axis in core.FOUR:
                        row["axes"][axis]["holm_family_fraction"]=adjusted[row["unit"]+"|"+axis]
    # Co-shift pair co-observation is formed BEFORE any taxon averaging.
    for a,b in PAIRS:
        if a not in registry["groups"] or b not in registry["groups"]:
            raise ValueError("Source multitraits changed under pair hypothesis")
        cols_a=registry["groups"][a];cols_b=registry["groups"][b]
        common_cols=tuple(dict.fromkeys((*cols_a,*cols_b)))
        matched=observed_tip_means(frame,common_cols,reference_keys,n_min=1)
        for panel_name,info in trees.items():
            for gate in (1,10):
                sub=matched.loc[matched.n_original_observations>=gate]
                for scope in ("all_edges","internal_edges"):
                    row=_test_coscore(panel_name,info,sub,cols_a,cols_b,nperm,scope,gate)
                    row["group_a"]=a;row["group_b"]=b
                    row["strict_same_original_observations_for_both_modules"]=True
                    pairs.append(row)
    # Pairwise family is separate from the endpoint x climate family, but
    # all proposed cross-module pairs are jointly corrected by Holm.
    for panel_name in trees:
        for gate in (1,10):
            for scope in ("all_edges","internal_edges"):
                selected=[r for r in pairs if r["panel"]==panel_name and
                         r["min_original_observations"]==gate and r["scope"]==scope
                         and r["status"]=="COMPUTED_SAME_PHOTO_MODULE_COSHIFT"]
                adjusted=core.holm_family([(r["unit"],r["conditional_tip_shuffle_fraction"])
                                           for r in selected])
                for r in selected:r["holm_pair_family_fraction"]=adjusted[r["unit"]]
    return {"version":"s1s2_multitrait_true_nuclear_branch_four_20261025_v2",
       "status":"SOURCE_VERIFIED_S1_PHOTO_BROWNIAN_MULTITRAIT_NOT_BOTANICAL_CTMC",
       "source_scope":"Original 313350 geolocated native observations / 306 accepted keys, filtered to frozen 251 original S1 source-name concept universe; exact 4 site climate raster pixels",
       "axes":list(core.FOUR),"n_registered_image_endpoints":27,
       "n_measured_image_endpoints":22,"n_unmeasured":5,
       "n_technical_only_not_ecological":4,
       "no_hue_double_testing":True,
       "no_photo_angle_thresholded_as_botanical_U_D":True,
       "no_gloss_interpreted_as_stickiness":True,
       "no_bract_contour_called_individual_phyllary_posture":True,
       "n_authenticated_genuine_multilocus_species_trees":2,
       "n_authentic_matched_tip_union":54,
       "NA_and_Japan_never_joined_into_invented_global_tree":True,
       "source_tree_provenance":treesrc,
       "hierarchy":"same-original-observation endpoint-specific taxon means -> nuclear Brownian branch deltas; no independent biological origin count",
       "unmeasured_endpoints":registry["unmeasured"],
       "technical_only_endpoints":registry["technical_only"],
       "S2_introduced_extension":"NOT_TESTED_UNTIL_SAME_MEASURED_ENDPOINTS_AND_FOUR_EXPOSURES_REJOINED_TO_INTRODUCED_OBSERVATIONS",
       "null":"conditioned original continuous taxon trait tip-shuffles, not fixed-count botanical-state maps",
       "multiple_testing":"Holm across all admitted endpoint/block x four climates, separately per genuine tree, edge scope and photo gate",
       "claim_ceiling":["No transition-time historical climate", "No true U-to-D from photo vertical",
          "No species-level ancestral pigment gain from JPEG L/C/hue", "No demonstrated genetic modularity, adaptation or independent origins"],
       "n_permutations":nperm,"endpoint_coverage":coverage,
       "nuclear_single_and_multivariate_tests":results,
       "joint_photo_module_co_shifts":pairs}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("decoded-native","s1-source251","s1-four251","registry","raster-root","out"):
        p.add_argument("--"+name,type=Path,required=True)
    p.add_argument("--permutations",type=int,default=9999)
    args=p.parse_args()
    frame,reference,photo,four,reg=load_sources(args.decoded_native,args.s1_source251,
                 args.s1_four251,args.registry)
    frame,source=add_physical_four_environment(frame,args.raster_root)
    pixel=verify_original_angle_four(frame,photo,four)
    result=compute(frame,reference,four,reg,reference=args.s1_source251,nperm=args.permutations)
    result["S1_S2_same_pixel_angle_reproduction"]=pixel
    result["physical_raster_sha256"]={k:v["source_raster_sha256"] for k,v in source.items()}
    result["frozen_original_source_sha256"]={"S1_original251":sha(args.s1_source251),
                                             "four_environment_species251":sha(args.s1_four251)}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({"status":result["status"],"n_measured":22,
          "n_run_models":sum(x["status"].startswith("COMPUTED") for x in result["nuclear_single_and_multivariate_tests"]),
          "n_joint_models":sum(x["status"].startswith("COMPUTED") for x in result["joint_photo_module_co_shifts"]),
          "n_blocked_or_unidentifiable":sum(not x["status"].startswith("COMPUTED") for x in result["nuclear_single_and_multivariate_tests"]),
          "same_169020_pixel_source_checked":True},sort_keys=True))

if __name__=="__main__":main()
