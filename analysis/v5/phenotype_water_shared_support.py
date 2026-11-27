"""Restrict the already fitted water-response model to its frozen joint support.

Same paired observations and coefficient definitions, no new response,
climate threshold or model. Convex hull support is necessary but not
sufficient for comparability; not an experiment or fitness inference.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from analysis.v5.phenotype_water_interaction import (
    ROOT,SPECIES,DEFINITIONS,ENDPOINTS,VARIABLES,RAIN,SEASON,P,
    authority_source,angular_cells,background_mesh,apply_original_five_CHELSA,
    load_three_axis,_background_species_audit,_pca_env_2d,run_species,
    common_scale,prepare,in_hull,supported_design,fit_water,
    summarize_exposure,stack_model,contrast_vectors,
)
OUT=ROOT/"phenotype_water_shared_support"
EXPECTED_SUPPORT={
    (SPECIES[0],"native_core95"):(18,7),
    (SPECIES[1],"native_core95"):(70,13),
    (SPECIES[0],"density_surplus"):(147,13),
    (SPECIES[1],"density_surplus"):(226,15),
}
# From downloaded artifact 15164060490, SHA256
# f7f3659218def525cbd68241efde668cde72901d00ef061d445d9932268f4c03.
# Original row order: each water term, CA/CV/CV-CA.
EXPECTED_FULL={
    ("native_core95","existing_high20"):[-.08570125,.0810082,.16670946,-.33825906,.04512507,.38338412,-.01738041,-.05289619,-.03551578],
    ("native_core95","existing_low20"):[-.02617434,-.04343486,-.01726052,-.24585085,-.04355132,.20229953,.03998656,.02699783,-.01298874],
    ("density_surplus","existing_high20"):[.0048952,.04555906,.04066386,-.11403371,.03186101,.14589472,-.0357998,-.00521266,.03058713],
    ("density_surplus","existing_low20"):[-.07738989,-.09658893,-.01919904,.00369526,-.02117512,-.02487039,.01155506,.03915782,.02760276],
}
PAIR_KEYS=("left","right","blocks","tail_left","tail_right","design")

def shared_water_mask(d,other):
    jj=[RAIN+1,SEASON+1,6]
    kk=[RAIN,SEASON]
    gd=in_hull(other["design"][:,jj],d["design"][:,jj])
    gl=in_hull(other["left"][:,kk],d["left"][:,kk])
    gr=in_hull(other["right"][:,kk],d["right"][:,kk])
    if any(g is None for g in (gd,gl,gr)):
        raise ValueError("Previously nondegenerate water hull is now degenerate")
    return gd&gl&gr

def subset_pairs(d,mask):
    mask=np.asarray(mask,bool)
    assert mask.shape==(len(d["blocks"]),)
    out=d.copy()
    for key in PAIR_KEYS:
        out[key]=d[key][mask]
    return out

def loo_diagnostic(prepared,endpoint):
    # No new hypothesis tests or refitting hyperparameters. Keep all folds.
    contrasts=contrast_vectors()
    if not all(supported_design(d) for d in prepared):
        return [dict(term=name,comparison=sp,status="insufficient_support",
                     valid_folds=0) for name,sp,c in contrasts]
    z,y,blocks=stack_model(prepared,endpoint)
    base=np.linalg.lstsq(z,y,rcond=None)[0]
    folds=[]
    ids=np.unique(blocks)
    for g in ids:
        sub=[subset_pairs(d,d["blocks"]!=g) for d in prepared]
        if not all(supported_design(d) for d in sub):
            continue
        xx,yy,_=stack_model(sub,endpoint)
        folds.append(np.linalg.lstsq(xx,yy,rcond=None)[0])
    out=[]
    for name,sp,c in contrasts:
        row=dict(term=name,comparison=sp,full_coefficient=float(c@base),
                 planned_folds=len(ids),valid_folds=len(folds),
                 status="estimated" if len(folds)==len(ids) else "incomplete_folds")
        if len(folds)==len(ids):
            v=np.asarray(folds)@c
            row.update(loo_min=float(v.min()),loo_max=float(v.max()),
                       same_sign_folds=int((v*(c@base)>0).sum()))
        out.append(row)
    return out

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    auth=authority_source((ROOT/"private/authority350.zip").read_bytes())
    wild=pd.read_csv(ROOT/"private/focal_wild.csv.gz",dtype={"obs_id":str},low_memory=False)
    assert len(wild)==297622
    wild["_cell"]=angular_cells(wild.latitude,wild.longitude)
    wild=wild.sort_values("obs_id").drop_duplicates(["taxon_name","range_state","_cell","user_id"]).copy()
    land=background_mesh(auth)
    wild,land,_=apply_original_five_CHELSA(wild,land,source_meta=Path("analysis/v4/environment_source_metadata.json"),work=ROOT/"private/rasters")
    measured=load_three_axis(ROOT/"private/original-measured.csv.gz")
    wild=wild.drop(columns=[c for c in measured if c!="obs_id"],errors="ignore").merge(measured,on="obs_id",how="left",validate="one_to_one")
    frozen=json.loads(Path("analysis/v4/results/20261021_unified_abcde/S2_TWO_SPECIES_FULL_WCVP_NICHE_VERIFIED_20261021.json").read_text())
    out={}
    for sp,rows in wild.groupby("taxon_name"):
        env,_=_background_species_audit(auth,land.dropna(subset=list(VARIABLES)),sp,land)
        complete=rows.dropna(subset=list(VARIABLES))
        sn=complete.loc[complete.range_state=="native"];si=complete.loc[complete.range_state=="introduced"]
        transform,edges,nat,inv=_pca_env_2d(env[env._range=="native"],env[env._range=="introduced"],sn,si,return_grid=True)
        assert abs(1-abs(nat-inv).sum()/2-frozen["species"][sp]["PCA_env_background_corrected_Schoener_D"])<1e-5
        run_species(sp,sn,si,transform,edges,nat,inv,out)
    coef=[];profiles=[];support=[];stability=[]
    for definition in DEFINITIONS:
        ca=out[(SPECIES[0],definition)]
        cv=out[(SPECIES[1],definition)]
        center,scale=common_scale(ca["native_climate"],cv["native_climate"])
        prepared=[prepare(d,center,scale) for d in (ca,cv)]
        restricted=[]
        for i,d in enumerate(prepared):
            mask=shared_water_mask(d,prepared[1-i])
            expected=EXPECTED_SUPPORT[(SPECIES[i],definition)]
            assert (int(mask.sum()),len(np.unique(d["blocks"][mask])))==expected
            rd=subset_pairs(d,mask)
            restricted.append(rd)
            support.append(dict(species=SPECIES[i],definition=definition,
                all_pairs=len(mask),supported_pairs=int(mask.sum()),
                supported_blocks=expected[1],supported_fraction=float(mask.mean()),
                supported_design=bool(supported_design(rd))))
            profiles.extend(summarize_exposure(SPECIES[i],definition,rd))
        shared=len(set(restricted[0]["blocks"]) & set(restricted[1]["blocks"]))
        for r in support[-2:]:
            r["shared_supported_blocks"]=shared
        for endpoint in ENDPOINTS:
            full=fit_water(prepared,endpoint)
            assert all(r["status"]=="estimated" for r in full)
            np.testing.assert_allclose([r["coefficient"] for r in full],
                EXPECTED_FULL[(definition,endpoint)],rtol=0,atol=1e-6)
            rows=fit_water(restricted,endpoint)
            for original,r in zip(full,rows):
                assert (original["term"],original["comparison"])==(r["term"],r["comparison"])
                coef.append(dict(definition=definition,endpoint=endpoint,
                    full_coefficient=original["coefficient"],
                    full_p_t=original["p_t"],full_p_score=original["p_score"],**r))
            stability.extend(dict(definition=definition,endpoint=endpoint,**r)
                for r in loo_diagnostic(restricted,endpoint))
    d=pd.DataFrame(coef)
    assert len(d)==36
    for p in ("p_t","p_score"):
        d["holm_all36_"+p]=multipletests(d[p].fillna(1),method="holm")[1]
        d["holm_definition18_"+p]=d.groupby("definition")[p].transform(lambda p:multipletests(p.fillna(1),method="holm")[1])
    d.to_csv(OUT/"shared_water_coefficients.csv",index=False)
    pd.DataFrame(profiles).to_csv(OUT/"shared_water_profiles.csv",index=False)
    pd.DataFrame(support).to_csv(OUT/"shared_water_cohort.csv",index=False)
    pd.DataFrame(stability).to_csv(OUT/"shared_water_loo.csv",index=False)
    print("Complete; aggregate outputs only.",flush=True)

if __name__=="__main__":
    main()
