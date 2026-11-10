"""Why did the BIO10 swap not restore cross-taxon capitulum resemblance?

This source-matched audit checks whether the set of climatic-analogue
taxon pairs and cell pairs changes when BIO05 becomes BIO10, and tests
BIO10 minus BIO05 excess resemblance using *paired* within-taxon
phenotype permutations. No original coordinates, species names,
photographer names, image measurements or pair identifiers are exported.
This is a post-result diagnostic, not a prespecified primary test.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from analysis.v4.s1_climate_trait_hypervolume import geo_aggregate, standardise_cells
from analysis.v4 import s1_h3_climate_analogue_alignment as h
from analysis.v5.bio10_s1 import (
    EXTRA, CLIMATE, N_DRAWS, N_NULL, THRESH, _fixed_taxa,
    verified_photo_source
)

COMPARE=("annual_plus_seasonality","warm_plus_seasonality","warm_bio10_plus_seasonality")


def pair_details(p):
    """All temporary identities stay private; return map for matched comparison."""
    if p is None:
        return {}
    i=p["i"];j=p["j"]
    all_means=p["d"].mean(axis=(1,2))
    analog_means=(p["d"]*p["mask"]).sum(axis=(1,2))/p["analog_n"]
    return {(int(a),int(b)):{
        "mask":m.copy(),
        "adv":float(aa-ab),
        "n_analog":int(nn)
    } for a,b,m,aa,ab,nn in zip(
        i,j,p["mask"],all_means,analog_means,p["analog_n"])}


def selection_summary(a,b):
    """Identity overlap of taxon pairs and of analogous cell pairs, aggregate."""
    aset=set(a);bset=set(b);both=aset&bset;union=aset|bset
    common_cells=0;either_cells=0
    for key in both:
        m1=a[key]["mask"];m2=b[key]["mask"]
        common_cells+=int((m1&m2).sum())
        either_cells+=int((m1|m2).sum())
    return {
        "eligible_taxon_pairs_first":len(aset),
        "eligible_taxon_pairs_second":len(bset),
        "same_taxon_pair_intersection":len(both),
        "taxon_pair_jaccard":float(len(both)/len(union)) if union else None,
        "cell_pair_jaccard_within_shared_taxon_pairs":
            float(common_cells/either_cells) if either_cells else None,
        "same_pair_mean_phenotype_advantage_first":
            float(np.mean([a[k]["adv"] for k in both])) if both else None,
        "same_pair_mean_phenotype_advantage_second":
            float(np.mean([b[k]["adv"] for k in both])) if both else None,
    }


def run(native, decoded, receipt):
    d,proof=verified_photo_source(native,decoded,receipt)
    rows=[];comparisons=[];support_by_degree={}
    for deg in (1,5):
        cells,_=geo_aggregate(d,deg,EXTRA)
        z,_=standardise_cells(cells,EXTRA)
        taxa=_fixed_taxa(z)
        regions=np.array([t["sector"] for t in taxa])
        rng=np.random.default_rng(20261109+deg)
        draws={k:[] for k in COMPARE}
        pair_refs={k:[] for k in COMPARE}
        for rep in range(N_DRAWS):
            picks=[rng.choice(len(t["y"]),8,replace=False) for t in taxa]
            y=np.stack([t["y"][ids] for t,ids in zip(taxa,picks)])
            source={}
            for model in COMPARE:
                x=np.stack([t["x"][model][ids] for t,ids in zip(taxa,picks)])
                pairs,support=h.matched_pairs(x,y,regions,THRESH)
                if support["status"]!="qualified":
                    # Original V5 drops unsupported rarefactions; never
                    # manufacture extra climate-analogue cell pairs.
                    draws[model].append(None)
                    pair_refs[model].append(None)
                    continue
                obs,null=h.observed_and_permutations(
                    pairs,len(taxa),np.random.default_rng(20261119+deg*1000+rep),N_NULL)
                draws[model].append((obs,null))
                pair_refs[model].append(pair_details(pairs))
            for k in range(len(COMPARE)):
                for j in range(k+1,len(COMPARE)):
                    model1,model2=COMPARE[k],COMPARE[j]
                    if pair_refs[model1][-1] is None or pair_refs[model2][-1] is None:
                        continue
                    q=selection_summary(pair_refs[model1][-1],pair_refs[model2][-1])
                    rows.append({
                        "degree":deg,"model_a":model1,"model_b":model2,
                        "draw_number_for_diagnostic_only":rep,
                        **q
                    })
        support_by_degree[str(deg)]={
            model:int(sum(x is not None for x in draws[model])) for model in COMPARE
        }
        for a,b in (
            ("warm_bio10_plus_seasonality","warm_plus_seasonality"),
            ("warm_bio10_plus_seasonality","annual_plus_seasonality")
        ):
            paired=[(x,y) for x,y in zip(draws[a],draws[b])
                    if x is not None and y is not None]
            if len(paired)<12:
                comparisons.append({
                    "degree":deg,"model_new":a,"model_reference":b,
                    "status":"NOT_ESTIMABLE_TOO_FEW_COMMON_QUALIFIED_DRAWS",
                    "n_paired_qualified_draws":len(paired),
                    "min_required":12
                })
                continue
            av=np.array([x[0][0] for x in paired],float)
            bv=np.array([x[1][0] for x in paired],float)
            an=np.mean([x[0][1] for x in paired],axis=0)
            bn=np.mean([x[1][1] for x in paired],axis=0)
            # Under the joint within-taxon null, the same phenotype
            # permutation stream is used for each climate definition.
            diff_obs=float(av.mean()-bv.mean())
            diff_null=an-bn
            excess=float(diff_obs-diff_null.mean())
            p=float((1+np.sum(np.abs(diff_null-diff_null.mean())>=abs(excess)-1e-12))/(N_NULL+1))
            comparisons.append({
                "degree":deg,"model_new":a,"model_reference":b,
                "status":"ESTIMATED_PAIRED_NULL_TEST",
                "n_paired_qualified_draws":len(paired),
                "observed_model_difference":diff_obs,
                "paired_difference_excess_over_null":excess,
                "paired_taxon_permutation_p_exploratory":p,
                "n_paired_permutations":N_NULL,
                "identical_taxon_cell_rare_factions":True
            })
    table=pd.DataFrame(rows)
    agg=[]
    if table.empty:
        return {
            "status":"NO_SHARED_QUALIFIED_RAREFACTION_DRAWS",
            "support_counts":support_by_degree,
            "no_taxon_identity_or_original_points_exported":True,
        }
    for (deg,a,b),g in table.groupby(["degree","model_a","model_b"],sort=True):
        rec={"degree":int(deg),"model_a":a,"model_b":b,"draws":len(g)}
        for field in (
            "eligible_taxon_pairs_first","eligible_taxon_pairs_second",
            "same_taxon_pair_intersection","taxon_pair_jaccard",
            "cell_pair_jaccard_within_shared_taxon_pairs",
            "same_pair_mean_phenotype_advantage_first",
            "same_pair_mean_phenotype_advantage_second"
        ):
            rec[field+"_mean"]=float(g[field].mean())
            rec[field+"_min"]=float(g[field].min())
            rec[field+"_max"]=float(g[field].max())
        agg.append(rec)
    return {
        "status":"POSTHOC_PAIR_SELECTION_DIAGNOSTIC_SUCCESS",
        "original_native_matched_photos":int(len(d)),
        "source_receipt_sha":proof["original_decoded_source_sha"],
        "all_original_photo_cells_and_taxa_shared":True,
        "climate_analogue_threshold":float(THRESH),
        "resolutions":[1,5],
        "selected_taxon_pair_overlap":agg,
        "bio10_minus_other_model_paired_null_tests":comparisons,
        "support_counts":support_by_degree,
        "negative_trait_advantage_not_evidence_for_adaptive_divergence":True,
        "no_taxon_identity_or_original_points_exported":True,
    }
