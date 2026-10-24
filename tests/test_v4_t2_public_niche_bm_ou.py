"""BM/OU inference cannot bypass original taxon/phylogeny evidence gates."""
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from analysis.v4.t2_public_niche_bm_ou import (
    analyse,fit_BM_OU,read_rooted_trees,shared_paths,
    SOURCE_COUNTS,AXES
)

ROOT=Path("analysis/v4/public_niche_inputs")
DATA=ROOT/"EAzami_JP_TW_CHELSA_species_centroids_source_20261022.csv"
PRIMARY=ROOT/"EAzami_full20_comp1061_primary_tree.nwk"
ALT=ROOT/"EAzami_comp1061_nine_AU_tree_scenarios.nwk"

def test_public_published_source_phylogenetic_provenance_exact():
    frame=pd.read_csv(DATA)
    assert dict(zip(frame.taxon,frame.n_observations))==SOURCE_COUNTS
    assert len(frame)==12
    trees=read_rooted_trees(PRIMARY,ALT)
    assert len(trees)==7
    for name,tree in trees:
        assert tree.rooted
        tip={x.name for x in tree.get_terminals()}
        assert "OUTGROUP_saff" in tip
        assert len(tip)>=20
        # Exact original tree has independently estimated substitution lengths.
        d,m=shared_paths(tree,[
            "Cirsium_brevicaule","Cirsium_irumtiense",
            "Cirsium_japonicum_var_japonicum","Cirsium_kamtschaticum",
            "Cirsium_yezoense","Cirsium_kawakamii","Cirsium_tatakaense"
        ])
        assert np.linalg.eigvalsh(m).min()>0
        assert np.median(d)>0

def test_ou_is_nested_gaussian_extension_of_brownian():
    lengths=np.array([.05,.07,.08,.09,.055,.11,.095,.065])
    cov=np.diag(lengths)
    # Each independent tip originates from an ingroup root; this is the
    # mathematically simplest valid shared-path matrix for ML checks.
    y=np.array([-1,.1,.5,.9,-.6,.8,-.2,-.8])
    for nugget in (0,.05,.25):
        fit=fit_BM_OU(y,lengths,cov,nugget,grid_points=70)
        assert fit["OU_logLik"]+1e-8>=fit["BM_logLik"]
        assert fit["n_taxa"]==8
        assert np.isfinite(fit["OU_minus_BM_AICc"])
        assert fit["AICc_preferred"] in {"BM","OU"}
        assert np.isfinite(fit["WN_independent_tip_AICc"])
        assert fit["AICc_preferred_WN_BM_OU"] in {"WN","BM","OU"}
        assert np.isclose(
            fit["WN_minus_BM_AICc"],
            fit["WN_independent_tip_AICc"]-fit["BM_AICc"],atol=1e-8)
        if fit["AICc_preferred_WN_BM_OU"]=="BM":
            assert fit["BM_AICc"]<=min(
                fit["WN_independent_tip_AICc"],fit["OU_AICc"])+1e-8


def test_real_fdt4_panel_is_not_all_global_taxa():
    result=analyse(DATA,PRIMARY,ALT)
    assert result["status"]=="EXECUTED_PUBLIC_CIRSIUM_EAST_ASIA_BM_VS_OU_LIMITED_NICHE_SCREEN"
    assert result["claim_stop"]["not_global_niche_radiation"]
    assert result["inputs"]["original_Japan_actions_artifact"]==9516784077
    assert result["inputs"]["original_Taiwan_actions_artifact"]==9706696624
    assert "six_AU_alternative_tree_file_sha256" in result["inputs"]
    assert result["claim_stop"]["not_final_Moreyra_2025_350_locus_tree"]
    for panel_name,n in (("japan7",7),("east_asia9",9),("east_asia12",12)):
        case=result["panels"][panel_name]
        assert case["n_taxa"]==n
        assert case["tree_scenarios"]==7
        assert len(case["fits"])==7*len(AXES)*3
        assert case["nAICc_preferred_BM"]+case["nAICc_preferred_OU"]==len(case["fits"])
        assert set(x["axis"] for x in case["fits"])==set(AXES)
        tally=case["three_model_AICc_preference_counts"]
        assert set(tally)=={"WN","BM","OU"}
        assert sum(tally.values())==len(case["fits"])
        assert tally["OU"]==0
        assert all(row["AICc_preferred_WN_BM_OU"] in {"WN","BM","OU"}
                   for row in case["fits"])


def test_taxon_membership_is_stopped_if_centroids_change(tmp_path):
    frame=pd.read_csv(DATA)
    frame.loc[0,"n_observations"]=2
    fake=tmp_path/"tampered.csv"
    frame.to_csv(fake,index=False)
    with pytest.raises(ValueError,match="occurrence counts"):
        analyse(fake,PRIMARY,ALT)
