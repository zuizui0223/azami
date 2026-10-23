"""Source strictness and permutation invariance for author-scored Cirsium U/D.

Synthetic fixture validates computation; numerical CI uses 251 real original
source-key climate observations and independent botanical authority.
"""
import pytest
from analysis.v4 import t2_full251_independent_ud_hydric_permutation as m


def original251():
    photo=[];trait=[]
    def metadata(i):
        if i<6 or 111<=i<117: return "D"
        if 6<=i<38 or 117<=i<139: return "U"
        return "?"
    for i in range(251):
        key=f"key{(i if i<237 else 237+(i-237)//2):03d}"
        name=f"Cirsium_test{i:03d}"
        n=677 if i<250 else 816
        state=metadata(i)
        old={
            "accepted_key":key,"original_source_name":name,
            "original_photo_count":str(n),
            "mean_image_angle_deg":str(72+i*.003),
            "chelsa_bio12":str(1100+((i*13)%107)+(400 if state=="D" else 0)),
            "chelsa_bio15":str(43+((i*7)%29)+(8 if state=="D" else 0)),
            "chelsa_bio01":str(12+((i*11)%13)-(2 if state=="D" else 0)),
        }
        source=("Siniscalchi2023_NA63_HeadPosition" if i%3==0
                else "nmns_japan38_seed")
        independent={
            "accepted_key":key,"original_source_name":name.replace("_"," "),
            "original_S1_photo_count":str(n),
            "botanical_U_D_state":state,
            "botanical_source_names":source if state!="?" else "",
            "nuclear_real120k_sample_in_111":str(i<111),
        }
        photo.append(old);trait.append(independent)
    return photo,trait


def test_251_real_source_scope_and_111_nuclear_authored_scope():
    photo,trait=original251()
    rows=m.merge(photo,trait)
    assert len(rows)==251
    assert sum(x["original_photo_count"] for x in rows)==170066
    assert len({x["accepted_key"] for x in rows})==244
    assert sum(x["botanical_U_D_state"]=="D" for x in rows)==12
    assert sum(x["botanical_U_D_state"]=="U" for x in rows)==54
    assert sum(x["botanical_U_D_state"]=="?" for x in rows)==185
    assert sum(x["nuclear_real120k_BioSample"] for x in rows)==111
    assert sum(x["botanical_U_D_state"]=="D" and x["nuclear_real120k_BioSample"]
               for x in rows)==6
    assert sum(x["WCVP_key_contains_multiple_original_source_names"] for x in rows)==14


def test_original_authority_stratified_rain_association_and_maxT_adjustment():
    photo,trait=original251()
    rows=m.merge(photo,trait)
    v=m.compare(rows,iterations=199,seed=20261023)
    assert v["n_authority_scored"]==66
    assert v["n_D_nodding"]==12
    assert v["n_U_erect"]==54
    assert v["authority_lanes_with_both_U_D"]==2
    key=m.LABELS["chelsa_bio12"]
    assert v["axes"][key]["pooled_descriptive_difference_D_minus_U"]>300
    assert 0<v["axes"][key]["permutation_two_sided_p"]<=1
    assert v["axes"][key]["maxT_three_climate_axes_adjusted_p"]>=v["axes"][key]["permutation_two_sided_p"]
    assert v["axes"][key]["maxT_three_climate_axes_adjusted_p"]<=1
    assert v["OU_BM_and_adaptive_radiation_not_inferred"]
    assert v["other_185_original_source_name_phenotypes_not_fabricated"]


def test_authority_pooled_names_fail_closed():
    photo,trait=original251()
    trait[0]["original_source_name"]="Cirsium imaginary"
    with pytest.raises(ValueError,match="EXACT source name"):
        m.merge(photo,trait)


def test_unscored_taxa_are_not_zero_posture():
    photo,trait=original251()
    trait[0]["botanical_U_D_state"]="?"
    with pytest.raises(ValueError,match="Botanical scored"):
        m.merge(photo,trait)


def test_all251_and_photo_threshold_cohorts_cannot_invent_nodal_history():
    photo,trait=original251()
    _,report=m.analyze(photo,trait,draws=99)
    assert report["all_original_source_name_units"]==251
    assert report["all_original_accepted_WCVP_keys"]==244
    assert report["known_independent_botanical_U_D_source_names"]==66
    assert report["unscored_or_ambiguous_original_taxa"]==185
    assert not report["phylogenetic_regimes_and_orientation_gain_count_fitted"]
    assert report["cases"]["same_real111_public_SRA_taxa"]["n_D_nodding"]==6
    assert report["cases"]["single_source_WCVP_key_only"]["n_original_taxa_in_screen"]==237


def test_insufficient_photographic_depth_never_produces_hidden_state_imputation():
    photo,trait=original251()
    rows=m.merge(photo,trait)
    d=[x for x in rows if x["botanical_U_D_state"]=="D"][:2]
    u=[x for x in rows if x["botanical_U_D_state"]=="U"][:8]
    result=m.compare(d+u,iterations=99)
    assert result["status"]=="NO_GROUNDED_U_D_COMPARISON_TOO_FEW_STATES"
