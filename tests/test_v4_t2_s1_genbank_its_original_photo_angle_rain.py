"""Actual S1 photos, five original climate axes and own nuclear ITS tree contracts."""
from io import StringIO
from copy import deepcopy
import numpy as np
from Bio import Phylo
import pytest
from analysis.v4 import t2_s1_genbank_its_original_photo_angle_rain as m


def source(n=12):
    out={}
    for k in range(n):
        out[f"Cirsium_species{k:02}"]={
            "original_source_name":f"Cirsium_species{k:02}",
            "accepted_key":str(206000000+k),
            "original_photo_count":str(4+k),
            "WCVP_key_collides_source_names":"False",
            "mean_image_angle_deg":str(40+2*k),
            "median_image_angle_deg":str(41+3*k),
            "botanical_gravity_posture_inferred_from_photo_angle":"False",
            "chelsa_bio01":str(8+3*k+3*np.sin(k)),
            "chelsa_bio04":str(400+37*k+19*np.cos(k*2)),
            "chelsa_bio12":str(750+87*k+29*np.sin(k*3)),
            "chelsa_bio15":str(31+4*k+3*np.cos(k*5)),
            "chelsa_rsds_mean":str(17+2*np.sin(k*4)+k/7),
        }
    return out


def test_all_five_original_S1_axes_return_photo_orientation_not_external_flora():
    rows=source()
    tips=list(rows)
    X,environment=m.source_design(rows,tips)
    assert X.shape==(12,6)
    assert environment.shape==(12,5)
    assert np.allclose(X[:,0],1.)
    assert np.isclose(X[:,3],(environment[:,2]-environment[:,2].mean())/100).all()
    y=np.array([float(rows[t]["mean_image_angle_deg"]) for t in tips])
    C=np.eye(12)
    r=m.fit_models(y,X,C)
    assert set(m.EXPOSURES)=={
        "chelsa_bio01","chelsa_bio04","chelsa_bio12","chelsa_bio15","chelsa_rsds_mean"
    }
    assert set(k for k in r if isinstance(r[k],dict))=={
        "IID_equal_species_weight","Pagel_lambda0_diagonal_nonultrametric",
        "fixed_Brownian_lambda1","estimated_Pagel_lambda"
    }
    beta=r["IID_equal_species_weight"]["source_yearly_rain_BIO12_beta_image_degrees_per100mm"]
    assert np.isfinite(beta)
    assert not r["photos_do_not_measure_absolute_gravity_direction"] is False


def test_IID_anglegls_matches_closed_form_equal_species_weight_OLS():
    rows=source()
    tips=list(rows)
    X,_=m.source_design(rows,tips)
    y=np.array([float(rows[t]["mean_image_angle_deg"]) for t in tips])
    expected=np.linalg.lstsq(X,y,rcond=None)[0][3]
    result=m.fit_gls(y,X,np.eye(len(tips)),0,structure="IID_equal_species_weight")
    assert result["photo_derived_angle_is_dependent_variable"]
    assert result["source_yearly_rain_BIO12_beta_image_degrees_per100mm"]==pytest.approx(expected)
    assert result["BIO12_95pct_CI_image_degrees_per100mm"][0] < result["BIO12_95pct_CI_image_degrees_per100mm"][1]


def test_synthetic_real_phylo_covariance_changes_without_photographic_nodding_labels(monkeypatch):
    tree=Phylo.read(StringIO(
      "(((Cirsium_species00:.3,Cirsium_species01:.4):.1,(Cirsium_species02:.3,Cirsium_species03:.5):.1):.7,"
      "((Cirsium_species04:.1,Cirsium_species05:.4):.4,((Cirsium_species06:.2,Cirsium_species07:.25):.25,"
      "((Cirsium_species08:.2,Cirsium_species09:.25):.3,(Cirsium_species10:.2,Cirsium_species11:.4):.3):.1):.1):.2);"),
      "newick")
    monkeypatch.setattr(m,"EXPECTED_ITS_TIPS",12)
    rows=source()
    analysis=m.evaluate_tree(tree,rows,rooting="midpoint_ITS_unrooted")
    assert analysis["taxon_count"]==12
    assert analysis["model_DOES_use_original_170066_image_derived_continuous_angle"]
    assert analysis["model_does_not_use_botanical_flora_U_D_or_species_label_classifier"]
    assert not analysis["min_covariance_eigenvalue"] <= 0
    assert set(analysis["angle_models"])==set(m.RESPONSES)
    original=deepcopy(rows)
    for key in original: original[key]["botanical_U_D_state"]="D"
    twin=m.evaluate_tree(tree,original,rooting="midpoint_ITS_unrooted")
    # Explicit regression test: changing EVERY flora D/U CANNOT change OWN photo outcome.
    assert twin["angle_models"]==analysis["angle_models"]


def test_unrooted_ITS_109_is_exploratory_not_an_authenticated_species_origin():
    assert m.EXPECTED_ITS_TIPS==109
    assert m.EXPECTED_MASTER==251 and m.EXPECTED_WCVP==244
    assert m.EXPECTED_PHOTOS==170066
    assert m.EXPECTED_ITS_REF_SITES==645
    assert m.EXPECTED_ITS_INFORMATIVE==103
    assert "S1 photographic orientation" in m.__doc__
    assert "gravity-calibrated tilt" in m.__doc__
    assert m.ITS_ML_TREE_SHA256=="8fd10d299b4d520ce08f2552993b02d3cc1b10350f0fb68ea24874dbb9ccfe92"


def test_five_axis_exposure_cannot_be_relabelled_or_omitted():
    r=source()
    tips=list(r)
    r[tips[3]]["chelsa_bio01"]="nan"
    with pytest.raises(ValueError,match="Non-finite"):
        m.source_design(r,tips)
    r=source()
    for x in r.values():x["chelsa_bio12"]="200"
    with pytest.raises(ValueError,match="collinear|invariant"):
        m.source_design(r,tips)
