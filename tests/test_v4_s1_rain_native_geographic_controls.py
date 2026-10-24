"""Geographic non-causality audit on the original continuous S1 observation source."""
import numpy as np
import pandas as pd
import pytest

from analysis.v4 import s1_rain_native_geographic_controls as sp
from analysis.v4.joint_s1_t2_phylogeny import EXPOSURES


def sample_taxa(n=135):
    rng=np.random.default_rng(321)
    pieces=[]
    for i in range(n):
        sector=i%3
        lon=[-104.,18.,124.][sector]+rng.normal(0,17)
        lat=38+rng.normal(0,13)
        c=rng.normal(size=5)*[5,190,200,20,3]+[12,450,1200,47,11]
        y=57+0.6*(c[2]/100)+rng.normal(0,3.0)+2*(sector==0)
        for k in range(12):
            lat_r=lat+rng.normal(0,1)
            lon_r=lon+rng.normal(0,1)
            row={"accepted_key":f"k{i:03d}","taxon_name":f"Cirsium_sp{i:03d}",
                 "latitude":lat_r,"longitude":lon_r,
                 sp.ENDPOINT:float(y+rng.normal(0,1))}
            row.update(zip(EXPOSURES,c))
            pieces.append(row)
    work=pd.DataFrame(pieces)
    b=sp.geo_features(work.latitude.to_numpy(),work.longitude.to_numpy())
    for j,name in enumerate(sp.GEO_COLS):
        work[name]=b[:,j]
    regions=sp.longitude_sector(work.longitude.to_numpy())
    for i,name in enumerate(sp.SECTOR_NAMES):
        work[f"sector_{name}"]=(regions==i).astype(float)
    work["_cell_lat"]=np.floor((work.latitude+90)/5).astype(int)
    work["_cell_lon"]=np.floor((work.longitude+180)/5).astype(int)
    return work


def test_geo_basis_on_sphere_no_dateline_artifact():
    lat=np.array([0.,0.,90.])
    lon=np.array([-180.,180.,179.])
    z=sp.geo_features(lat,lon)
    assert z.shape==(3,8)
    assert np.allclose(z[0],z[1],atol=1e-12)
    assert sp.longitude_sector([-130.,5.,110.]).tolist()==[0,1,2]
    with pytest.raises(ValueError,match="outside valid degrees"):
        sp.geo_features([91.],[0.])


def test_each_species_gets_one_equal_weight_and_cell_equilibration_is_a_real_sensitivity():
    work=sample_taxa(135)
    photo=sp.taxon_rosters(work,"photo_weighted")
    cell=sp.taxon_rosters(work,"equal_5deg_spatial_cell")
    assert len(photo)==len(cell)==135
    assert photo.n_source_images.sum()==len(work)
    assert cell.n_source_images.sum()==len(work)
    assert all(photo.source_concept_collision==False)
    assert all(cell.source_concept_collision==False)
    assert (cell.n_geo_cells>=1).all()
    assert np.allclose(
      photo[[f"sector_{s}" for s in sp.SECTOR_NAMES]].sum(axis=1),1.)
    assert np.allclose(
      cell[[f"sector_{s}" for s in sp.SECTOR_NAMES]].sum(axis=1),1.)


def test_spatial_controls_five_original_climate_axes_and_hc3_positive_beta():
    work=sample_taxa()
    roster=sp.taxon_rosters(work,"photo_weighted")
    for geographic in sp.GEO_ADJUSTMENTS:
        model=sp.model_report(roster,"image_mean",geographic)
        assert model["status"]=="estimated",model
        bio=model["exposures"]["chelsa_bio12"]
        assert bio["increment"]==100.
        assert bio["estimate_image_degree_per_registered_increment"]>0
        assert bio["HC3_taxon_interval95"][0]>0
        assert model["n_source_taxa"]==135
        assert model["n_source_name_collisions"]==0
    holdouts=sp.directional_sector_holdouts(roster,"image_mean")
    assert set(holdouts["fixed_predefined_sector_holdouts"])==set(sp.SECTOR_NAMES)
    assert all(v["status"]=="estimated"
          for v in holdouts["fixed_predefined_sector_holdouts"].values())


def test_reject_nonidentifiable_spatial_controls_not_fabricate_a_pvalue():
    work=sample_taxa()
    roster=sp.taxon_rosters(work,"photo_weighted")
    roster["chelsa_bio12"]=roster["chelsa_bio01"]*10
    bad=sp.model_report(roster,"image_mean","spherical_second_order")
    assert bad["status"]=="not_identifiable"
    assert "reason" in bad
    with pytest.raises(ValueError,match="Unregistered spatial control"):
        sp.design(roster,"arbitrary_background_map")


def test_unresolved_WCVP_taxon_collision_stays_excluded():
    work=sample_taxa()
    work.loc[work.accepted_key=="k000","taxon_name"]="Cirsium_alias0"
    work.loc[work.accepted_key=="k000","taxon_name"]=np.where(
        (work.accepted_key=="k000").loc[work.accepted_key=="k000"].index%2==0,
        "Cirsium_alias1","Cirsium_alias2")
    roster=sp.taxon_rosters(work,"photo_weighted")
    assert int(roster.source_concept_collision.sum())==1
    assert len(roster[~roster.source_concept_collision])==134


def test_within_each_fixed_longitude_sector_uses_original_five_climates_only():
    native=sample_taxa(135)
    roster=sp.taxon_rosters(native,"photo_weighted")
    result=sp.regional_within_sector_models(roster,"image_mean")
    assert result["sector_selection_predeclared_without_response"] is True
    assert set(result["regional_results"])==set(sp.SECTOR_NAMES)
    for name,segment in result["regional_results"].items():
        assert segment["n_species"]==45
        assert len(segment["BIO12_taxon_mean_yearly_mm_q10_q50_q90"])==3
        assert segment["within_sector_5_CHELSA_HC3_result"]["status"]=="estimated"
        assert segment["within_sector_5_CHELSA_HC3_result"]["n_source_taxa"]==45
        assert segment["geographic_qualification"].startswith("fixed longitude sector")
