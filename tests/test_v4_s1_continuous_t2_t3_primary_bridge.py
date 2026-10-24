"""S1 photo orientation is continuous for every active V4 analysis."""
import pytest
from analysis.v4 import s1_continuous_t2_t3_primary_bridge as p


def t3_record(taxon,tribe,n=12):
    return {
        "taxon_name":taxon,"tribe":tribe,"n_original_photo_observations":str(n),
        "measurement_endpoint":p.ENDPOINT,"camera_vertical_EXIF_aligned":"True",
        "mean_image_angle_deg":"73.52",
        "source_photo_digest":"sha256:sourcephoto","source_climate_digest":"sha256:sourceclimate",
        "original_phenotype_observed":"True",
    }


def test_t3_blocks_without_photos_even_if_flora_and_ITS_tree_exist():
    a=p.validate_t3_continuous_source_manifest([])
    assert a["status"]=="BLOCKED_NO_MEASURED_EXTERNAL_TRIBE_PHOTOS"


def test_t3_rejects_floristic_binary_disguised_as_s1_angle():
    rows=[t3_record("a","Senecioneae"),t3_record("b","Cardueae")]
    rows[0]["measurement_endpoint"]="flora_nodding_U_D"
    with pytest.raises(ValueError,match="categorical"):
        p.validate_t3_continuous_source_manifest(rows)


def test_t3_requires_replication_across_tribes_and_nonbinary_measures():
    one=[t3_record("a","Senecioneae"),t3_record("b","Senecioneae")]
    assert p.validate_t3_continuous_source_manifest(one)["status"].startswith("BLOCKED")
    both=one+[t3_record("c","Cardueae")]
    assert p.validate_t3_continuous_source_manifest(both)["status"]=="ELIGIBLE_FOR_T3_TAXON_COVERAGE_REVIEW"
    both[-1]["mean_image_angle_deg"]="200"
    with pytest.raises(ValueError,match="Invalid"):
        p.validate_t3_continuous_source_manifest(both)


def test_mainline_taxa_climate_and_camera_semantics_are_fixed():
    assert p.ENDPOINT=="orientation_image_vertical_angle"
    assert len(p.FIVE)==5
    assert p.S2_FOCALS==("Cirsium arvense","Cirsium vulgare")
    assert p.EXPECTED_251_SHA256=="47b32972b20564ade163907ef6fdf3eb2a69f8b42c8d87b8cf4ea75ca1682dc1"
