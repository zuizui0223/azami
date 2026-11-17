import numpy as np
from analysis.v4.rf_support_atlas import support_lattice,GRAINS,MIN_PHOTOS_PER_CELL
from analysis.v4.rf_multitrait_predictability import aggregate_endpoint
from analysis.v4.rf_archived_atmosphere import ATMOSPHERE23
from analysis.v4.rf_extended_environment import BIO9
from test_v4_rf_multitrait_predictability import tiny_synthetic


def test_photographic_coverage_grid_is_cross_product_of_predefined_conditions():
    d=tiny_synthetic(n_taxa=30,n_cells=8,reps=6)
    rng=np.random.default_rng(2)
    for name in BIO9:
        d[name] = rng.normal(size=len(d))
    d['chelsa_vpd_mean']=rng.normal(size=len(d))
    d['chelsa_sfcwind_mean']=rng.normal(size=len(d))
    report=support_lattice(d,constructs=('angle','chroma'))
    assert len(report['records'])==2*len(GRAINS)*len(MIN_PHOTOS_PER_CELL)
    assert report['status']=='source_eligibility_counts_only_no_selected_grid'
    assert all('taxon_0' not in str(r) for r in report['records'])
    five=[x for x in report['records'] if x['cell_degrees']==5 and x['minimum_photos_per_cell']==3]
    assert {x['n_taxa'] for x in five}=={30}
    assert {x['n_taxon_cells'] for x in five}=={240}
    assert all(x['eligible_for_grouped_region_holdout'] for x in five)


def test_strict_five_degree_remains_byte_equivalent_in_sampling_structure():
    d=tiny_synthetic(n_taxa=30,n_cells=8,reps=6)
    for name in BIO9: d[name]=2.
    d['chelsa_vpd_mean']=2.;d['chelsa_sfcwind_mean']=3.
    a=aggregate_endpoint(d,'angle',predictors=ATMOSPHERE23)
    b=aggregate_endpoint(d,'angle',predictors=ATMOSPHERE23,degrees=5)
    assert a.equals(b)
    with np.testing.assert_raises_regex(ValueError,'UNDECLARED_SPATIAL_GRAIN'):
        aggregate_endpoint(d,'angle',predictors=ATMOSPHERE23,degrees=3)