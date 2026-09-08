import numpy as np
import pytest
from analysis.v3.dependence_resampling import source_partition,cohort_partition,crossed_draw


def test_full_source_bridges_survive_missing_module_rows():
    # A--B--C taxon closure runs through two observations later absent from the module.
    source = source_partition(list('abcdef'),list('ABBCDD'),['u','u','v','v','w','z'],
                              [1,11,11,21,31,41],[1,11,11,21,31,41],grid_degrees=2)
    positions,t,s = cohort_partition(source,np.array([0,3,4,5]))
    assert t[0] == t[1]
    assert s[0] == s[1]
    assert t[0] != t[2]
    assert source.observation_ids.flags.writeable is False


@pytest.mark.parametrize('degrees',[2,5])
def test_crossed_draw_preserves_whole_components_across_taxa_and_space(degrees):
    taxa = np.repeat(np.arange(6),5)
    components = np.arange(30).astype(str)
    components[7] = components[0]
    lat = np.tile(np.array([1,11,21,31,41]),6)
    lon = np.tile(np.array([2,12,22,32,42]),6)
    source = source_partition(np.arange(30),taxa,components,lat,lon,grid_degrees=degrees)
    for replicate in range(25):
        indices,labels = crossed_draw(source,np.arange(30),seed=4431,replicate=replicate)
        weights = np.bincount(indices,minlength=30)
        assert weights[0] == weights[7]
        assert len(labels)==len(indices)
        # One copied taxon label can never merge different original taxa.
        for label in np.unique(labels):
            assert len(np.unique(taxa[indices[labels==label]])) == 1
        repeated,second_labels = crossed_draw(source,np.arange(30),seed=4431,replicate=replicate)
        np.testing.assert_array_equal(indices,repeated)
        np.testing.assert_array_equal(labels,second_labels)


def test_coordinate_wrap_and_poles_share_the_same_geographic_cell():
    source = source_partition(range(6),range(6),range(6),[0,0,90,90,-90,-90],
                              [-180,180,-170,170,-170,170],grid_degrees=2)
    assert source.spatial_blocks[0]==source.spatial_blocks[1]
    assert source.spatial_blocks[2]==source.spatial_blocks[3]
    assert source.spatial_blocks[4]==source.spatial_blocks[5]


def test_collapsed_factor_is_not_replaced_with_iid_rows():
    source = source_partition(range(4),range(4),['same']*4,[0,10,20,30],[0,10,20,30],grid_degrees=5)
    with pytest.raises(ValueError,match='collapsed'):
        crossed_draw(source,np.arange(4),seed=1,replicate=0)


@pytest.mark.parametrize('positions',[[0,0],[0,4],[-1,2],[.5,2],[]])
def test_module_positions_cannot_duplicate_or_invent_source_records(positions):
    source = source_partition(range(4),range(4),range(4),[0,10,20,30],[0,10,20,30],grid_degrees=2)
    with pytest.raises(ValueError):
        cohort_partition(source,np.array(positions))


@pytest.mark.parametrize('bad',[[0,np.nan],[0,91]])
def test_nonfinite_or_impossible_source_coordinates_stop(bad):
    with pytest.raises(ValueError):
        source_partition(['a','b'],['A','B'],['u','v'],bad,[1,3],grid_degrees=2)
