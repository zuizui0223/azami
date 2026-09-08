import json
import numpy as np
import pytest
from analysis.v3.simulate_crossed_bootstrap import generate,SPEC
from analysis.v3.dependence_resampling import source_partition
from analysis.v3.workflow import ROOT


@pytest.mark.parametrize('scenario',['iid_null','crossed_spatial_null','heterogeneous_slopes','scale_difference'])
def test_synthetic_design_has_explicit_joint_nulls_and_both_spatial_scales(scenario):
    y,x,g,b,c,lat,lon,truth = generate(scenario,2026090841)
    assert y.shape==(len(g),2) and x.shape==(len(g),9) and b.shape==(len(g),10)
    assert np.isfinite(y).all() and np.isfinite(truth).all()
    assert len(np.unique(g))==64
    assert np.all(truth[:,:,:4]==0)
    if scenario=='scale_difference':
        assert np.all(truth[0,:,4]==0) and np.all(truth[2,:,4]!=0)
    else:
        assert np.all(truth[:,:,4]==0) and np.all(truth[2]==0)
    partitions = [source_partition(np.arange(len(g)),g,c,lat,lon,grid_degrees=d) for d in (2,5)]
    assert [len(np.unique(p.spatial_blocks)) for p in partitions]==[32,16]
    for component in set(c):
        i = np.flatnonzero(c==component)
        if len(i)>1:
            np.testing.assert_array_equal(x[i[0]],x[i[1]])
            np.testing.assert_array_equal(y[i[0]],y[i[1]])
    repeated = generate(scenario,2026090841)
    np.testing.assert_array_equal(repeated[0],y)


def test_calibration_slices_cover_every_planned_case_without_admitting_primary_family():
    spec = json.loads(SPEC.read_text())
    assert spec['bootstrap_replicates']==199
    assert spec['replicates_per_scenario']==40
    assert spec['grid_degrees']==[2,5]
    assert spec['taxa']==64 and spec['predictors']==9 and spec['response_coordinates']==2
    assert spec['ecological_fitting_authorized'] is False
    assert '36slot' in spec['qualification_boundary']
    workflow = (ROOT/'.github/workflows/ch1-v3-crossed-bootstrap-calibration.yml').read_text()
    assert 'max-parallel: 4' in workflow and 'cancel-in-progress: false' in workflow
    assert 'start: [0, 5, 10, 15, 20, 25, 30, 35]' in workflow
    assert 'contents: read' in workflow and 'retention-days: 30' in workflow
    assert 'shapely==2.1.2' in workflow  # protected numerical-output helper imports it
    assert '--end "$(( ${{ matrix.start }} + 5 ))"' in workflow
    assert sorted(i for start in range(0,40,5) for i in range(start,start+5))==list(range(40))
