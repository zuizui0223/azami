import copy
import json
import re
from pathlib import Path

import pytest

from analysis.v3 import measurement_chunks as chunks
from analysis.v3.reconciled_stream_input import pilot_input
from analysis.v3 import reconciled_photo_schedule as schedule
from analysis.v3.workflow import digest
from test_v3_reconciled_photo_schedule import fixture as source_fixture


def test_cloud_workflow_references_existing_tests():
    root=Path(__file__).resolve().parents[1]
    workflow=(root/'.github/workflows/ch1-v3-native-raw-chunks.yml').read_text()
    files=re.findall(r'tests/test_[a-z0-9_]+\.py',workflow)
    assert files and all((root/p).is_file() for p in files)


def test_consecutive_chunks_cover_all_components_without_splitting_or_repeating_pilot():
    groups=[dict(component_id=str(i),observations=2,requests=i%3,photo_jobs=4) for i in range(12)]
    actual=chunks.group_chunks(groups,{'0','1'},max_obs=6,max_requests=3)
    assert [g for c in actual for g in c['component_ids']]==[str(i) for i in range(2,12)]
    assert all(c['observations']<=6 and c['requests']<=3 for c in actual)
    assert sum(c['observations'] for c in actual)==20
    assert sum(c['photo_jobs'] for c in actual)==40
    assert chunks.group_chunks(groups,{'0','1'},6,3)==actual


@pytest.mark.parametrize('bad',[
    [dict(component_id='x',observations=129,requests=1,photo_jobs=1)],
    [dict(component_id='x',observations=1,requests=65,photo_jobs=65)],
    [dict(component_id='x',observations=0,requests=0,photo_jobs=1)]])
def test_oversized_components_stop_instead_of_being_split_or_lost(bad):
    with pytest.raises(ValueError,match='complete component'):
        chunks.group_chunks(bad,set())


def test_all_blocked_photo_groups_are_retained():
    groups=[dict(component_id=str(i),observations=1,requests=0,photo_jobs=3) for i in range(200)]
    actual=chunks.group_chunks(groups,set())
    assert [c['observations'] for c in actual]==[128,72]
    assert sum(c['photo_jobs'] for c in actual)==600


def test_explicit_component_packet_preserves_shared_photo_links_and_blocked_states(tmp_path):
    args=source_fixture(tmp_path)
    schedule.build(**args)
    path=args['out']/'reconciled_photo_schedule_private.sqlite'
    a=pilot_input(path,digest(path),128,component_ids=['component-a'])
    b=pilot_input(path,digest(path),128,component_ids=['component-b'])
    assert a['selected']==['1','2'] or a['selected']==['2','1']
    assert len(a['queue'])==1 and a['queue'][0]['known_metadata_obs_ids']==['2','99']
    assert b['selected']==['3'] and b['queue']==[]
    assert {r['photo_id'] for r in a['links']}.isdisjoint(r['photo_id'] for r in b['links'])
    a.update(plan_id='a'*64,chunk_id='c000000')
    unit=chunks.photo_unit_packet(a,a['queue'][0])
    assert unit['selected']==['2'] and unit['queue'][0]['component_id']=='component-a'
    assert unit['report']['mode']=='reconciled_photo_unit'
    assert 'SECRET_LOCATION' not in json.dumps(unit)
    with pytest.raises(ValueError,match='absent or exceed'):
        pilot_input(path,digest(path),128,component_ids=['absent'])
    with pytest.raises(ValueError,match='unique'):
        pilot_input(path,digest(path),128,component_ids=['component-a','component-a'])


def test_current_measurement_functions_match_completed_pilot():
    chunks.validate_algorithm()
