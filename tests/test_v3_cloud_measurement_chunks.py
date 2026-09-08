import copy
import json
from pathlib import Path

import pytest

from analysis.v3 import cloud_measurement_chunks as cloud
from analysis.v3 import measurement_chunks as chunks
from analysis.v3 import stream_original_traits as worker
from analysis.v3.workflow import digest
from test_v3_original_stream_traits import offline_run, photo_row, write_metadata


@pytest.fixture
def prepared(offline_run,monkeypatch):
    write_metadata(offline_run.metadata,[photo_row(photo='11'),photo_row(photo='12')])
    links=[]
    queue,_=worker.photo_schedule(offline_run.metadata,{'1'},links)
    for row in queue: row['component_id']='component-one'
    packet={'selected':['1'],'selection_scores':{'1':'score'},'queue':queue,'links':links,
            'input_sha256':{'enriched':'e'*64},'schedule_sha256':chunks.contract()['source_schedule_sha256'],
            'plan_id':'a'*64,'chunk_id':'c000000',
            'report':{'mode':'reconciled_whole_component_chunk','selected_observations':1,
                      'selected_photo_links':2,'selected_photo_jobs':2,'request_candidates':2}}
    for row in packet['links']: row['status']='request_candidate_not_authorized'
    monkeypatch.setattr(cloud,'validate_algorithm',lambda:None)
    monkeypatch.setattr(chunks,'validate_algorithm',lambda:None)
    monkeypatch.setenv('GITHUB_RUN_ID','123')
    monkeypatch.setenv('GITHUB_RUN_ATTEMPT','1')
    class Store:
        def __init__(self): self.assets=[]; self.blobs={}
        def check(self): return {'draft':True,'assets':self.assets}
        def upload(self,path,name):
            assert not any(a['name']==name for a in self.assets)
            entry={'id':len(self.assets)+1,'name':name,'state':'uploaded','size':path.stat().st_size,'digest':'sha256:'+digest(path)}
            self.assets.append(entry); self.blobs[entry['id']]=path.read_bytes()
            return {'asset_id':entry['id'],'asset_name':name}
        def download(self,asset,out):
            out.parent.mkdir(parents=True,exist_ok=True)
            assert not out.exists()
            out.write_bytes(self.blobs[asset['asset_id']])
    return packet,Store(),offline_run


def test_cloud_chunk_all27_and_complete_resume_never_refetches(prepared,tmp_path,monkeypatch,capsys):
    packet,store,fixture=prepared
    first=cloud.execute_packet(packet,store,tmp_path/'first',fixture.weights)
    assert first['all27_raw_slots']==108 and first['bbox_slots']==100
    assert first['new_photo_units_executed']==2 and first['units_restored_without_requests']==0
    assert len(store.assets)==1
    def forbidden(*a,**k): raise AssertionError('Completed photo must not be requested again')
    monkeypatch.setattr(worker,'_download',forbidden)
    second=cloud.execute_packet(packet,store,tmp_path/'second',fixture.weights)
    assert second['new_photo_units_executed']==0 and second['units_restored_without_requests']==2
    assert second['all27_raw_slots']==first['all27_raw_slots'] and len(store.assets)==1
    public=json.dumps(second)
    assert 'original.jpg' not in public and 'photo_id' not in public and 'obs_id' not in public
    logs=capsys.readouterr().out
    assert all(secret not in logs for secret in ('original.jpg','photo_id','obs_id','component-one'))


def test_interrupted_chunk_restores_only_completed_photo_then_continues(prepared,tmp_path,monkeypatch):
    packet,store,fixture=prepared
    original=worker.run_photo_unit
    calls=[]
    def interrupted(unit,*a,**k):
        calls.append(unit['queue'][0]['photo_id'])
        if len(calls)==2: raise KeyboardInterrupt('simulated runner interruption')
        return original(unit,*a,**k)
    monkeypatch.setattr(worker,'run_photo_unit',interrupted)
    with pytest.raises(KeyboardInterrupt):
        cloud.execute_packet(packet,store,tmp_path/'first',fixture.weights)
    assert len(store.assets)==1 and '-incomplete-' in store.assets[0]['name']
    new_calls=[]
    def resumed(unit,*a,**k):
        new_calls.append(unit['queue'][0]['photo_id']); return original(unit,*a,**k)
    monkeypatch.setattr(worker,'run_photo_unit',resumed)
    result=cloud.execute_packet(packet,store,tmp_path/'resume',fixture.weights)
    assert new_calls==[calls[1]]
    assert result['units_restored_without_requests']==1 and result['new_photo_units_executed']==1
    assert result['all27_raw_slots']==108 and len(store.assets)==2


def test_terminal_download_failure_is_not_retried_by_resume(prepared,tmp_path,monkeypatch):
    packet,store,fixture=prepared
    def fail(*a,**k): raise ConnectionError('offline source failure')
    monkeypatch.setattr(worker,'_download',fail)
    first=cloud.execute_packet(packet,store,tmp_path/'first',fixture.weights)
    assert first['all27_raw_slots']==0 and first['new_photo_units_executed']==2
    second=cloud.execute_packet(packet,store,tmp_path/'second',fixture.weights)
    assert second['new_photo_units_executed']==0 and second['units_restored_without_requests']==2


def test_changed_packet_cannot_reuse_existing_chunk(prepared,tmp_path):
    packet,store,fixture=prepared
    cloud.execute_packet(packet,store,tmp_path/'first',fixture.weights)
    changed=copy.deepcopy(packet); changed['selection_scores']['1']='different'
    with pytest.raises(ValueError,match='checkpoint input differs'):
        cloud.execute_packet(changed,store,tmp_path/'bad',fixture.weights)


def test_invalid_host_or_duplicate_request_is_rejected_before_network(prepared):
    packet,_,_=prepared
    bad=copy.deepcopy(packet); bad['queue'][0]['original_url']='https://example.org/private'
    with pytest.raises(ValueError,match='source identity'):
        cloud.validate_packet(bad)
    bad=copy.deepcopy(packet); bad['queue'][1]=bad['queue'][0]
    with pytest.raises(ValueError,match='conservation'):
        cloud.validate_packet(bad)


def test_runtime_contract_checks_versions_and_fixed_geometry(monkeypatch):
    spec=json.loads(cloud.RUNTIME.read_text())
    monkeypatch.setattr(worker,'software_runtime',lambda:{'python':spec['python'],'system':spec['system'],
                                                          'distribution_versions':spec['distribution_versions']})
    cloud.runtime_guard()
    monkeypatch.setattr(worker,'software_runtime',lambda:{'python':'wrong','system':spec['system'],
                                                          'distribution_versions':spec['distribution_versions']})
    with pytest.raises(ValueError,match='runtime differs'):
        cloud.runtime_guard()


def test_interrupted_partial_numerics_are_retained_but_not_promoted(prepared,tmp_path,monkeypatch):
    packet,store,fixture=prepared
    original=worker.run_photo_unit
    def partial(unit,weights,decision,out,**kwargs):
        out.mkdir(parents=True)
        (out/'photo_detection_private.jsonl').write_text('{"incomplete":true}\n')
        raise KeyboardInterrupt()
    monkeypatch.setattr(worker,'run_photo_unit',partial)
    with pytest.raises(KeyboardInterrupt):
        cloud.execute_packet(packet,store,tmp_path/'interrupted',fixture.weights)
    recovered,_=cloud.download_result(store,store.assets[0],tmp_path/'inspection')
    assert (recovered/'interrupted_unverified/u0000/photo_detection_private.jsonl').is_file()
    assert json.loads((recovered/'chunk_state.json').read_text())['completed_unit_indices']==[]
    monkeypatch.setattr(worker,'run_photo_unit',original)
    result=cloud.execute_packet(packet,store,tmp_path/'resume',fixture.weights)
    assert result['units_restored_without_requests']==0 and result['new_photo_units_executed']==2


@pytest.mark.parametrize('change',['receipt_hash','overlapping_chunk'])
def test_next_wave_rejects_changed_predecessor_before_runtime_or_network(tmp_path,monkeypatch,change):
    batch=json.loads((cloud.ROOT/'analysis/v3/native_measurement_wave_20260908_b.json').read_text())
    if change=='receipt_hash':
        batch['previous_completion_receipt_canonical_sha256']='0'*64
    else:
        batch['chunks']['c000000']=batch['chunks']['c000002']
    path=tmp_path/'bad_batch.json'; path.write_text(json.dumps(batch))
    def forbidden(): raise AssertionError('Invalid predecessor must fail before runtime or network')
    monkeypatch.setattr(cloud,'runtime_guard',forbidden)
    with pytest.raises(ValueError,match='Previous batch completion'):
        cloud.run_cloud(path,'c000002',tmp_path/'out')
    assert not (tmp_path/'out').exists()
