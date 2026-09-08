import copy
import json

import pytest

from analysis.v3 import archived_measurement as archived
from analysis.v3 import cloud_measurement_chunks as cloud
from analysis.v3 import stream_original_traits as worker
from analysis.v3.protected_artifacts import new_json
from analysis.v3.workflow import canonical_digest, digest
from test_v3_cloud_measurement_chunks import prepared
from test_v3_original_stream_traits import offline_run


def completed(prepared,tmp_path):
    packet,store,fixture=prepared
    report=cloud.execute_packet(packet,store,tmp_path/'cloud',fixture.weights)
    packet_path=tmp_path/'pinned_packet.json'; new_json(packet_path,packet)
    batch={'plan_id':packet['plan_id'],'chunks':{packet['chunk_id']:{'packet_sha256':digest(packet_path)}}}
    return packet,store,report,packet_path,batch,tmp_path/'cloud/completed/result.zip'


def test_archive_restores_one_verified_photo_at_a_time(prepared,tmp_path,monkeypatch):
    packet,store,report,pinned,batch,bundle=completed(prepared,tmp_path)
    def forbidden(*a,**k): raise AssertionError('Local archive verification must not request source images')
    monkeypatch.setattr(worker,'_download',forbidden)
    seen=[]
    def collect(unit,original,checked):
        assert unit.is_dir() and len(list(unit.iterdir()))==11
        assert canonical_digest(original)==canonical_digest(packet) and checked['raw_endpoint_slots']==54
        assert not seen or not seen[-1].exists()
        seen.append(unit)
    result=archived.verify_chunk(bundle,report['output_asset'],batch,packet['chunk_id'],tmp_path/'verify',
                                 packet_path=pinned,on_unit=collect)
    assert result['raw_endpoint_slots']==108 and result['bbox_slots']==100
    assert result['verified_photo_units']==2 and result['persistent_uncompressed_unit_files']==0
    assert len(seen)==2 and all(not p.exists() for p in seen) and bundle.exists()
    assert result['local_image_requests']==0 and result['ecological_fitting_authorized'] is False
    assert not list((tmp_path/'verify').rglob('*private*'))


def test_pinned_input_bytes_and_reserialized_cloud_packet_are_distinct(prepared,tmp_path):
    packet,store,report,pinned,batch,bundle=completed(prepared,tmp_path)
    # The original packet is independently byte-pinned; whitespace is allowed
    # only between that checked input and the cloud's reserialization of it.
    pinned.write_bytes(json.dumps(packet,separators=(',',':')).encode())
    batch['chunks'][packet['chunk_id']]['packet_sha256']=digest(pinned)
    result=archived.verify_chunk(bundle,report['output_asset'],batch,packet['chunk_id'],tmp_path/'verify',packet_path=pinned)
    assert result['input_packet_sha256']==digest(pinned)
    changed=copy.deepcopy(packet); changed['selection_scores']['1']='changed'
    pinned.write_text(json.dumps(changed)); batch['chunks'][packet['chunk_id']]['packet_sha256']=digest(pinned)
    with pytest.raises(ValueError,match='Restored packet differs'):
        archived.verify_chunk(bundle,report['output_asset'],batch,packet['chunk_id'],tmp_path/'bad',packet_path=pinned)


@pytest.mark.parametrize('key',['bundle_sha256','manifest_sha256','restored_bytes','files'])
def test_archive_rejects_changed_transport_evidence(prepared,tmp_path,key):
    packet,store,report,pinned,batch,bundle=completed(prepared,tmp_path)
    bad=copy.deepcopy(report['output_asset'])
    bad[key]='0'*64 if key.endswith('sha256') else bad[key]+1
    with pytest.raises(ValueError):
        archived.verify_chunk(bundle,bad,batch,packet['chunk_id'],tmp_path/'bad',packet_path=pinned)
    assert not (tmp_path/'bad/public_report.json').exists()
    assert (tmp_path/'bad/incomplete_verification.json').exists() and bundle.exists()


def test_callback_failure_cannot_promote_partial_view(prepared,tmp_path):
    packet,store,report,pinned,batch,bundle=completed(prepared,tmp_path)
    def interrupted(*a): raise RuntimeError('interrupted derived view')
    with pytest.raises(RuntimeError):
        archived.verify_chunk(bundle,report['output_asset'],batch,packet['chunk_id'],tmp_path/'bad',
                              packet_path=pinned,on_unit=interrupted)
    assert not list((tmp_path/'bad').glob('verified-unit-*')) and bundle.exists()
    assert not (tmp_path/'bad/public_report.json').exists()


def test_protected_recovery_uses_same_saved_zip_without_refetch(prepared,tmp_path,monkeypatch):
    packet,store,report,pinned,batch,bundle=completed(prepared,tmp_path)
    batch_path=tmp_path/'batch.json'; new_json(batch_path,batch)
    out=tmp_path/'saved'; out.mkdir(); (out/'bundle.zip').write_bytes(bundle.read_bytes())
    def forbidden(*a,**k): raise AssertionError('Exact saved ZIP must be reused')
    monkeypatch.setattr(store,'download',forbidden)
    result=archived.recover_chunk(batch_path,packet['chunk_id'],pinned,out,store=store)
    assert result['raw_endpoint_slots']==108 and len(store.assets)==1


def test_zero_request_chunk_has_explicit_zero_counts(prepared,tmp_path):
    packet,store,fixture=prepared
    packet['queue']=[]
    packet['report']['request_candidates']=0
    for row in packet['links']: row['status']='rights_restricted'
    report=cloud.execute_packet(packet,store,tmp_path/'cloud',fixture.weights)
    pinned=tmp_path/'packet.json'; new_json(pinned,packet)
    batch={'plan_id':packet['plan_id'],'chunks':{packet['chunk_id']:{'packet_sha256':digest(pinned)}}}
    result=archived.verify_chunk(tmp_path/'cloud/completed/result.zip',report['output_asset'],batch,packet['chunk_id'],
                                 tmp_path/'verify',packet_path=pinned)
    assert result['raw_endpoint_slots']==result['detected_heads']==result['bbox_slots']==0
