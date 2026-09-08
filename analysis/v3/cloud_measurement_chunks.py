"""Protected, resumable native raw-image measurement on GitHub Actions.

No ecological join. A completed photo is reusable only after its complete all27
numerical transaction is verified. Partial runs preserve verified photo units;
terminal download failures remain outcomes, not automatic retry candidates.
"""
from __future__ import annotations

import argparse
import hashlib
import inspect
import json
import os
from pathlib import Path
import re
import signal
import zipfile

from .measurement_chunks import contract, photo_unit_packet, validate_algorithm
from .private_replay import snapshot
from .protected_artifacts import DraftStore, MANIFEST, manifest_entries, new_json, pack, require, unpack
from .verify_original_stream import FILES, verify
from .workflow import ROOT, canonical_digest, digest, text_digest

DECISION = ROOT/'analysis/v3/measurement_qualification_decision_20260908.json'
RUNTIME = ROOT/'analysis/v3/measurement_runtime_contract.json'


def runtime_guard():
    from . import stream_original_traits as worker
    from .detect_cached_images import PARAMETERS
    spec=json.loads(RUNTIME.read_text())
    actual=worker.software_runtime()
    require(actual['python']==spec['python'] and actual['system']==spec['system']
            and all(actual['distribution_versions'][k]==v for k,v in spec['distribution_versions'].items()), 'Pinned measurement runtime differs')
    payload={'bbox_conditions':worker.BBOX_CONDITIONS,'identity':worker.IDENTITY,
             'oriented_pixel_function':inspect.getsource(worker.oriented_pixel_identity),'detector_parameters':PARAMETERS}
    require(canonical_digest(payload)==spec['bbox_identity_pixel_and_detector_canonical_sha256']
            and all(text_digest(ROOT/'analysis/v3'/k)==v for k,v in spec['helper_sha256_text_lf'].items()),'Measurement geometry or helper implementation differs')


def validate_packet(packet):
    spec=contract()
    require(packet['schedule_sha256']==spec['source_schedule_sha256']
            and re.fullmatch('[0-9a-f]{64}',packet['plan_id'])
            and re.fullmatch('c[0-9]{6}',packet['chunk_id'])
            and packet['report']['mode']=='reconciled_whole_component_chunk', 'Wrong native chunk identity')
    require(0 < len(packet['selected']) <= spec['maximum_observations_per_chunk']
            and len(packet['selected']) == len(set(packet['selected']))
            and len(packet['queue']) <= spec['maximum_requests_per_chunk'], 'Chunk exceeds fixed budget or duplicates observations')
    queue_ids={r['photo_id'] for r in packet['queue']}
    links=packet['links']
    require(len(queue_ids)==len(packet['queue'])
            and {r['obs_id'] for r in links}==set(packet['selected'])
            and queue_ids=={r['photo_id'] for r in links if r['status']=='request_candidate_not_authorized'}
            and len(links)==packet['report']['selected_photo_links']
            and len(queue_ids)==packet['report']['request_candidates'], 'Chunk request/link conservation differs')
    from .reconciled_photo_schedule import original_identity
    allowed=json.loads((ROOT/'analysis/v3/reconciled_photo_schedule_contract.json').read_text())['license_codes']
    for item in packet['queue']:
        require(item['license_code'] in allowed and original_identity(item['photo_id'],{'urls':[item['original_url']]})==('valid',item['original_url']), 'Invalid photo rights or source identity')
        photo_unit_packet(packet,item)


def download_result(store, metadata, out):
    """Authenticate GitHub's immutable asset digest before inspecting manifest."""
    sha=metadata.get('digest','')
    require(re.fullmatch('sha256:[0-9a-f]{64}',sha) is not None, 'Existing result lacks server SHA-256')
    asset={'asset_id':metadata['id'],'asset_name':metadata['name'],
           'bundle_bytes':metadata['size'],'bundle_sha256':sha[7:]}
    store.download(asset,out/'bundle.zip')
    with zipfile.ZipFile(out/'bundle.zip') as zipped:
        require(zipped.getinfo(MANIFEST).file_size <= 10_000_000,'Oversized result manifest')
        raw=zipped.read(MANIFEST)
    entries,_=manifest_entries(raw)
    asset.update(manifest_sha256=hashlib.sha256(raw).hexdigest(), files=len(entries),
                 restored_bytes=sum(r['bytes'] for r in entries),source_images_included=False)
    unpack(out/'bundle.zip',out/'verified',asset)
    return out/'verified/restored',asset


def verify_unit(path, expected):
    receipt=verify(path,DECISION)
    execution=json.loads((path/'execution_contract.json').read_text())
    require(execution['selection']==expected['report'],'Completed unit belongs to a different packet')
    from . import stream_original_traits as worker
    from .detect_cached_images import MODEL_SHA, PARAMETERS
    require(execution['feature_specification']==worker.features.specification()
            and execution['detector_parameters']==PARAMETERS and execution['model_sha256']==MODEL_SHA,
            'Completed unit measurement algorithm differs')
    detection=json.loads((path/'photo_detection_private.jsonl').read_text())
    require(all(canonical_digest(detection[k])==canonical_digest(v) for k,v in expected['queue'][0].items()),'Completed unit source versions differ')
    return receipt


def persist(store, packet, units, out, name, complete, interrupted=None):
    out.mkdir(parents=True,exist_ok=False)
    packet_path=out/'chunk_packet_private.json'
    new_json(packet_path,packet)
    state={'plan_id':packet['plan_id'],'chunk_id':packet['chunk_id'],
           'packet_canonical_sha256':canonical_digest(packet),'complete':complete,
           'completed_unit_indices':sorted(units),'interrupted_unit_index':interrupted,
           'runtime_contract_canonical_sha256':canonical_digest(json.loads(RUNTIME.read_text())),
           'ecological_models_executed':0,'source_images_persisted':0}
    state_path=out/'chunk_state.json'; new_json(state_path,state)
    entries=[{'name':p.name,'path':str(p.resolve()),'sha256':digest(p)} for p in (packet_path,state_path)]
    for index,path in sorted(units.items()):
        expected=photo_unit_packet(packet,packet['queue'][index])
        verify_unit(path,expected)
        for f in sorted(FILES|{'original_stream_report.json','independent_verification.json'}):
            source=path/f
            entries.append({'name':f'units/u{index:04d}/{f}','path':str(source.resolve()),'sha256':digest(source)})
    if interrupted is not None:
        partial=out.parent/f'new_units/u{interrupted:04d}'
        for f in sorted(FILES|{'original_stream_report.json'}):
            source=partial/f
            if source.is_file():
                entries.append({'name':f'interrupted_unverified/u{interrupted:04d}/{f}',
                                'path':str(source.resolve()),'sha256':digest(source)})
    selection=out/'selection.json'; new_json(selection,{'schema_version':1,'files':entries})
    snapshot(selection,out/'snapshot')
    asset=pack(out/'snapshot',out/'result.zip')
    asset.update(store.upload(out/'result.zip',name))
    new_json(out/'upload_receipt.json',asset)
    # A receipt is not complete until a fresh protected download also restores.
    store.download(asset,out/'returned.zip')
    unpack(out/'returned.zip',out/'return',asset)
    require(json.loads((out/'return/restored/chunk_state.json').read_text())==state,'Returned checkpoint state differs')
    return asset


def execute_packet(packet, store, out, weights):
    """Resume committed units, measure only missing units, checkpoint on interruption."""
    validate_packet(packet)
    validate_algorithm()
    require(not out.exists(),'Preserve previous cloud execution directory')
    out.mkdir(parents=True)
    prefix=f"v3-raw-{packet['plan_id'][:16]}-{packet['chunk_id']}"
    existing=[r for r in store.check()['assets'] if r['name']==prefix+'.zip' or r['name'].startswith(prefix+'-incomplete-')]
    units, original_hashes, output_asset={}, {}, None
    for meta in sorted(existing,key=lambda r:r['id']):
        restored, asset=download_result(store,meta,out/f"recovered-{meta['id']}")
        saved=json.loads((restored/'chunk_packet_private.json').read_text())
        state=json.loads((restored/'chunk_state.json').read_text())
        require(canonical_digest(saved)==canonical_digest(packet) and state['packet_canonical_sha256']==canonical_digest(packet)
                and state['plan_id']==packet['plan_id'] and state['chunk_id']==packet['chunk_id']
                and state['runtime_contract_canonical_sha256']==canonical_digest(json.loads(RUNTIME.read_text())), 'Stored checkpoint input differs')
        indices=state['completed_unit_indices']
        require(len(indices)==len(set(indices)) and all(type(i) is int and 0 <= i < len(packet['queue']) for i in indices),'Invalid checkpoint unit indices')
        for i in indices:
            path=restored/f'units/u{i:04d}'
            verified=verify_unit(path,photo_unit_packet(packet,packet['queue'][i]))
            require(i not in original_hashes or original_hashes[i]==verified['verified_file_sha256'],'Conflicting completed numerical units; do not choose a result')
            units[i],original_hashes[i]=path,verified['verified_file_sha256']
        if meta['name']==prefix+'.zip':
            require(state['complete'] is True and set(indices)==set(range(len(packet['queue']))),'Final asset is not complete')
            output_asset=asset
        else:
            require(state['complete'] is False,'Incomplete asset relabelled as final')
    resumed=len(units)
    image_requests=0
    if output_asset is None:
        from ultralytics import YOLO
        from .stream_original_traits import run_photo_unit
        model=YOLO(str(weights)) if len(units)<len(packet['queue']) else None
        interrupted=None
        try:
            for i,item in enumerate(packet['queue']):
                if i in units:
                    continue
                interrupted=i
                path=out/f'new_units/u{i:04d}'
                expected=photo_unit_packet(packet,item)
                run_photo_unit(expected,weights,DECISION,path,model=model)
                checked=verify_unit(path,expected)
                new_json(path/'independent_verification.json',checked)
                units[i]=path
                image_requests+=1
                print(json.dumps({'chunk':packet['chunk_id'],'photo_units_committed_locally':len(units),'request_slots':len(packet['queue'])}),flush=True)
            interrupted=None
        except BaseException:
            run_id=os.environ.get('GITHUB_RUN_ID','local')
            attempt=os.environ.get('GITHUB_RUN_ATTEMPT','1')
            require(run_id.isdigit() and attempt.isdigit(),'Interrupted cloud checkpoint requires run identity')
            # Includes all verified transactions. An interrupted photo without a
            # complete all27 report is not promoted to a completed measurement.
            persist(store,packet,units,out/'checkpoint',f'{prefix}-incomplete-{run_id}-{attempt}.zip',False,interrupted)
            raise
        output_asset=persist(store,packet,units,out/'completed',prefix+'.zip',True)
    summaries=[verify_unit(units[i],photo_unit_packet(packet,item)) for i,item in enumerate(packet['queue'])]
    report={'status':'NATIVE_RAW_MEASUREMENT_CHUNK_PROTECTED_AND_VERIFIED_NO_ECOLOGY',
            'plan_id':packet['plan_id'],'chunk_id':packet['chunk_id'],'packet_canonical_sha256':canonical_digest(packet),
            'selected_observations':len(packet['selected']),'source_photo_jobs':packet['report']['selected_photo_jobs'],
            'source_photo_links':len(packet['links']),'request_slots':len(packet['queue']),
            'units_restored_without_requests':resumed,'new_photo_units_executed':image_requests,
            'detected_heads':sum(r['detected_heads'] for r in summaries),
            'all27_raw_slots':sum(r['raw_endpoint_slots'] for r in summaries),
            'bbox_slots':sum(r['bbox_slots'] for r in summaries),'output_asset':output_asset,
            'restored_checkpoint_asset_ids':[r['id'] for r in existing],
            'run_id':os.environ.get('GITHUB_RUN_ID'),'run_attempt':os.environ.get('GITHUB_RUN_ATTEMPT'),
            'commit':os.environ.get('GITHUB_SHA'),'source_images_persisted':0,'environment_values_read':0,
            'ecological_models_executed':0,'ecological_fitting_authorized':False,
            'limits':['This is raw image measurement, not full-native completion, independent physical accuracy or ecological admission.',
                      'Reuse includes failed downloads and no-detection outcomes; no success-conditioned retry or selection.',
                      'Source dependence components remain together in the chunk; per-photo transactions are computational checkpoints, not independent biological replicates.']}
    new_json(out/'public_report.json',report)
    return report


def run_cloud(batch_path,chunk_id,out):
    batch=json.loads(batch_path.read_text())
    require(batch['status']=='bounded_native_raw_measurement_batch_no_ecology'
            and batch['measurement_contract_canonical_sha256']==canonical_digest(contract())
            and batch['runtime_contract_canonical_sha256']==canonical_digest(json.loads(RUNTIME.read_text()))
            and chunk_id in batch['chunks'],'Batch does not authorize this chunk')
    require(not out.exists(),'Preserve previous batch output')
    out.mkdir(parents=True)
    runtime_guard()
    store=DraftStore(batch)
    try:
        store.download(batch['input_asset'],out/'input.zip')
        unpack(out/'input.zip',out/'input',batch['input_asset'])
        packet_path=out/'input/restored'/f'{chunk_id}_packet_private.json'
        require(digest(packet_path)==batch['chunks'][chunk_id]['packet_sha256'],'Batch packet bytes differ')
        packet=json.loads(packet_path.read_text())
        require(packet['plan_id']==batch['plan_id'] and packet['chunk_id']==chunk_id,'Batch plan differs')
        from .recover_revision_inputs import download
        download(8076736948,out/'detector.zip','2bc7cc49c2f213d4a5c5dab96e42fcdda7ceda23ffc3130520a5dd44b8c76c05')
        with zipfile.ZipFile(out/'detector.zip') as zipped:
            weights=zipped.read('recovery/model/weights/best.pt')
        require(hashlib.sha256(weights).hexdigest()==contract()['algorithm']['model_sha256'],'Detector model differs')
        path=out/'best.pt'
        with path.open('xb') as handle: handle.write(weights)
        return execute_packet(packet,store,out/'measurement',path)
    finally:
        store.close()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch',type=Path,required=True)
    parser.add_argument('--chunk',required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    def interrupted(signum,frame):
        raise KeyboardInterrupt('Runner requested interruption; preserve completed numerical units')
    signal.signal(signal.SIGTERM,interrupted)
    try:
        print(json.dumps(run_cloud(args.batch,args.chunk,args.out)))
    except (Exception, KeyboardInterrupt) as error:
        # Public Actions logs must not receive private URLs, signed download
        # addresses, photo identifiers or traceback locals from an exception.
        print(json.dumps({'status':'RAW_MEASUREMENT_RUN_INCOMPLETE','error_type':type(error).__name__}),flush=True)
        raise SystemExit(1) from None


if __name__=='__main__':
    main()
