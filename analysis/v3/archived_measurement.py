"""Verify locally restored numerical units without persisting millions of copies.

The complete protected ZIP stays on disk. Every blob is decompressed and hashed,
then each photo's files are restored, independently checked and optionally read
by an observation-view builder. Only the temporary verification copies are removed.
No source images, traits/environment join or measurement admission is involved.
"""
from __future__ import annotations

from collections import Counter
import argparse
import hashlib
import json
from pathlib import Path
import stat
import tempfile
import zipfile

from .cloud_measurement_chunks import DECISION, RUNTIME, validate_packet, verify_unit
from .measurement_chunks import photo_unit_packet
from .protected_artifacts import DraftStore, MANIFEST, manifest_entries, new_json, require
from .recover_native_source_authority import private_directory, recover_lf_source
from .verify_original_stream import FILES, read_rows
from .workflow import canonical_digest, digest

UNIT_FILES=FILES|{'original_stream_report.json','independent_verification.json'}


class NumericalSnapshot:
    """Exact, bounded local decompression of the manifest's numerical bytes."""
    def __init__(self,bundle: Path,expected: dict):
        require(bundle.stat().st_size==expected['bundle_bytes']
                and digest(bundle)==expected['bundle_sha256'],'Saved numerical bundle differs')
        self.archive=zipfile.ZipFile(bundle)
        try:
            infos=self.archive.infolist()
            names=[i.filename for i in infos]
            require(len(names)==len(set(names)) and names.count(MANIFEST)==1,'Duplicate or missing snapshot member')
            require(self.archive.getinfo(MANIFEST).file_size<=10_000_000,'Manifest exceeds bound')
            raw=self.archive.read(MANIFEST)
            require(hashlib.sha256(raw).hexdigest()==expected['manifest_sha256'],'Snapshot manifest differs')
            entries,sizes=manifest_entries(raw)
            require(len(entries)==expected['files'] and sum(r['bytes'] for r in entries)==expected['restored_bytes'],
                    'Restoration denominator differs')
            require(set(names)=={MANIFEST,*('blobs/'+sha for sha in sizes)},'Unexpected archive member')
            for info in infos:
                require(not info.is_dir() and not stat.S_ISLNK(info.external_attr>>16),'Directory or link in snapshot')
                require(info.filename==MANIFEST or info.file_size==sizes[info.filename[6:]],'Blob size differs')
            for sha,size in sizes.items():
                h=hashlib.sha256(); restored=0
                with self.archive.open('blobs/'+sha) as source:
                    for part in iter(lambda:source.read(1024*1024),b''):
                        restored+=len(part); require(restored<=size,'Decompressed blob exceeds bound'); h.update(part)
                require(restored==size and h.hexdigest()==sha,'Locally decompressed bytes differ')
            self.entries={r['name']:r for r in entries}
        except BaseException:
            self.archive.close()
            raise

    def read(self,name):
        row=self.entries[name]
        data=self.archive.read('blobs/'+row['sha256'])
        require(len(data)==row['bytes'] and hashlib.sha256(data).hexdigest()==row['sha256'],'Restored member differs')
        return data

    def json(self,name):
        return json.loads(self.read(name))

    def close(self):
        self.archive.close()


def exact_decision_view(out: Path,expected_sha: str):
    if digest(DECISION)==expected_sha:
        return DECISION
    target=out/'measurement_decision_exact_lf.json'
    if not target.exists():
        recover_lf_source(DECISION,target,expected_sha)
    require(digest(target)==expected_sha,'Exact measurement decision view differs')
    return target


def verify_chunk(bundle: Path,asset: dict,batch: dict,chunk_id: str,out: Path,*,packet_path: Path,on_unit=None):
    """Verify all saved units; an optional local callback sees one verified unit.

    The callback must not keep temporary paths. It may build a separate numerical
    view; it cannot change the raw bundle or promote a missing/unverified unit.
    """
    out=private_directory(out)
    require(not out.exists(),'Preserve earlier archive verification')
    require(chunk_id in batch['chunks'],'Chunk outside batch')
    require(digest(packet_path)==batch['chunks'][chunk_id]['packet_sha256'],'Pinned input packet bytes differ')
    expected_packet=json.loads(packet_path.read_text(encoding='utf-8'))
    out.mkdir(parents=True)
    snapshot=None
    try:
        snapshot=NumericalSnapshot(bundle,asset)
        packet=snapshot.json('chunk_packet_private.json')
        state=snapshot.json('chunk_state.json')
        validate_packet(packet)
        # The original batch bytes remain exact-pinned. The cloud serializes a
        # new packet file on Linux, so its JSON identity (not Windows line
        # endings) must equal that separately verified original input.
        require(canonical_digest(packet)==canonical_digest(expected_packet),
                'Restored packet differs from pinned batch input')
        require(packet['chunk_id']==chunk_id and packet['plan_id']==batch['plan_id']
                and state['chunk_id']==chunk_id and state['plan_id']==batch['plan_id']
                and state['packet_canonical_sha256']==canonical_digest(packet)
                and state['complete'] is True
                and state['runtime_contract_canonical_sha256']==canonical_digest(json.loads(RUNTIME.read_text()))
                and state['completed_unit_indices']==list(range(len(packet['queue']))),'Final checkpoint state differs')
        expected_names={'chunk_packet_private.json','chunk_state.json'}|{
            f'units/u{i:04d}/{name}' for i in range(len(packet['queue'])) for name in UNIT_FILES}
        require(set(snapshot.entries)==expected_names,'Final snapshot unit inventory differs')
        counts=Counter(dict.fromkeys(('detected_heads','raw_endpoint_slots','bbox_slots'),0))
        transfer=Counter(); eligible=Counter(); elapsed=0.0; total_bytes=0; no_detection=0
        for i,item in enumerate(packet['queue']):
            # The disposable directory is created below the validated private
            # output path. Only these reconstructible copies are ever removed.
            with tempfile.TemporaryDirectory(prefix='verified-unit-',dir=out) as tmp:
                unit=Path(tmp).resolve()
                require(unit.is_relative_to(out) and unit!=out,'Temporary restore left its private parent')
                for name in sorted(UNIT_FILES):
                    with (unit/name).open('xb') as handle:
                        handle.write(snapshot.read(f'units/u{i:04d}/{name}'))
                original=json.loads((unit/'original_stream_report.json').read_text())
                decision=exact_decision_view(out,original['measurement_decision_sha256'])
                checked=verify_unit(unit,photo_unit_packet(packet,item),decision_path=decision)
                for key in ('detected_heads','raw_endpoint_slots','bbox_slots'): counts[key]+=checked[key]
                eligible.update(checked['measurement_eligible_heads'])
                row=read_rows(unit/'transfer_private.csv')[0]
                transfer[row['status']]+=1
                elapsed+=float(row['elapsed_seconds']); total_bytes+=int(row['bytes'] or 0)
                no_detection+=int(row['status']=='success' and int(row['detections'])==0)
                if on_unit is not None:
                    on_unit(unit,packet,checked)
        report={'status':'ARCHIVED_NATIVE_RAW_CHUNK_LOCALLY_RESTORED_AND_VERIFIED_NO_ECOLOGY',
                'chunk_id':chunk_id,'plan_id':batch['plan_id'],'asset':asset,
                'source_packet_canonical_sha256':canonical_digest(packet),
                'input_packet_sha256':digest(packet_path),
                'selected_observations':len(packet['selected']),'verified_photo_units':len(packet['queue']),
                'source_photo_jobs':packet['report']['selected_photo_jobs'],'source_photo_links':len(packet['links']),
                **dict(counts),'transfer_status_counts':dict(transfer),'photos_without_detection':no_detection,
                'measurement_eligible_heads':dict(eligible),'request_elapsed_seconds_sum':elapsed,
                'downloaded_source_bytes':total_bytes,'numerical_files_restored_and_verified':asset['files'],
                'numerical_bytes_restored_and_verified':asset['restored_bytes'],
                'persistent_uncompressed_unit_files':0,'source_images_persisted':0,
                'local_image_requests':0,'environment_values_read':0,'ecological_models_executed':0,
                'ecological_fitting_authorized':False,
                'limits':['Complete original numerical bytes remain in the saved protected ZIP; temporary verification copies are reconstructed from it.',
                          'This verifies data integrity and operational support, not botanical accuracy or full-cohort ecological admission.']}
        new_json(out/'public_report.json',report)
        return report
    except BaseException as error:
        new_json(out/'incomplete_verification.json',{'status':'INCOMPLETE_DO_NOT_USE','error_type':type(error).__name__})
        raise
    finally:
        if snapshot is not None:
            snapshot.close()


def recover_chunk(batch_path: Path,chunk_id: str,packet_path: Path,out: Path,*,store=None):
    """Recover an already uploaded final numerical asset; never request images.

    An interrupted download remains on disk. It is never overwritten or treated
    as complete. A complete saved ZIP may be reused after exact remote checks.
    """
    batch=json.loads(batch_path.read_text(encoding='utf-8'))
    require(chunk_id in batch['chunks'],'Chunk outside batch')
    require(digest(packet_path)==batch['chunks'][chunk_id]['packet_sha256'],'Pinned input packet bytes differ')
    out=private_directory(out); out.mkdir(parents=True,exist_ok=True)
    require(not (out/'verified').exists(),'Preserve earlier local verification')
    own_store=store is None
    store=store or DraftStore(batch)
    try:
        name=f"v3-raw-{batch['plan_id'][:16]}-{chunk_id}.zip"
        matches=[r for r in store.check()['assets'] if r['name']==name and r['state']=='uploaded']
        require(len(matches)==1,'One completed protected numerical asset required')
        meta=matches[0]
        require(meta.get('digest','').startswith('sha256:') and len(meta['digest'])==71,'Remote asset digest unavailable')
        asset={'asset_id':meta['id'],'asset_name':name,'bundle_bytes':meta['size'],
               'bundle_sha256':meta['digest'][7:]}
        bundle=out/'bundle.zip'
        if not bundle.exists():
            store.download(asset,bundle)
        require(bundle.stat().st_size==asset['bundle_bytes'] and digest(bundle)==asset['bundle_sha256'],
                'Saved bundle incomplete or changed; preserve it, do not overwrite')
        # These manifest fields are extracted only from an exact server-digest
        # verified ZIP. NumericalSnapshot independently checks every member.
        with zipfile.ZipFile(bundle) as zipped:
            require(zipped.getinfo(MANIFEST).file_size<=10_000_000,'Manifest exceeds bound')
            raw=zipped.read(MANIFEST)
        entries,_=manifest_entries(raw)
        asset.update(manifest_sha256=hashlib.sha256(raw).hexdigest(),files=len(entries),
                     restored_bytes=sum(r['bytes'] for r in entries),source_images_included=False)
        return verify_chunk(bundle,asset,batch,chunk_id,out/'verified',packet_path=packet_path)
    finally:
        if own_store:
            store.close()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch',type=Path,required=True)
    parser.add_argument('--chunk',required=True)
    parser.add_argument('--packet',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    try:
        print(json.dumps(recover_chunk(args.batch,args.chunk,args.packet,args.out)))
    except Exception as error:
        print(json.dumps({'status':'ARCHIVE_RESTORATION_INCOMPLETE','error_type':type(error).__name__}))
        raise SystemExit(1) from None


if __name__=='__main__':
    main()
