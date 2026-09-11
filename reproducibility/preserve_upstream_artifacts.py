"""Draft-only preservation of upstream image-processing provenance; no publication."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import urllib.request
import zipfile

EXTRA = {
    8066675131: ('8acd73dcd227aec49d60c268fd97e8145b4fc9479226dc6227863b83464dc6bd', 'original_screening_queue'),
    8068122589: ('fa755db8de0288de0ed09b6858a535bfc8cf0a4eb583a132e8671bc07889aaa5', 'original_screening_images_qc_not_independent_validation'),
}

def hash_file(path, algorithm='sha256'):
    h=hashlib.new(algorithm)
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

def preserve(draft, request, base, out):
    source=json.loads(Path('reproducibility/actions_artifact_catalog.json').read_text(encoding='utf-8'))
    plan={r['artifact_id']:(r.get('verified_download_zip_sha256') or r.get('local_archive_sha256') or r.get('github_digest','').removeprefix('sha256:'),r['role'])
          for r in source['artifacts'] if r['artifact_id'] not in [9612943217,8227254443,8983877726,9632715852]}
    plan.update(EXTRA)
    # Large merged tables go last; preserve the acquisition/model packages first.
    order=sorted(plan,key=lambda aid:(aid==8269246732,aid))
    root=Path('work/zenodo-upstream');root.mkdir(parents=True,exist_ok=True)
    report={'scope':'Upstream acquisition, detector and measurement artifact preservation',
            'zero_from_scratch_reproduction_verified':False,'published':False,'artifacts':[],
            'remaining_limits':[
                'Full original photographic collection not permanently preserved or redownload-tested.',
                'Model package contains dataset manifests, not all training image and label files.',
                'External Grounding DINO / CLIP weights and revisions require a separate completeness check.',
                'Annotation packet, private audit mapping and hidden predictions must remain non-public until blinding and image rights are reviewed.',
                'Artifact preservation does not establish independent detector or biological accuracy.']}
    bucket=draft['links']['bucket']
    assert bucket.startswith('https://zenodo.org/api/files/')
    def upload(path):
        state=request(base)
        assert not state['submitted'] and state['state']=='unsubmitted'
        url=bucket+'/'+path.name
        digest=hash_file(path);md5=hash_file(path,'md5')
        existing=[f for f in state.get('files',[]) if f['filename']==path.name]
        if existing:
            assert existing[0]['checksum'].removeprefix('md5:')==md5, 'Refusing to overwrite different bytes'
        else:
            with path.open('rb') as f:
                req=urllib.request.Request(url,data=f,method='PUT',headers={
                    'Authorization':'Bearer '+os.environ['ZENODO_TOKEN'],
                    'Content-Type':'application/octet-stream','Content-Length':str(path.stat().st_size)})
                with urllib.request.urlopen(req,timeout=1800) as r: uploaded=json.load(r)
            assert uploaded['checksum']=='md5:'+md5
        req=urllib.request.Request(url,headers={'Authorization':'Bearer '+os.environ['ZENODO_TOKEN']})
        h=hashlib.sha256()
        with urllib.request.urlopen(req,timeout=1800) as r:
            for b in iter(lambda:r.read(1024*1024),b''):h.update(b)
        assert h.hexdigest()==digest, 'Readback mismatch'
        return {'filename':path.name,'sha256':digest,'bytes':path.stat().st_size,'authenticated_readback_verified':True}
    for aid in order:
        expected,role=plan[aid]
        assert len(expected)==64
        meta=json.loads(subprocess.check_output(['gh','api',f'repos/zuizui0223/azami/actions/artifacts/{aid}']))
        assert not meta['expired'] and meta['digest']=='sha256:'+expected
        path=root/f'upstream-artifact-{aid}.zip'
        with path.open('wb') as f:
            subprocess.run(['gh','api',f'repos/zuizui0223/azami/actions/artifacts/{aid}/zip'],stdout=f,check=True)
        assert hash_file(path)==expected
        with zipfile.ZipFile(path) as z:
            assert z.testzip() is None
            members=[]
            for info in z.infolist():
                if info.is_dir():continue
                h=hashlib.sha256()
                with z.open(info) as f:
                    for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
                members.append({'path':info.filename,'bytes':info.file_size,'sha256':h.hexdigest()})
            if aid==8076736948:
                weight=z.read('recovery/model/weights/best.pt')
                assert hashlib.sha256(weight).hexdigest()=='4078e0510532852681b65ee529cd82237b649ec99b17c4ca5f1da460a62d2bed'
        result=upload(path)
        report['artifacts'].append({**result,'artifact_id':aid,'role':role,'run':meta['workflow_run'],'members':members})
        (out/'upstream_progress.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
        print(f'Upstream artifact {aid}: upload and readback verified',flush=True)
    catalog=root/'UPSTREAM_ARTIFACT_CATALOG.json'
    catalog.write_text(json.dumps(report,indent=2),encoding='utf-8')
    upload(catalog)
    state=request(base);metadata=state['metadata'].copy()
    metadata['description'] += ('<p>Upstream processing artifacts are additionally preserved: acquisition metadata and screening queue/images, '
        'the frozen detector package with best/last weights, dataset manifest and training diagnostics, historical trait outputs, '
        'exhaustive merged continuous measurements, and the uncompleted independent-audit materials. '
        'See UPSTREAM_ARTIFACT_CATALOG.json for exact source runs, member hashes and remaining gaps. '
        'This is NOT certified zero-from-scratch reproduction: full original photographs and all training inputs or external model weights have not been verified complete. '
        'Audit mappings, hidden predictions and image-containing packages require blinding and rights review before any publication.</p>')
    request(base,'PUT',{'metadata':metadata})
    final=request(base)
    assert not final['submitted'] and final['state']=='unsubmitted'
    report['status']='UPSTREAM_DRAFT_UPLOADS_READBACK_VERIFIED'
    (out/'upstream_final.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
