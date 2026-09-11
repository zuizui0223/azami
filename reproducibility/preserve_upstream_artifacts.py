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

def preserve(draft, request, base, out, extra_only=False):
    source=json.loads(Path('reproducibility/actions_artifact_catalog.json').read_text(encoding='utf-8'))
    plan={r['artifact_id']:(r.get('verified_download_zip_sha256') or r.get('local_archive_sha256') or r.get('github_digest','').removeprefix('sha256:'),r['role'])
          for r in source['artifacts'] if r['artifact_id'] not in [9612943217,8227254443,8983877726,9632715852]}
    plan.update(EXTRA)
    # User scope: current manuscript dependencies only, not the whole legacy tree.
    paper_ids={8066010557,8066675131,8068122589,8076736948,8099953404,8225059018,8269246732}
    plan={aid:value for aid,value in plan.items() if aid in paper_ids}
    assert not extra_only, 'Historical ML expansion excluded by current manuscript-only scope'
    if extra_only:
        plan={
            8071529579:('ebcaab40fecc49a8515004a724c587717130b9531a89fa1f88f8e7c84f3953e4','pseudo_label_training_source'),
            8069610715:('d6a706658c8735727caa25506a6d14828381457cb23c2d5e6a9b67d9352a2f39','early_grounding_dino_proposals'),
            8077189280:('8aa211f5907dcb22fbdae287a892abeb3e84af6068f63ae8e44ae8d1c40c14b9','historical_zero_shot_CLIP_outputs_and_ontology')}
    # Large merged tables go last; preserve the acquisition/model packages first.
    order=sorted(plan,key=lambda aid:(aid==8269246732,aid))
    root=Path('work/zenodo-upstream');root.mkdir(parents=True,exist_ok=True)
    report={'scope':'Upstream acquisition, detector and measurement artifact preservation',
            'zero_from_scratch_reproduction_verified':False,'published':False,'artifacts':[],
            'remaining_limits':[
                'Full original photographic collection not permanently preserved or redownload-tested.',
                'Model package contains dataset manifests, not all training image and label files.',
                'External Grounding DINO / CLIP weights and revisions require a separate completeness check.',
                'Unused categorical ML development artifacts and incomplete private independent-audit packets are excluded from this manuscript deposit.',
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
            payload=path.read_bytes()
            req=urllib.request.Request(url,data=payload,method='PUT',headers={
                    'Authorization':'Bearer '+os.environ['ZENODO_TOKEN'],
                    'Content-Type':'application/octet-stream','Content-Length':str(path.stat().st_size)})
            with urllib.request.urlopen(req,timeout=1800) as r: uploaded=json.load(r)
            del payload
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
    catalog=root/('ML_HISTORY_ARTIFACT_CATALOG.json' if extra_only else 'UPSTREAM_ARTIFACT_CATALOG.json')
    catalog.write_text(json.dumps(report,indent=2),encoding='utf-8')
    upload(catalog)
    note=root/'SOURCE_TO_ANALYSIS_README.txt'
    note.write_bytes(Path('reproducibility/zenodo_upstream_readme.txt').read_bytes())
    upload(note)
    state=request(base);metadata=state['metadata'].copy()
    metadata['description'] += ('<p>Earlier Grounding DINO pseudo-label and CLIP zero-shot artifacts are preserved separately in ML_HISTORY_ARTIFACT_CATALOG.json. '
        'They document development history, not current continuous-trait inference. The external pretrained weights themselves are not included.</p>' if extra_only else
        '<p>Upstream processing artifacts are additionally preserved: acquisition metadata and screening queue/images, '
        'the frozen detector package with best/last weights, dataset manifest and training diagnostics, historical trait outputs, '
        'and exhaustive merged continuous measurements. Unused categorical ML development and incomplete private independent-audit packets are excluded. '
        'See UPSTREAM_ARTIFACT_CATALOG.json for exact source runs, member hashes and remaining gaps. '
        'This is NOT certified zero-from-scratch reproduction: full original photographs and all training inputs or external model weights have not been verified complete. '
        'Image-containing packages and upstream data terms require rights review before any publication.</p>')
    request(base,'PUT',{'metadata':metadata})
    final=request(base)
    assert not final['submitted'] and final['state']=='unsubmitted'
    report['status']='UPSTREAM_DRAFT_UPLOADS_READBACK_VERIFIED'
    (out/'upstream_final.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
