"""One-record draft-only maintenance. No publication endpoint is implemented."""
import json
import os
import urllib.request
import urllib.parse
import hashlib
from pathlib import Path

BASE = 'https://zenodo.org/api/deposit/depositions/22703215'
OUT = Path('work/zenodo-receipt')
OUT.mkdir(parents=True, exist_ok=True)

def request(url, method='GET', data=None):
    assert url.startswith('https://zenodo.org/api/')
    headers = {'Authorization': 'Bearer ' + os.environ['ZENODO_TOKEN']}
    if data is not None:
        data = json.dumps(data).encode()
        headers['Content-Type'] = 'application/json'
    req = urllib.request.Request(url, data=data, headers=headers, method=method)
    with urllib.request.urlopen(req, timeout=120) as response:
        return json.load(response)

draft = request(BASE)
assert str(draft['id']) == '22703215'
assert not draft.get('submitted'), 'Refusing to modify a submitted record'
assert draft.get('state') == 'unsubmitted', 'Not an unpublished draft'
assert str(draft.get('conceptrecid')) == '22295790', 'Unexpected version family'
receipt = {k: draft.get(k) for k in ('id', 'state', 'submitted', 'conceptrecid', 'metadata', 'files')}
(OUT / 'before.json').write_text(json.dumps(receipt, indent=2), encoding='utf-8')
print(json.dumps(receipt, indent=2))
if os.environ.get('UPDATE_MODE') == 'update':
    from reproducibility.build_zenodo_data_archive import build, REV
    files = build()
    allowed = {p.name for p in files}
    assert all(f['filename'] in allowed for f in draft.get('files', [])), 'Unknown existing files; no deletion allowed'
    bucket = draft['links']['bucket']
    assert bucket.startswith('https://zenodo.org/api/files/')
    uploads = []
    for p in files:
        payload = p.read_bytes()
        digest = hashlib.sha256(payload).hexdigest()
        url = bucket + '/' + urllib.parse.quote(p.name)
        req = urllib.request.Request(url, data=payload, method='PUT', headers={
            'Authorization': 'Bearer ' + os.environ['ZENODO_TOKEN'],
            'Content-Type': 'application/octet-stream'})
        with urllib.request.urlopen(req, timeout=600) as r:
            uploaded = json.load(r)
        assert uploaded['checksum'] == 'md5:' + hashlib.md5(payload).hexdigest()
        req = urllib.request.Request(url, headers={'Authorization': 'Bearer ' + os.environ['ZENODO_TOKEN']})
        with urllib.request.urlopen(req, timeout=600) as r:
            readback = hashlib.sha256(r.read()).hexdigest()
        assert readback == digest, 'Readback differs'
        uploads.append({'filename':p.name,'bytes':len(payload),'sha256':digest,'authenticated_readback_verified':True})
        print('Uploaded and read back: ' + p.name, flush=True)
    metadata = draft['metadata'].copy()
    metadata['title'] = 'Azami Chapter 1 continuous capitulum phenotypes: analysis data and GitHub Actions artifacts'
    metadata['version'] = 'ch1-analysis-data-20260911-52c207a'
    metadata['description'] = (
        '<p>Permanent-archive preparation for the numerical data and GitHub Actions artifacts supporting the current Azami Chapter 1 manuscript. '
        'Nine original artifact ZIPs are retained byte-for-byte: four numerical input artifacts (continuous measurements, nine-predictor environment, '
        'broad-region lookup and placement trees), four current construct-analysis/result artifacts, and one earlier complete-18 upgrade provenance artifact. '
        'The frozen native-status table is supplied separately. ARTIFACT_CATALOG.json records artifact IDs, runs, source commits, file membership and SHA-256 hashes.</p>'
        '<p>Analysis code remains on <a href="https://github.com/zuizui0223/azami/tree/' + REV + '">GitHub at commit ' + REV + '</a>. '
        'No manuscript, private reviewer correspondence, code snapshot or original photograph collection is added. '
        'This preserves the inputs used by the seven-stage current numerical workflow; it is not a new image-measurement run or an independent biological validation.</p>'
        '<p>Current and historical content are distinguished in README.txt. Original artifact bytes can include superseded outputs and retired fields; '
        'these are retained as provenance, not promoted to current findings. The earlier published v2 record 10.5281/zenodo.22295791 is unchanged.</p>'
        '<p>Third-party data retain their original terms; the software MIT license does not relicense all data. '
        'This draft is NOT approved for publication: backbone-data redistribution terms and final record-level license scope remain to be confirmed. '
        'The inherited CC BY 4.0 metadata is retained pending that decision and is not blanket clearance for third-party material.</p>')
    metadata.setdefault('custom', {})['code:codeRepository'] = 'https://github.com/zuizui0223/azami/tree/' + REV
    related = metadata.get('related_identifiers', [])
    related = [r for r in related if r.get('relation') != 'isSupplementTo' or 'github.com/zuizui0223/azami' not in r.get('identifier','')]
    related.append({'identifier':'https://github.com/zuizui0223/azami/tree/'+REV,'relation':'isSupplementTo','scheme':'url'})
    metadata['related_identifiers'] = related
    request(BASE, 'PUT', {'metadata':metadata})
    after = request(BASE)
    assert not after['submitted'] and after['state'] == 'unsubmitted'
    assert {f['filename'] for f in after['files']} == allowed
    assert after['metadata']['version'] == metadata['version']
    assert after['metadata']['creators'] == draft['metadata']['creators']
    (OUT/'after.json').write_text(json.dumps({k:after.get(k) for k in ('id','state','submitted','metadata','files')},indent=2),encoding='utf-8')
    (OUT/'verification.json').write_text(json.dumps({'status':'DRAFT_UPDATED_READBACK_VERIFIED','published':False,'code_commit':REV,
        'files':uploads,'artifact_count':9,'license_scope_pending':True},indent=2),encoding='utf-8')
    print('DRAFT_UPDATED_READBACK_VERIFIED; NOT PUBLISHED')
