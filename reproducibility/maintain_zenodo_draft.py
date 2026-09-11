"""One-record draft-only maintenance. No publication endpoint is implemented."""
import json
import os
import urllib.request
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
    raise SystemExit('Update payload not yet admitted; inspection only')
