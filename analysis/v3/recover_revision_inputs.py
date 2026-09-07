#!/usr/bin/env python3
"""Recover exact frozen inputs and per-head records into an external work directory.

Read-only remote access. No manuscript, reviewer text, images or fitted results
are committed. Recovery is not a completed measurement-stability analysis.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
import os
import shutil
import tempfile
import urllib.error
import urllib.request
import zipfile
from pathlib import Path

EXPECTED = {
    9612943217: '101e996b638996a0c5ae79d358bf51293c3585f0e84c4a961b91dcbedf96211e',
    9632715852: '51e7a26b5bd09e030b67b9342586699abaaf46e630f45b6bb4ee7bfc9152ced6',
    8269246732: '5f18b42d18cfcb81691c38ce0f04bcef754e6a67382025ea90110dbc50ae194b',
    8225059018: '9af18f9a2595f6c43dae4b933833f1721d153dcacef5b7c80df543f16235b580',
}
FULL_ENV_SHA = 'e242aa7ce69d12b11937c1335e84b9638799c50b42ef36b95725e77190df98e7'

class NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return None

def api(path: str):
    repo = os.environ.get('GH_REPOSITORY', 'zuizui0223/azami')
    if repo != 'zuizui0223/azami':
        raise ValueError('This recovery is pinned to zuizui0223/azami')
    req = urllib.request.Request('https://api.github.com/repos/' + repo + path,
        headers={'Authorization': 'Bearer ' + os.environ['GH_TOKEN'],
                 'Accept': 'application/vnd.github+json', 'User-Agent': 'azami-v3-input-recovery'})
    return urllib.request.build_opener(NoRedirect).open(req, timeout=120)

def get_json(path: str):
    with api(path) as response:
        return json.load(response)

def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def download(artifact_id: int, destination: Path, expected: str | None = None):
    meta = get_json(f'/actions/artifacts/{artifact_id}')
    if meta.get('expired'):
        raise ValueError(f'Artifact {artifact_id} expired')
    try:
        response = api(f'/actions/artifacts/{artifact_id}/zip')
    except urllib.error.HTTPError as error:
        if error.code not in (301, 302, 303, 307, 308):
            raise
        target = error.headers['Location']
        if not target.startswith('https://'):
            raise ValueError('Insecure archive redirect')
        # No GitHub credential is forwarded to the signed archive host.
        response = urllib.request.urlopen(target, timeout=180)
    with response, destination.open('wb') as f:
        shutil.copyfileobj(response, f, 1024 * 1024)
    digest = sha(destination)
    required = expected or str(meta.get('digest', '')).removeprefix('sha256:')
    if not required or digest != required:
        raise ValueError(f'Artifact {artifact_id} hash mismatch')
    return {'id': artifact_id, 'name': meta['name'], 'sha256': digest,
            'size_in_bytes': destination.stat().st_size, 'workflow_run': meta.get('workflow_run')}

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--out-dir', type=Path, required=True)
    args = p.parse_args()
    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    report = {'status': 'started', 'archives': [], 'inventories': {}, 'errors': [],
              'measurement_stability_analysis_completed': False}
    def checkpoint():
        (out / 'recovery_report.json').write_text(json.dumps(report, indent=2) + '\n')
    with tempfile.TemporaryDirectory() as tmp_name:
        tmp = Path(tmp_name)
        primary = tmp / 'primary.zip'
        report['archives'].append(download(9612943217, primary, EXPECTED[9612943217]))
        with zipfile.ZipFile(primary) as z:
            core = z.read('environment/strict_spatial_chelsa.csv')
            (out / 'strict_spatial_chelsa.csv').write_bytes(core)
            ids = {r['obs_id'] for r in csv.DictReader(io.StringIO(core.decode('utf-8-sig')))}
            assert len(ids) == 46276
        report['primary_observation_ids'] = len(ids)
        checkpoint()
        for aid in (8269246732, 8225059018):
            archive = tmp / f'{aid}.zip'
            try:
                report['archives'].append(download(aid, archive, EXPECTED[aid]))
                inventory = []
                with zipfile.ZipFile(archive) as z:
                    for member in z.infolist():
                        item = {'name': member.filename, 'bytes': member.file_size}
                        if member.filename.endswith('.csv'):
                            with z.open(member) as raw:
                                reader = csv.DictReader(io.TextIOWrapper(raw, encoding='utf-8-sig'))
                                cols = reader.fieldnames or []
                                item['columns'] = cols
                                if 'orientation_angle_degrees' in cols:
                                    keep = [c for c in cols if c not in {'user_login','user_id','photo_attribution'}]
                                    target = out / f'{aid}_{Path(member.filename).name}'
                                    n = 0
                                    with target.open('w', encoding='utf-8', newline='') as f:
                                        writer = csv.DictWriter(f, fieldnames=keep)
                                        writer.writeheader()
                                        for row in reader:
                                            if 'obs_id' in cols and row['obs_id'] not in ids:
                                                continue
                                            writer.writerow({k: row.get(k,'') for k in keep})
                                            n += 1
                                    item['recovered_rows'] = n
                                    item['recovered_file'] = target.name
                        inventory.append(item)
                report['inventories'][str(aid)] = inventory
            except Exception as e:
                report['errors'].append({'artifact': aid, 'error': type(e).__name__ + ': ' + str(e)})
            finally:
                archive.unlink(missing_ok=True)
                checkpoint()
        # Locate an exact all-nine environment copy by identity, never by outcome.
        artifacts = []
        for page in range(1, 21):
            chunk = get_json(f'/actions/artifacts?per_page=100&page={page}').get('artifacts', [])
            artifacts.extend(chunk)
            if len(chunk) < 100:
                break
        (out / 'artifact_metadata_index.json').write_text(json.dumps(artifacts, indent=2) + '\n')
        candidates = [a for a in artifacts if not a['expired'] and a['size_in_bytes'] < 200_000_000
                      and any(w in a['name'].lower() for w in ('full27','full-27','frozen-analysis','numerical-rebuild'))
                      and not a['name'].startswith('ch1-v3-')]
        report['environment_archive_candidates'] = [{'id': a['id'], 'name': a['name']} for a in candidates]
        report['full_environment_recovered'] = False
        for meta in candidates[:12]:
            archive = tmp / 'candidate.zip'
            try:
                receipt = download(meta['id'], archive)
                with zipfile.ZipFile(archive) as z:
                    for name in z.namelist():
                        if name.endswith('strict_spatial_chelsa_process.csv'):
                            data = z.read(name)
                            if hashlib.sha256(data).hexdigest() == FULL_ENV_SHA:
                                (out / 'strict_spatial_chelsa_process.csv').write_bytes(data)
                                report['full_environment_recovered'] = True
                                report['full_environment_source'] = dict(receipt, member=name)
                                break
            except Exception as e:
                report['errors'].append({'environment_candidate': meta['id'], 'error': type(e).__name__ + ': ' + str(e)})
            finally:
                archive.unlink(missing_ok=True)
                checkpoint()
            if report['full_environment_recovered']:
                break
    report['status'] = 'recovery_complete_with_reported_gaps'
    checkpoint()
    print(json.dumps({k: v for k,v in report.items() if k != 'inventories'}, indent=2))

if __name__ == '__main__':
    main()
