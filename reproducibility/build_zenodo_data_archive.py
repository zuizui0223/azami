"""Preserve selected current-paper Actions artifacts byte-for-byte, without code snapshot."""
import csv
import hashlib
import io
import json
import subprocess
import zipfile
from pathlib import Path
from reproducibility.run_current_analysis import INPUTS, NATIVE_REF, NATIVE_SHA, native_bytes, verified

ROOT = Path(__file__).resolve().parents[1]
REV = '52c207a64a09b243e021dc6e7b20598ab4ab7f02'
REF = ROOT / 'reproducibility/current_reference/manifest.json'

def gh(path):
    return subprocess.check_output(['gh', 'api', 'repos/zuizui0223/azami/' + path])

def build():
    out = ROOT / 'work/zenodo-data'
    out.mkdir(parents=True, exist_ok=True)
    expected = {aid: (sha, 'input_' + role) for role, (aid, sha, _) in INPUTS.items()}
    reference = json.loads(REF.read_text(encoding='utf-8'))
    for row in reference['files']:
        expected[row['artifact']] = (row['archive_sha256'], 'current_results')
    expected[10135139012] = ('65438082976af8aa0135edc964a59f531e8763fde25f92739ae31c5056429f7e', 'pre_scale_contrast_provenance')
    entries, catalog = {}, []
    for aid, (sha, role) in sorted(expected.items()):
        meta = json.loads(gh(f'actions/artifacts/{aid}'))
        assert not meta['expired']
        assert meta['digest'] == 'sha256:' + sha
        payload = verified(gh(f'actions/artifacts/{aid}/zip'), sha, str(aid))
        name = f'artifacts/artifact-{aid}-{role}.zip'
        with zipfile.ZipFile(io.BytesIO(payload)) as z:
            assert z.testzip() is None
            members = [{'path': n, 'bytes': z.getinfo(n).file_size,
                        'sha256': hashlib.sha256(z.read(n)).hexdigest()}
                       for n in z.namelist() if not n.endswith('/')]
            for row in reference['files']:
                if row['artifact'] == aid:
                    verified(z.read(row['member']), row['sha256'], row['member'])
            for _, (input_id, _, input_members) in INPUTS.items():
                if input_id == aid:
                    for member, (_, digest) in input_members.items():
                        verified(z.read(member), digest, member)
        entries[name] = payload
        catalog.append({'artifact_id': aid, 'role': role, 'file': name,
                        'sha256': sha, 'bytes': len(payload), 'name': meta['name'],
                        'run': meta['workflow_run'], 'created_at': meta['created_at'],
                        'run_url': f"https://github.com/zuizui0223/azami/actions/runs/{meta['workflow_run']['id']}",
                        'members': members})
        print(f'Verified artifact {aid}: {len(payload)} bytes', flush=True)
    subprocess.run(['git','fetch','--no-tags','origin',
                    'refs/tags/azami-ch1-v2-2026-08-27:refs/tags/azami-ch1-v2-2026-08-27'],cwd=ROOT,check=True)
    entries['supplemental_inputs/native_status.csv'] = native_bytes(subprocess.check_output(['git','show',NATIVE_REF],cwd=ROOT))
    documents = ['analysis/ch1/chelsa_process_environment_sources.json',
                 'analysis/ch1/capitulum_environment_blocks_contract.json',
                 'analysis/ch1/native_range_sensitivity_contract.json',
                 'reproducibility/current_reference/manifest.json']
    for name in documents:
        entries['source_metadata/' + Path(name).name] = (ROOT/name).read_bytes()
    entries['ARTIFACT_CATALOG.json'] = json.dumps({'code_commit':REV,
        'scientific_reference_commit':reference['scientific_reference_commit'],
        'artifacts':catalog,'supplemental_native_status':{'git_ref':NATIVE_REF,'sha256':NATIVE_SHA},
        'retired_fields_not_current_evidence': {
            'construct_scale_integration_report.json':['environment_signature_alignment'],
            'construct_scale_upgrade_report.json':['integration_environment_coupling']},
        'scope':'Current numerical input and construct-analysis artifact preservation; not all historical CI runs or original images.'},indent=2).encode()
    entries['README.txt'] = (ROOT/'reproducibility/zenodo_data_readme.txt').read_bytes()
    manifest = {n:{'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()} for n,b in sorted(entries.items())}
    entries['MANIFEST.json'] = json.dumps(manifest,indent=2).encode()
    path = out/'Azami_Chapter1_analysis_data_artifacts_20260911.zip'
    with zipfile.ZipFile(path,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for n,b in sorted(entries.items()):
            info=zipfile.ZipInfo(n,date_time=(2026,9,11,0,0,0))
            info.compress_type=zipfile.ZIP_DEFLATED
            z.writestr(info,b)
    with zipfile.ZipFile(path) as z:
        for n,m in manifest.items(): verified(z.read(n),m['sha256'],n)
    checksum = out/'SHA256SUMS.txt'
    checksum.write_text(hashlib.sha256(path.read_bytes()).hexdigest()+'  '+path.name+'\n',encoding='utf-8')
    (out/'ARTIFACT_CATALOG.json').write_bytes(entries['ARTIFACT_CATALOG.json'])
    (out/'README.txt').write_bytes(entries['README.txt'])
    return [path,checksum,out/'ARTIFACT_CATALOG.json',out/'README.txt']

if __name__ == '__main__': build()
