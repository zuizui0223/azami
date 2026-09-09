"""One protected all27 point-estimate run; never requests images or publishes data."""
from __future__ import annotations

import argparse
from importlib.metadata import version
import json
import os
from pathlib import Path
import platform

from .all27_exploration import definition, run
from .private_replay import snapshot
from .protected_artifacts import DraftStore, new_json, pack, require, unpack
from .recover_native_source_authority import private_directory
from .workflow import ROOT, canonical_digest, digest, text_digest


def validate_batch(batch):
    spec, _ = definition()
    require(batch['status'] == 'single_partial_checkpoint_all27_point_exploration', 'Wrong execution scope')
    require(batch['exploratory_contract_canonical_sha256'] == canonical_digest(spec), 'All27 contract changed')
    require(batch['image_requests_authorized'] is False and batch['primary_inference_authorized'] is False,
            'Exploratory execution cannot authorize images or primary inference')
    require(batch['source_asset']['source_images_included'] is False, 'Numerical input required')
    for name, sha in batch['implementation_sha256_text_lf'].items():
        path = ROOT / name
        require(path.resolve().is_relative_to(ROOT) and text_digest(path) == sha, 'Pinned implementation changed')
    required = {'analysis/v3/all27_exploration.py', 'analysis/v3/cloud_all27_exploration.py',
                'analysis/v3/module_ecology.py', 'analysis/v3/joint_partial_pooling.py',
                'analysis/v3/environment_model.py', 'analysis/v3/nuisance_design.py'}
    require(required <= set(batch['implementation_sha256_text_lf']), 'Missing core execution identities')
    return batch


def execute(batch_path: Path, out: Path):
    batch = validate_batch(json.loads(batch_path.read_text(encoding='utf-8')))
    require(platform.python_version() == batch['runtime']['python'], 'Python runtime differs')
    require(all(version(k) == v for k, v in batch['runtime']['packages'].items()), 'Numerical runtime differs')
    run_id, attempt = os.environ.get('GITHUB_RUN_ID', ''), os.environ.get('GITHUB_RUN_ATTEMPT', '')
    require(run_id.isdigit() and attempt.isdigit(), 'This entrypoint requires GitHub Actions')
    require(os.environ.get('GITHUB_REPOSITORY') == batch['repository']
            and os.environ.get('GITHUB_REF') == 'refs/heads/' + batch['branch'], 'Wrong repository or branch')
    out = private_directory(out)
    require(not out.exists(), 'Preserve previous cloud result')
    out.mkdir(parents=True)
    store = DraftStore(batch)
    try:
        print(json.dumps({'stage': 'restore_protected_input'}), flush=True)
        store.download(batch['source_asset'], out / 'source.zip')
        unpack(out / 'source.zip', out / 'source', batch['source_asset'])
        prepared = out / 'source/restored/prepared'
        require(digest(prepared / 'preparation_report.json') == batch['preparation_report_sha256'], 'Prepared source receipt differs')
        print(json.dumps({'stage': 'all27_point_estimation'}), flush=True)
        summary = run(prepared, out / 'result')
        require(summary['endpoints_attempted'] == 27 and summary['coefficient_slots_retained'] == 729,
                'Incomplete all27 attempt inventory')
        paths = sorted((out / 'result').iterdir())
        selection = out / 'selection_private.json'
        new_json(selection, {'schema_version': 1, 'files': [
            {'name': 'result/' + p.name, 'path': str(p.resolve()), 'sha256': digest(p)} for p in paths if p.is_file()]})
        snapshot(selection, out / 'snapshot')
        asset = pack(out / 'snapshot', out / 'output.zip')
        print(json.dumps({'stage': 'protect_numerical_results'}), flush=True)
        asset.update(store.upload(out / 'output.zip', f'v3-all27-point-{run_id}-{attempt}.zip'))
        # Keep this durable identity if return confirmation subsequently fails.
        new_json(out / 'protected_output_receipt.json', asset)
        print(json.dumps({'stage': 'verify_protected_return'}), flush=True)
        store.download(asset, out / 'returned.zip')
        returned = unpack(out / 'returned.zip', out / 'returned', asset)
        require(digest(out / 'returned/restored/result/public_report.json') == digest(out / 'result/public_report.json'),
                'Protected returned result differs')
        public = {'status': 'ALL27_POINT_EXPLORATION_PROTECTED_ROUNDTRIP_VERIFIED',
                  'run_id': run_id, 'run_attempt': attempt, 'execution_commit': os.environ.get('GITHUB_SHA'),
                  'batch_sha256': digest(batch_path), 'output_asset': asset,
                  'returned_files': returned['files'], 'returned_bytes': returned['restored_bytes'],
                  'draft_verified_unpublished': True, 'anonymous_release_and_asset_requests': '404',
                  'analysis_summary': summary, 'primary_inference_authorized': False,
                  'limit': 'Partial source-ordered checkpoint. Point estimates only, not biological validation or formal inference.'}
        new_json(out / 'public_report.json', public)
        return public
    finally:
        store.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    try:
        result = execute(args.batch, args.out)
        print(json.dumps({k: result[k] for k in ['status', 'run_id', 'run_attempt']}))
    except BaseException as error:
        # No response values, taxon identities or private source paths in public logs.
        args.out.mkdir(parents=True, exist_ok=True)
        target = args.out / 'public_failure.json'
        if not target.exists():
            new_json(target, {'status': 'ALL27_EXPLORATION_INCOMPLETE', 'error_type': type(error).__name__,
                              'full_attempt_completion_claimed': False, 'primary_inference_authorized': False})
        print(json.dumps({'status': 'ALL27_EXPLORATION_INCOMPLETE', 'error_type': type(error).__name__}))
        raise SystemExit(1) from None


if __name__ == '__main__':
    main()
