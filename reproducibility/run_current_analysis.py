"""Single numerical replay of the current analysis; no image acquisition.

Only path orchestration and input verification live here. Statistical methods,
seeds, cohort choices and permutation counts remain in the original modules.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import urllib.request
import zipfile

ROOT = Path(__file__).resolve().parents[1]
REPOSITORY = 'zuizui0223/azami'
INPUTS = {
    'continuous': (9612943217, '101e996b638996a0c5ae79d358bf51293c3585f0e84c4a961b91dcbedf96211e', {
        'universe/continuous_trait_universe_observation_long.csv': ('traits.csv', 'd775794f2bce2dfd0c1f63c5c8e01778c518f6eeb327bf0d9944045143a02344')}),
    'environment': (9633419268, 'd7c0c466f55b67695d06ae46c21a6452dbe6cfd92a52db8042caa200429e97f4', {
        'process_environment/strict_spatial_chelsa_process.csv': ('environment.csv', 'f86a3418e9b21453026bba1aaf350b061f03048949547e2092602620af98cbf6')}),
    'spatial': (8983877726, '151161f926dbc92d35832a0bd71622e4f1f3018c9047a7ba8b4ad94926273dca', {
        'spatial_regions/broad_region_lookup.csv': ('regions.csv', '085c4e8d45ceb34d32c6c961675ce74a4f0a33580f6cdd8ecd2ff1800a6364ff')}),
    'historical': (8227254443, '499061e7a49f9455cf8c367fe26e313b7e0e33b2280d2354717e61a90ea8c6bc', {
        'historical_trees/gbotb_lcvp_scenario1.tre': ('trees/gbotb_lcvp_scenario1.tre', '8ef5d5ea5f4e0c2f166071244a838cae77a0fe582817d729bba0b36f6b5ccd92'),
        'historical_trees/gbotb_lcvp_scenario3.tre': ('trees/gbotb_lcvp_scenario3.tre', '8ef5d5ea5f4e0c2f166071244a838cae77a0fe582817d729bba0b36f6b5ccd92'),
        'historical_trees/gbotb_lcvp_scenario2_randomized.trees': ('trees/gbotb_lcvp_scenario2_randomized.trees', '82655f79297e44a6630a599d8b0a1dc6f85e792812b8b01f2234e567a478e3af')})}
NATIVE_SHA = 'c01eeb9ff245d7f73da1a12fa4eede904dd9770467655f20e3d85de2ac8dd84a'
NATIVE_REF = 'azami-ch1-v2-2026-08-27:analysis_outputs/native_range_sensitivity_v1/observation_native_status.csv'
ZENODO_URL = 'https://zenodo.org/api/records/22703215/files/azami_ch1_v3_analysis_inputs.zip/content'
ZENODO_SHA = 'be57e9ba80773c97afdc200f9a3982eaa85044cfc6791b62187ccfcd41f54438'


def verified(data: bytes, expected: str, name: str) -> bytes:
    actual = hashlib.sha256(data).hexdigest()
    if actual != expected:
        raise ValueError(f'{name}: SHA-256 mismatch {actual}; expected {expected}')
    return data


def native_bytes(data: bytes) -> bytes:
    # Exactly the transport-normalization alternatives used in the frozen CI.
    lf = data.replace(b'\r\n', b'\n')
    for candidate in (data, lf, lf.replace(b'\n', b'\r\n')):
        if hashlib.sha256(candidate).hexdigest() == NATIVE_SHA:
            return candidate
    raise ValueError('Native-status hash mismatch after permitted newline normalization')


def restore_zenodo(archives: Path, inputs: Path) -> None:
    """Restore only checksum-locked inputs from the published v3 record."""
    bundle = archives / 'azami_ch1_v3_analysis_inputs.zip'
    if bundle.exists():
        data = verified(bundle.read_bytes(), ZENODO_SHA, bundle.name)
    else:
        print('Downloading published Zenodo v3 numerical inputs (22703215)', flush=True)
        with urllib.request.urlopen(ZENODO_URL, timeout=120) as response:
            data = verified(response.read(), ZENODO_SHA, 'Zenodo v3 bundle')
        bundle.write_bytes(data)
    payloads = {}
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        for role, (artifact, digest, _) in INPUTS.items():
            name = f'artifact-{artifact}-{role}.zip'
            payloads[archives / name] = verified(z.read(f'inputs/{name}'), digest, name)
        payloads[inputs / 'native_status.csv'] = native_bytes(
            z.read('inputs/observation_native_status.csv'))
    # Validate every required input before writing any extracted files.
    for path, payload in payloads.items():
        path.write_bytes(payload)


def prepare(archives: Path, inputs: Path, download: bool) -> dict:
    archives.mkdir(parents=True, exist_ok=True)
    inputs.mkdir(parents=True, exist_ok=True)
    missing = any(not (archives / f'artifact-{artifact}-{role}.zip').exists()
                  for role, (artifact, _, _) in INPUTS.items())
    if download and (missing or not (inputs / 'native_status.csv').exists()):
        restore_zenodo(archives, inputs)
    receipt = {}
    for role, (artifact, digest, members) in INPUTS.items():
        archive = archives / f'artifact-{artifact}-{role}.zip'
        if not archive.exists():
            if not download:
                raise FileNotFoundError(f'{archive}; supply verified archive or use --download')
            raise FileNotFoundError(f'Published bundle did not restore {archive}')
        data = verified(archive.read_bytes(), digest, archive.name)
        with zipfile.ZipFile(io.BytesIO(data)) as z:
            for member, (relative, sha) in members.items():
                payload = verified(z.read(member), sha, member)
                out = inputs / relative  # only fixed paths above, never archive-provided paths
                out.parent.mkdir(parents=True, exist_ok=True)
                out.write_bytes(payload)
                receipt[relative] = sha
    native = inputs / 'native_status.csv'
    if native.exists():
        native.write_bytes(native_bytes(native.read_bytes()))
    else:
        check = subprocess.run(['git', 'cat-file', '-e', NATIVE_REF], cwd=ROOT, capture_output=True)
        if check.returncode:
            if not download:
                raise FileNotFoundError('Native status missing; fetch the immutable tag or use --download')
            subprocess.run(['git', 'fetch', '--no-tags', 'origin',
                            'refs/tags/azami-ch1-v2-2026-08-27:refs/tags/azami-ch1-v2-2026-08-27'], cwd=ROOT, check=True)
        native.write_bytes(native_bytes(subprocess.check_output(['git','show',NATIVE_REF],cwd=ROOT)))
    receipt['native_status.csv'] = NATIVE_SHA
    (inputs/'verified_inputs.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    return receipt


def commands(inputs: Path, out: Path) -> list[list[str]]:
    common = ['--traits',str(inputs/'traits.csv'),'--environment',str(inputs/'environment.csv')]
    axes = ['--axis-among',str(out/'axes/biological_axes_among_min5.csv'),
            '--axis-within',str(out/'axes/biological_axes_within.csv')]
    seed = ['--seed','20260910']
    return [
        ['analysis.v3.run_biological_axis_reanalysis',*common,'--out-dir',str(out/'axes'),'--permutations','9999',*seed],
        ['analysis.v3.run_biological_axis_sensitivity_chain',*common,*axes,'--regions',str(inputs/'regions.csv'),
         '--native-status',str(inputs/'native_status.csv'),'--tree-dir',str(inputs/'trees'),
         '--out-dir',str(out/'sensitivity'),'--spatial-permutations','999','--moran-permutations','999',
         '--minimum-taxa-historical','30',*seed],
        ['analysis.v3.run_construct_scale_integration_entrypoint',*common,*axes,'--out-dir',str(out/'integration'),
         '--minimum-paired-observations-per-taxon','5','--minimum-taxa','20','--qap-permutations','9999',*seed],
        ['analysis.v3.run_construct_scale_upgrade',*common,*axes,'--out-dir',str(out/'upgrade'),
         '--minimum-complete-observations-per-taxon','5','--bootstrap-replicates','1000','--permutations','9999',*seed],
        ['analysis.v3.run_construct_scale_contrast_summary','--pairwise',str(out/'upgrade/complete18_construct_pairwise.csv'),
         '--bootstrap',str(out/'upgrade/complete18_taxon_bootstrap.csv'),'--out',str(out/'upgrade/construct_scale_contrast_summary.json')],
        ['analysis.v3.run_assessability_selection_audit',*common,'--out-dir',str(out/'assessability')],
        ['analysis.v3.run_frozen_technical_error_stress',*common,'--technical-audit-summary',
         str(ROOT/'analysis/ch1/image_to_trait_automated_technical_audit_summary.json'),
         '--out-dir',str(out/'technical_stress'),'--replicates','2000',*seed],
        ['analysis.v3.run_rv_estimator_validity',*common,'--out-dir',str(out/'estimator_validity'),
         '--minimum-complete-observations-per-taxon','5','--equal-n-replicates','1000',
         '--null-permutations','499','--qap-permutations','9999','--seed','20260915'],
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir',type=Path,default=Path('work/current'))
    parser.add_argument('--archive-dir',type=Path)
    parser.add_argument('--download',action='store_true',help='Fetch missing numerical inputs from pinned public Zenodo v3; never images')
    parser.add_argument('--prepare-only',action='store_true')
    args = parser.parse_args()
    work = args.work_dir.resolve()
    inputs, out = work/'inputs', work/'results'
    prepare(args.archive_dir.resolve() if args.archive_dir else work/'archives',inputs,args.download)
    if args.prepare_only:
        print('Input verification PASS; numerical reproduction not yet run.')
        return
    out.mkdir(parents=True,exist_ok=True)
    commands_to_run = commands(inputs,out)
    total = len(commands_to_run)
    for i, cmd in enumerate(commands_to_run,1):
        print(f'[{i}/{total}] {cmd[0]}',flush=True)
        with (out/f'step-{i}.log').open('w',encoding='utf-8') as log:
            subprocess.run([sys.executable,'-m',*cmd],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,check=True)
    subprocess.run([sys.executable,'-m','reproducibility.validate_current_analysis','--results',str(out)],cwd=ROOT,check=True)


if __name__ == '__main__':
    main()
