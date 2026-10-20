from pathlib import Path
import json
import hashlib
import importlib
import pytest

from reproducibility.run_current_analysis import commands, verified, native_bytes
from reproducibility.validate_current_analysis import REFERENCE, compare, validate, active_scope
from reproducibility import run_current_analysis as runner

ROOT=Path(__file__).resolve().parents[1]


@pytest.mark.parametrize('corruption', [None, 'outer', 'embedded', 'native'])
def test_zenodo_restore_checks_all_bytes_before_extraction(tmp_path, monkeypatch, corruption):
    import io
    import zipfile
    archives, inputs = tmp_path / 'archives', tmp_path / 'inputs'
    archives.mkdir()
    inputs.mkdir()
    original, native = b'frozen nested archive', b'native\n'
    sha = lambda data: hashlib.sha256(data).hexdigest()
    monkeypatch.setattr(runner, 'INPUTS', {'continuous': (123, sha(original), {})})
    monkeypatch.setattr(runner, 'NATIVE_SHA', sha(native))
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, 'w') as z:
        z.writestr('inputs/artifact-123-continuous.zip',
                   b'changed' if corruption == 'embedded' else original)
        z.writestr('inputs/observation_native_status.csv',
                   b'changed' if corruption == 'native' else native)
        z.writestr('../unexpected.txt', b'never extracted')
    data = buffer.getvalue()
    monkeypatch.setattr(runner, 'ZENODO_SHA', '0' * 64 if corruption == 'outer' else sha(data))
    calls = []
    def fetch(url, timeout):
        calls.append(url)
        return io.BytesIO(data)
    monkeypatch.setattr(runner.urllib.request, 'urlopen', fetch)
    if corruption:
        with pytest.raises(ValueError):
            runner.restore_zenodo(archives, inputs)
        assert not list(inputs.iterdir())
        assert not (archives / 'artifact-123-continuous.zip').exists()
    else:
        runner.restore_zenodo(archives, inputs)
        runner.restore_zenodo(archives, inputs)
        assert calls == [runner.ZENODO_URL]  # Verified cache is reusable offline.
        assert (archives / 'artifact-123-continuous.zip').read_bytes() == original
        assert (inputs / 'native_status.csv').read_bytes() == native
    assert not (tmp_path / 'unexpected.txt').exists()


def test_current_entrypoints_import_and_counts_stay_fixed():
    cmds=commands(Path('inputs'),Path('results'))
    assert len(cmds)==8
    for cmd in cmds:
        assert importlib.import_module(cmd[0])
    assert cmds[0][-4:]==['--permutations','9999','--seed','20260910']
    assert cmds[-1][0]=='analysis.v3.run_rv_estimator_validity'
    assert cmds[-1][-2:]==['--seed','20260915']
    assert '--native-only' not in str(cmds)


def test_reference_files_match_original_artifacts():
    manifest=json.loads((REFERENCE/'manifest.json').read_text())
    assert len(manifest['files'])==16
    for row in manifest['files']:
        assert hashlib.sha256((REFERENCE/row['path']).read_bytes()).hexdigest()==row['sha256']
    assert validate(REFERENCE)['aggregate_files_compared']==16


def test_verifier_fails_closed_and_detects_changed_results(tmp_path):
    with pytest.raises(ValueError): verified(b'wrong','0'*64,'input')
    with pytest.raises(ValueError): native_bytes(b'wrong')
    with pytest.raises(FileNotFoundError): validate(tmp_path)
    with pytest.raises(AssertionError): compare({'beta':-.345},{'beta':.345})
    with pytest.raises(AssertionError): compare({'pass':False},{'pass':True})
    with pytest.raises(AssertionError): compare([1,2],[1])
    result=active_scope({'median_within_rv':.1,'integration_environment_coupling':{}},
                        'upgrade/construct_scale_upgrade_report.json')
    assert result=={'median_within_rv':.1}


def test_legacy_paths_are_explicit_and_complete():
    mapping=json.loads((ROOT/'reproducibility/legacy_path_map.json').read_text())['moves']
    for old,new in mapping.items():
        assert not (ROOT/old).exists(),old
        assert new.startswith('legacy/')
        assert (ROOT/new).is_file(),new
    assert not (ROOT/'ch1_global').exists() or not list((ROOT/'ch1_global').rglob('*.py'))
    assert not list((ROOT/'analysis').glob('*.py'))


def test_current_runbook_does_not_inherit_v2_publication_status():
    text=(ROOT/'reproducibility/README.md').read_text(encoding='utf-8')
    assert 'run_current_analysis' in text
    assert 'do not certify current v3 public availability' in text
    assert 'C:\\Users' not in text
    assert (ROOT/'LICENSE').is_file()

def test_current_reproducibility_docs_match_eight_stage_boundary():
    runbook = (ROOT / "reproducibility/CURRENT_ANALYSIS.md").read_text(encoding="utf-8")
    assert "runs all **eight** current numerical stages" in runbook
    assert "compares **16** aggregate outputs" in runbook
    boundary = (ROOT / "reproducibility/UPSTREAM_DATA_BOUNDARY_20260926.md").read_text(encoding="utf-8")
    assert "Starting point of the current manuscript replay" in boundary
    assert "current v3 runner does not redownload photographs" in boundary
    assert "8 numerical stages" in boundary
    assert "16 checksum-locked files" in boundary
    zenodo = (ROOT / "reproducibility/ZENODO_UPDATE_AUDIT.md").read_text(encoding="utf-8")
    assert "Zenodo is the **durable analysis-input store**" in zenodo
    assert "8-stage numerical replay" in zenodo
    assert "16-file reference comparison" in zenodo

def test_current_requirements_match_completed_replay_environment():
    req = (ROOT / "reproducibility/requirements-current.txt").read_text(encoding="utf-8")
    assert "numpy==2.4.6" in req
    assert "pandas==3.0.5" in req
    assert "scipy==1.17.1" in req
    assert "statsmodels==0.14.6" in req
    assert "biopython==1.88" in req
    assert "biopython>=" not in req


def test_native_first_27_endpoint_ledger_matches_frozen_ontology_and_coverage():
    """Guard the published aggregate source gate; never replace missing by zero traits."""
    import csv
    ledger_path = ROOT / 'reproducibility/native_first_27_endpoint_coverage_20261020.json'
    data = json.loads(ledger_path.read_text(encoding='utf-8'))
    with (ROOT / 'analysis_outputs/v2_full27_environment_atlas_2026-08-27/v2_full27_endpoint_inventory.csv').open(encoding='utf-8-sig', newline='') as handle:
        v2 = {r['endpoint_id']: r for r in csv.DictReader(handle)}
    items = data['endpoints']
    assert data['schema_version'] == 'native_first_27_endpoint_coverage_20261020_v2'
    assert len(items) == len(v2) == 27
    assert len({r['endpoint_id'] for r in items}) == 27
    assert {r['endpoint_id'] for r in items} == set(v2)
    d = data['denominators']
    assert d['original_pre_detector_observations'] == 665115
    assert d['native_primary_observations'] == 313350
    assert d['native_wcvp_accepted_keys'] == 306
    assert d['native_with_inherited_measurement_row'] + d['native_without_inherited_measurement_row'] == d['native_primary_observations']
    assert d['native_climate_all_nine_complete'] <= d['native_climate_four_anchor_complete'] <= d['native_primary_observations']
    p = data['provenance']
    assert p['artifact_id'] == 12257893990
    assert p['receipt_sha256'] == '246a5ea9643617c6dc766e0b8d27a444a5496afe150d7768032740cf6a53f0a6'
    assert p['existing_joint_model_run_id'] == 39416125522
    assert p['existing_joint_model_response_sha256'] == '707e3cf0041d127e085bfa47ef83a6ff89a1ad84740cd4ec6c09a893bbe1c517'
    statuses = {'finite': 0, 'eligible': 0, 'unexecuted': 0, 'finite_unusable': 0}
    for row in items:
        source = v2[row['endpoint_id']]
        finite = row['n_finite_native']
        eligible = row['n_eligible_native']
        taxa = row['n_eligible_native_taxa']
        climate = row['n_eligible_with_four_climate']
        assert 0 <= climate <= eligible <= finite <= d['native_primary_observations']
        assert 0 <= taxa <= d['native_wcvp_accepted_keys']
        if eligible == 0:
            assert taxa == 0
        if finite:
            statuses['finite'] += 1
        if eligible:
            statuses['eligible'] += 1
        if source['measurement_status'] == 'unexecuted_no_measurement':
            statuses['unexecuted'] += 1
            assert finite == eligible == taxa == 0
        if finite and not eligible:
            statuses['finite_unusable'] += 1
            assert row['status'] == 'technical_validation_only_not_analysis_eligible'
    assert statuses == {'finite': 22, 'eligible': 18, 'unexecuted': 5, 'finite_unusable': 4}
    assert data['status_summary'] == {
        'registered_endpoints': 27, 'have_any_numeric_native': 22,
        'currently_computationally_eligible_native': 18,
        'finite_but_zero_eligible': 4, 'unexecuted_no_measurement': 5,
    }
    assert next(x for x in items if x['endpoint_id']=='corolla_lab_chroma')['n_eligible_native_taxa']==244
    assert next(x for x in items if x['endpoint_id']=='capitulum_outline_solidity')['n_eligible_native_taxa']==243
    assert next(x for x in items if x['endpoint_id']=='bract_projection_p95')['n_eligible_native_taxa']==111
