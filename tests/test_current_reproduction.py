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


def test_native_space_gate_requires_eligible_joint_original_observations():
    """A real native-space input must be the frozen private source; only synthetic
    in-memory data are used here. Never promote a synthetic ellipsoid to evidence."""
    import numpy as np
    import pandas as pd
    from reproducibility.native_space_gate import analyze_native_space, CORE9

    rows = []
    for taxon in range(4):
        for j in range(32):
            v = {
                'obs_id': f'{taxon}-{j}', 'accepted_key': str(100+taxon),
                'latitude': -30+taxon + j*0.6, 'longitude': 30+taxon*.15 + j*.61,
                'chelsa_bio12': 950+taxon*160 + j*3,
                'chelsa_rsds_mean': 120+taxon*14+j*.5,
                'orientation_image_vertical_angle': 12+taxon*.5 + j*.11 + np.sin(j*.7)*3,
                'corolla_lab_lightness': 52+taxon + np.cos(j*.37),
                'corolla_lab_chroma': 20+taxon*3 + np.cos(j*.42)*6 + j*.09,
                'corolla_hue_sin': np.sin(j*.08),
                'corolla_hue_cos': np.cos(j*.08),
                'capitulum_outline_aspect_ratio': 1+taxon*.1 + np.sin(j*.31)*.13,
                'capitulum_outline_circularity': .6+np.cos(j*.19)*.03,
                'capitulum_outline_solidity': .85+np.sin(j*.18)*.02,
                'capitulum_width_profile_cv': .12+np.cos(j*.21)*.02,
            }
            for col in CORE9:
                v[col+'__eligible'] = True
            rows.append(v)
    frame = pd.DataFrame(rows)
    a=analyze_native_space(frame, rarefaction_n=20, replicates=12)
    b=analyze_native_space(frame.sample(frac=1, random_state=17), rarefaction_n=20, replicates=12)
    assert a == b
    assert a['complete_core3_observations'] == len(frame)
    assert a['complete_core9_observations'] == len(frame)
    assert a['core3_taxa_with_enough_distinct_cells'] == 4
    assert a['occupancy_status'] == 'descriptive_space_proxy_available'
    assert a['rarefied_gaussian_95pct_ellipsoid_volume_proxy']['median'] > 0
    assert 'obs_id' not in str(a) and 'latitude' not in str(a) and 'longitude' not in str(a)
    assert a['core3_scale_specific_covariance']['identity_max_abs_residual'] < 1e-10
    assert len(a['core3_scale_specific_covariance']['within_taxon_pairwise_correlations']) == 3
    assert len(a['core3_scale_specific_covariance']['among_taxon_pairwise_correlations']) == 3
    assert set(a['core3_among_covariance_by_climate_halves']) == {'chelsa_bio12', 'chelsa_rsds_mean'}
    assert all(v['status'] == 'insufficient_taxa_without_inferential_interpretation'
               for v in a['core3_among_covariance_by_climate_halves'].values())

    # A source-generated eligibility failure must remove only the affected module.
    frame.loc[0,'capitulum_outline_aspect_ratio__eligible'] = False
    c=analyze_native_space(frame, rarefaction_n=20, replicates=12)
    assert c['complete_core3_observations'] == len(frame)-1
    assert c['complete_core9_observations'] == len(frame)-1
    frame['capitulum_outline_aspect_ratio__eligible'] = frame['capitulum_outline_aspect_ratio__eligible'].astype(object)
    frame.loc[0,'capitulum_outline_aspect_ratio__eligible'] = 'unknown'
    with pytest.raises(ValueError, match="eligibility"):
        analyze_native_space(frame, rarefaction_n=20, replicates=12)
    with pytest.raises(ValueError, match="Observation IDs"):
        analyze_native_space(pd.concat([frame.iloc[:1],frame.iloc[:1]],ignore_index=True),rarefaction_n=6,replicates=2)


def test_native_space_scale_specific_covariance_sign_reversal():
    """A pooled summary hides independent within- and among-taxon directions."""
    import numpy as np
    import pandas as pd
    from reproducibility.native_space_gate import analyze_native_space, CORE9, CORE3

    rows = []
    for taxon in range(6):
        for k in range(28):
            # Within taxon, presentation and colour increase together. Across
            # taxa, presentation increases but colour decreases. Real nonlinearity
            # in all three variables keeps the ellipsoid nondegenerate.
            wave = np.sin(k * .35)
            angle = 10 + 4 * taxon + wave * 3 + k * .03
            row = {
                'obs_id': f'o-{taxon}-{k}', 'accepted_key': str(taxon+101),
                'latitude': 10+taxon*.5+k*.6, 'longitude': 35+taxon*.7+k*.72,
                'chelsa_bio12': 300+taxon*90+k,
                'chelsa_rsds_mean': 130+taxon*20+k*.3,
                'orientation_image_vertical_angle': angle,
                'corolla_lab_chroma': 70-7*taxon+wave*3.6+np.cos(k*.7)*.6,
                'capitulum_outline_aspect_ratio': 1.2+taxon*.04+np.cos(k*.25)*.08,
                'corolla_lab_lightness': 45+taxon*.1+np.sin(k*.23),
                'corolla_hue_sin': np.sin(k*.11), 'corolla_hue_cos': np.cos(k*.11),
                'capitulum_outline_circularity': .5+np.sin(k*.13)*.02,
                'capitulum_outline_solidity': .8+np.cos(k*.17)*.01,
                'capitulum_width_profile_cv': .15+np.sin(k*.19)*.02,
            }
            for name in CORE9:
                row[name+'__eligible'] = True
            rows.append(row)
    frame = pd.DataFrame(rows)
    result = analyze_native_space(frame, rarefaction_n=20, replicates=8)
    key = CORE3[0] + '__' + CORE3[1]
    cov = result['core3_scale_specific_covariance']
    assert cov['within_taxon_pairwise_correlations'][key] > .6
    assert cov['among_taxon_pairwise_correlations'][key] < -.8
    assert cov['among_minus_within_pairwise_correlations'][key] < -1.4
    assert cov['identity_max_abs_residual'] < 1e-10
    assert result['occupancy_status'] == 'descriptive_space_proxy_available'
