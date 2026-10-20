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
    assert data['schema_version'] == 'native_first_27_endpoint_coverage_20261020_v3'
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


def test_native_space_gate_rejects_superseded_and_unknown_environment_bytes(tmp_path, monkeypatch):
    from reproducibility import native_space_gate as gate
    import hashlib
    source = tmp_path / "native.csv.gz"
    source.write_bytes(b"raw-uncorrected-test")
    monkeypatch.setattr(gate, "RAW_UNCORRECTED_SHA256", hashlib.sha256(source.read_bytes()).hexdigest())
    with pytest.raises(ValueError, match="Uncorrected environment"):
        gate.verify_corrected_native_input(source)
    source.write_bytes(b"unrecognized-native-test")
    with pytest.raises(ValueError, match="Unknown native input bytes"):
        gate.verify_corrected_native_input(source)
    corrected = hashlib.sha256(source.read_bytes()).hexdigest()
    monkeypatch.setattr(gate, "EXPECTED_SHA256", corrected)
    assert gate.verify_corrected_native_input(source) == corrected


def test_native_space_gate_source_identity_is_corrected_20261020():
    from reproducibility import native_space_gate as gate
    assert gate.EXPECTED_N == 313350
    assert gate.RAW_UNCORRECTED_SHA256 == "8ceed49e4ac73b31f9a5c726ca84f4198bfc3c51e9ee43303e29524e915a5c29"
    assert gate.EXPECTED_SHA256 == "342d78e9a8a960e8d9bdcf0bc2cdacd13bee48214a8b962419d5a71d5aef46ed"
    evidence=json.loads((ROOT/"reproducibility/native_first_27_endpoint_coverage_20261020.json").read_text())
    d=evidence["denominators"]
    assert d["native_climate_all_nine_complete"] == 313240  # legacy numeric-only!
    assert d["native_climate_all_nine_complete_corrected"] == 312001
    assert evidence["environment_correction"]["corrected_eligible_nine_environment"] == 312001
    assert evidence["environment_correction"]["invalid_gsp_uint32_max"] == 1239


def test_gate_a_fréchet_bounds_are_source_backed_and_not_exact_overlaps():
    """Real source aggregates give strict inclusion-exclusion limits, not fabricated joins."""
    from reproducibility.native_gate_a import guaranteed_bounds, load_ledger
    path = ROOT/'reproducibility/native_first_27_endpoint_coverage_20261020.json'
    ledger = load_ledger(path)
    result = guaranteed_bounds(ledger)
    assert result['shared_inherited_measurement_universe'] == 217104
    for name, limits in {
        'orientation_x_chroma': (139401, 170066),
        'orientation_x_outline_aspect': (129010, 170066),
        'chroma_x_outline_aspect': (145383, 176048),
        'matched_core3': (98345, 170066),
    }.items():
        assert (result['intersections'][name]['lower'], result['intersections'][name]['upper']) == limits
        assert result['intersections'][name]['status'] == 'aggregate_frechet_bound_not_exact_join'
    committed = json.loads((ROOT/'reproducibility/native_gate_a_bounds_20261020.json').read_text())
    assert committed == result
    assert result['five_historical_fields']['native_obs_recovery_status'] == 'NOT_YET_VERIFIED'


def test_gate_a_exact_intersection_is_different_from_marginals_and_private_keys():
    import pandas as pd
    import numpy as np
    from reproducibility.native_gate_a import (
        audit_gate_a, CORE3, CORE9, RECOVERED_FIVE,
    )
    fields = list(dict.fromkeys((*CORE9, *RECOVERED_FIVE)))
    rows = []
    for i in range(40):
        r = {
            'obs_id': str(i+100), 'accepted_key': str(1+i//10),
            'latitude': -10. + i/2, 'longitude': 12. + i/3,
            'observed_on': '2026-05-20' if i%3 else None,
        }
        for f in fields:
            r[f] = float(i/5)
            r[f+'__eligible'] = True
        for f in RECOVERED_FIVE:
            r[f] = np.nan
            r[f+'__eligible'] = False
        if i<5:
            r[CORE3[0]+'__eligible'] = False
        if 5 <= i < 10:
            r[CORE3[1]+'__eligible'] = False
        if 10 <= i < 15:
            r[CORE3[2]+'__eligible'] = False
        rows.append(r)
    data = pd.DataFrame(rows)
    agg = audit_gate_a(data)
    assert agg['modules']['primary_core3']['n_observations_exact'] == 25
    assert agg['modules']['complete_primary_core9']['n_observations_exact'] == 25
    assert agg['intersections']['all_core3'] == 25
    assert agg['intersections']['__'.join(CORE3[:2])] == 30
    assert agg['modules']['colour']['n_observations_exact'] == 35
    assert agg['observer_audit'] == 'NOT_EVALUABLE_NO_ORIGINAL_OBSERVER_JOIN'
    assert not any(k in json.dumps(agg) for k in ('"obs_id"', '"user_id"', '"accepted_key"'))
    shuffled = audit_gate_a(data.sample(frac=1, random_state=9).reset_index(drop=True))
    assert agg == shuffled
    with pytest.raises(ValueError, match="Duplicated|duplicated|observation ID"):
        audit_gate_a(pd.concat([data.iloc[:1],data.iloc[:1]], ignore_index=True))
    bad = data.copy()
    bad[CORE3[0]+'__eligible'] = bad[CORE3[0]+'__eligible'].astype(object)
    bad.loc[0, CORE3[0]+'__eligible'] = 'unknown'
    with pytest.raises(ValueError, match="eligibility"):
        audit_gate_a(bad)
    bad2 = data.copy()
    bad2.loc[10, CORE3[1]] = np.nan
    with pytest.raises(ValueError, match="nonfinite"):
        audit_gate_a(bad2)
    observers = pd.DataFrame({'obs_id':data.obs_id,'user_id':['user'+str(i%4) for i in range(40)]})
    with_observers = audit_gate_a(data, observers=observers)
    assert with_observers['modules']['primary_core3']['n_distinct_original_observers'] == 4
    with pytest.raises(ValueError, match="Duplicate"):
        audit_gate_a(data, observers=pd.concat([observers, observers.iloc[:1]]))


def test_gate_a_recovered_five_requires_authenticated_head_aggregation(tmp_path):
    import pandas as pd
    import numpy as np
    from reproducibility.native_gate_a import audit_recovered_five, sha256, RECOVERED_FIVE, HEAD_ARCHIVE_SHA256
    native = pd.DataFrame({'obs_id':['10','20','30'], 'accepted_key':['a','a','b']})
    data = pd.DataFrame({
        'obs_id':['10','20','40'],
        'n_visible_photos':[1,0,1],
        'n_composition_photos':[1,0,1],
        'corolla_visible_fraction':[.6,np.nan,.3],
        'corolla_white_fraction':[.1,np.nan,.4],
        'corolla_redmagenta_fraction':[.2,np.nan,.3],
        'corolla_purple_fraction':[.4,np.nan,.2],
        'corolla_yellow_fraction':[.3,np.nan,.1],
    })
    csv=tmp_path/'recovered.csv'
    data.to_csv(csv,index=False)
    receipt={'status':'FIVE_PREVIOUSLY_MEASURED_FIELDS_RECOVERED_TO_VERSIONED_OBSERVATION_VIEW',
             'source_head_sha256':HEAD_ARCHIVE_SHA256, 'counts':{'source_heads':1255791,'head_records_deleted':0},
             'new_image_operations':False, 'observation_values_sha256':sha256(csv)}
    rp=tmp_path/'receipt.json'
    rp.write_text(json.dumps(receipt))
    result = audit_recovered_five(native,csv,rp)
    assert result['historical_observations_matching_native_source'] == 2
    assert result['native_observations_with_valid_visible_fraction'] == 1
    assert result['native_observations_with_closed_four_colour_fractions'] == 1
    assert result['n_native_accepted_taxa_with_recovered_colour_composition'] == 1
    data.loc[0,'corolla_yellow_fraction'] = .9
    data.to_csv(csv,index=False)
    receipt['observation_values_sha256'] = sha256(csv)
    rp.write_text(json.dumps(receipt))
    with pytest.raises(ValueError, match="not closed"):
        audit_recovered_five(native,csv,rp)



def test_existing_native_gate_a_trait_overlap_receipt_is_verified():
    """Aggregate original source audit already established the core nine overlap."""
    from reproducibility.native_gate_a import load_ledger, validate_existing_overlap_receipt
    ledger = load_ledger(ROOT/'reproducibility/native_first_27_endpoint_coverage_20261020.json')
    receipt = json.loads((ROOT/'reproducibility/native_gate_a_existing_overlap_20261020.json').read_text())
    result = validate_existing_overlap_receipt(ledger, receipt)
    assert result['status'] == 'VERIFIED_EXISTING_AGGREGATE_NOT_REPLAY'
    assert (result['complete_core9_n'], result['complete_core9_taxa']) == (157779, 242)
    assert result['core3_observation_range'] == {'lower': 157779, 'upper': 157828}
    assert result['photographer_join_required']
    assert result['no_new_trait_values_estimated']
    assert receipt['module_overlap']['orientation/colour'] == {'n': 164167, 'taxa': 243}
    assert receipt['module_overlap']['orientation/outline'] == {'n': 159967, 'taxa': 242}
    assert receipt['module_overlap']['colour/outline'] == {'n': 170227, 'taxa': 243}
    tampered = json.loads(json.dumps(receipt))
    tampered['complete_core9']['n'] += 1
    with pytest.raises(ValueError, match="complete-nine intersection"):
        validate_existing_overlap_receipt(ledger, tampered)
    tampered2 = json.loads(json.dumps(receipt))
    tampered2['module_census']['outline']['n'] += 2
    with pytest.raises(ValueError, match="module marginal"):
        validate_existing_overlap_receipt(ledger, tampered2)


def test_gate_a_historical_sampling_scope_and_head_metadata_dont_certify_core9():
    from reproducibility.native_gate_a import validate_existing_sampling_context
    import copy
    receipt = json.loads((ROOT/'reproducibility/native_gate_a_existing_overlap_20261020.json').read_text())
    r = validate_existing_sampling_context(receipt)
    assert r['status'] == 'HISTORICAL_SOURCE_EFFORT_AND_SPATIAL_SUPPORT_VERIFIED_ONLY'
    assert r['source_observers'] == 92839
    assert r['core9_measured_observations'] == 157779
    assert (r['core9_fitted_previous_observations'],r['core9_fitted_previous_taxa']) == (157654,174)
    assert r['core9_original_observer_strata'] == 'NOT_VERIFIED'
    assert r['legacy_head_parts_with_metadata'] == 20
    assert r['legacy_five_native_values'] == 'NOT_VERIFIED'
    s = receipt['source_sampling_context']['original_observation_effort']
    assert 0.45 < s['top2_taxa_record_fraction'] < 0.451
    assert 0.60 < s['same_taxon_cell_observer_repeated_record_fraction'] < 0.61
    bad = copy.deepcopy(receipt)
    bad['source_sampling_context']['original_observation_effort']['observers'] += 1
    with pytest.raises(ValueError, match='Native source effort'):
        validate_existing_sampling_context(bad)
    bad = copy.deepcopy(receipt)
    bad['legacy_five_archive_metadata']['archive_raw_bytes_reassembled_and_verified'] = True
    with pytest.raises(ValueError, match='Cannot certify'):
        validate_existing_sampling_context(bad)
    bad = copy.deepcopy(receipt)
    bad['source_sampling_context']['pending_exact_core9_sampling']['native_matched_observer_n'] = 92839
    with pytest.raises(ValueError, match='cannot close'):
        validate_existing_sampling_context(bad)



def test_historical_native_five_restored_measurements_do_not_replace_new_ledger():
    """PR92 executed real head-level aggregation in a different native universe."""
    from reproducibility.native_gate_a import (
        validate_existing_sampling_context, coarse_cell_retention_ceiling,
    )
    original = json.loads(
        (ROOT/"reproducibility/native_gate_a_existing_overlap_20261020.json").read_text()
    )
    out = validate_existing_sampling_context(original)
    assert out["legacy_five_restored_observations_per_endpoint_on_historical_27066_frame"] == 23881
    assert out["legacy_five_native_values"] == "NOT_VERIFIED"
    ceiling = coarse_cell_retention_ceiling(original)
    assert ceiling["upper_bound_matched_nine_taxon_cells"] == 1419
    assert ceiling["maximum_matched_nine_taxa_with_at_least_k_cells"] == {
        "5":242,"10":130,"20":61,"50":24
    }
    assert ceiling["status"] == "MATHEMATICAL_UPPER_BOUNDS_NOT_ACTUAL_CELL_COVERAGE"
    counterfeit = json.loads(json.dumps(original))
    counterfeit["legacy_five_prior_native_recovery"]["new_313350_native_restored_observations_per_endpoint"] = 23881
    with pytest.raises(ValueError, match="Do not treat old 27k"):
        validate_existing_sampling_context(counterfeit)
    counterfeit = json.loads(json.dumps(original))
    counterfeit["source_sampling_context"]["individual_endpoint_coarse_spatial_support"]["endpoints"]["orientation_image_vertical_angle"]["taxa"] = 241
    with pytest.raises(ValueError, match="Impossible prior spatial"):
        coarse_cell_retention_ceiling(counterfeit)


def test_gate_a_corrected_source_contract_matches_native_space_and_rejects_legacy():
    from reproducibility import native_gate_a, native_space_gate
    assert native_gate_a.PREPARED_SHA256 == native_space_gate.RAW_UNCORRECTED_SHA256
    assert native_gate_a.CORRECTED_PREPARED_SHA256 == native_space_gate.EXPECTED_SHA256
    assert native_gate_a.CORRECTED_PREPARED_SHA256 != native_gate_a.PREPARED_SHA256


def test_gate_a_exact_core9_matched_effort_is_safe_and_order_invariant():
    import numpy as np
    import pandas as pd
    from reproducibility.native_gate_a import audit_gate_a, CORE9, RECOVERED_FIVE
    rows = []
    for i in range(60):
        taxon = str(200 + i//20)
        r = {'obs_id':str(10000+i), 'accepted_key':taxon,
             'latitude': 10+(i//10)*.1, 'longitude': 30+(i//5)*.1,
             'observed_on': '2026-05-14' if i < 30 else '2026-08-20',
             'source_photo_count': 1}
        for name in CORE9:
            r[name] = .5+i*.01
            r[name+'__eligible'] = True
        if i % 4 == 0:
            r[CORE9[0]+'__eligible'] = False
        for name in RECOVERED_FIVE:
            r[name] = np.nan
            r[name+'__eligible'] = False
        rows.append(r)
    frame = pd.DataFrame(rows)
    obs = pd.DataFrame({
        'obs_id': frame.obs_id, 'user_id': ['private_user_'+str(i//10) for i in range(60)]
    })
    out = audit_gate_a(frame, observers=obs)
    d = out['matched_core9_sampling_diagnostics']
    assert d['source_observations'] == 60 and d['matched_observations'] == 45
    assert d['observer_status'] == 'EXACT_MATCHED_ORIGINAL_OBSERVER_SUPPORT'
    assert d['matched_observer_effort']['n_original_observers'] == 6
    assert 0 < d['matched_observer_effort']['matched_top_1pct_observer_fraction'] <= 1
    assert d['matched_observer_effort']['n_taxon_cell_observer_strata'] <= d['matched_taxon_cells']*6
    assert d['matched_observer_effort']['matched_observations_in_repeated_taxon_cell_observer_strata'] > 0
    assert len(d['source_vs_matched_monthly']) == 12
    assert sum(m['source_n'] for m in d['source_vs_matched_monthly']) == 60
    assert sum(m['matched_n'] for m in d['source_vs_matched_monthly']) == 45
    assert d['photographic_quality_status'].startswith('NOT_EVALUABLE')
    assert d['endpoint_failure_reason_status'].startswith('NOT_DISAGGREGATED')
    # Public output should contain only counts and fixed labels.
    exported = json.dumps(out, ensure_ascii=False)
    assert 'private_user_' not in exported
    assert '"obs_id"' not in exported and '"accepted_key"' not in exported
    permuted = frame.sample(frac=1, random_state=29).reset_index(drop=True)
    assert out == audit_gate_a(permuted, observers=obs.sample(frac=1, random_state=13))


def test_gate_a_no_original_observers_never_claims_sampling_gate_complete():
    import pandas as pd
    from reproducibility.native_gate_a import audit_gate_a, CORE9, RECOVERED_FIVE
    base = {'obs_id':['1','2'], 'accepted_key':['A','B'],
            'latitude':[0.,20.], 'longitude':[0.,40.]}
    for name in CORE9:
        base[name] = [.5,.6]
        base[name+'__eligible'] = [True,True]
    for name in RECOVERED_FIVE:
        base[name] = [float('nan'),float('nan')]
        base[name+'__eligible'] = [False,False]
    out = audit_gate_a(pd.DataFrame(base))
    d = out['matched_core9_sampling_diagnostics']
    assert d['matched_observer_effort'] is None
    assert d['observer_status'] == 'NOT_EVALUABLE_NO_ORIGINAL_OBSERVER_JOIN'
    assert d['season_status'] == 'NOT_EVALUABLE_NO_OBSERVATION_DATE'
    assert out['observer_audit'] == 'NOT_EVALUABLE_NO_ORIGINAL_OBSERVER_JOIN'


def test_native_gate_b_fixed_aggregate_contract_and_module_census():
    from reproducibility.native_gate_b_existing_rv import (
        CONSTRUCT_GROUP, CORRECTED_INPUT_SHA256, EXPECTED_N, EXPECTED_TAXA,
        EXPECTED_PAIR_NAMES, SOURCE_RUN, ARTIFACT_ID, MEMBER_SHA256,
    )
    from reproducibility.native_space_gate import EXPECTED_SHA256
    assert (EXPECTED_N, EXPECTED_TAXA, len(CONSTRUCT_GROUP),
            len(EXPECTED_PAIR_NAMES)) == (157654, 174, 8, 28)
    assert (SOURCE_RUN, ARTIFACT_ID) == (39428065666, 12266812186)
    assert CORRECTED_INPUT_SHA256 == EXPECTED_SHA256
    assert MEMBER_SHA256 == (
        '1545037041b62b2494934d022c7b21e7a7109f54f3582579c6557dfadf998752'
    )
    assert sum(group == 'colour' for group in CONSTRUCT_GROUP.values()) == 3
    assert sum(group == 'outline' for group in CONSTRUCT_GROUP.values()) == 4


def test_native_gate_b_grouping_distinguishes_cross_module_from_within_module():
    import itertools
    from reproducibility.native_gate_b_existing_rv import (
        CONSTRUCT_GROUP, EXPECTED_N, EXPECTED_TAXA, summarize_existing_aggregate,
    )
    constructions = list(CONSTRUCT_GROUP)
    items = []
    for i, (a, b) in enumerate(itertools.combinations(constructions, 2)):
        # Distinguish stronger among for precisely three cross-domain pairs.
        cross = CONSTRUCT_GROUP[a] != CONSTRUCT_GROUP[b]
        positive = cross and (a, b) in (
            ('angle', 'chroma'), ('angle', 'hue'), ('hue', 'width_profile')
        )
        within, between = .05, .10 if positive else .035
        difference = between - within
        record = {'rv': within, 'rv_interval': {'lower': .025, 'upper': .08}}
        alt = {'rv': between, 'rv_interval': {'lower': .015, 'upper': .13}}
        ci = {'estimate': difference,
              'interval': {'lower': .02, 'upper': .09} if positive else
                          {'lower': -.05, 'upper': .01}}
        items.append({
            'left': a, 'right': b,
            'estimates': {'raw/within':record, 'raw/between':alt,
                          'raw/between_minus_within_rv':ci,
                          'adjusted/within':record, 'adjusted/between':alt,
                          'adjusted/between_minus_within_rv':ci},
        })
    example = {'common_cohort': {
        'status':'estimated','n':EXPECTED_N,'species':EXPECTED_TAXA,
        'pairs':items,'equal_n_diagnostic':{'raw':{'scope':'synthetic'}}
    }}
    out = summarize_existing_aggregate(example)
    assert out['all_pairs']['pairs'] == 28
    assert out['cross_module_only']['pairs'] == 19
    assert out['cross_module_only']['adjusted_pointwise_interval_entirely_positive'] == 3
    assert out['all_pairs']['adjusted_pointwise_interval_entirely_positive'] == 3
    assert out['by_module_pair']['colour/outline']['pairs'] == 12
    assert out['by_module_pair']['colour/colour']['pairs'] == 3
    assert out['by_module_pair']['outline/outline']['pairs'] == 6
    assert out['by_module_pair']['outline/presentation']['pairs'] == 4
    assert out['by_module_pair']['colour/presentation']['pairs'] == 3
    assert out['pointwise_bootstrap_95_intervals_not_multiplicity_adjusted']
    assert all('genetic' not in str(i).lower() or 'not genetic' in str(i).lower()
               for i in out['limitations'])
    # A changed coefficient must fail, never become a spurious quantitative result.
    import copy
    altered = copy.deepcopy(example)
    altered['common_cohort']['pairs'][0]['estimates']['adjusted/between_minus_within_rv']['estimate'] += .1
    with pytest.raises(ValueError, match='identity broken'):
        summarize_existing_aggregate(altered)
    duplicate = copy.deepcopy(example)
    duplicate['common_cohort']['pairs'][0]['left'] = duplicate['common_cohort']['pairs'][1]['left']
    duplicate['common_cohort']['pairs'][0]['right'] = duplicate['common_cohort']['pairs'][1]['right']
    with pytest.raises(ValueError, match='Missing, duplicated'):
        summarize_existing_aggregate(duplicate)


def test_native_gate_b_fails_closed_for_changed_archive_member_bytes(tmp_path):
    from reproducibility.native_gate_b_existing_rv import verify_aggregate_file
    f = tmp_path / 'fake.json'
    f.write_text('{"run_id": 39428065666}')
    with pytest.raises(ValueError, match='Wrong native integration aggregate bytes'):
        verify_aggregate_file(f)


def test_native_gate_b_source_derived_receipt_separates_internal_from_cross_module():
    receipt=json.loads((
        ROOT/'reproducibility/native_gate_b_existing_rv_receipt_20261021.json'
    ).read_text(encoding='utf-8'))
    from reproducibility.native_gate_b_existing_rv import MEMBER_SHA256, CORRECTED_INPUT_SHA256
    assert receipt['source']['member_sha256'] == MEMBER_SHA256
    assert receipt['source']['corrected_prepared_table_sha256'] == CORRECTED_INPUT_SHA256
    assert receipt['cohort']['native_measurement_9_outputs_n'] == 157779
    assert receipt['cohort']['model_eligible_n'] == 157654
    assert receipt['cohort']['native_measurement_9_outputs_taxa'] == 242
    assert receipt['cohort']['model_eligible_taxa'] == 174
    groups=receipt['analysis']['groups']
    assert sum(g['pairs'] for g in groups.values()) == 28
    assert sum(g['adjusted_positive_deltas'] for g in groups.values()) == 24
    assert sum(g['adjusted_pointwise_positive'] for g in groups.values()) == 7
    cross=[g for k,g in groups.items() if k.split('/')[0] != k.split('/')[1]]
    assert sum(g['pairs'] for g in cross) == 19
    assert sum(g['adjusted_positive_deltas'] for g in cross) == 16
    assert sum(g['adjusted_pointwise_positive'] for g in cross) == 3
    assert receipt['analysis']['pointwise_not_multiplicity_adjusted']
    assert len(receipt['adjusted_pointwise_positive_pairs']) == 7
    assert len({p['pair'] for p in receipt['adjusted_pointwise_positive_pairs']}) == 7
    assert all(p['low'] < p['delta'] < p['high']
               for p in receipt['adjusted_pointwise_positive_pairs'])
    assert receipt['empirical_new_result'] is False
    assert receipt['no_decrypted_native_rows_in_receipt'] is True
