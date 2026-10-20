"""Fail-closed tests for full-source WCVP eligibility and environment cache pins."""
import pandas as pd
import pytest
from analysis.v4.native_665115_pipeline import (
    add_explicit_non_species_placeholders, CLIMATE, ORIGINAL_SHA, AUTHORITY_SHA, measurement_frame,
)

def test_all_species_names_queried_and_non_species_explicitly_unresolved():
    frame=pd.DataFrame({'taxon_name':['Cirsium vulgare','Cirsium','Cirsium vulgare vulgare','Cirsium × affine'],
                        'taxon_rank':['species','genus','variety','hybrid']})
    src={'resolutions':[{'input_name':'Cirsium vulgare','resolution_status':'resolved_unique_accepted_key',
                         'accepted_key':123,'accepted_name':'Cirsium vulgare'}]}
    updated,number,missing=add_explicit_non_species_placeholders(frame,src)
    assert (number,missing)==(1,3)
    assert len(updated['resolutions'])==4
    assert all(r['accepted_key'] is None and r['resolution_status']=='unresolved_taxon'
               for r in updated['resolutions'][1:])
    assert len(src['resolutions'])==1

def test_missing_species_name_fails_before_native_classification():
    frame=pd.DataFrame({'taxon_name':['Cirsium vulgare','Cirsium discolor'],
                        'taxon_rank':['species','species']})
    src={'resolutions':[{'input_name':'Cirsium vulgare','accepted_key':1}]}
    with pytest.raises(ValueError,match='Missing exact species'):
        add_explicit_non_species_placeholders(frame,src)

def test_nine_climate_archives_and_two_source_hashes_pinned():
    assert len(CLIMATE)==9 and len({v[0] for v in CLIMATE.values()})==9
    assert all(len(v[1])==64 for v in CLIMATE.values())
    assert len(ORIGINAL_SHA)==len(AUTHORITY_SHA)==64


def join_fixture():
    roster = pd.DataFrame({'obs_id':['1','2','3','4'], 'taxon_name':['A']*4,
                           'accepted_key':[7]*4})
    env = roster.copy()
    for name, _ in CLIMATE.values():
        env[name] = [1., 2., 3., 4.]
    measured = pd.DataFrame({'obs_id':['1','2'], 'taxon_name':['A','A'],
                            'x_median':[1.,9.], 'x_n_usable_heads':[1,0]})
    archived = pd.DataFrame({'obs_id':['3'], 'taxon_name':['A'], 'endpoint_id':['bract'],
                             'value':[5.], 'analysis_eligible':[True]})
    contract = pd.DataFrame({'endpoint_id':['colour','bract','unmeasured'],
                            'observation_variable':['x_median','b_median','u_median'],
                            'measurement_variable':['x','b','u']})
    return roster, env, measured, archived, contract


def test_join_preserves_native_population_partial_traits_and_raw_ineligible_values():
    inputs = join_fixture()
    frame, report = measurement_frame(*inputs)
    assert list(frame.obs_id) == ['1','2','3','4']
    assert frame.colour.iloc[1] == 9. and not frame.colour__eligible.iloc[1]
    assert frame.bract__eligible.tolist() == [False,False,True,False]
    assert frame.unmeasured.isna().all()
    assert frame.chelsa_bio01__between.tolist() == [2.5]*4
    assert report['without_archived_measurement_row'] == 2
    assert not report['all_endpoint_complete_case_required']
    assert report['endpoints']['bract']['eligible_with_four_climate'] == 1


@pytest.mark.parametrize('problem', ['duplicate','taxon','missing_environment','invalid_qc'])
def test_measurement_join_rejects_ambiguous_identity_and_eligibility(problem):
    roster, env, measured, archived, contract = join_fixture()
    if problem == 'duplicate':
        measured = pd.concat([measured,measured.iloc[:1]])
    elif problem == 'taxon':
        archived.loc[0,'taxon_name'] = 'Other'
    elif problem == 'missing_environment':
        env = env.iloc[:3]
    else:
        archived['analysis_eligible'] = 'unknown'
    with pytest.raises(ValueError):
        measurement_frame(roster,env,measured,archived,contract)


@pytest.mark.parametrize('case', ['missing', 'wrong_key', 'tampered', 'changed_bytes'])
def test_private_source_authentication_stops_before_export(tmp_path, monkeypatch, case):
    import os
    import textwrap
    from pathlib import Path
    from cryptography.hazmat.primitives.ciphers.aead import AESGCM
    workflow = (Path(__file__).resolve().parents[1] /
                '.github/workflows/v4-native-665115-replay.yml').read_text()
    # Exercise the actual runner's acquisition code, with no external network.
    script = textwrap.dedent(workflow.split("python - <<'PY'\n", 1)[1].split("\n          PY", 1)[0])
    key = AESGCM.generate_key(bit_length=256)
    nonce = os.urandom(12)
    payload = nonce + AESGCM(key).encrypt(nonce, b'not the frozen source',
                                         b'azami-original-metadata-v1')
    if case == 'tampered':
        payload = payload[:-1] + bytes([payload[-1] ^ 1])
    supplied_key = os.urandom(32) if case == 'wrong_key' else key
    monkeypatch.setenv('AZAMI_METADATA_KEY', '' if case == 'missing' else supplied_key.hex())
    monkeypatch.chdir(tmp_path)
    source = Path('work/native-665115/private/original-metadata-665115-dfa98fa162a8029c.aesgcm')
    source.parent.mkdir(parents=True)
    source.write_bytes(payload)
    with pytest.raises(SystemExit) as exc:
        exec(compile(script, '<private-source-acquisition>', 'exec'), {})
    assert supplied_key.hex() not in str(exc.value)
    if case == 'changed_bytes':
        assert 'checksum mismatch' in str(exc.value)
    else:
        assert 'authentication failed' in str(exc.value)
    assert not list(tmp_path.rglob('*.zip'))
