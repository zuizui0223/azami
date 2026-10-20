"""Fail-closed tests for full-source WCVP eligibility and environment cache pins."""
import pandas as pd
import pytest
from analysis.v4.native_665115_pipeline import (
    add_explicit_non_species_placeholders, CLIMATE, ORIGINAL_SHA, AUTHORITY_SHA,
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
