from pathlib import Path
import json
import hashlib
import importlib
import pytest

from reproducibility.run_current_analysis import commands, verified, native_bytes
from reproducibility.validate_current_analysis import REFERENCE, compare, validate

ROOT=Path(__file__).resolve().parents[1]


def test_current_entrypoints_import_and_counts_stay_fixed():
    cmds=commands(Path('inputs'),Path('results'))
    assert len(cmds)==7
    for cmd in cmds:
        assert importlib.import_module(cmd[0])
    assert cmds[0][-4:]==['--permutations','9999','--seed','20260910']
    assert '--native-only' not in str(cmds)


def test_reference_files_match_original_artifacts():
    manifest=json.loads((REFERENCE/'manifest.json').read_text())
    assert len(manifest['files'])==15
    for row in manifest['files']:
        assert hashlib.sha256((REFERENCE/row['path']).read_bytes()).hexdigest()==row['sha256']
    assert validate(REFERENCE)['aggregate_files_compared']==15


def test_verifier_fails_closed_and_detects_changed_results(tmp_path):
    with pytest.raises(ValueError): verified(b'wrong','0'*64,'input')
    with pytest.raises(ValueError): native_bytes(b'wrong')
    with pytest.raises(FileNotFoundError): validate(tmp_path)
    with pytest.raises(AssertionError): compare({'beta':-.345},{'beta':.345})
    with pytest.raises(AssertionError): compare({'pass':False},{'pass':True})
    with pytest.raises(AssertionError): compare([1,2],[1])


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
