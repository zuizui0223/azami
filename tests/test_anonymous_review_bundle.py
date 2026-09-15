from pathlib import Path
import json

from reproducibility import build_anonymous_review_bundle as anon
from reproducibility.run_current_analysis import commands

ROOT = Path(__file__).resolve().parents[1]


def test_anonymous_runner_tracks_current_eight_stage_surface(tmp_path):
    current_commands = commands(tmp_path / 'inputs', tmp_path / 'results')
    current = [row[0] for row in current_commands]
    for module in current:
        assert module in anon.RUNNER
    assert len(current) == 8
    assert 'run_rv_estimator_validity' in anon.RUNNER
    assert '--equal-n-replicates' in anon.RUNNER
    # The first seven current stages use the 20260910 seed contract, while the
    # post-hoc estimator-validity stage is frozen separately at 20260915.
    assert current_commands[-1][-2:] == ['--seed', '20260915']
    assert '"--seed", "20260915"' in anon.RUNNER


def test_anonymous_bundle_source_surface_excludes_direct_identity_tokens():
    paths = list(anon.STATIC_SOURCE_PATHS) + sorted(
        p.relative_to(ROOT) for p in (ROOT / 'analysis/v3').glob('*.py')
    )
    for rel in paths:
        text = (ROOT / rel).read_text(encoding='utf-8', errors='ignore').lower()
        for token in anon.BANNED_TEXT:
            assert token.lower() not in text, (rel, token)


def test_anonymous_reference_manifest_can_be_sanitized(tmp_path):
    package = tmp_path / 'package'
    package.mkdir()
    rows = anon.sanitized_reference(package)
    assert len(rows) == 16
    manifest = json.loads((package / 'current_reference/manifest.json').read_text())
    assert manifest['scientific_reference_commit'] == 'anonymous-review'
    assert all(set(row) == {'path', 'sha256'} for row in manifest['files'])
    text = (package / 'current_reference/manifest.json').read_text().lower()
    assert 'github' not in text
    assert 'artifact' not in text
    assert 'zenodo' not in text


def test_generated_review_text_is_identity_clean(tmp_path):
    package = tmp_path / 'package'
    package.mkdir()
    anon.copy_source(package)
    assert anon.identity_hits(package) == []
    assert 'public repository url' in (package / 'README.md').read_text().lower()
