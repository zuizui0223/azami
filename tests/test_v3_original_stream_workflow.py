"""Keep cloud streaming stopped until private numerical persistence exists."""
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
WORKFLOW = ROOT / ".github" / "workflows" / "ch1-v3-original-stream-pilot.yml"


def test_private_archive_stop_precedes_checkout_and_all_retrieval():
    source = WORKFLOW.read_text(encoding="utf-8")
    stop = source.index("- name: Require verified durable private numerical archive")
    checkout = source.index("- uses: actions/checkout@v4")
    retrieve = source.index("- name: Reconstruct exact private source cohort")
    assert stop < checkout < retrieve
    preflight = source[stop:checkout]
    assert "id: private_archive" in preflight
    assert "::error::STOP:" in preflight
    assert "exit 1" in preflight


def test_archive_stop_cannot_emit_a_measurement_receipt_or_delete_outputs():
    source = WORKFLOW.read_text(encoding="utf-8")
    assert "if: always()" not in source
    assert "rm -rf" not in source
    assert "workflow_incomplete.json" not in source
    guard = "if: ${{ success() && steps.private_archive.outputs.verified == 'true' }}"
    assert source.count(guard) == 2


def test_integrity_runs_hierarchical_and_stream_guard_tests():
    integrity = (ROOT / ".github" / "workflows" / "reproducibility-integrity.yml").read_text(encoding="utf-8")
    assert "tests/test_v3_*.py" in integrity
    assert (ROOT / "tests/test_v3_hierarchical_ecology.py").is_file()
    assert (ROOT / "tests/test_v3_original_stream_workflow.py").is_file()
