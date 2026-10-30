"""Offline tests for the read-only GitHub branch reference inventory."""
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from io import BytesIO
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/repository_ref_audit.py"
spec = importlib.util.spec_from_file_location("repository_ref_audit", SCRIPT)
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def branch(name, sha):
    return {"name": name, "commit": {"sha": sha}}


def pr(number, head, base, state="closed"):
    return {"number": number, "title": f"PR {number}", "state": state,
            "head": {"ref": head, "sha": f"head-{number}"},
            "base": {"ref": base, "sha": f"base-{number}"}}


class BranchAuditTests(unittest.TestCase):
    def test_groups_historical_and_open_dependencies(self):
        branches = [branch("main", "sha-main"), branch("paper", "sha-paper"),
                    branch("old-a", "sha-old"), branch("old-b", "sha-old"),
                    branch("old-c", "sha-old")]
        pulls = [pr(14, "old-a", "main"), pr(159, "paper", "main", "open")]
        report = audit.build_report(branches, pulls, [], [], "main")
        group = report["shared_tip_groups"][0]
        self.assertEqual(group["count"], 3)
        self.assertEqual(group["no_matching_pr_head_names"], ["old-b", "old-c"])
        self.assertEqual(group["historical_pr_head_match_count"], 1)
        self.assertEqual(group["open_pr_or_default_refs"], [])
        self.assertEqual(report["open_pr_dependencies"][0]["head"], "paper")
        self.assertEqual(report["open_pr_dependencies"][0]["base"], "main")
        self.assertEqual(group["assessment"], "REVIEW_ONLY_NEVER_AUTODELETE")

    def test_protects_open_pr_base_even_with_duplicate_tip(self):
        report = audit.build_report(
            [branch("main", "m"), branch("paper", "same"), branch("other", "same")],
            [pr(160, "other", "paper", "open")], [], [], "main",
        )
        self.assertEqual(report["shared_tip_groups"][0]["open_pr_or_default_refs"], ["other", "paper"])

    def test_annotated_tags_are_not_treated_as_commit_proof(self):
        tag = {"ref": "refs/tags/v1", "object": {"sha": "tag-object", "type": "tag"}}
        report = audit.build_report([], [], [tag], [{"draft": True}], "main")
        self.assertEqual(report["tag_refs"][0]["target_type"], "tag")
        self.assertEqual(report["draft_release_count"], 1)

    def test_duplicate_pagination_fails_closed(self):
        with self.assertRaisesRegex(ValueError, "Duplicate branch"):
            audit.build_report([branch("main", "a"), branch("main", "a")], [], [], [], "main")

    def test_query_pagination_preserves_state_filter(self):
        requested = []
        class Response(BytesIO):
            pass
        def fake_open(request, timeout):
            requested.append(request.full_url)
            return Response(b"[]")
        with patch.object(audit, "urlopen", side_effect=fake_open):
            self.assertEqual(audit.get_all("zuizui0223/azami", "pulls?state=all"), [])
        self.assertEqual(requested, ["https://api.github.com/repos/zuizui0223/azami/pulls?state=all&per_page=100&page=1"])

    def test_offline_cli_and_no_delete_commands(self):
        inputs = {"branches": [branch("main", "a")], "pulls": [], "tags": [], "releases": []}
        with tempfile.TemporaryDirectory() as temp:
            source = Path(temp) / "input.json"
            output = Path(temp) / "output.json"
            source.write_text(json.dumps(inputs))
            self.assertEqual(audit.main(["--input-json", str(source), "--output-json", str(output)]), 0)
            self.assertEqual(json.loads(output.read_text())["branch_count"], 1)
        source_code = SCRIPT.read_text()
        self.assertNotIn('method="DELETE"', source_code)
        self.assertNotIn('method="PATCH"', source_code)
        self.assertNotIn('method="POST"', source_code)

    def test_frozen_branch_snapshot(self):
        source = ROOT / "docs/BRANCH_REF_INPUT_20261030.json"
        inputs = json.loads(source.read_text(encoding="utf-8"))
        report = audit.build_report(
            inputs["branches"], inputs["pulls"], inputs["tags"], inputs["releases"],
            inputs["default_branch"],
        )
        self.assertEqual(report["branch_count"], 137)
        self.assertEqual(report["pull_request_count"], 157)
        self.assertEqual(report["open_pull_request_count"], 6)
        shared = [group for group in report["shared_tip_groups"]
                  if group["sha"] == "7d28f2340e1d1b7fed5bd7a5611f52853e09a388"]
        self.assertEqual(len(shared), 1)
        self.assertEqual(shared[0]["count"], 49)
        self.assertEqual(len(shared[0]["no_matching_pr_head_names"]), 7)
        self.assertEqual(shared[0]["open_pr_or_default_refs"], [])


if __name__ == "__main__":
    unittest.main()
