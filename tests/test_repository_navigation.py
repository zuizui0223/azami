"""Repository navigation and canonical manuscript checks (no network access)."""
from __future__ import annotations

import re
import unittest
from pathlib import Path
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[1]
INDEXES = (
    "README.md",
    "analysis/README.md",
    "manuscripts/README.md",
    "legacy/README.md",
    "docs/REPOSITORY_MAP.md",
    "docs/REPRODUCIBILITY_AND_SUBMISSION_RUNBOOK.md",
    "docs/BRANCH_ARCHIVE_TRIAGE_20261030.md",
)
LINK = re.compile(r"\[[^\]]*\]\(([^)]+)\)")


class RepositoryNavigationTests(unittest.TestCase):
    def test_index_files_exist(self):
        for name in INDEXES:
            with self.subTest(index=name):
                self.assertTrue((ROOT / name).is_file())

    def test_internal_index_links_resolve(self):
        for name in INDEXES:
            src = ROOT / name
            content = src.read_text(encoding="utf-8")
            for url in LINK.findall(content):
                if "://" in url or url.startswith(("#", "mailto:")):
                    continue
                target = unquote(url.split("#", 1)[0])
                if not target:
                    continue
                resolved = (src.parent / target).resolve()
                with self.subTest(index=name, target=target):
                    self.assertTrue(resolved.is_relative_to(ROOT.resolve()))
                    self.assertTrue(resolved.exists(), f"Broken link: {name} -> {target}")

    def test_current_results_not_obsolete_index(self):
        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        paper = (ROOT / "manuscripts/v4_s1_s2_frozen_20261025/GEB_RESEARCH_ARTICLE.md").read_text(encoding="utf-8")
        for identifier in ("3.1", "3.2", "3.3"):
            self.assertIn(identifier, readme)
            self.assertIn(f"### {identifier}.", paper)
        self.assertIn("H3c", readme)
        self.assertIn("168,907", readme)
        self.assertIn("7,117", readme)
        self.assertIn("128-test", readme)
        self.assertNotIn("Part 1: Phenotypic integration", readme)


if __name__ == "__main__":
    unittest.main()
