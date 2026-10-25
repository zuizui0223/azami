"""Read-only integrity checks for the frozen V4 GEB manuscript and V3 archive."""
from __future__ import annotations
import csv
from pathlib import Path
import re
import unittest

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "manuscripts" / "v4_s1_s2_frozen_20261025"
V4 = ROOT / "analysis" / "v4"

class V4GEBMainContract(unittest.TestCase):
    def test_geb_structure(self):
        src = (PAPER / "GEB_RESEARCH_ARTICLE.md").read_text(encoding="utf-8")
        abstract = src.split("## Abstract", 1)[1].split("## 1. Introduction", 1)[0]
        for heading in ("Aim:", "Location:", "Time period:", "Major taxa studied:", "Methods:", "Results:", "Main conclusions:"):
            self.assertIn("**"+heading+"**", abstract)
        words = re.findall(r"\b[\w]+(?:[–-][\w]+)*\b", abstract.split("**Keywords:**")[0])
        self.assertLessEqual(len(words), 300, "GEB abstract >300 words")
        kw = abstract.split("**Keywords:**", 1)[1].split("\n",1)[0].strip().split(",")
        self.assertEqual(len(kw), 8)
        self.assertEqual([a.strip() for a in kw], sorted((a.strip() for a in kw), key=str.lower))
        blind = (PAPER / "GEB_RESEARCH_ARTICLE_BLINDED.md").read_text(encoding="utf-8")
        self.assertIn("## Data and code availability", blind)
        self.assertNotIn("zuizui0223", blind)

    def test_frozen_holm_four(self):
        with (PAPER / "table1_holm_associations.csv").open(newline="",encoding="utf-8") as f:
            rows=list(csv.DictReader(f))
        self.assertEqual(len(rows), 4)
        expected=[("S1",168907,0.045447,0.00920635),
                  ("S2",26171,-0.173732,0.00131233),
                  ("S2",27200,0.134504,0.01456813),
                  ("S2",27200,-0.142801,0.01925414)]
        for row,(study,n,b,p) in zip(rows,expected):
            self.assertEqual(row["stratum"],study)
            self.assertEqual(int(row["n"]),n)
            self.assertAlmostEqual(float(row["beta"]),b,places=6)
            self.assertAlmostEqual(float(row["holm_p"]),p,places=7)

    def test_v4_rv_and_source(self):
        with (PAPER / "module_summary.csv").open(newline="",encoding="utf-8") as f:
            pairs=list(csv.DictReader(f))
        self.assertEqual(len(pairs), 4)
        self.assertEqual({(r["n_photographs"],r["n_taxa"]) for r in pairs},{("156683","174")})
        self.assertAlmostEqual(float(pairs[0]["median_RV"]),0.085245,places=6)
        self.assertAlmostEqual(float(pairs[-1]["median_RV"]),0.015453,places=6)
        freeze = (V4 / "FINAL_V4_BIO18_HOLM_LEGO_FREEZE_20261024_JA.md").read_text(encoding="utf-8")
        for value in ("168,907","156,683","BIO18","Holm","24/28","12547876704"):
            self.assertIn(value, freeze)
        self.assertTrue((V4 / "final_v4_bio18_holm_freeze.py").is_file())

    def test_v3_archived(self):
        self.assertFalse((ROOT / "analysis/v3").exists())
        self.assertFalse((ROOT / ".github/workflows/publish-ch1-zenodo-v3-data.yml").exists())
        for path in ("legacy/v3/analysis/v3/README.md",
                     "legacy/v3/reproducibility/CURRENT_ANALYSIS.md",
                     "legacy/v3/workflows/active_as_of_20260928/publish-ch1-zenodo-v3-data.yml",
                     "legacy/v3/workflows/previously_archived/ch1-v3-scale-integration.yml"):
            self.assertTrue((ROOT/path).is_file(), path)
        archive=(ROOT/"legacy/v3/README.md").read_text(encoding="utf-8")
        self.assertIn("b463dbc938e454baa54b2e480f7f1d8badd7795d",archive)

if __name__=="__main__":
    unittest.main()
