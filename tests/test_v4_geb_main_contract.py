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

    def test_frozen_geb_scope_excludes_s2_niche(self):
        article=(PAPER/"GEB_RESEARCH_ARTICLE.md").read_text(encoding="utf-8")
        blind=(PAPER/"GEB_RESEARCH_ARTICLE_BLINDED.md").read_text(encoding="utf-8")
        self.assertEqual(article,blind)
        self.assertIn("not a comparison of climatic niche occupancy",article)
        self.assertIn("no native-versus-introduced climatic-niche analysis",article)
        self.assertIn("photo",article.lower())
        self.assertIn("24/28",article)
        self.assertIn("0.0487/0.00532",article)
        self.assertIn("**Figure 4.** Four original S2",article)
        forbidden=("+48.2","270.4 versus 222.2","−0.172","0.448",
                   "1.791","23/91","−1.03°","40166532404")
        for item in forbidden:
            self.assertNotIn(item,article)
        for path in ("analysis/v4/s2_four_axis_joint_niche.py",
                     "analysis/v4/s2_four_axis_source_recheck.py",
                     "analysis/v4/S2_NATURAL_CONTRAST_BIO18_FOUR_AXIS_AUDIT_20261025.md",
                     "analysis/v4/results/20261025_S2_FOUR_AXIS_SOURCE_VERIFIED_SUMMARY.json"):
            self.assertFalse((ROOT/path).exists(),path)
            self.assertTrue((ROOT/"legacy/v4_s2_precipitation_niche_20261025"/path).exists(),path)
        self.assertFalse((PAPER/"s2_bio18_range_occupancy_summary.csv").exists())
        self.assertTrue((ROOT/"legacy/v4_s2_precipitation_niche_20261025/manuscripts/v4_s1_s2_frozen_20261025/s2_bio18_range_occupancy_summary.csv").exists())
        policy=(PAPER/"SCOPE_LOCK_20261025_JA.md").read_text(encoding="utf-8")
        self.assertIn("2026-10-24",policy)
        self.assertIn("降水ニッチ",policy)
        plotting=(PAPER/"plot_frozen_summary.py").read_text(encoding="utf-8")
        self.assertIn("s2_phenotype_domains",plotting)
        self.assertNotIn("s2_wet_quarter_occupancy",plotting)

    def test_s2_colour_native_counterparts_are_inference_bounded(self):
        article=(PAPER/"GEB_RESEARCH_ARTICLE.md").read_text(encoding="utf-8")
        blind=(PAPER/"GEB_RESEARCH_ARTICLE_BLINDED.md").read_text(encoding="utf-8")
        self.assertEqual(article,blind)
        self.assertIn("native coefficient (**+0.136883**)",article)
        self.assertIn("Holm P of 0.102398",article)
        self.assertIn("introduced coefficient (**+0.134504**)",article)
        self.assertIn("a direct interaction test",article)
        self.assertIn("without** comparing native–introduced rainfall niches",article)
        si=(PAPER/"SI_S2_COLOUR_NATIVE_INTRO_FIXED_HOLM_CONTEXT.md").read_text(encoding="utf-8")
        self.assertIn("22,936",si)
        self.assertIn("17,335",si)
        self.assertIn("0.017760",si)
        self.assertIn("128 biological tests",si)
        self.assertIn("not a new model",si)
        self.assertTrue((PAPER/"verify_s2_frozen_colour_context.py").is_file())

    def test_s2_native_introduced_evolutionary_theory_is_bounded(self):
        manuscript=(PAPER/"GEB_RESEARCH_ARTICLE.md").read_text(encoding="utf-8")
        blinded=(PAPER/"GEB_RESEARCH_ARTICLE_BLINDED.md").read_text(encoding="utf-8")
        self.assertEqual(manuscript,blinded)
        for name in ("Hulme and Barrett (2013)",
                     "Colautti, Maron & Barrett, 2009",
                     "Guggisberg et al. (2013)",
                     "rather than assuming that introduced plants have undergone a recent *species radiation*",
                     "did not fit a range-status interaction",
                     "does **not** compare the climatic niches"):
            self.assertIn(name,manuscript)
        note=(PAPER/"S2_NATIVE_INTRO_EVOLUTIONARY_RATIONALE_REFERENCES_JA.md").read_text(encoding="utf-8")
        self.assertIn("急速な種分化・適応放散",note)
        self.assertIn("1600年代",note)
        self.assertIn("10.1111/nph.12258",note)
        self.assertIn("降水ニッチ",note)

    def test_s2_focal_species_full_350_source_eligibility(self):
        import json
        source=(V4/"results/20261021_unified_abcde/S2_ALL_350_SPECIES_SOURCE_SCREEN_VERIFIED_20261021.json")
        d=json.loads(source.read_text(encoding="utf-8"))
        self.assertEqual(d["evidence_census"]["source_species_rank_names"],350)
        self.assertEqual(d["evidence_census"]["wild_native_and_explicit_introduced_species_names"],7)
        self.assertEqual(d["evidence_census"]["both_range_stronger_gate_species"],2)
        self.assertEqual(d["evidence_census"]["WCVP_explicit_introduced_wild_observations"] if "WCVP_explicit_introduced_wild_observations" in d["evidence_census"] else d["source"]["WCVP_explicit_introduced_wild_observations"],156750)
        strong=d["predeclared_gates"]["stronger"]
        self.assertEqual((strong["orientation_eligible_per_range"],
                          strong["5_degree_equal_sin_cells_per_range"],
                          strong["original_observers_per_range"]),(100,15,10))
        seven=d["all_both_range_species"]
        self.assertEqual(len(seven),7)
        expected={"Cirsium arvense","Cirsium vulgare"}
        for case in seven:
            computed=all(case[f"{range}_orientation_eligible"]>=100 and
                         case[f"{range}_cells"]>=15 and
                         case[f"{range}_observers"]>=10
                         for range in ("native","introduced"))
            self.assertEqual(computed,case["stronger_pass"])
        self.assertEqual({s["source_name"] for s in seven if s["stronger_pass"]},expected)
        self.assertTrue(d["selection_and_interpretation"]["two_focal_species_selected_posthoc_before_global_source_screen"])
        self.assertTrue(d["source"]["all_source_taxa_checked_without_rain_or_orientation_values"])
        doc=(PAPER/"SI_S2_FOCAL_SPECIES_SELECTION_AUDIT.md").read_text(encoding="utf-8")
        paper=(PAPER/"GEB_RESEARCH_ARTICLE.md").read_text(encoding="utf-8")
        for species in seven:
            self.assertIn(species["source_name"].split()[1],doc)
        for token in ("156,608","156,750","99.91%","350","seven","100","15","10"):
            self.assertIn(token,paper)
        self.assertIn("Table S3",paper)
        self.assertIn("not random",paper)
        self.assertNotIn("+48.2 mm/month",paper)

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
