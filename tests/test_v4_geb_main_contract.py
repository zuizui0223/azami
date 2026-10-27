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
        self.assertIn("not an independent source cohort",article)
        self.assertIn("no native-versus-introduced climatic-niche analysis",article)
        self.assertIn("photo",article.lower())
        si_cov=(PAPER/"SI_V4_FULL_COVARIATION_AND_CONFIGURATION_CONTEXT.md").read_text(encoding="utf-8")
        self.assertIn("24/28",si_cov)
        self.assertIn("0.0487/0.00532",si_cov)
        self.assertIn("**123/135**",si_cov)
        self.assertIn("**Figure 2 (Part 1).** Source-verified descriptive signed correlations",article)
        self.assertNotIn("**Figure 5 (Part 1).**",article)
        for biological_pair in ("angle", "chroma", "aspect"):
            self.assertIn(biological_pair,(PAPER/"plot_cross_functional_pairs.py").read_text(encoding="utf-8"))
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
        self.assertIn("SI_S2_trait_domain_RV",plotting)
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
                     "do **not** compare native and introduced climatic niches"):
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

    def test_three_part_manuscript_has_no_formal_lego_hypothesis(self):
        article=(PAPER/"GEB_RESEARCH_ARTICLE.md").read_text(encoding="utf-8")
        blinded=(PAPER/"GEB_RESEARCH_ARTICLE_BLINDED.md").read_text(encoding="utf-8")
        self.assertEqual(article,blinded)
        self.assertNotIn("LEGO",article)
        headings=re.findall(r"^### 3\.(\d+)\. (.+)$",article,re.MULTILINE)
        self.assertEqual([int(i) for i,_ in headings],[1,2,3])
        for prefix in ("### 2.4. Part 1 — Visible phenotypic integration",
                       "### 2.5. Part 2 — Native-taxon capitulum angle",
                       "### 2.6. Part 3 — Native/introduced",
                       "### 3.1. Part 1 — Partial integration",
                       "### 3.2. Part 2 — Among-taxon image angle",
                       "### 3.3. Part 3 — Native/introduced colour",
                       "### 4.1. Part 1 — Integration",
                       "### 4.2. Part 2 — Climatic precipitation",
                       "### 4.3. Part 3 — Colour gradients"):
            self.assertIn(prefix,article)
        self.assertIn("not an independent source cohort",article)
        self.assertIn("does **not** create a source S3",article)
        self.assertIn("phenotypic integration",article.lower())
        self.assertIn("156,683 shared native photographs",article)
        outline=(PAPER/"GEB_THREE_PART_STRUCTURE_JA.md").read_text(encoding="utf-8")
        for val in ("## 1. 頭花形質の定量化、表現型統合と組合せ",
                    "## 2. 頭花の画像角度とBIO18",
                    "## 3. 在来域と導入域の花色―環境関係",
                    "新しいS3ソース"):
            self.assertIn(val,outline)
        self.assertNotIn("降水ニッチの比較を行った",article)

    def test_v3_to_v4_minimal_story_and_ecological_figure_evidence(self):
        paper=(PAPER/"GEB_RESEARCH_ARTICLE.md").read_text(encoding="utf-8")
        blind=(PAPER/"GEB_RESEARCH_ARTICLE_BLINDED.md").read_text(encoding="utf-8")
        self.assertEqual(paper,blind)
        for key in ("Dalrymple et al. (2020)","Berardi et al., 2016",
                    "greater image chroma and lower image lightness",
                    "**Figure 3 (Part 2: featured ecological result).**",
                    "**Figure 4 (Part 3).**",
                    "https://doi.org/10.1111/nph.16737",
                    "https://doi.org/10.3732/ajb.1600106"):
            self.assertIn(key,paper)
        for key in ("Table Part 2","V4 the three-part study"):
            self.assertNotIn(key,paper)
        plan=(PAPER/"GEB_V3_V4_MINIMAL_DELTA_STORY_AND_FIGURES_JA.md").read_text(encoding="utf-8")
        for key in ("46,276","1,734","33/36","24/28","156,683",
                    "V3 chroma","native–introduced","Fig. 3","Fig. 4",
                    "https://doi.org/10.1111/nph.16737"):
            self.assertIn(key,plan)
        self.assertTrue((PAPER/"GEB_COVER_LETTER_CONTRIBUTION_PITCH.md").is_file())
        figure=(PAPER/"plot_s2_native_introduced_colour.py").read_text(encoding="utf-8")
        for field in ("read_source","EXPECTED","holm_p","ci_lower","ci_upper","native","introduced"):
            self.assertIn(field,figure)

    def test_featured_bio18_main_and_full_uncorrected_si(self):
        paper=(PAPER/"GEB_RESEARCH_ARTICLE.md").read_text(encoding="utf-8")
        blinded=(PAPER/"GEB_RESEARCH_ARTICLE_BLINDED.md").read_text(encoding="utf-8")
        self.assertEqual(paper,blinded)
        self.assertEqual([int(v) for v in re.findall(r"^### 3\.(\d+)\.",paper,re.MULTILINE)],[1,2,3])
        for term in ("**Part 1** identifies weak among-taxon associations",
                     "**Part 2** identifies the strongest",
                     "**Part 3** identifies three",
                     "**Figure 3 (Part 2: featured ecological result).**",
                     "main figure includes only this Holm-supported term",
                     "all 80 native-taxon tests with raw, BH and Holm P values"):
            self.assertIn(term,paper)
        main_table=paper.split("**Table 2.",1)[1].split("**Out of scope",1)[0]
        self.assertNotIn("| Raw P |",main_table)
        self.assertIn("| Holm P |",main_table)
        self.assertIn("0.00920635",main_table)
        note=(PAPER/"SI_S1_BIO18_ANGLE_RAW_BH_HOLM_CONTEXT.md").read_text(encoding="utf-8")
        for value in ("0.039865108","0.395836573","1.000000000",
                      "+0.058700384","+0.045446843","80 biological tests",
                      "numerically larger","168,907","0.009206351"):
            self.assertIn(value,note)
        script=(PAPER/"plot_bio18_angle_part2.py").read_text(encoding="utf-8")
        for value in ("SOURCE_SHA256","original_frozen_rows",
                      "SI_S1_all_80_original_raw_BH_Holm.csv",
                      "SI_S1_BIO18_image_angle_both_scales.csv",
                      "figure3_BIO18_among_taxa_image_angle",
                      "SI_Figure_BIO18_image_angle_both_scales",
                      "p_original","p_BH","p_Holm"):
            self.assertIn(value,script)
        ci=(ROOT/".github/workflows/v4-geb-main-contract.yml").read_text(encoding="utf-8")
        for value in ("SI_S1_all_80_original_raw_BH_Holm.csv",
                      "figure3_BIO18_among_taxa_image_angle.png",
                      "figure2_cross_functional_trait_pairs.png",
                      "figure4_S2_native_introduced_colour_coefficients.png"):
            self.assertIn(value,ci)

    def test_photo_cross_domain_interaction_discussion_has_primary_sources(self):
        article=(PAPER/"GEB_RESEARCH_ARTICLE.md").read_text(encoding="utf-8")
        blind=(PAPER/"GEB_RESEARCH_ARTICLE_BLINDED.md").read_text(encoding="utf-8")
        self.assertEqual(article,blind)
        audit=(PAPER/"GEB_INTERACTION_FILTER_LITERATURE_AUDIT_20261027_JA.md").read_text(encoding="utf-8")
        for value in ("−0.2251","−0.2469","−0.0281","−0.0694","r=+0.0199",
                      "P(landing | approach)", "photographic CIELAB chroma",
                      "neither a fixed three-trait syndrome", "not pollinator choice",
                      "capitulum receptacle width", "Neither measure is the photograph-derived height-to-width ratio"):
            self.assertIn(value,article)
        for value in ("## 2. 実験・比較文献のエビデンス階層",
                      "### H1.", "### H2.", "### H3.", "### H4.", "### H5.",
                      "分類群間", "分類群内", "視覚", "食害"):
            self.assertIn(value,audit)
        for doi in ("10.1007/s00114-023-01845-w",
                    "10.1111/ele.13204",
                    "10.1002/ece3.73221",
                    "10.1007/s004420100773",
                    "10.1139/b99-182",
                    "10.1111/plb.13414",
                    "10.1073/pnas.071053098"):
            self.assertIn(doi,article)
            self.assertIn(doi,audit)
        self.assertIn("a larger **number of flowering heads per plant**",article)
        self.assertIn("lower aspect ratio",audit)
        self.assertNotIn("LEGO",article)

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
