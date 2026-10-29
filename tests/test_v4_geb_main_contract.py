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
        self.assertIn("**Figure 2 (planned, NOT YET SOURCE-GENERATED).**",article)
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
        self.assertIn("direct interaction test",article)
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
        for name in ("Hulme & Barrett, 2013",
                     "Colautti, Maron & Barrett, 2009",
                     "Guggisberg et al., 2013",
                     "not a randomised climate intervention",
                     "**direct interaction test**",
                     "no native-versus-introduced climatic-niche analysis"):
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
        for prefix in ("### 2.4. Supplementary matched-photograph trait associations",
                       "### 2.5. Climatic associations within and among native taxa",
                       "### 2.6. Within-species native and introduced",
                       "### 3.1. Global image measurements and precipitation-associated",
                       "### 3.2. Radiation and elevation gradients",
                                              "### 4.1. Precipitation, head presentation",
                       "### 4.2. Floral colour and climate",
                       "### 4.4. What public photographs measure"):
            self.assertIn(prefix,article)
        self.assertIn("not an independent source cohort",article)
        self.assertIn("does **not** create a source S3",article)
        self.assertIn("SI_V4_FULL_COVARIATION_AND_CONFIGURATION_CONTEXT.md",article)
        self.assertIn("156,683 photographs across 174 taxa",article)
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
        for key in ("Dalrymple et al., 2020","Berardi et al., 2016",
                    "Greater chroma and lower lightness",
                    "**Figure 3.**",
                    "**Figure 4.**",
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
        for term in ("Globally distributed public photographs",
                     "a larger downward-looking image angle",
                     "The conspecific native–introduced comparison provides a second ecological perspective.",
                     "**Figure 3.**",
                     "The native within-taxon coefficient and the complete 80-test family are in Supporting Information",
                     "80-test family"):
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

    def test_v4_manuscript_results_before_mechanistic_speculation(self):
        article=(PAPER/"GEB_RESEARCH_ARTICLE.md").read_text(encoding="utf-8")
        blind=(PAPER/"GEB_RESEARCH_ARTICLE_BLINDED.md").read_text(encoding="utf-8")
        self.assertEqual(article,blind)
        self.assertTrue(article.startswith("# Global photographic phenomics reveals"))
        self.assertIn("regional radiations",article)
        self.assertIn("Moreyra et al., 2025",article)
        self.assertIn("Siniscalchi, Ackerfield & Folk, 2023",article)
        self.assertIn("Hulme & Barrett, 2013",article)
        self.assertIn("avoidance of pollen wetting and wash-off",article)
        self.assertIn("only Holm-surviving association",article)
        self.assertIn("1981–2010 climatic normal",article)
        self.assertIn("**141,014 native, 156,608 introduced**",article)
        self.assertIn("three Holm-supported comparisons from a single 128-test family",article)
        self.assertIn("native coefficient (**+0.136883**)",article)
        self.assertIn("Holm P of 0.102398",article)
        self.assertIn("**direct interaction test**",article)
        self.assertIn("not evidence of rainfall-mediated selection",article)
        self.assertIn("no native-versus-introduced climatic-niche analysis",article)
        self.assertIn("**without** comparing native–introduced rainfall niches",article)
        for doi in ("10.1016/j.ympev.2025.108285",
                    "10.1086/724310",
                    "10.1080/17550874.2012.702793",
                    "10.1111/1365-2745.12009",
                    "10.1111/nph.16737"):
            self.assertIn(doi,article)
        si=(PAPER/"SI_V4_FULL_COVARIATION_AND_CONFIGURATION_CONTEXT.md").read_text(encoding="utf-8")
        for pair in ("−0.225128","−0.246946","+0.019940"):
            self.assertIn(pair,si)
        for pair in ("r=−0.2251","r=−0.2469","r=+0.0199"):
            self.assertNotIn(pair,article)
        audit=(PAPER/"GEB_INTERACTION_FILTER_LITERATURE_AUDIT_20261027_JA.md").read_text(encoding="utf-8")
        for value in ("## 2. 実験・比較文献のエビデンス階層",
                      "### H1.", "### H2.", "### H3.", "### H4.", "### H5.",
                      "分類群間", "分類群内", "視覚", "食害"):
            self.assertIn(value,audit)
        self.assertNotIn("LEGO",article)

    def test_geb_published_style_main_findings_and_display_priorities(self):
        article=(PAPER/"GEB_RESEARCH_ARTICLE.md").read_text(encoding="utf-8")
        blinded=(PAPER/"GEB_RESEARCH_ARTICLE_BLINDED.md").read_text(encoding="utf-8")
        self.assertEqual(article,blinded)
        for term in ("Kuppler et al., 2020","Grenié et al., 2025",
                     "https://doi.org/10.1111/geb.13077","https://doi.org/10.1111/geb.70131",
                     "**Table 1.** Source observations and image-angle eligibility",
                     "| *Cirsium arvense* | 78,646 | 70,964 | 49,072 | 43,195 |",
                     "| *Cirsium vulgare* | 62,368 | 85,644 | 33,637 | 47,833 |",
                     "**Total** | **141,014** | **156,608**",
                     "Figure 2 (planned, NOT YET SOURCE-GENERATED)",
                     "The full within- and among-taxon image-metric covariance analysis"):
            self.assertIn(term,article)
        self.assertLess(article.index("### 3.1. Global image measurements"),
                        article.index("### 3.2. Radiation and elevation gradients"))
        self.assertLess(article.index("### 4.1. Precipitation, head presentation"),
                        article.index("### 4.4. What public photographs measure"))
        benchmark=(PAPER/"GEB_PUBLISHED_RESEARCH_STYLE_AUDIT_20261029_JA.md").read_text(encoding="utf-8")
        for doi in ("10.1111/geb.13077","10.1111/geb.13507",
                    "10.1111/geb.70131","10.1111/geb.70245"):
            self.assertIn(doi,benchmark)

    def test_bio18_geographic_sensitivity_is_si_only_and_not_causal_partition(self):
        import json
        path=V4/"results/20261029_S1_BIO18_GEOGRAPHY_SOURCE_VERIFIED.json"
        data=json.loads(path.read_text(encoding="utf-8"))
        self.assertEqual((data["n"],data["taxa"]),(168907,178))
        self.assertEqual(data["source"]["original_secure_workflow_run"],40837107465)
        self.assertEqual(data["source"]["public_artifact_id"],12881655588)
        a=data["geographic_fits"]["photo_weighted"]["global_geographic_adjustments"]
        self.assertAlmostEqual(a["none"]["standardised_beta"],0.04544684307686103,places=7)
        self.assertAlmostEqual(a["longitude_sector"]["standardised_beta"],0.02003485434695006,places=9)
        self.assertAlmostEqual(a["spherical_degree2"]["standardised_beta"],-0.003294794324869813,places=9)
        paper=(PAPER/"GEB_RESEARCH_ARTICLE.md").read_text(encoding="utf-8")
        si=(PAPER/"SI_S1_BIO18_SPATIAL_SENSITIVITY_20261029.md").read_text(encoding="utf-8")
        self.assertIn("SI_S1_BIO18_SPATIAL_SENSITIVITY_20261029.md",paper)
        self.assertIn("Spatial predictors share information with climate", si)
        self.assertIn("no causal partition", si)
        self.assertNotIn("The **exact-cohort geographic sensitivity**",paper)
        self.assertNotIn("geographically structured precipitation–capitulum associations",paper)
        self.assertTrue((PAPER/"plot_bio18_geographic_sensitivity.py").is_file())

    def test_new_third_ecological_question_not_old_trait_r_claim(self):
        article=(PAPER/"GEB_RESEARCH_ARTICLE.md").read_text(encoding="utf-8")
        self.assertIn("### 3.1. Global image measurements",article)
        self.assertIn("### 3.2. Radiation and elevation gradients",article)
        self.assertNotIn("r=−0.2251",article)
        self.assertNotIn("r=−0.2469",article)
        self.assertTrue((V4/"s1_climate_trait_hypervolume.py").is_file())
        self.assertTrue((V4/"HYPERVOLUME_CLIMATE_MORPHOSPACE_H3_PREDECLARED_20261029_JA.md").is_file())

    def test_h3_climatic_breadth_and_flower_morphospace_original_source(self):
        import json
        h3=(V4/"results/20261029_H3_CLIMATE_TRAIT_HYPERVOLUME_SOURCE_VERIFIED.json")
        source=json.loads(h3.read_text(encoding="utf-8"))
        self.assertEqual(source["verified_execution"]["run_id"],40841007115)
        self.assertEqual(source["verified_execution"]["artifact_id"],12883021095)
        self.assertEqual(source["source"]["quality_eligible_matched_observations"],157822)
        self.assertEqual(source["source"]["matched_source_taxa"],242)
        a=source["analyses"]["1degree"]
        b=source["analyses"]["5degree"]
        self.assertEqual((a["taxa"],b["taxa"]),(63,31))
        self.assertAlmostEqual(a["breadth"]["trait_logvolume_slope_per_climate_logvolume"],
                               0.06689477516770817,places=12)
        self.assertAlmostEqual(a["breadth"]["sector_conditioned_residual_permutation_p_exploratory"],
                               0.3875,places=6)
        self.assertAlmostEqual(a["similarity"]["climate_distance_vs_trait_distance_spearman"],
                               -0.008144762186869107,places=12)
        self.assertTrue(a["similarity"]["phylogenetic_evolutionary_convergence_NOT_tested"])
        manuscript=(PAPER/"GEB_RESEARCH_ARTICLE.md").read_text(encoding="utf-8")
        for term in ("### 3.3. Climate-conditioned floral resemblance beyond niche breadth",
                     "### 4.3. Local climatic correspondence without a single global floral syndrome",
                     "**157,822 records across 242 taxa**","**63 taxa**",
                     "permutation P=0.3875","**ρ=−0.0081**",
                     "not a thresholded hypervolume intersection",
                     "SI_H3_CLIMATE_ANALOGUE_CONDITIONAL_OVERLAP_20261029.md"):
            self.assertIn(term,manuscript)
        self.assertIn("The standardized among-taxon coefficient was **+0.045447**",manuscript)
        self.assertIn("β=**+0.134504**",manuscript)
        self.assertIn("β=**−0.173732**",manuscript)
        self.assertIn("β=**−0.142801**",manuscript)
        self.assertTrue((PAPER/"plot_h3_climate_morphospace.py").is_file())
        self.assertTrue((PAPER/"SI_H3_CLIMATE_ANALOGUE_CONDITIONAL_OVERLAP_20261029.md").is_file())

    def test_h3c_climate_analogue_and_h3d_covariance_sources_frozen(self):
        import json
        c=json.loads((V4/"results/20261029_H3C_CLIMATE_ANALOGUE_TRAIT_OVERLAP_SOURCE_VERIFIED.json").read_text(encoding="utf-8"))
        d=json.loads((V4/"results/20261029_H3D_ENVIRONMENT_MODERATED_TRAIT_COVARIANCE_SOURCE_VERIFIED.json").read_text(encoding="utf-8"))
        for payload in (c,d):
            self.assertEqual(payload["verified_execution"]["github_actions_run"],40847696537)
            self.assertEqual(payload["verified_execution"]["aggregate_artifact_id"],12885860903)
            self.assertEqual(payload["source"]["quality_eligible_matched_observations"],157822)
        for key,n in (("1degree",63),("5degree",31)):
            r=c["results"][key]
            self.assertEqual(r["n_taxa"],n)
            self.assertEqual(r["successful_rarefaction_draws"],24)
            self.assertAlmostEqual(r["holm2_distance_vs_overlap_p"],0.01,places=6)
            self.assertEqual(d["summary"][key]["taxa"],n)
            self.assertFalse(any(t["exploratory_Holm_6_moderations"]<.05
                                 for t in d["summary"][key]["moderation"]["results"]))
            self.assertFalse(any(t["exploratory_Holm_3_traits"]<.05
                                 for t in d["summary"][key]["breadth"]["results"]))
        self.assertAlmostEqual(c["results"]["1degree"]["primary_trait_distance_advantage_of_climate_matched_cells"],
                               .040963661993538125,places=12)
        self.assertAlmostEqual(c["results"]["5degree"]["primary_trait_distance_advantage_of_climate_matched_cells"],
                               .0943310892554725,places=12)
        paper=(PAPER/"GEB_RESEARCH_ARTICLE.md").read_text(encoding="utf-8")
        for claim in ("**+0.04096**", "**+0.09433**",
                      "Holm-2 P=0.010",
                      "conditional climate",
                      "SI_H3_CLIMATE_ANALOGUE_CONDITIONAL_OVERLAP_20261029.md",
                      "SI_H3_ENVIRONMENT_MODERATED_TRAIT_COVARIANCE_20261029.md"):
            self.assertIn(claim if claim!="conditional climate" else "local climate-analogue analysis",paper)
        self.assertTrue((PAPER/"plot_h3_climate_morphospace.py").is_file())

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
