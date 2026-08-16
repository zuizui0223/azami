import unittest

import pandas as pd

from analysis.export_eazami_japan38_exhaustive_traits import ENDPOINTS, binomial, build


class TestJapan38ExhaustiveTraits(unittest.TestCase):
    def observations(self):
        rows = []
        for index, taxon in enumerate(["Cirsium alpha", "Cirsium alpha", "Cirsium beta", "Cirsium other"]):
            row = {
                "taxon_name": taxon,
                "obs_id": str(index),
                "spatial_block_5deg": f"b{index % 2}",
                "coordinate_usable_bool": True,
            }
            for endpoint in ENDPOINTS:
                row[endpoint] = float(index + 1)
            rows.append(row)
        return pd.DataFrame(rows)

    def membership(self):
        return pd.DataFrame([
            {"paper_japan_member_id": "J1", "paper_taxon_concept": "Cirsium alpha var. one"},
            {"paper_japan_member_id": "J2", "paper_taxon_concept": "Cirsium alpha var. two"},
            {"paper_japan_member_id": "J3", "paper_taxon_concept": "Cirsium beta"},
            {"paper_japan_member_id": "J4", "paper_taxon_concept": "Cirsium missing"},
        ])

    def test_binomial_reduction(self):
        self.assertEqual(binomial("C. alpha var. one"), "Cirsium alpha")

    def test_build_preserves_infraspecific_boundary(self):
        traits, coverage, summary = build(self.observations(), self.membership())
        self.assertEqual(set(traits["taxon_name"]), {"Cirsium alpha", "Cirsium beta"})
        alpha = traits.loc[traits["taxon_name"].eq("Cirsium alpha")].iloc[0]
        self.assertEqual(alpha["n_paper_concepts"], 2)
        self.assertIn("no_infraspecific_assignment", alpha["trait_identity_scope"])
        self.assertEqual(summary["n_paper_concepts_represented_at_binomial_level"], 3)
        missing = coverage.loc[coverage["species_binomial"].eq("Cirsium missing")].iloc[0]
        self.assertFalse(bool(missing["azami_detector_positive_trait_taxon_present"]))


if __name__ == "__main__":
    unittest.main()
