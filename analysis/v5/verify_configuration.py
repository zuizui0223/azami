"""Check the aggregate tables against the previously exported V5 comparisons."""
from pathlib import Path
import numpy as np
import pandas as pd


def main():
    results = Path(__file__).parent / "results"
    root = results / "configuration"
    convergence = pd.read_csv(root / "v5_configuration_convergence.csv")
    previous = pd.read_csv(results / "directional_traits/v5_trait_convergence.csv").query("metric == 'joint'")
    np.testing.assert_allclose(convergence.query("dimensions == 3").set_index("degrees").observed,
                               previous.set_index("degrees").observed)
    head = pd.read_csv(results / "directional_traits/novel_head_overlap.csv")
    paired = head.query("design == 'geographic_pairs' and native95_mass == .95 and bandwidth == .5")
    control = pd.read_csv(root / "configuration_overlap.csv").query("within_block")
    for row in control.itertuples():
        traits = "angle+chroma" if row.dimensions == 2 else "angle+chroma+outline"
        original = paired[(paired.traits == traits) & (paired.species == row.species)].iloc[0]
        np.testing.assert_allclose([row.observed_D_sim, row.observed_D_new, row.observed_drop],
                                   [original.D_native_analogous, original.D_native_novel, original.drop_in_overlap], atol=1e-3)
    vectors = pd.read_csv(root / "matched_centroids.csv")
    assert len(vectors) == 6 and not vectors.duplicated(["species", "trait"]).any()
    np.testing.assert_allclose(vectors.novel - vectors.analogous, vectors.novel_minus_analogous, atol=1e-10)
    assert vectors.p_block_flip.between(0,1).all() and (vectors.holm_six >= vectors.p_block_flip-1e-12).all()
    assert (vectors.groupby("species").pairs.first() == control.query("dimensions == 3").set_index("species").pairs).all()
    slopes = pd.read_csv(root / "s1_slope_heterogeneity.csv").query("spec == 'warm'")
    assert (slopes.observed_positive + slopes.observed_negative == slopes.n_taxa).all()
    assert set(slopes.loc[slopes.trait=="chroma", "axis"]) == {"bio05", "rsds_mean"}
    print("Reproduced distance and overlap contrasts, signed centroid differences, and taxon counts.")


if __name__ == "__main__":
    main()

