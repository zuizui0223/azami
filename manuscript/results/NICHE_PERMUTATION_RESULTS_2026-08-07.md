# Completed environmental-niche permutation/null tests

Workflow run: `31152400487`  
Permutation count: 10,000 per trait and scope  
Primary scope: 148 species complete for all nine primary traits and available environmental variables  
Direct-backbone sensitivity: 49 complete species directly represented in the dated backbone

The null shuffled each trait among the same species while keeping the environmental PCA, environmental availability and trait distribution fixed. Benjamini-Hochberg correction was applied separately across the nine primary traits for centroid-distance and Bhattacharyya-overlap tests within each scope.

## All complete species

Significant centroid-distance contrasts after BH correction were:

- orientation: centroid distance 2.398, q = 0.00090;
- chroma: 1.632, q = 0.0185;
- hue sine: 1.807, q = 0.00990;
- hue cosine: 2.271, q = 0.00090;
- width-profile CV: 1.663, q = 0.0185.

Significantly lower-than-null Gaussian niche overlap was supported for:

- orientation: overlap 0.746, q = 0.00060;
- chroma: 0.694, q = 0.00060;
- hue sine: 0.817, q = 0.0227;
- hue cosine: 0.641, q = 0.00060.

Lightness, aspect ratio, circularity and solidity were not BH-supported in either metric. Width-profile CV was supported for centroid separation but not overlap after BH correction (q = 0.0581).

## Direct-backbone sensitivity

Within the 49 complete taxa directly represented in the dated backbone, support was narrower:

- hue cosine retained both centroid-distance support (q = 0.00180) and overlap support (q = 0.00090);
- lightness retained overlap support (q = 0.00765);
- hue sine retained overlap support (q = 0.0480).

Orientation, chroma and the four outline traits did not retain BH support in the direct-backbone subset. This sensitivity is expected to have lower power and a changed taxonomic/geographic composition; it is not a replacement for the all-species descriptive contrast.

## Submission interpretation

The permutation test confirms that the strongest all-species environmental sorting of capitulum phenotype is concentrated in orientation and visible colour rather than being a generic consequence of quartile splitting. Historical/direct-backbone restriction weakens several contrasts, so the manuscript should distinguish broad among-species environmental sorting from historical robustness and should not infer causal adaptation from these contrasts.
