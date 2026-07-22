# Methods

## Study design and analytical separation

We used public biodiversity photographs to quantify continuous capitulum traits across *Cirsium* and related thistles while retaining repeated observations within species. The design separated three structures that comparative trait studies often conflate: visible dispersion among records assigned to the same species, spatial environment–trait associations within species and environmental sorting among species. These layers were not combined into an inference about plasticity, adaptation or evolutionary rate.

## Cohort definitions

Two executed data streams served different purposes. The balanced image-comparison atlas contained 3,725 public observations, 6,626 detected capitula and 216 accepted image-analysis taxa. It was used for visible variance partitioning, species-level trait PCA, among-species summaries and historical-placement sensitivity.

The exhaustive stream ran the detector and continuous measurement workflow before trait-based thinning. It contained 406,582 observations with detected capitula from 286 taxa. Coordinate filtering retained 392,989 observations from 271 taxa; a positional-accuracy threshold of ≤10 km retained 297,293 observations from 259 taxa; and one observation per taxon × 0.25° cell yielded the **exhaustive spatially thinned primary cohort** of 46,276 observations from 259 taxa. The complete cohort flow, filenames and analysis permissions are fixed in `manuscript/COHORT_FLOW_AND_ANALYSIS_LEDGER.md`.

## Capitulum detection and continuous image measurements

Capitula were detected from source photographs using the frozen production detector. Each retained detection preserved the source observation identifier, image provenance, bounding box and contextual crop. Tight crops were used for colour and two-dimensional outline; context crops retained stem and display information needed for image-referenced orientation.

Nine primary endpoints represented three capitulum modules: orientation; visible corolla lightness, chroma and circular hue coordinates; and outline aspect ratio, circularity, solidity and width-profile variation. Trait measurements were deterministic production functions rather than human categories. A failed or unassessable measurement remained missing rather than being converted to zero or biological absence. Horizontal mirroring served as a technical repeatability check. Detector precision and continuous-trait accuracy require independent manual validation and are not inferred from quality-control retention.

## Environmental predictors and primary within-species coefficients

Four predeclared CHELSA v2.1 predictors represented mean annual temperature (BIO1), temperature seasonality (BIO4), annual precipitation (BIO12) and precipitation seasonality (BIO15). For each endpoint–predictor combination in the exhaustive spatially thinned cohort, outcome and predictor values were demeaned within species and standardized after demeaning. Ordinary least-squares coefficients were fitted without an intercept and used species-clustered standard errors. Benjamini–Hochberg correction was applied across the 36 endpoint-component models. Hue sine and cosine rows were retained for computation but were not interpreted as independent biological colour tests.

## Visible variance and species-level trait structure

For each endpoint in the image-comparison atlas, total sums of squares were separated into within-assigned-species and among-species-mean components. The within fraction describes uncontrolled image variance and can combine biological heterogeneity, developmental state, illumination, viewpoint and measurement error. Species medians were standardized and analysed by PCA to describe multivariate trait architecture.

## Reviewer-driven precision audit

The legacy lability analysis summarized each species by the RMS absolute value of separately fitted slopes. Because an absolute noisy estimate is positive even under a zero true effect, this score is expected to increase as slope uncertainty increases. We therefore audited its relation to median species-specific sample size and median slope standard error and calculated a rank-based partial correlation between the legacy axes while controlling median sample size.

The negative legacy relation and median-split quadrants were withdrawn after this audit. The revised precision-aware cohort required all seven linear endpoints, all four climate predictors for each endpoint and at least 10 observations for every species × endpoint × predictor slope. This yielded 101 taxa. Circular hue was excluded from this specific correction because the archived joint species hue vectors lack component standard errors.

## Equal-module visible variation and sampling-noise-adjusted association energy

Within-variation values were first averaged within orientation, colour and shape and then averaged across the three modules, so modules rather than endpoint counts received equal weight.

For slope estimate \(\hat\beta\) with standard error \(s\), we used \(\hat\beta^2-s^2\) as a sampling-noise-adjusted estimator of squared true association magnitude, because under an approximately normal slope estimator \(E(\hat\beta^2-s^2)=\beta^2\). Values were averaged across predictors within trait, across traits within module and equally across modules. Negative realized values indicate that the data do not resolve excess association beyond sampling noise; they are not negative biological responsiveness. Spearman correlation described the relation between equal-module visible variation and this association-energy score. A 5,000-replicate species bootstrap provided a confidence interval.

## Hierarchical variance meta-regression

As the primary uncertainty-aware test, all 2,828 archived linear slope estimates from the 101 complete taxa were modelled with their standard errors. For trait–predictor group \(g\) and species \(i\),

\[
\hat\beta_{ig} \sim \mathrm{Normal}(\mu_g,\; s_{ig}^2 + \tau_g^2\exp(bV_i)),
\]

where \(V_i\) is the standardized equal-module visible-variation index, \(\mu_g\) is a group-specific mean, \(\tau_g^2\) is group-specific latent slope variance at mean visible variation and \(b\) is a common log-variance change per standard deviation of visible variation. Group means were profiled analytically and variance parameters by maximum likelihood. The common coefficient was tested against \(b=0\) with a likelihood-ratio test, and a 95% interval was obtained by profile likelihood.

## Among-species environmental sorting and historical sensitivity

Grouped spatial models and trait-extreme environmental PCA contrasts were kept separate from within-species coefficient analyses. Their outputs describe among-species environmental sorting, not within-species change. Historical sensitivity was evaluated across direct-backbone and alternative grafted placements. Direct and grafted results remained separate, and opportunistic molecular records were summarized only as a coverage audit.

## Multiplicity, terminology and reproducibility

Every FDR count is reported with its cohort, endpoint family and number of tests. The 46,276-observation exhaustive primary cohort is not conflated with the earlier balanced-atlas ≤10 km sensitivity. The manuscript uses *visible dispersion*, *within-species spatial environment–trait association*, *among-species environmental sorting* and *historical-placement sensitivity*. It does not use *climate tracking* or *environmental responsiveness* as if temporal or experimental response had been measured.

The revised analysis is implemented in `analysis/reanalyze_lability_precision.py`, rerun from frozen artifact `8330350031` by `.github/workflows/ch1-reviewer-precision-reanalysis.yml`, and validated against `manuscript/final_claims.json`.
