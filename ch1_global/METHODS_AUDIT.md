# Chapter 1 — methods audit: peer-review resilience, robustness, and the foundations of the macro-analysis

Purpose: enumerate, in the order a reviewer will attack them, the ways the
current (v1) global orientation analysis could fail scrutiny, and state for
each the concrete remediation that the v2 redesign must implement. This is a
critical self-review, not a defence. Findings are graded:

- **BLOCKER** — the reported result is not defensible until this is fixed.
- **MAJOR** — a competent reviewer will demand this; likely a revision-blocker.
- **MINOR** — should be addressed for a strong paper but unlikely to sink it.

Each finding links to where the v2 plan already covers it
([`README.md`](./README.md), [`../data/ch1_image_traits/cirsium_ch1_image_trait_dictionary.csv`](../data/ch1_image_traits/cirsium_ch1_image_trait_dictionary.csv))
or flags it as **NEW** (not yet in the redesign and needing to be added).

---

## A. Label validity and leakage (the classifier itself)

### A1 — Species-level labels applied to individual images (BLOCKER)
`ch1_shared/make_training_data_from_labels.py` derives orientation from
**Kahaku specimen species descriptions** (`orientation_phrase` text) and joins
them to iNaturalist images **by species name**. Every image of a species
described as "nodding" is labeled nodding regardless of what that photo shows.
Consequences:
- The label is a **species attribute, not an image measurement.** The
  classifier can reach high accuracy by recognizing species gestalt (leaf
  outline, colour, habit) rather than orientation.
- Any environmental signal downstream is partly "where does this *species*
  grow," i.e. a handful of evolutionary events, not thousands of independent
  observations (see B2).

Remediation (v2): image-level labelling by ≥2 annotators with adjudication and
reported Cohen's/weighted kappa, per the trait dictionary
(`head_orientation_class` protocol). The species-name join is retired.
**Covered** in README §3.2, but the audit stresses this is the single most
important fix — without it, nothing else in Chapter 1 is trustworthy.

### A2 — Random train/test split → species and observation leakage (BLOCKER)
`03_train_head_direction_classifier.py` uses
`train_test_split(..., stratify=label_id)`, so crops from the same observation
and the same species appear in both train and test. Combined with A1, the
reported Acc 0.896 is almost certainly **optimistically biased**: the model is
tested largely on species it has already seen.
Remediation (v2): grouped CV by source photo/observation, by species, and by
geographic block, plus **leave-one-species-out** with the full per-species
distribution reported (not just the mean). **Covered** in README §3.3. The
random-split number is retained only as a leakage-quantifying reference, never
as the headline.

### A3 — No image-level ground truth / inter-annotator agreement (MAJOR)
There is currently no human double-labelling of the images themselves, so
there is no measurable label reliability floor. A reviewer cannot tell whether
model error or label noise dominates. Remediation: A1's two-annotator protocol
produces a kappa that bounds achievable accuracy. **Covered** (README §3.2).

### A4 — Decision thresholds and the discarded "uncertain" band (MINOR)
`NODDING_HIGH = 0.70 / UPWARD_HIGH = 0.30` are set a priori; the middle band
is dropped when forming "conservative" labels. Dropping uncertain cases is
**non-random** (harder/intermediate morphologies, certain species/angles) and
can bias both the accuracy estimate and the downstream environmental sample.
Remediation: calibrate thresholds on a held-out calibration set, report the
sensitivity of all downstream results to the threshold, and prefer carrying
predicted probabilities into the analysis over hard-thresholding (see D1).
**Partially covered**; make the threshold-sensitivity explicit in README §3.4.

---

## B. Macro-analysis foundations (the environmental models)

These are the load-bearing assumptions of any global macroecological claim.
If they fail, "trait X corresponds to climate Y" is not supported regardless
of how good the classifier is.

### B1 — Pseudoreplication from spatial sampling bias (BLOCKER)
`05_rf_glm_gam_env_analysis.R` treats each iNaturalist observation as an
independent data point. iNat sampling is heavily clustered by road access,
population centres, and observer density. Thousands of near-duplicate records
in a few well-surveyed regions inflate effective sample size and significance.
The script does include space terms (`z_lon`, `z_lat`, polynomials) and a
`space-only` vs `env` vs `env+space` comparison — a genuine strength — but a
polynomial trend surface is a **weak** control for fine-scale clustering.
Remediation (v2):
- Spatial thinning / grid aggregation of records before modelling.
- Report results at photo-level, species-mean level, **and** species×grid-cell
  level (aggregation-level sensitivity, README §3.4) — the single most
  convincing robustness check for a reviewer.
- Consider a spatial random effect or spatial block CV rather than only a
  trend-surface polynomial. **NEW: add explicit spatial CV / mixed-model
  option to README §3.4; currently only aggregation-level sensitivity is
  named.**

### B2 — Phylogenetic / taxonomic non-independence (MAJOR)
Because v1 labels are species-level (A1), the response is effectively
one-value-per-species smeared over many photos. Even after fixing A1 to
image-level labels, species remain non-independent (shared ancestry). Chapter
1 correctly **defers formal phylogenetic correction to Ch.2/3**, but must not
present a raw per-observation GLM as if observations were independent draws.
Remediation: species-level and species×grid aggregation (README §3.4) as the
honest macro unit; state plainly that Chapter 1's inference is descriptive and
that phylogenetic non-independence is addressed in Ch.2/3. **Covered** as a
boundary (README §5) — but the audit adds: **the species-level model, not the
per-observation model, should be the headline** for any climate claim.

### B3 — Environmental extraction not reproducible in-repo (MAJOR)
`05` starts from a pre-computed `..._CHELSA_topography_soilgrids.csv` that no
script in the repository produces. The download/extraction (CHELSA, SRTM-
derived topography, SoilGrids) — including CRS handling, resolution, and how
point values are sampled — is invisible and therefore unreviewable.
Remediation (v2 / **NEW**): add a committed extraction script (or a documented,
version-pinned data-provenance record: variable list, source version,
resolution, extraction date, buffer/aggregation rule). Add to the README §4
roadmap as an explicit phase.

### B4 — Coordinate quality: obscured/imprecise records (MAJOR)
The collection pipeline sets `EXCLUDE_OBSCURED = False` and does not filter on
`positional_accuracy`. iNat obscured coordinates are randomized within a ~0.2°
box; low-accuracy records add location error that propagates into every
extracted environmental value (attenuating climate–trait relationships).
Remediation (v2 / **NEW**): exclude obscured records and filter by
`positional_accuracy` for the environmental analysis (a coordinate-quality
sensitivity tier), while they may still be usable for the trait-only syndrome
description. Add to README §3.4 sensitivity tiers.

### B5 — Class imbalance handled by downsampling (MINOR)
RF uses `balance_binary_df` (downsampling the majority class); GLM/GAM use
class-balanced case weights. Downsampling discards data and adds run-to-run
variance; the script mitigates with 20 repeats, which is good. Keep the
repeats, report base rates alongside balanced metrics, and prefer AUC/PR over
raw accuracy for the imbalanced classes. **Mostly fine**; document the choice.

---

## C. Detection / measurement chain

### C1 — Unreported detection performance (MAJOR)
Everything downstream is conditional on YOLO detecting the right region, yet no
detection precision/recall/mAP is reported, and detection failures are a
**non-random** filter on which images (species, backgrounds, image quality)
enter the analysis. Remediation (v2): report detector precision/recall on a
held-out manually boxed test set; characterize what gets dropped. Add to
README §3.2 explicitly for the detector (currently framed mainly around trait
classifiers).

### C2 — Non-portable, non-turnkey pipeline (MINOR, but reproducibility-relevant)
Hardcoded `C:\Users\...` paths and reliance on uncommitted intermediate files
mean the analysis cannot be re-run by a reviewer or a collaborator. Remediation
(v2): config-driven paths, a manifest of expected inputs/outputs, and the
missing extraction step from B3. See v1 README "Known non-portability."

### C3 — Colour, camera tilt, and viewing angle (MAJOR for v2 traits)
For the new colour/cover/shape traits: citizen-science photos have no colour
calibration, uncontrolled white balance, arbitrary camera tilt (contaminating
the gravity-frame orientation), and arbitrary viewing angle (a head looks more
"covered" from the front than the side). Remediation: per-image white balance
with relative/ordinal colour only; tilt screening for orientation; a
lateral-vs-frontal view-angle qualifier for cover — all already specified in
the trait dictionary and README §2. **Covered**, listed here so the reviewer-
facing limitations section is complete.

---

## D. Inference hygiene

### D1 — Errors-in-variables: predicted traits used as if observed (MAJOR)
v1 (and, if not designed carefully, v2) feeds **hard-thresholded predicted
labels** into the RF/GLM/GAM as though they were ground truth. Classification
error in a predictor causes regression dilution (attenuated, sometimes
sign-distorted coefficients) and understated uncertainty. Remediation (v2 /
**NEW — strengthen README §3.4**):
- carry predicted **probabilities/continuous estimates** into the models
  rather than hard labels where possible;
- restrict to high-confidence predictions as one sensitivity tier and report
  how conclusions move;
- ideally propagate classifier uncertainty (e.g. multiple imputation over
  predicted labels, or a measurement-error model).
This is distinct from A2 (which is about honest accuracy estimation); D1 is
about not letting known measurement error silently bias the ecological
coefficients.

### D2 — Multiple testing across traits and variables (MINOR)
v2 runs many models (5+ traits × many environmental predictors × several
aggregation levels). Without a stated inference philosophy this invites
cherry-picking. Remediation: pre-declare the primary trait(s)/hypotheses,
separate confirmatory from exploratory analyses, and use effect sizes with CIs
plus consistency-across-sensitivity-tiers as the evidence standard rather than
a wall of p-values. **NEW: add an "inference philosophy" note to README §3.4.**

### D3 — Redundant syndrome axes presented as independent (MINOR)
`involucral_cover_ratio` and `corolla_involucre_display_ratio` are derived from
the same masks and may be near-collinear; treating them as two independent
"findings" would overstate structure. Remediation: report their correlation
explicitly and collapse if r is high (already H1-Ch1d in README §1.2, and the
trait dictionary note). **Covered.**

---

## E. Foundations summary — what must hold for the macro-claim to stand

A global "image-derived trait ↔ environment" claim in Chapter 1 rests on five
pillars. The table states each pillar, the v1 status, and the v2 gate.

| Pillar | v1 status | Gate that must pass in v2 |
|---|---|---|
| **Measurement validity** — the image trait measures the biological trait | Fails (A1: label is species attribute) | Image-level labels, kappa reported, extraction models validated vs humans |
| **Generalization** — the model reads the trait, not the species | Untested (A2: leakage) | LOSO + grouped CV, per-species distribution reported |
| **Independent replication** — data points are not pseudoreplicates | Fails (B1/B2) | Species and species×grid aggregation as the headline units; spatial control |
| **Environmental data quality** — predictors are accurate and reproducible | Weak (B3/B4) | Committed/provenanced extraction; obscured & low-accuracy record handling |
| **Honest uncertainty** — measurement error doesn't masquerade as signal | Fails (D1) | Probabilities not hard labels; confidence-tier sensitivity; error propagation |

Chapter 1's scientific contribution is only as strong as its weakest pillar.
As of v1, three of the five **fail** and two are **weak**; the v2 redesign is
built specifically to move all five to "pass," and it should not be written up
as a finished macroecological result until they do.

---

## F. New items to fold into the v2 plan

These are findings above that are **not yet** explicit in `README.md` and
should be added when the plan is next revised:

- B1: spatial CV / spatial random effect option (not only aggregation-level).
- B3: committed environmental-extraction / data-provenance phase in the §4 roadmap.
- B4: obscured-coordinate and positional-accuracy sensitivity tier in §3.4.
- C1: detector precision/recall reporting alongside the classifier validation in §3.2.
- D1: errors-in-variables handling (probabilities, confidence tiers, propagation) in §3.4.
- D2: a pre-declared inference philosophy (confirmatory vs exploratory) in §3.4.
