# Chapter 1 v3 offline numerical audit

`audit_revision_consistency.py` reads existing local CSV files only. It does not
retrieve or transform images, rerun the detector, or modify frozen v2 inputs,
models, multiplicity families or result tables. Manuscripts and reviewer responses
remain outside this repository.

## Verified scope

The executed report is `reproducibility/v3_offline_audit_20260907.json`.

- All 37,251 original orientation observation medians were reconstructed from
  102,982 usable heads; maximum absolute discrepancy was 2.84e-14 degrees.
- The declared hash-ordered subset was reconstructed: 142 taxa, 3,829
  observations, 10,431 heads and 6,561 photographs. This verifies selection, not
  image remeasurement or measurement stability.
- The audit recomputes 22 conditional model rows and ten minimum-replication rows.
  The orientation `climate_core` and `plus_chelsa_bio01` rows intentionally have
  the same design and must not be counted as independent evidence.
- Two local executions produced byte-identical copies of all six audit outputs.
  Twenty-nine new synthetic unit tests and 25 existing integrity tests passed.
  These checks establish numerical/code consistency, not independent biological
  validation or a new confirmatory analysis.

The all-nine conditional orientation coefficient has an exploratory HC3 interval
that includes zero. Chroma retains a negative coefficient in the listed models.
These results do not identify independent causal environmental effects. Strongly
correlated predictors and unmodelled spatial/phylogenetic dependencies remain.

## Required local inputs

The command verifies the exact SHA-256 values defined in `HASHES` before fitting.

| Argument | Original member / provenance |
|---|---|
| `--traits` | `universe/continuous_trait_universe_observation_long.csv` in continuous artifact `9612943217` |
| `--heads` | `8269246732_exhaustive_continuous_head_level.csv` in recovered artifact `10004125659` |
| `--environment` | Exact reconstructed `strict_spatial_chelsa_process.csv`, SHA-256 `e242aa7ce69d12b11937c1335e84b9638799c50b42ef36b95725e77190df98e7` |

Use the frozen public v2 runbook for original numerical-input recovery. The head
recovery artifact is a separate v3 input; its identity does not make this a new
public archival release. No unverified replacement CSV is accepted.

From the repository root, in an environment containing numpy, pandas and scipy:

```bash
python analysis/v3/audit_revision_consistency.py \
  --traits /path/to/continuous_trait_universe_observation_long.csv \
  --heads /path/to/8269246732_exhaustive_continuous_head_level.csv \
  --environment /path/to/strict_spatial_chelsa_process.csv \
  --out-dir /path/to/external-v3-audit

python -m pytest -q tests/test_v3_revision_consistency.py
```

Conditional coefficients use separately standardized predictors and response,
with HC3 sandwich covariance and a normal-reference 95% interval. Environmental
medians are from all original primary observations, not recalculated after
endpoint missingness or minimum-replication restrictions. These exploratory
intervals are not spatial, phylogenetic or multiplicity-adjusted. Output hashes
record the executed software environment; numerical agreement should also be
checked when software versions or CSV serialization differ.

## Linked records remain a separate verification step

The executed audit had no linked bounding-box perturbation CSV. Its recorded
status is therefore `NOT_VERIFIED_RECORDS_NOT_SUPPLIED`, not failed remeasurement
and not a claim that remeasurement was never conducted elsewhere. Aggregate
reported values cannot substitute for per-head linked records.

An optional `--bbox-records` argument accepts a normalized CSV with:

```text
annotation_unit_id,obs_id,photo_id,condition,angle_deg,usable,head_clipped,source_sha256
```

Exactly five rows are required per scheduled head, with conditions `baseline`,
`x_minus_5pct`, `x_plus_5pct`, `y_minus_5pct` and `y_plus_5pct`. Failed executions
must remain explicit rows with `usable=False`, rather than disappearing.

The validator checks IDs, condition completeness, image-hash consistency,
booleans, angle bounds, original-baseline agreement and nested 5/10/20-degree
stability flags. It does not verify image authenticity, perform image operations,
or fit stable-subset environmental models. Adaptation to another record schema
must be explicit and provenance-preserving; missing fields must not be invented.
