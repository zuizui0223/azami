# Methods source ledger

This ledger maps each subsection of `02_methods.md` to the repository implementation or control file that supports it. It also records which statements remain provisional until external audit products are supplied.

| Methods subsection | Repository source | Status |
|---|---|---|
| Study design and analytical separation | `manuscript/FINAL_MANUSCRIPT_STRATEGY.md`; `manuscript/final_claims.json` | Frozen |
| Public-image observations and taxonomic scope | `manuscript/final_claims.json`; `analysis/audit_taxonomic_freeze.py`; `manuscript/TAXONOMIC_FREEZE_PROTOCOL.md` | Counts frozen; decision table pending |
| Capitulum detection and crops | frozen production scripts in `ch1_global/v2`; `ch1_global/v2/79_evaluate_capitulum_detector.py`; Figure 1 audit workflow | Production workflow frozen; independent metrics pending |
| Continuous trait measurements | frozen production scripts in `ch1_global/v2`; `ch1_global/v2/78_build_figure1_measurement_audit.py`; `manuscript/FIGURE1_MEASUREMENT_AUDIT.md` | Measurement definitions frozen; demonstration panels pending |
| Quality control and observation retention | `analysis/ch1/pipeline.json`; `manuscript/FIGURE_TABLE_MAP.md`; workflow manifests | Frozen |
| Environmental data and spatial cohorts | `ch1_global/v2/75_run_exhaustive_within_species_climate_analysis.py`; frozen climate artifact; submission manifest | Core analysis frozen; exact source metadata to copy into final text |
| Visible variance and PCA | continuous-analysis outputs mapped in `manuscript/FIGURE_TABLE_MAP.md` | Frozen |
| Within-species variation and responsiveness | `ch1_global/v2/77_build_lability_axes_and_supplement.py`; `species_lability_axes.csv`; `sensitivity_minimum_sample_size.csv` | Frozen |
| Among-species environmental sorting | grouped SPDE-INLA scripts and outputs; trait-extreme niche outputs; `manuscript/FINAL_MANUSCRIPT_STRATEGY.md` | Frozen |
| Residual spatial and broad-region robustness | `analysis/prepare_spatial_diagnostic_input.py`; `analysis/audit_spatial_robustness.py`; associated workflows | Infrastructure frozen; diagnostic values pending |
| Historical-placement sensitivity | phylogenetic evidence workflows; `manuscript/final_claims.json`; molecular database audit | Frozen |
| Multiplicity and claim control | `manuscript/final_claims.json`; `analysis/validate_final_claims.py`; `manuscript/RUNBOOK.md` | Frozen |
| Reproducibility and archive | `analysis/ch1/pipeline.json`; `ch1_global/v2/69_validate_submission_outputs.py`; `ch1_global/v2/70_build_submission_manifest.py` | Archive pending |

## Mandatory replacement fields before submission

The draft intentionally avoids inventing values that are not yet present in frozen outputs. Before submission, insert and verify:

1. public-data source name, query date, filters and API/version details;
2. detector architecture, weight identifier, input resolution, confidence threshold and independent audit metrics;
3. exact mathematical definitions and units for all nine endpoints;
4. environmental dataset names, versions, resolutions, variables and extraction dates;
5. spatial-thinning distance and implementation details;
6. full SPDE-INLA formulas, priors, mesh construction and model-comparison rules;
7. PCA preprocessing and missing-data rules;
8. exact construction and scaling of the two lability indices;
9. multiplicity families and BH correction implementation;
10. phylogenetic backbone identity, placement rules and software versions;
11. residual Moran and leave-one-region-out diagnostic values;
12. software versions, workflow run IDs, artifact digests and permanent archive DOI.

## Freeze rule

A missing reporting detail may be filled from the exact frozen implementation or artifact. A change to taxa, trait definitions, thresholds, environmental inputs, model formulas, priors, tree scenarios or accepted output values requires a new analysis version and cannot be treated as manuscript editing.
