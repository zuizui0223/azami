# Chapter 1 v3 workflow-only revision scope

Status: active from 2026-09-09

## Purpose

The current v3 revision is restricted to **workflow improvement only**. Do not add a new scientific analysis, endpoint, predictor, hypothesis, model family, multiplicity family, synthesis, or post hoc rescue.

The goal is a short, reproducible route from the recovered source to the existing ecological analysis.

## Canonical path

There is one canonical execution order:

1. `source`
2. `measurement_qc`
3. `ecology`
4. `output`

`python -m analysis.v3.integrated_preflight` is the single readiness entry point.

Assessability is not a separate scientific stage. Hypervolume/breadth synthesis is not a canonical stage.

## Minimal Measurement / QC gate

Measurement/QC is complete when all of the following are true:

1. the existing measurement definitions and QC route decisions remain frozen;
2. the already-scheduled raw measurement outputs are recovered and integrated;
3. all 27 raw endpoints remain available in the retained numerical record;
4. the predeclared primary module representation is built on that reconciled cohort;
5. a simple completion/missingness ledger records complete, pending and source-blocked observations plus endpoint/module missingness.

Nothing else is required to enter Ecology.

In particular, the following are **not canonical gates**:

- a replacement/alternative chroma definition;
- taxon-, region- or environment-specific assessability models;
- a separate model-diagnostics research stage;
- additional tail/repair/calibration simulation branches;
- all-27 exploratory ecological coefficients;
- hypervolume or breadth synthesis.

Existing outputs from those activities may remain as technical or historical evidence, but they cannot block the canonical path or change endpoint admission.

## Ecology gate

Once Measurement/QC is complete, run only the already-declared ecological model and nuisance design on the fixed module cohorts. Do not add a second model because a result is weak, null, or difficult to interpret.

## Output

After the canonical ecology run, generate figures, tables, manuscript-facing summaries and reproducibility records. Output generation must not trigger new scientific analyses.

## Allowed workflow changes

- clarify stage ownership and ordering;
- remove duplicated execution paths;
- improve checkpoint/restart and artifact recovery;
- route status checks through the canonical preflight;
- separate canonical, technical, exploratory and historical surfaces;
- consolidate documentation and runbooks;
- add integrity tests that verify this workflow without changing calculations.

## Prohibited changes

Do **not**:

- add another ecological analysis to rescue a weak or null result;
- open additional endpoints beyond already-declared routes;
- add predictors, interactions, response transformations or model families;
- add a new multiple-testing family;
- promote a held measurement from exploratory results;
- turn measurement success/failure into a separate ecological research program;
- extend tail/repair/calibration branches after the current bounded work;
- reinterpret technical repeatability as physical trait accuracy;
- change frozen v2 outputs.

## Stop rule

A workflow issue is closed when:

1. the canonical stage can run, or is blocked by one named missing input;
2. restart/replay preserves the same scientific inputs and calculations; and
3. the status is visible from the integrated preflight.

Once those conditions hold, stop adding diagnostics for that issue.

## Definition of progress

Progress means fewer gates and a shorter path from source to canonical result. It does **not** mean more models, endpoints, coefficient tables, calibration branches, or exploratory findings.
