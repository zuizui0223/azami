# Chapter 1 v3 workflow-only revision scope

Status: active from 2026-09-09

## Purpose

The current v3 revision is restricted to **workflow improvement only**.
It must not add a new scientific analysis, endpoint, predictor, hypothesis,
model family, multiplicity family, ecological synthesis, or post hoc rescue.
The scientific scope and estimands remain those already declared in
`integrated_workflow_contract.json`.

The purpose of this phase is to make the existing design easier to execute,
audit, and explain after incorporating advisor/reviewer comments upstream in
the source-to-inference workflow rather than as downstream manuscript defence.

## Canonical path

There is one canonical execution order:

1. `source`
2. `measurement`
3. `assessability`
4. `ecology`
5. `synthesis` only if its already-declared optional qualification is met

`python -m analysis.v3.integrated_preflight` remains the single readiness
entry point. Earlier contracts and receipts remain historical evidence unless
the integrated contract explicitly supersedes them.

## What may change in this phase

Allowed changes are limited to workflow mechanics:

- clarify stage ownership and ordering;
- make prerequisites and stop conditions explicit;
- remove duplicated execution paths;
- route existing checks through the canonical preflight;
- separate canonical, sensitivity, exploratory, and historical surfaces;
- improve checkpoint/restart behaviour and artifact provenance;
- rename ambiguous workflow surfaces without changing calculations;
- consolidate documentation and runbooks;
- add integrity tests that verify the existing contract is followed.

## What must not change in this phase

Do **not**:

- add another ecological analysis because an existing result is weak or null;
- open additional endpoints beyond already-declared routes;
- add new environmental predictors or interactions;
- add new response transformations or alternative model families;
- add a new multiple-testing family;
- extend all-27 point exploration into confirmatory inference;
- promote a held measurement because its exploratory coefficient looks useful;
- create further tail/repair/aggregate analysis branches to explain an existing
  sensitivity result;
- reinterpret technical repeatability as physical trait accuracy;
- change frozen v2 outputs.

If an existing stage fails, record the failure and either repair the workflow
mechanics or lower/withhold the affected claim. Do not create a new analysis to
rescue it during this workflow-only phase.

## Surface classification

### Canonical

The five stages in `integrated_workflow_contract.json` and the integrated
preflight are the only canonical scientific workflow surface.

### Sensitivity / technical audit

Existing calibration, geometry, resolution, dependence, and measurement
checks may be retained when they test a declared assumption of the canonical
workflow. They are supporting evidence, not additional scientific questions.
No new sensitivity family should be created unless a canonical stage cannot be
executed correctly without it.

### Exploratory

The existing all-27 point-estimate exploration and other explicitly labelled
exploratory checkpoints remain non-gating. They cannot change route admission,
primary hypotheses, predictor choice, multiplicity, or manuscript headline
claims in this phase.

### Historical

Superseded contracts, pilot receipts, repair receipts, and previous execution
waves remain for provenance. They are not alternative active workflows.

## Stop rule

A workflow issue is considered closed when all three are true:

1. the canonical stage can be executed or is explicitly blocked with a named
   unmet prerequisite;
2. restart/replay preserves the same scientific inputs and calculations; and
3. the resulting state is visible from the integrated preflight or its pinned
   evidence index.

Once these conditions hold, do not continue adding diagnostics for that issue
unless they are required to correct a newly observed workflow failure.

## Definition of progress

During this phase, progress means fewer ambiguous paths and a shorter route from
source inventory to a reproducible canonical result. It does **not** mean more
models, more endpoints, more coefficient tables, or more exploratory findings.
