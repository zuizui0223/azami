#!/usr/bin/env python3
"""Build and execute the local PR #69 × continuous-trait audit notebook."""
from __future__ import annotations

import ast
from contextlib import redirect_stderr, redirect_stdout
import io
import os
from pathlib import Path
import traceback

import nbformat as nbf
from nbclient import NotebookClient


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "analysis" / "notebooks" / "ch1_pr69_continuous_reanalysis.ipynb"


def markdown(text: str):
    return nbf.v4.new_markdown_cell(text.strip())


def code(text: str):
    return nbf.v4.new_code_cell(text.strip())


def _result_output(value, execution_count: int):
    data = {"text/plain": repr(value)}
    html = getattr(value, "_repr_html_", None)
    if callable(html):
        rendered = html()
        if rendered:
            data["text/html"] = rendered
    return nbf.v4.new_output(
        "execute_result",
        data=data,
        metadata={},
        execution_count=execution_count,
    )


def execute_in_process(notebook):
    """Execute cells sequentially when sandbox policy blocks a Jupyter TCP kernel."""
    namespace = {"__name__": "__main__"}
    previous = Path.cwd()
    os.chdir(ROOT)
    try:
        execution_count = 0
        for cell in notebook["cells"]:
            if cell["cell_type"] != "code":
                continue
            execution_count += 1
            cell["execution_count"] = execution_count
            outputs = []
            stdout = io.StringIO()
            stderr = io.StringIO()
            try:
                tree = ast.parse(cell["source"], filename=f"notebook-cell-{execution_count}")
                final_expression = None
                if tree.body and isinstance(tree.body[-1], ast.Expr):
                    final_expression = ast.Expression(tree.body.pop().value)
                with redirect_stdout(stdout), redirect_stderr(stderr):
                    if tree.body:
                        exec(compile(tree, f"notebook-cell-{execution_count}", "exec"), namespace)
                    result = (
                        eval(compile(final_expression, f"notebook-cell-{execution_count}", "eval"), namespace)
                        if final_expression is not None
                        else None
                    )
                if stdout.getvalue():
                    outputs.append(nbf.v4.new_output("stream", name="stdout", text=stdout.getvalue()))
                if stderr.getvalue():
                    outputs.append(nbf.v4.new_output("stream", name="stderr", text=stderr.getvalue()))
                if final_expression is not None and result is not None:
                    outputs.append(_result_output(result, execution_count))
            except Exception as error:
                outputs.append(
                    nbf.v4.new_output(
                        "error",
                        ename=type(error).__name__,
                        evalue=str(error),
                        traceback=traceback.format_exc().splitlines(),
                    )
                )
                cell["outputs"] = outputs
                raise
            cell["outputs"] = outputs
    finally:
        os.chdir(previous)
    notebook["metadata"]["execution"] = {
        "status": "executed_top_to_bottom",
        "executor": "in_process_fallback_after_sandbox_blocked_jupyter_tcp",
    }
    return notebook


def main() -> None:
    notebook = nbf.v4.new_notebook()
    notebook["metadata"] = {
        "kernelspec": {
            "display_name": "Python 3",
            "language": "python",
            "name": "python3",
        },
        "language_info": {"name": "python", "version": "3.12"},
    }
    notebook["cells"] = [
        markdown(
            """
# Chapter 1 PR #69 × continuous-trait local reanalysis

## tl;dr

- 27 continuous endpoints are explicitly defined; 17 have frozen numerical measurements locally and 10 remain genuinely unexecuted.
- The local ≤10 km high-resolution diagnostic screen ran 48 endpoint–climate models after a CSV round-trip bug was fixed. Four unadjusted candidate rows passed the tier-wise BH screen, but all are screening-only and none survives the locked PR #69 image-quality/stratum gate.
- All 72 locked PR #69 involucre model rows were reproduced to numerical tolerance, and 0/3 canonical endpoints passed the frozen resolution/sharpness retention rule.
- The submission-facing biological conclusion therefore remains the two small primary rows retained by PR #69: orientation–BIO1 and outline aspect ratio–BIO4. The local 17-endpoint run is a diagnostic bridge, not a replacement FDR family.
"""
        ),
        markdown(
            """
## Context & Methods

This notebook joins two frozen streams without changing their scientific roles. The 3,725-observation balanced image atlas supplies already measured continuous orientation, colour/display and outline values. The PR #69 high-resolution table supplies three executed contour formulas, mapped to their sole canonical endpoint identifiers. Missing endpoints remain missing.

For diagnostics, all available primary and candidate endpoints are joined to the CHELSA values present for the 1,292 high-resolution observations. Linear responses use taxon-demeaned standardized slopes with taxon-clustered uncertainty; hue is one joint circular test per predictor. BH correction is separate for primary and candidate tiers. This diagnostic subset lacks the full 46,276-observation seasonal, dominant-taxon and native-range control structure, so it cannot create or replace submission claims.
"""
        ),
        code(
            """
from pathlib import Path
import json
import pandas as pd
import numpy as np

ROOT = Path.cwd()
assert (ROOT / "manuscript" / "final_claims.json").exists(), "Run from the repository root"
LOCAL = ROOT / "analysis_outputs" / "local_pr69_continuous_reanalysis"

with (ROOT / "manuscript" / "final_claims.json").open(encoding="utf-8") as handle:
    final_claims = json.load(handle)
with (LOCAL / "universe_17" / "continuous_trait_universe_report.json").open(encoding="utf-8") as handle:
    universe_report = json.load(handle)
with (LOCAL / "involucre_bridge" / "bridge_report.json").open(encoding="utf-8") as handle:
    bridge_report = json.load(handle)
with (LOCAL / "diagnostic_climate_screen_17_le10km" / "continuous_trait_universe_climate_report.json").open(encoding="utf-8") as handle:
    screen_report = json.load(handle)

traits = pd.read_csv(LOCAL / "universe_17" / "continuous_trait_universe_observation_long.csv", low_memory=False)
coverage = pd.read_csv(LOCAL / "diagnostic_climate_screen_17_le10km" / "continuous_trait_universe_coverage.csv")
coefficients = pd.read_csv(LOCAL / "diagnostic_climate_screen_17_le10km" / "continuous_trait_universe_climate_coefficients.csv")
gate = pd.read_csv(LOCAL / "involucre_bridge" / "canonical_involucre_headline_audit.csv")
gate_comparison = pd.read_csv(LOCAL / "reviewer_gate_comparison" / "candidate_screen_vs_pr69_gates.csv")
contract = pd.read_csv(ROOT / "ch1_global" / "v2" / "ontology" / "ch1_continuous_trait_contract.csv", keep_default_na=False)
"""
        ),
        markdown(
            """
## Data

The hundreds of thousands of photographs were not reduced through one arbitrary convenience sample. There are two analysis streams. The exhaustive stream supports the primary within-taxon climate analysis after coordinate and spatial thinning; the smaller balanced atlas supports repeated head-level measurement, variance decomposition and development/audit. The local bridge uses the latter because the 923.7 MB exhaustive artifact is not available inside this local connector session.
"""
        ),
        code(
            """
datasets = final_claims["datasets"]
attrition = pd.DataFrame([
    {"stream": "exhaustive", "stage": "detector-positive observations", "n_observations": 406_582, "n_taxa": 286, "role": "post-detection source"},
    {"stream": "exhaustive", "stage": "finite coordinates", "n_observations": 392_989, "n_taxa": 271, "role": "coordinate QC"},
    {"stream": "exhaustive", "stage": "positional accuracy ≤10 km", "n_observations": 297_293, "n_taxa": 259, "role": "geographic accuracy QC"},
    {"stream": "exhaustive", "stage": "taxon × 0.25° spatial thinning", "n_observations": 46_276, "n_taxa": 259, "role": "primary climate cohort"},
    {"stream": "balanced atlas", "stage": "one photograph per observation", "n_observations": 3_725, "n_taxa": 216, "role": "measurement/variance/historical audit"},
    {"stream": "high-resolution", "stage": "usable involucre measurements", "n_observations": 1_292, "n_taxa": 210, "role": "candidate contour traits"},
    {"stream": "high-resolution", "stage": "≤10 km complete cases", "n_observations": 904, "n_taxa": 165, "role": "locked PR69 bias-control family"},
])
attrition
"""
        ),
        code(
            """
assert traits.duplicated(["obs_id", "endpoint_id"]).sum() == 0
assert set(traits["endpoint_id"]).issubset(set(contract["endpoint_id"]))
assert universe_report["category_columns_in_output"] == []
assert bridge_report["n_involucre_outside_primary"] == 0
assert bridge_report["n_taxon_name_mismatches"] == 0
assert bridge_report["model_reproduction"]["status"] == "exact_within_tolerance"

endpoint_quality = (
    traits.groupby(["endpoint_id", "module", "analysis_tier"], dropna=False)
    .agg(
        n_rows=("obs_id", "size"),
        n_observations_measured=("measurement_available", "sum"),
        n_taxa_measured=("taxon_name", lambda s: s[traits.loc[s.index, "measurement_available"].astype(bool)].nunique()),
    )
    .reset_index()
)
endpoint_quality["measurement_coverage"] = endpoint_quality["n_observations_measured"] / endpoint_quality["n_rows"]
endpoint_quality.sort_values(["analysis_tier", "module", "endpoint_id"])
"""
        ),
        code(
            """
quality_summary = pd.DataFrame([
    {"check": "contract endpoints", "value": len(contract), "status": "pass"},
    {"check": "executed endpoints in local universe", "value": universe_report["n_available_endpoints"], "status": "partial by design"},
    {"check": "unexecuted endpoints", "value": universe_report["n_missing_unexecuted_endpoints"], "status": "requires image remeasurement"},
    {"check": "duplicate obs_id × endpoint_id", "value": int(traits.duplicated(["obs_id", "endpoint_id"]).sum()), "status": "pass"},
    {"check": "category columns leaked", "value": len(universe_report["category_columns_in_output"]), "status": "pass"},
    {"check": "PR69 coefficient rows reproduced", "value": bridge_report["model_reproduction"]["n_compared_rows"], "status": "pass"},
    {"check": "successful diagnostic models", "value": screen_report["n_successful_endpoint_predictor_models"], "status": "diagnostic only"},
])
quality_summary
"""
        ),
        markdown(
            """
## Results

The diagnostic 17-endpoint universe contains all currently available numerical measurements and no biological categories. Four descriptive floral-pixel composition fractions are retained as numerical description but are not fitted as four independent ecological outcomes. Ten endpoints remain absent because their required high-resolution image functions have not been executed on the frozen exhaustive cohort.
"""
        ),
        code(
            """
model_overview = pd.DataFrame([
    {"metric": "linear endpoints modelled", "value": screen_report["n_linear_endpoints_modelled"]},
    {"metric": "joint circular traits modelled", "value": screen_report["n_circular_traits_modelled"]},
    {"metric": "successful endpoint–predictor models", "value": screen_report["n_successful_endpoint_predictor_models"]},
    {"metric": "primary BH signals", "value": screen_report["n_primary_fdr_signals"]},
    {"metric": "candidate BH signals", "value": screen_report["n_candidate_fdr_signals"]},
])
model_overview
"""
        ),
        code(
            """
display_columns = [
    "endpoint_id", "analysis_tier", "predictor", "inferential_unit",
    "beta_std_within", "effect_magnitude", "p_value", "q_fdr_bh_within_tier",
    "n_observations", "n_taxa",
]
top_models = coefficients.sort_values("q_fdr_bh_within_tier").head(12)
top_models[[column for column in display_columns if column in top_models.columns]]
"""
        ),
        code(
            """
gate_comparison[[
    "endpoint_id", "predictor", "screen_beta", "screen_q",
    "quality_adjusted_beta", "quality_adjusted_q",
    "positive_in_all_successful_strata", "failure_reasons",
    "submission_claim_eligible",
]]
"""
        ),
        code(
            """
gate_display = gate[[
    "endpoint", "adjusted_bio4_beta", "adjusted_bio4_q",
    "positive_in_all_successful_strata", "strict_resolution_control_retained",
    "claim_status", "submission_eligible",
]].copy()
gate_display
"""
        ),
        code(
            """
native_rows = pd.DataFrame(
    final_claims["bias_control_reanalysis"]["native_range_sensitivity"]["rows"]
)
native_rows[[
    "endpoint", "predictor", "native_only_beta", "native_only_q",
    "n_observations", "n_taxa", "native_range_robust",
]]
"""
        ),
        markdown(
            """
## Takeaways

1. **The current submission conclusion does not change.** PR #69 retains orientation–BIO1 and outline aspect ratio–BIO4 only; it withdraws all three high-resolution involucre/BIO4 rows.
2. **The continuous-trait redesign is now testable but incomplete.** Seventeen endpoints are in one category-free long table. Six candidate architecture/armature endpoints and four validation-only surface endpoints still require image execution; missingness is explicit.
3. **The “hundreds of thousands versus a few thousand” issue is a cohort-ledger issue, not simple data waste.** The 46,276-observation exhaustive primary cohort is the main climate-analysis stream. The 3,725-observation atlas and 1,292-observation involucre subset serve different measurement and validation purposes and must not share FDR labels.
4. **A small q value in the generic screen is not a result to publish.** Four candidate rows have q < 0.05 before endpoint-specific quality controls; the locked PR #69 reanalysis reduces the corresponding involucre headline to 0/3 retained endpoints. The diagnostic subset also lacks the full seasonal, dominant-taxon, native-range and developmental-stage controls.
5. **Next executable gate:** obtain the frozen exhaustive artifact in an environment without the 512 MB connector cap, run all 27 endpoint functions, and then pass every promoted endpoint through its matching PR #69 bias/validation gate.
"""
        ),
    ]

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    client = NotebookClient(notebook, timeout=180, kernel_name="python3", resources={"metadata": {"path": str(ROOT)}})
    try:
        executed = client.execute()
        executed["metadata"]["execution"] = {
            "status": "executed_top_to_bottom",
            "executor": "nbclient",
        }
    except RuntimeError as error:
        if "Kernel died before replying to kernel_info" not in str(error):
            raise
        executed = execute_in_process(notebook)
    nbf.write(executed, OUTPUT)
    print(OUTPUT)


if __name__ == "__main__":
    main()
