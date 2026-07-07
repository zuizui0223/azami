#!/usr/bin/env python3
"""Report-safe entry point for the hardened Ch.1 trait–environment atlas.

The original report helper attempted to read ``eligibility_status`` from the
model-status table, although eligibility is stored in a separate feasibility
table. This wrapper keeps all numerical fits and grouped validation unchanged
and replaces only report generation with a status-table-compatible summary.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd


HARDENED_PATH = Path(__file__).with_name("44_run_trait_environment_association_atlas_hardened.py")
spec = importlib.util.spec_from_file_location("ch1_trait_environment_hardened", HARDENED_PATH)
if spec is None or spec.loader is None:
    raise RuntimeError(f"Unable to load hardened association atlas: {HARDENED_PATH}")
hardened = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hardened)
base = hardened.base


def outcome_summary_markdown(
    statuses: pd.DataFrame,
    coefficients: pd.DataFrame,
    validation: pd.DataFrame,
    predictors: list[str],
    plan: dict,
) -> str:
    """Generate a report using columns actually present in model-status output."""
    lines = [
        "# Ch.1 trait–environment association atlas",
        "",
        str(plan["semantic_status"]),
        "",
        "## Inputs",
        f"- Predictors: {', '.join(predictors)}",
        "- Measurement modes: all-ensemble winner and conservative strict-consensus.",
        "- Validation: grouped 10-degree spatial blocks, plus held-out species for climate-only models.",
        "- No random-split headline metrics, p-values, phylogenetic claims, or causal interpretation are included.",
        "",
        "## Model-state outcomes with successful climate-only fit",
    ]
    if statuses.empty:
        lines.append("No model-status rows were generated.")
    else:
        successful = statuses.loc[
            statuses["model_family"].eq("climate_only") & statuses["status"].eq("success")
        ].copy()
        if successful.empty:
            lines.append("No climate-only state outcome completed a full-data fit.")
        else:
            display_columns = [column for column in ["measurement_mode", "trait_id", "trait_state", "trait_tier", "status"] if column in successful.columns]
            lines.append(successful[display_columns].drop_duplicates().sort_values(display_columns[:-1]).to_markdown(index=False))

    lines.extend(["", "## Grouped validation summary", ""])
    if validation.empty:
        lines.append("No validation rows were generated.")
    else:
        successful_validation = validation.loc[validation["status"].eq("success")].copy()
        if successful_validation.empty:
            lines.append("No grouped-validation fold was jointly scorable for any fitted state outcome.")
        else:
            grouping = ["measurement_mode", "trait_id", "trait_state", "trait_tier", "model_family", "validation_scheme"]
            grouping = [column for column in grouping if column in successful_validation.columns]
            aggregate = successful_validation.groupby(grouping, as_index=False).agg(
                n_scorable_folds=("fold", "count"),
                mean_roc_auc=("roc_auc", "mean"),
                mean_balanced_accuracy=("balanced_accuracy", "mean"),
            )
            lines.append(aggregate.to_markdown(index=False, floatfmt=".3f"))

    lines.extend(["", "## Coefficient output", ""])
    lines.append(f"- Standardized penalized coefficients are provided for {len(coefficients)} predictor/outcome/model rows.")
    lines.append("- Coefficient sign and magnitude are descriptive associations conditional on the stated model family; they are not adaptive-effect estimates.")
    return "\n".join(lines) + "\n"


base.outcome_summary_markdown = outcome_summary_markdown

if __name__ == "__main__":
    base.main()
