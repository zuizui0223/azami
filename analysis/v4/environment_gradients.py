"""Conditional changes over observed exposure ranges; not absolute predictions."""
import numpy as np


def exposure_support(rows, predictors):
    support = {}
    for variable in predictors:
        for level in ("between", "within"):
            term = variable + "__" + level
            values = (rows.drop_duplicates("accepted_key")[term]
                      if level == "between" else rows[term])
            values = values[np.isfinite(values)]
            if len(values) == 0:
                continue
            q = np.quantile(values, [.05, .25, .5, .75, .95])
            support[term] = {
                "q05": float(q[0]), "q25": float(q[1]), "median": float(q[2]),
                "q75": float(q[3]), "q95": float(q[4]),
                "n": len(values),
                "weighting": "one_per_species" if level == "between" else "observations",
                "scope": "marginal exposure support, not joint covariate overlap",
            }
            if level == "within":
                ranges = rows.groupby("accepted_key")[term].agg(["min", "max"])
                support[term]["species_spanning_q05_q95"] = int(
                    ((ranges["min"] <= q[0]) & (ranges["max"] >= q[4])).sum())
                support[term]["species_n"] = len(ranges)
    return support


def gradient_points(effect, support, environment_sd, count=51):
    """Slope contrast relative to median; reverse interval endpoints for x<median."""
    if not np.isfinite(environment_sd) or environment_sd <= 0:
        raise ValueError("Positive finite environmental SD required")
    x = np.linspace(support["q05"], support["q95"], count)
    delta = (x - support["median"]) / environment_sd
    lower, upper = effect["raw_ci95"]
    limits = np.array([delta * lower, delta * upper])
    return {
        "x": x.tolist(),
        "change": (delta * effect["raw_units_per_environment_sd"]).tolist(),
        "lower": limits.min(axis=0).tolist(),
        "upper": limits.max(axis=0).tolist(),
        "reference": support["median"],
        "scope": "linear conditional change; pointwise slope uncertainty; not prediction interval",
    }
