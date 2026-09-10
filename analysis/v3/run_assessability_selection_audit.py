#!/usr/bin/env python3
"""Outcome-blind assessability/environment audit for Chapter 1 v3.

The response is whether a frozen v2 endpoint (or all members of a biological
construct) was measurable in an observation. Trait values are never used as the
response. Environment is taxon-centred, so the fitted slope asks whether
measurement opportunity changes along an environmental gradient within taxa.

This is a retrospective selection-bias diagnostic. It does not replace the
frozen v2 endpoint atlas, change its multiplicity family, or prove that a small
availability gradient cannot cause value-dependent selection.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

from analysis.v3 import run_biological_axis_reanalysis as axis


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--traits", type=Path, required=True)
    p.add_argument("--environment", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    return p.parse_args()


def as_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def bh(values: pd.Series) -> np.ndarray:
    p = np.asarray(values, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.minimum(q, 1.0)
    out = np.empty(n, dtype=float)
    out[order] = q
    return out


def centred_predictors(env: pd.DataFrame) -> pd.DataFrame:
    out = env[["obs_id", "taxon_name", *axis.PREDICTORS]].copy()
    for predictor in axis.PREDICTORS:
        x = pd.to_numeric(out[predictor], errors="coerce")
        xc = x - x.groupby(out.taxon_name).transform("mean")
        sd = float(xc.std(ddof=0))
        if not np.isfinite(sd) or sd <= 0:
            raise ValueError(f"no within-taxon environmental variation: {predictor}")
        out[f"zwithin__{predictor}"] = xc / sd
    return out


def cluster_lpm(frame: pd.DataFrame, response: str, predictor: str) -> dict[str, float | int | str]:
    y = pd.to_numeric(frame[response], errors="coerce")
    x = pd.to_numeric(frame[f"zwithin__{predictor}"], errors="coerce")
    keep = y.notna() & x.notna()
    y = y[keep].to_numpy(float)
    x = x[keep].to_numpy(float)
    taxa = frame.loc[keep, "taxon_name"].astype(str).to_numpy()
    if len(y) < 100 or len(np.unique(taxa)) < 10:
        return {"status": "insufficient_support"}
    yc = y - pd.Series(y).groupby(taxa).transform("mean").to_numpy(float)
    den = float(np.dot(x, x))
    if den <= 1e-15 or float(np.std(yc)) <= 1e-15:
        return {"status": "constant_within_taxon_response"}
    beta = float(np.dot(x, yc) / den)
    resid = yc - beta * x
    scores = pd.DataFrame({"taxon": taxa, "score": x * resid}).groupby("taxon").score.sum().to_numpy(float)
    g = len(scores)
    n = len(y)
    correction = (g / (g - 1)) * ((n - 1) / max(n - 1, 1)) if g > 1 else 1.0
    var = correction * float(np.dot(scores, scores)) / (den * den)
    se = float(np.sqrt(max(var, 0.0)))
    z = beta / se if se > 0 else np.inf
    p = float(2.0 * norm.sf(abs(z)))
    taxon_means = pd.DataFrame({"taxon": taxa, "available": y}).groupby("taxon").available.mean()
    return {
        "status": "ok",
        "n_observations": int(n),
        "n_taxa": int(g),
        "availability_fraction_observation_weighted": float(np.mean(y)),
        "availability_fraction_equal_taxon": float(taxon_means.mean()),
        "beta_probability_per_within_taxon_sd": beta,
        "beta_percentage_points_per_within_taxon_sd": beta * 100.0,
        "cluster_robust_se": se,
        "z_value": float(z),
        "p_value": p,
    }


def wide_availability(traits: pd.DataFrame) -> pd.DataFrame:
    t = traits[["obs_id", "taxon_name", "endpoint_id", "measurement_available"]].copy()
    t["measurement_available"] = as_bool(t.measurement_available)
    dup = t.duplicated(["obs_id", "taxon_name", "endpoint_id"]).any()
    if dup:
        raise ValueError("trait universe is not one row per observation x endpoint")
    wide = t.pivot(index=["obs_id", "taxon_name"], columns="endpoint_id", values="measurement_available").reset_index()
    wide.columns.name = None
    return wide


def run_family(base: pd.DataFrame, responses: dict[str, pd.Series], id_name: str) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for item_id, response in responses.items():
        frame = base[["obs_id", "taxon_name", *[f"zwithin__{p}" for p in axis.PREDICTORS]]].copy()
        frame["available"] = response.to_numpy(bool).astype(float)
        for predictor in axis.PREDICTORS:
            result = cluster_lpm(frame, "available", predictor)
            rows.append({id_name: item_id, "predictor": predictor, **result})
    out = pd.DataFrame(rows)
    ok = out.status.eq("ok")
    out["q_bh"] = np.nan
    if ok.any():
        out.loc[ok, "q_bh"] = bh(out.loc[ok, "p_value"])
    out["fdr_0_05"] = out.q_bh < 0.05
    return out


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    traits = pd.read_csv(args.traits, low_memory=False)
    env = pd.read_csv(args.environment, low_memory=False)
    traits["obs_id"] = traits.obs_id.astype(str)
    traits["taxon_name"] = traits.taxon_name.astype(str)
    env["obs_id"] = env.obs_id.astype(str)
    env["taxon_name"] = env.taxon_name.astype(str)
    traits = traits[~traits.endpoint_id.isin(axis.RESTORED5)].copy()
    endpoints = sorted(traits.endpoint_id.unique())
    expected = sorted({m for d in axis.CONSTRUCTS.values() for m in d["members"]})
    if endpoints != expected:
        raise ValueError(f"frozen 22-endpoint universe mismatch: {len(endpoints)} vs {len(expected)}")
    if traits.groupby("endpoint_id").size().nunique() != 1:
        raise ValueError("endpoint rows do not share the same observation denominator")

    envc = centred_predictors(env)
    wide = wide_availability(traits)
    base = envc.merge(wide, on=["obs_id", "taxon_name"], how="inner", validate="one_to_one")
    if len(base) != len(env):
        raise ValueError(f"observation denominator mismatch: {len(base)} vs {len(env)}")

    endpoint_responses = {endpoint: base[endpoint].astype(bool) for endpoint in endpoints}
    endpoint_results = run_family(base, endpoint_responses, "endpoint_id")

    construct_responses: dict[str, pd.Series] = {}
    for construct, definition in axis.CONSTRUCTS.items():
        members = definition["members"]
        construct_responses[construct] = base[members].astype(bool).all(axis=1)
    construct_results = run_family(base, construct_responses, "construct_id")

    reasons = traits.copy()
    reasons["measurement_available"] = as_bool(reasons.measurement_available)
    reasons["analysis_eligible"] = as_bool(reasons.analysis_eligible)
    reasons["exclusion_reason"] = reasons.exclusion_reason.fillna("").replace("", "included")
    exclusion = (
        reasons.groupby(["endpoint_id", "exclusion_reason"], dropna=False)
        .size().rename("n_rows").reset_index()
    )
    totals = exclusion.groupby("endpoint_id").n_rows.transform("sum")
    exclusion["fraction_of_endpoint_denominator"] = exclusion.n_rows / totals

    anchors = [
        ("corolla_lab_chroma", "floral_chroma", "chelsa_rsds_mean", -0.3453720170895292),
        ("orientation_image_vertical_angle", "presentation_angle", "chelsa_bio12", 0.30435928589775146),
    ]
    anchor_rows = []
    for endpoint, construct, predictor, trait_beta in anchors:
        er = endpoint_results[(endpoint_results.endpoint_id == endpoint) & (endpoint_results.predictor == predictor)].iloc[0]
        cr = construct_results[(construct_results.construct_id == construct) & (construct_results.predictor == predictor)].iloc[0]
        anchor_rows.append({
            "endpoint_id": endpoint,
            "construct_id": construct,
            "predictor": predictor,
            "frozen_trait_beta_reference": trait_beta,
            "endpoint_availability_fraction": float(er.availability_fraction_observation_weighted),
            "endpoint_assessability_beta_pp_per_within_taxon_sd": float(er.beta_percentage_points_per_within_taxon_sd),
            "endpoint_assessability_p": float(er.p_value),
            "endpoint_assessability_q": float(er.q_bh),
            "construct_assessability_beta_pp_per_within_taxon_sd": float(cr.beta_percentage_points_per_within_taxon_sd),
            "same_numeric_direction_as_trait_beta": bool(np.sign(er.beta_probability_per_within_taxon_sd) == np.sign(trait_beta)),
            "interpretation_boundary": "availability slope and trait-value slope are different estimands; same sign does not establish selection-induced mimicry",
        })
    anchor_df = pd.DataFrame(anchor_rows)

    endpoint_results.to_csv(args.out_dir / "endpoint_assessability_environment.csv", index=False)
    construct_results.to_csv(args.out_dir / "construct_assessability_environment.csv", index=False)
    exclusion.to_csv(args.out_dir / "endpoint_exclusion_reason_counts.csv", index=False)
    anchor_df.to_csv(args.out_dir / "headline_anchor_assessability.csv", index=False)

    ok_e = endpoint_results[endpoint_results.status.eq("ok")]
    ok_c = construct_results[construct_results.status.eq("ok")]
    report = {
        "analysis_id": "ch1_v3_assessability_selection_audit_20260910",
        "response": "measurement availability only; trait values are not used as outcomes",
        "model": "taxon-centred linear probability slope with taxon-cluster sandwich SE",
        "observations": int(len(base)),
        "taxa": int(base.taxon_name.nunique()),
        "endpoints": len(endpoints),
        "endpoint_tests": int(len(ok_e)),
        "endpoint_fdr_rows": int(ok_e.fdr_0_05.sum()),
        "constructs": len(construct_responses),
        "construct_tests": int(len(ok_c)),
        "construct_fdr_rows": int(ok_c.fdr_0_05.sum()),
        "largest_absolute_endpoint_gradient_pp_per_sd": float(ok_e.beta_percentage_points_per_within_taxon_sd.abs().max()),
        "largest_absolute_construct_gradient_pp_per_sd": float(ok_c.beta_percentage_points_per_within_taxon_sd.abs().max()),
        "headline_anchors": anchor_df.to_dict("records"),
        "claim_boundary": "outcome-blind post-hoc selection diagnostic; does not rule out trait-value-dependent missingness or replace frozen v2 inference",
    }
    (args.out_dir / "assessability_selection_report.json").write_text(json.dumps(report, indent=2, allow_nan=False), encoding="utf-8")
    print(json.dumps(report, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
