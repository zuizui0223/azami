import json
from pathlib import Path

import pandas as pd

from analysis.v3.prepare_authority_taxonomy_sensitivity import canonical_key
from analysis.v3.summarize_authority_taxonomy_sensitivity import summarize


def test_canonical_key_is_stable_and_rejects_missing():
    assert canonical_key(123) == "wcvp:123"
    assert canonical_key("123.0") == "wcvp:123"
    try:
        canonical_key(None)
    except ValueError:
        pass
    else:
        raise AssertionError("missing accepted key must fail closed")


def _write_inputs(tmp_path: Path, second_q: float = 0.02):
    tmp_path.mkdir(parents=True, exist_ok=True)
    prep = tmp_path / "prep.json"
    prep.write_text(json.dumps({"analysis_id": "prep", "source_taxa_start": 259}), encoding="utf-8")

    authority = pd.DataFrame([
        {"construct_id": "floral_chroma", "predictor": "chelsa_rsds_mean", "kind": "scalar", "n_taxa": 100, "beta_std": -0.4, "p_value": 0.001, "q_bh": 0.01},
        {"construct_id": "presentation_angle", "predictor": "chelsa_bio12", "kind": "scalar", "n_taxa": 99, "beta_std": 0.3, "p_value": 0.003, "q_bh": second_q},
    ])
    frozen = pd.DataFrame([
        {"construct_id": "floral_chroma", "predictor": "chelsa_rsds_mean", "kind": "scalar", "n_taxa": 143, "beta_std": -0.345, "p_value": 0.001, "q_bh": 0.002},
        {"construct_id": "presentation_angle", "predictor": "chelsa_bio12", "kind": "scalar", "n_taxa": 142, "beta_std": 0.304, "p_value": 0.002, "q_bh": 0.004},
    ])
    authority_path = tmp_path / "authority.csv"
    frozen_path = tmp_path / "frozen.csv"
    authority.to_csv(authority_path, index=False)
    frozen.to_csv(frozen_path, index=False)

    upgrade = tmp_path / "upgrade.json"
    upgrade.write_text(json.dumps({
        "common_cohort": {"observations": 1500, "taxa": 39, "minimum_complete_observations_per_taxon": 5},
        "common_cohort_matrix_alignment": {"rho": 0.4, "qap_p_one_sided": 0.01},
        "taxon_bootstrap": {"rho_median": 0.3, "rho_low95": 0.02, "rho_high95": 0.5},
    }), encoding="utf-8")

    contrast = tmp_path / "contrast.json"
    contrast.write_text(json.dumps({
        "observed": {"relations": 36, "median_within_rv": 0.01, "median_among_rv": 0.05, "relations_stronger_among": 30},
        "taxon_bootstrap": {
            "difference_of_median_rv_median": 0.05,
            "difference_of_median_rv_low95": 0.02,
            "difference_of_median_rv_high95": 0.10,
            "probability_median_among_exceeds_within": 1.0,
        },
    }), encoding="utf-8")
    return prep, authority_path, upgrade, contrast, frozen_path


def test_predeclared_taxonomy_gates_require_geometry_strength_and_both_anchors(tmp_path):
    paths = _write_inputs(tmp_path)
    out = tmp_path / "summary.json"
    result = summarize(*paths, out)
    assert result["authority_geometry"]["gate_pass"] is True
    assert result["authority_strength"]["gate_pass"] is True
    assert result["anchor_gate_pass"] is True
    assert result["headline_taxonomy_robust"] is True

    paths = _write_inputs(tmp_path / "fail", second_q=0.20)
    out = tmp_path / "fail" / "summary.json"
    result = summarize(*paths, out)
    assert result["anchor_gate_pass"] is False
    assert result["headline_taxonomy_robust"] is False
