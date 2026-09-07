from analysis.v3.resolution_size_followup import exact_two_sided_binomial_half, summarize_endpoint


def row(bstatus, cstatus, b=None, c=None):
    return {
        "baseline_status": bstatus,
        "comparison_status": cstatus,
        "baseline_value": b,
        "comparison_value": c,
    }


def test_binomial_discordance_is_two_sided():
    assert exact_two_sided_binomial_half(15, 0) <= 0.01
    assert exact_two_sided_binomial_half(8, 7) > 0.01
    assert exact_two_sided_binomial_half(0, 0) is None


def test_original_eligibility_gain_requires_predeclared_discordance_evidence():
    rows = [row("not_usable", "usable") for _ in range(15)]
    result = summarize_endpoint(rows, min_common=30, min_discordant=15)
    assert result["eligibility_decidable"] is True
    assert result["numeric_decidable"] is False
    assert result["status"] == "ORIGINAL_MAY_ADD_INFORMATION_BEYOND_LARGE"
    assert "original_eligibility_gain_over_large" in result["flags"]


def test_stable_common_pairs_can_clear_numeric_gate():
    rows = [row("usable", "usable", float(i), float(i) + 0.01) for i in range(40)]
    result = summarize_endpoint(rows, min_common=30, min_discordant=15)
    assert result["numeric_decidable"] is True
    assert result["rank_correlation"] == 1.0
    assert result["status"] == "NO_FOLLOWUP_EVIDENCE_ORIGINAL_REQUIRED_BEYOND_LARGE"


def test_insufficient_followup_information_is_not_promoted_to_stability():
    rows = [row("usable", "usable", float(i), float(i)) for i in range(8)]
    result = summarize_endpoint(rows, min_common=30, min_discordant=15)
    assert result["numeric_decidable"] is False
    assert result["eligibility_decidable"] is False
    assert result["status"] == "INSUFFICIENT_FOLLOWUP_INFORMATION"
