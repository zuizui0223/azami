"""Independent score, phylogeny and inference-stop tests for V4 T2/T3."""
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from analysis.v4.rain_history import (
    EvidenceGateError, audited_tips, conditional_poisson_events, execute,
    pruning_likelihood, rate_fit, score_botanical_records, stochastic_maps,
    validate_tree,
)


def row(record, taxon, concept, who, state, role="primary", stage="anthesis"):
    return dict(record_id=record, accepted_taxon=taxon, concept_id=concept,
                source_taxon=taxon, source_id="flora_2020",
                source_locator="page 4 line 2",
                evidence_quote="Head nodding at flowering.",
                head_stage=stage, annotator_id=who, annotation_role=role,
                blind_climate="yes", blind_phylogeny="yes", state=state)


def test_two_blind_scoring_and_adjudication():
    records = [
        row("r1", "Cirsium a", "a", "botanist1", "erect"),
        row("r1", "Cirsium a", "a", "botanist2", "erect"),
        row("r2", "Cirsium b", "b", "botanist1", "erect"),
        row("r2", "Cirsium b", "b", "botanist2", "nodding"),
        row("r3", "Cirsium c", "c", "botanist1", "variable"),
        row("r3", "Cirsium c", "c", "botanist2", "variable"),
    ]
    decisions, audits = score_botanical_records(records)
    assert decisions["a"]["state"] == "erect"
    assert decisions["b"]["state"] is None
    assert decisions["c"]["status"] == "variable_nonbinary"
    assert audits[1]["status"] == "awaiting_adjudication"
    records.append(row("r2", "Cirsium b", "b", "botanist3", "nodding", "adjudicator"))
    decisions, audits = score_botanical_records(records)
    assert decisions["b"]["state"] == "nodding"
    assert audits[1]["status"] == "third_adjudication"


def test_conflicting_sources_not_fabricated_as_variable():
    rows = [row("r1", "Cirsium a", "a", "1", "erect"),
            row("r1", "Cirsium a", "a", "2", "erect"),
            row("r2", "Cirsium a", "a", "1", "nodding"),
            row("r2", "Cirsium a", "a", "2", "nodding")]
    for r in rows[2:]:
        r["source_id"] = "flora_2022"
    states, _ = score_botanical_records(rows)
    assert states["a"]["state"] is None
    assert states["a"]["status"] == "conflicting_botanical_sources"


@pytest.mark.parametrize("change", [
    {"blind_climate": "no"}, {"blind_phylogeny": "no"},
    {"source_locator": ""}, {"evidence_quote": ""},
    {"state": "down"}, {"head_stage": "not_declared"},
])
def test_unverifiable_botanical_score_is_rejected(change):
    records = [row("r1", "Cirsium a", "a", "1", "erect"),
               row("r1", "Cirsium a", "a", "2", "erect")]
    records[0].update(change)
    with pytest.raises(EvidenceGateError):
        score_botanical_records(records)


def test_scoring_does_not_accept_climate_column():
    records = [row("r1", "Cirsium a", "a", "1", "erect"),
               row("r1", "Cirsium a", "a", "2", "erect")]
    for r in records:
        r["chelsa_bio12"] = "2000"
    with pytest.raises(EvidenceGateError, match="must not contain climate"):
        score_botanical_records(records)


TREE = "((a:0.2,b:0.2):0.3,((c:0.3,d:0.3):0.2,(e:0.2,f:0.2):0.3):0.1);\n"


def write_tree(tmp, kind="independent_nuclear_species_tree", unit="substitutions_per_site"):
    tree = tmp / "tree.nwk"
    tree.write_text(TREE, encoding="utf-8")
    manifest = tmp / "tree.json"
    data = {"tree_sha256": hashlib.sha256(tree.read_bytes()).hexdigest(),
            "source_doi": "10.1016/j.ympev.2025.108285",
            "tree_type": kind,
            "root_method": "externally justified outgroup separation",
            "branch_unit": unit,
            "review_status": "independently_reviewed"}
    manifest.write_text(json.dumps(data), encoding="utf-8")
    return tree, manifest, data


def test_unvalidated_grafting_and_fake_dates_fail_closed(tmp_path):
    tree, manifest, data = write_tree(tmp_path, kind="grafted_placement_sensitivity")
    with pytest.raises(EvidenceGateError, match="independently reviewed"):
        validate_tree(tree, manifest)
    data["tree_type"] = "independent_nuclear_species_tree"
    data["branch_unit"] = "million_years"
    manifest.write_text(json.dumps(data), encoding="utf-8")
    with pytest.raises(EvidenceGateError, match="calibration"):
        validate_tree(tree, manifest)
    data["branch_unit"] = "substitutions_per_site"
    data["tree_sha256"] = "0" * 64
    manifest.write_text(json.dumps(data), encoding="utf-8")
    with pytest.raises(EvidenceGateError, match="SHA256"):
        validate_tree(tree, manifest)


def test_symmetric_Mk_scaled_posteriors_and_bridge():
    from Bio import Phylo
    from io import StringIO
    tree = Phylo.read(StringIO(TREE), "newick")
    observed = {"a": 0, "b": 0, "c": 1, "d": 1, "e": 0, "f": 1}
    rate, ll, at_boundary = rate_fit(tree, observed)
    assert np.isfinite(rate) and rate > 0
    assert np.isfinite(ll)
    assert pruning_likelihood(tree, observed, rate)[0] == pytest.approx(ll)
    samples = stochastic_maps(tree, observed, rate, draws=200, seed=7)
    assert samples["erect_to_nodding_q025_median_q975"][0] >= 0
    assert samples["nodding_to_erect_q025_median_q975"][-1] >= 0
    assert 0 <= samples["root_nodding_p"] <= 1
    assert conditional_poisson_events(1e-12, True, np.random.default_rng(1)) == 1
    assert conditional_poisson_events(1e-12, False, np.random.default_rng(1)) == 0
    for parity in (True, False):
        for _ in range(15):
            n = conditional_poisson_events(0.5, parity, np.random.default_rng(_))
            assert (n % 2 == 1) == parity


def _write_csv(path, rows, fields):
    with path.open("w", newline="", encoding="utf-8") as fd:
        writer = csv.DictWriter(fd, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def inputs(tmp):
    tree, manifest, _ = write_tree(tmp)
    states = {"a": "erect", "b": "erect", "c": "nodding",
              "d": "nodding", "e": "erect", "f": "nodding"}
    rows = []
    for concept, state in states.items():
        for annotator in ("botanist1", "botanist2"):
            rows.append(row(concept, "Cirsium " + concept, concept, annotator, state))
    scores = tmp / "scores.csv"
    _write_csv(scores, rows, list(rows[0]))
    crosswalk_rows = []
    for concept in states:
        crosswalk_rows.append({
            "tip": concept, "concept_id": concept,
            "accepted_taxon": "Cirsium " + concept, "decision": "verified_exact",
            "tribe": "Cardueae" if concept in "abc" else "Senecioneae",
            "phenotype_reference": "botanical_head_posture", "head_stage": "anthesis",
            "homology_reviewed": "yes",
        })
    crosswalk = tmp / "crosswalk.csv"
    _write_csv(crosswalk, crosswalk_rows, list(crosswalk_rows[0]))
    climate_rows = []
    for concept in states:
        climate_rows.append({"concept_id": concept, "chelsa_bio12": "900",
                             "chelsa_bio01": "10", "elevation_m": "500",
                             "region": "temperate", "sampling_frame": "source_mean"})
    climate = tmp / "climate.csv"
    _write_csv(climate, climate_rows, list(climate_rows[0]))
    return scores, crosswalk, tree, manifest, climate


def test_end_to_end_T2_and_T3_are_conditional_not_adaptation(tmp_path):
    scores, crosswalk, tree, manifest, climate = inputs(tmp_path)
    for stage in ("T2", "T3"):
        output = tmp_path / f"{stage}.json"
        result = execute(scores, crosswalk, tree, manifest, output, stage=stage,
                         climate_path=climate, draws=100, min_scored=4)
        assert result["status"] == "conditional_trait_only_Mk_pilot"
        assert result["scored_taxa"] == 6
        assert result["climate_coverage"]["complete_climate_and_elevation"] == 6
        assert result["convergent_rain_adaptation_established"] is False
        assert result["transition_ages_estimated"] is False
        assert "rain_vs_temperature_evolutionary_association_fitted" in result
        assert json.loads(output.read_text()) == result
        assert result["scored_tip_fraction"] == 1
        assert result["unscored_tips_excluded_from_Mk_pilot"] == 0
        if stage == "T3":
            assert set(result["transitions"]["within_tribe_gain_distribution"]) == {"Cardueae", "Senecioneae"}
            p = result["transitions"]["p_gains_in_two_or_more_tribes_conditional_Mk"]
            assert 0 <= p <= 1


def test_T3_requires_more_than_one_reviewed_tribe_and_homology(tmp_path):
    scores, crosswalk, tree, manifest, _ = inputs(tmp_path)
    rows = list(csv.DictReader(crosswalk.open()))
    for row in rows:
        row["tribe"] = "Cardueae"
    _write_csv(crosswalk, rows, list(rows[0]))
    result = execute(scores, crosswalk, tree, manifest, tmp_path / "stopped.json",
                     stage="T3", draws=10, min_scored=4)
    assert result["status"] == "blocked_no_historical_inference"
    assert "two tribes" in result["stop_reason"]
    rows[0]["homology_reviewed"] = "no"
    _write_csv(crosswalk, rows, list(rows[0]))
    result = execute(scores, crosswalk, tree, manifest, tmp_path / "blocked.json",
                     stage="T3", draws=10, min_scored=4)
    assert result["mapping_exclusions"]["unreviewed_cross_tribe_homology"] == 1


def test_map_collisions_cannot_be_silently_pooled(tmp_path):
    scores, crosswalk, tree_path, manifest, _ = inputs(tmp_path)
    from Bio import Phylo
    tree, _ = validate_tree(tree_path, manifest)
    src = list(csv.DictReader(scores.open()))
    states, _ = score_botanical_records(src)
    rows = list(csv.DictReader(crosswalk.open()))
    rows[1]["concept_id"] = "a"
    with pytest.raises(EvidenceGateError, match="one-to-one"):
        audited_tips(tree, rows, states, "T2")

def test_unscored_tips_are_excluded_not_given_free_stochastic_transitions(tmp_path):
    scores, crosswalk, tree, manifest, _ = inputs(tmp_path)
    rows = list(csv.DictReader(scores.open()))
    for row in rows:
        if row["concept_id"] == "e":
            row["state"] = "unknown"
    _write_csv(scores, rows, list(rows[0]))
    result = execute(scores, crosswalk, tree, manifest, tmp_path / "reduced.json",
                     stage="T2", min_scored=4, draws=100)
    assert result["status"] == "conditional_trait_only_Mk_pilot"
    assert result["scored_taxa"] == 5
    assert result["unscored_tips_excluded_from_Mk_pilot"] == 1
    assert result["scored_tip_fraction"] == pytest.approx(5/6)
