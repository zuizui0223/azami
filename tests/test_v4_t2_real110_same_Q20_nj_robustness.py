"""Finite DNA source checks for independent 110-tip Q20/N neighbour joining."""
import json

import pytest
from Bio import Phylo

from analysis.v4 import t2_real110_same_Q20_nj_robustness as mod


def four():
    return {
        "WCVP100000001":"AAAAAACCCCCC",
        "WCVP100000002":"AAAAAACCCCCN",
        "WCVP100000003":"GGGGGGTTTTTT",
        "WCVP100000004":"GGGGGGTTTTTN",
    }


def test_unobserved_q20_nucleotides_never_filled_with_reference():
    ids,matrix,pairs=mod.pairwise_observed_pdist(four(),threshold=10)
    assert ids==sorted(four())
    assert len(pairs)==6 and len(matrix)==4
    target=next(x for x in pairs if x["tip_A"]=="WCVP100000001" and
                x["tip_B"]=="WCVP100000002")
    assert target["n_shared_Q20_positions"]==11
    assert target["n_differing_Q20_positions"]==0
    assert matrix[-1][-1]==0.0


def test_fail_closed_when_pairwise_overlap_has_only_missing_N():
    seq={"WCVP1":"AAAAAANNNNNN","WCVP2":"NNNNNNCCCCCC"}
    with pytest.raises(ValueError,match="Too few joint observed"):
        mod.pairwise_observed_pdist(seq,threshold=1)


def test_true_four_source_NJ_has_no_added_source_tip(monkeypatch):
    monkeypatch.setattr(mod,"N_TIPS",4)
    # The real-world default 500 jointly observed Q20 bases is fixed at
    # function-definition time, so make this synthetic example genuinely long.
    tree,dist,summary=mod.audit_nj({k:v*50 for k,v in four().items()},{})
    assert {x.name for x in tree.get_terminals()}==set(four())
    assert summary["n_original_public_BioSample_pairs"]==6
    assert not summary["nuclear_copy_number_single_orthology_verified"]
    assert not summary["nodding_ancestral_gains_or_rain_adaptation_inferred"]


def test_original111_failed_gate_source_cannot_be_promoted(tmp_path):
    path=tmp_path/"preflight.json"
    path.write_text(json.dumps({
        "posthoc_sensitivity_taxa":110,
        "original_all111_Q20_successful_real_SRA_BioSamples":111,
        "original111_technical_gate_passed":True,
    }))
    with pytest.raises(ValueError,match="failed gate"):
        mod.check_manifest(path)
