"""Synthetic accepted-key/full251/source-read depth accounting and false-phylogeny guards."""
from copy import deepcopy
import pytest
from analysis.v4 import t2_s1_all251_real120k_q20_bridge as m

def fixtures():
    units=[]
    distinct=[str(206745753)] + [str(206750000+i) for i in range(1,244)]
    for i in range(251):
        key=distinct[i] if i<244 else distinct[i-244]
        if i==0: name="Cirsium_brevicaule"
        elif i==244: name="Cirsium_irumtiense"
        else:name=f"Cirsium_test_{i}"
        source_matched=7<=i<120
        units.append({
            "accepted_key":key,"S1_original_source_name":name,
            "original_S1_photo_count":"678" if i<250 else "566",
            "S1_CHELSA_BIO12_mm":str(2000+i*.5),
            "original_WCVP_collapsed_pair":str(i<7 or i>=244),
            "real_241_site_depth_available":str(source_matched),
            "Q20_dp3_ge50pct_loci_best_single_run":"17" if source_matched else "",
            "n_EAzami294_name_candidate_samples":"0",
        })
    assert sum(int(x["original_S1_photo_count"]) for x in units)==170066
    taxa=[]
    runs=[]
    for i in range(7,120):
        source=units[i]
        key=source["accepted_key"]
        shorter=i in (7,8)
        nrun=2 if 9<=i<=20 else 1
        metric=105 if 9<=i<=15 else 71
        taxa.append({
            "accepted_key":key,
            "source_name":source["S1_original_source_name"],
            "best_30k_single_run_dp3_ge50pct":"17",
            "best_120k_single_run_dp3_ge50pct":"" if shorter else str(metric),
            "successfully_depth120k_mapped_runs":"0" if shorter else str(nrun),
        })
        for j in range(nrun):
            status="failed" if shorter else "real_ENA_120k_quality20_depth"
            runs.append({
                "run":f"SRR{i:04d}{j}","biosample":f"SAMN{i:04d}{j}",
                "accepted_key":key,
                "source_taxon":source["S1_original_source_name"],
                "original_first30k_read_pairs":"30000",
                "first30k_n_loci_dp3_ge50pct":"17",
                "full_120k_comparison_status":status,
                "source_30k_to_120k_valid":str(not shorter),
                "first120k_n_loci_dp3_ge50pct":"" if shorter else str(metric),
                "source_run_unavailability_reason":(
                    "ENA source full FASTQ has fewer than 120k paired reads" if shorter else ""),
            })
    assert len(runs)==125
    summary={
        "n_runs_with_recovered_120k_read_depth":123,
        "n_runs_without_complete_120k_prefix":2,
        "n_source_taxa_with_any_120k_test":111,
        "n_source_taxa_with_best_single_120k_dp3_ge50pct_ge100":7,
    }
    return units,taxa,runs,summary

def test_all_original_species_separate_and_7_high_depth_candidates():
    s,t,r,d=fixtures()
    lines,summary=m.merge(s,t,r,d)
    assert len(lines)==251
    assert len({x["accepted_key"] for x in lines})==244
    assert summary["S1_source_taxa_with_at_least100_loci_dp3_half_length"]==7
    assert summary["S1_source_taxa_with_any_120k"]==111
    assert summary["source_complete_exact_run_depth120k"]==123
    assert summary["short_public_SRA_read_archives"]==2
    assert summary["original_photo_depth_sensitivity"]["10"]["n_with_at_least100_loci_dp3_ge50pct"]==7
    assert sum(x["source_same_biosample_run_120k_Q20_successes"] for x in lines)==123
    pair=[x for x in lines if x["accepted_key"]=="206745753"]
    assert {x["S1_original_source_name"] for x in pair}==m.TROPICAL
    assert all(x["Q20_dp3_ge50pct_loci_best_120k_single_source_run"] is None for x in pair)
    # All 138 original names with no matched Moreyra SRA retain a NULL.
    assert sum(x["Q20_dp3_ge50pct_loci_best_120k_single_source_run"] is None for x in lines)==140
    assert all(x["no_nuclear_species_tree_estimated"] is True for x in lines)

def test_30k_source_baseline_read_pairs_cannot_be_substituted():
    s,t,r,d=fixtures()
    s[9]["Q20_dp3_ge50pct_loci_best_single_run"]="19"
    with pytest.raises(ValueError,match="Baseline 30k"):
        m.merge(s,t,r,d)

def test_falsely_imputed_nuclear_genes_on_wcvp_collided_pair_blocked():
    s,t,r,d=fixtures()
    s[0]["real_241_site_depth_available"]="True"
    with pytest.raises(ValueError,match="identity drift"):
        m.merge(s,t,r,d)

def test_same_voucher_nested_120k_q20_locus_calls_cannot_shrink():
    s,t,r,d=fixtures()
    bad=next(z for z in r if z["full_120k_comparison_status"].startswith("real_"))
    bad["first120k_n_loci_dp3_ge50pct"]="16"
    with pytest.raises(ValueError,match="cannot decrease"):
        m.merge(s,t,r,d)

def test_short_source_runs_are_not_zero_gene_coverage():
    s,t,r,d=fixtures()
    lines,_=m.merge(s,t,r,d)
    a=next(z for z in lines if z["S1_original_source_name"]==s[7]["S1_original_source_name"])
    assert a["source_same_biosample_run_120k_Q20_short_archive"]==1
    assert a["Q20_dp3_ge50pct_loci_best_120k_single_source_run"] is None
    assert a["original_30k_best_single_run_loci_Q20"]==17
    assert a["Q20_120k_read_state"]=="no120k_due_to_short_public_fastq"

def test_wrong_species_source_name_does_not_map_molecular_voucher():
    s,t,r,d=fixtures()
    t[3]["source_name"]="Cirsium_incorrect"
    with pytest.raises(ValueError,match="different S1 name"):
        m.merge(s,t,r,d)
