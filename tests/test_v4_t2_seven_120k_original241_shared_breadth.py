"""Same-source seven 120k high-depth shared-target numerical source gates."""
from copy import deepcopy
import pytest

from analysis.v4 import t2_seven_120k_original241_shared_breadth as m


def fixture():
    source=[];selected=[];paired=[];results=[]
    taxa=sorted(m.TAXA)
    genes=[f"L{j:03d}" for j in range(241)]
    for i,name in enumerate(taxa):
        run=f"SRR{20000+i}"
        biosample=f"SAMN{40000+i}"
        key=str(206700000+i)
        source.append(dict(run=run,biosample=biosample,accepted_key=key,
            S1_source_taxon=name,original_tree_code=name.replace("_"," ")))
        # shared 100 full half breadth plus five loci on six of seven;
        # a single non-consistent original ref copy for run index 0.
        locus=[]
        for j,gene in enumerate(genes):
            if j<100: frac=.65
            elif j<105 and i!=0: frac=.65
            else: frac=0.
            chosen_ref="saff" if not(i==0 and j==12) else "lett"
            locus.append({
              "locus":gene,"best_original_reference":f"{chosen_ref}-{gene}",
              "reference_choice_maximizes_coverage_not_orthology":True,
              "dp1_breadth":max(frac,.10),
              "dp3_breadth":frac,
              "dp5_breadth":0.2 if frac>.5 else 0.,
            })
        size=sum(l["dp3_breadth"]>=.5 for l in locus)
        selected.append(dict(run=run,biosample=biosample,accepted_key=key,
            source_taxon=name,same_run_120k_dp3_ge50pct=str(size)))
        paired.append(dict(run=run,biosample=biosample,accepted_key=key,
            source_30k_to_120k_valid="True",
            first120k_n_loci_dp3_ge50pct=str(size),
            first120k_n_loci_dp3_ge25pct=str(size)))
        results.append(dict(run=run,biosample=biosample,read_pairs=120000,
            n_loci_dp3_ge25pct=size,n_loci_dp3_ge50pct=size,loci=locus))
    return source,selected,paired,genes,results


def test_all_seven_source_vouchers_and_actual_joint_q20_breadth():
    src,chosen,old,genes,reads=fixture()
    approved=m.validate_sources(src,chosen,old)
    taxa,shared,report=m.combine(approved,genes,reads)
    assert len(approved)==7 and len(taxa)==7 and len(shared)==241
    assert report["n_original_target_loci_seven_taxa_ge25pct_Q20"]==100
    assert report["n_original_target_loci_seven_taxa_ge50pct_Q20"]==100
    assert report["n_original_target_loci_seven_taxa_ge25pct_Q20_same_single_ref"]==99
    assert report["n_original_target_loci_at_least6of7_ge50pct_Q20"]==105
    assert not report["actual_allelic_genotypes_called"]
    assert not report["nuclear_species_tree_constructed"]
    assert report["same_reference_is_not_single_copy_paralog_audit"]


def test_voucher_changed_between_two_read_depth_experiments_rejected():
    src,chosen,old,genes,reads=fixture()
    chosen[0]["biosample"]="SAMN_INVALID"
    with pytest.raises(ValueError,match="BioSample"):
        m.validate_sources(src,chosen,old)


def test_seven_source_120k_depth_threshold_cannot_be_relaxed():
    src,chosen,old,genes,reads=fixture()
    chosen[0]["same_run_120k_dp3_ge50pct"]="90"
    with pytest.raises(ValueError,match="selection"):
        m.validate_sources(src,chosen,old)


def test_sample_120k_original_241_index_totals_must_agree_with_independent_run():
    src,chosen,old,genes,reads=fixture()
    a=m.validate_sources(src,chosen,old)
    reads[0]["n_loci_dp3_ge50pct"]=0
    with pytest.raises(ValueError,match="differs"):
        m.combine(a,genes,reads)


def test_no_missing_taxa_and_no_different_reference_copy_declared_same():
    src,chosen,old,genes,reads=fixture()
    a=m.validate_sources(src,chosen,old)
    with pytest.raises(ValueError,match="seven"):
        m.combine(a,genes,reads[:-1])
    # Never conflate same target identifier with identical original ref copy.
    _,matrix,_=m.combine(a,genes,reads)
    assert matrix[12]["best_reference_copy_identical_all7"] is False


def test_actual_run_depth_must_have_all_241_named_loci():
    src,chosen,old,genes,reads=fixture()
    a=m.validate_sources(src,chosen,old)
    reads[0]["loci"].pop()
    with pytest.raises(ValueError,match="241"):
        m.combine(a,genes,reads)
