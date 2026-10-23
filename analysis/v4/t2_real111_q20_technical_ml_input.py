"""Technical 111-tip source Q20 FASTA validation for optional IQ-TREE ML.

Only starts from successful REAL 111-source, 118 nuclear reference-copy site
calls with observed nucleotide alleles (>=3Q20, >=90% majority, no reference
imputation). This is a reference-anchored, candidate-locus technical ML
screen; no copy-number or independent orthology, no accepted genus-wide
species phylogeny, head posture transitions or rainfall adaptation.
"""
from __future__ import annotations
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path

SOURCE_TAXA=111
SOURCE_MASTER_WCVP=244
SOURCE_MASTER_NAMES=251
SOURCE_TARGETS=118
N_QUORUM=89
MIN_GENE_LOCI=30
MIN_DNA_SITES=10000
MIN_INFORMATIVE_SITES=50
MIN_SOURCE_TAXON_CALLS=1000
NA63_COMMON=("Cirsium flodmanii","Cirsium scopulorum",
              "Cirsium tracyi","Cirsium vulgare")


def csv_read(path):
    with Path(path).open(newline="",encoding="utf-8-sig") as f:
        a=list(csv.DictReader(f))
    if not a:
        raise ValueError("Original public SRA source CSV missing")
    return a


def fasta_read(path):
    out={}
    current=None
    for line in Path(path).read_text(encoding="ascii").splitlines():
        line=line.strip()
        if not line:continue
        if line.startswith(">"):
            name=line[1:].split()[0]
            if not name.startswith("WCVP") or name in out:
                raise ValueError("Only one original WCVP source tip per SRA BioSample")
            out[name]=""
            current=name
        else:
            if current is None or set(line)-set("ACGTN"):
                raise ValueError("Invalid original Q20 site call, reference fill or FASTA")
            out[current]+=line
    if not out or len({len(v) for v in out.values()})!=1:
        raise ValueError("Real DNA source alignment columns have inconsistent lengths")
    return out


def validate(data,genes,samples,sequences):
    if (len(samples)!=SOURCE_TAXA or len({s["phylogenetic_tip"] for s in samples})!=SOURCE_TAXA or
        len(genes)!=SOURCE_TARGETS or len({g["locus"] for g in genes})!=SOURCE_TARGETS):
        raise ValueError("Source 111 genuine BioSample and 118 candidate original loci needed")
    if (data.get("original_S1_WCVP_keys_retained")!=SOURCE_MASTER_WCVP or
        data.get("original_S1_source_name_operational_units_retained")!=SOURCE_MASTER_NAMES or
        data.get("all_original_111_S1_biological_samples")!=SOURCE_TAXA or
        data.get("n_original_source_candidate_loci")!=SOURCE_TARGETS or
        data.get("predeclared_site_occupancy_taxa")!=N_QUORUM or
        not data.get("technical_ML_alignment_source_gate_passed") or
        data.get("strict_single_copy_orthology_gate_passed") or
        data.get("orientation_evolution_rainfall_test_permitted")):
        raise ValueError("Real source Q20/80pct preflight not sufficient for technical ML")
    if set(sequences)!={s["phylogenetic_tip"] for s in samples}:
        raise ValueError("Nucleotide source FASTA does not preserve all 111 specimen identities")
    selected=[g for g in sorted(genes,key=lambda x:x["locus"])
              if g["retain_min100_Q20_positions"].lower()=="true"]
    expected_len=sum(int(g["n_sites_same_reference_called_at_least80pct_taxa"]) for g in selected)
    if (len(selected)!=data["n_loci_with_atleast100_site_positions_in_80pct_taxa"] or
        expected_len!=data["n_full111_alignable_positions_called_by_80pct_taxa"] or
        len(selected)<MIN_GENE_LOCI or expected_len<MIN_DNA_SITES):
        raise ValueError("Technical ML alignment genomic position gate changed")
    if any(len(seq)!=expected_len for seq in sequences.values()):
        raise ValueError("Source Q20 callable FASTA coordinate count drift")
    if min(len(z)-z.count("N") for z in sequences.values())<MIN_SOURCE_TAXON_CALLS:
        raise ValueError("Some original SRA taxa have fewer than 1000 real callable bases")
    variable=informative=0
    for pos in range(expected_len):
        a=Counter(s[pos] for s in sequences.values() if s[pos]!="N")
        if sum(a.values())<N_QUORUM:
            raise ValueError("Source masked nucleotide position fails original 89-taxon gate")
        variable+=len(a)>1
        informative+=sum(n>=2 for n in a.values())>=2
    if informative<MIN_INFORMATIVE_SITES:
        raise ValueError("Actual SRA Q20 loci have insufficient informative positions for ML")
    parts=[]
    start=1
    for row in selected:
        length=int(row["n_sites_same_reference_called_at_least80pct_taxa"])
        parts.append({"locus":row["locus"],
                      "reference_copy":row["original_reference_copy"],
                      "start_1based":start,"end_1based":start+length-1,
                      "no_verified_orthology_or_phase":True})
        start+=length
    tip_to_taxon={s["taxon_name"]:s["phylogenetic_tip"] for s in samples}
    missing_na=sorted(set(NA63_COMMON)-set(tip_to_taxon))
    summary={
        "schema":"azami_v4_original111_Q20_80pct_alignment_preliminary_technical_ML_gate_v1",
        "source_111_independent_public_SRA_BioSamples":SOURCE_TAXA,
        "source_244_WCVP_and_251_original_names_retained":True,
        "original_Moreyra_low_warning_241_universe_unchanged":True,
        "all111_selected_original_reference_locus_partitions":len(parts),
        "all111_source_Q20_80pct_reference_positions":expected_len,
        "variable_observed_Q20_positions":variable,
        "parsimony_informative_Q20_positions":informative,
        "minimum_source_tip_observation_fraction":min(
            (len(x)-x.count("N"))/len(x) for x in sequences.values()),
        "real_q20_data_unobserved_are_N_no_reference_imputation":True,
        "source_author_botanical_nodding_names":["Cirsium conspicuum","Cirsium scopulorum"],
        "source_independent_NA63_quartet_exact_name_overlap_if_available":[
            {"taxon":name,"tip":tip_to_taxon[name]} for name in NA63_COMMON
            if name in tip_to_taxon],
        "source_NA63_four_name_overlap_missing":missing_na,
        "preliminary_GTRG4_IQTREE_permitted":True,
        "single_copy_and_homeolog_orthology_verified":False,
        "global_EAzami294_species_tree_accepted":False,
        "adaptive_radiation_or_rain_transition_supported":False,
    }
    return parts,summary


def write_csv(path,rows):
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),extrasaction="raise")
        w.writeheader()
        w.writerows(rows)


def prepare(source_summary,genes,samples,fasta,outdir):
    summary=json.loads(Path(source_summary).read_text(encoding="utf-8"))
    genes=csv_read(genes);samples=csv_read(samples)
    seqs=fasta_read(fasta)
    partitions,report=validate(summary,genes,samples,seqs)
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    (out/"REAL111_OBSERVED_Q20_TECHNICAL_IQTREE_INPUT.fasta").write_bytes(Path(fasta).read_bytes())
    write_csv(out/"REAL111_ORIGINAL_LOCUS_PARTITIONS.csv",partitions)
    report["source_input_SHA256"]={
        "all111_Q20_site_summary":hashlib.sha256(Path(source_summary).read_bytes()).hexdigest(),
        "Q20_real_111_source_Fasta":hashlib.sha256(Path(fasta).read_bytes()).hexdigest(),
    }
    (out/"REAL111_TECHNICAL_IQTREE_INPUT_RECEIPT.json").write_text(
        json.dumps(report,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    print("T2_SOURCE111_TECHNICAL_ML_INPUT="+json.dumps(report,sort_keys=True),flush=True)
    return report


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--site-summary",type=Path,required=True)
    ap.add_argument("--loci",type=Path,required=True)
    ap.add_argument("--samples",type=Path,required=True)
    ap.add_argument("--fasta",type=Path,required=True)
    ap.add_argument("--out-dir",type=Path,required=True)
    x=ap.parse_args()
    prepare(x.site_summary,x.loci,x.samples,x.fasta,x.out_dir)


if __name__=="__main__":
    main()
