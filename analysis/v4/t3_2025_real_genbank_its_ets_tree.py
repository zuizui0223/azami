"""Infer T3 Senecioneae original GenBank ITS+ETS linked-rDNA ML genealogy.

Real source: 2025 Liu & Jiang Taiwania Table S1, 163 public accession
records, 154 source taxon labels. Exactly matching modern NCBI organism
names/ITS accepted from original accession audit (never synonym-forced).
Each named species enters ONCE, after one original supplement row is
selected based on complete linked ITS+ETS sequence availability.
Genes are ITS / ETS spacers of SAME multicopy nuclear ribosomal repeat:
not two ortholog-independent nuclear gene trees, not a species tree.

A source-independent T3 direct species flora/experiment head posture
crosswalk maps known D/unknown developmental V but there are no validated
flowering-stage erect controls in this phylogenomic source set. NO independent
nodding origin count, coevolution with local rainfall, functional causal
comparison or species-tree posteriors are inferred from this topology.
"""
from __future__ import annotations
import argparse
import csv
from collections import Counter
import hashlib
import io
import json
from pathlib import Path

from Bio import SeqIO, Phylo

SOURCE_TABLE_PATH=Path(__file__).resolve().parent/"results"/"20261023_T3_TAIWANIA2025_ORIGINAL163_ITS_ETS_GENBANK_ACCESSIONS.csv"
FOCAL_TRAIT_PATH=Path(__file__).resolve().parent/"results"/"20261023_T3_REN2021_TRUE_SRA_ORIGINAL_POSTURE_SUPPLEMENT.csv"
ORIGINAL_TRAITS_PATH=Path(__file__).resolve().parent/"results"/"20261023_T3_SENECIONEAE_POSTURE_SOURCE_CANDIDATES.csv"
ORIGINAL_RECORDS=163
ORIGINAL_NAMES=154
MIN_ITS_TIPS=100
MIN_ETS_TIPS=50
MIN_ITS_INFORMATIVE=35
MIN_ETS_INFORMATIVE=15
OCCUPANCY_ITS=.57
OCCUPANCY_ETS=.55
MAX_NUCLEOTIDE_ALIGNMENT=10000

def rows(path):
    with Path(path).open(encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f))

def load_fasta(path,expected):
    source=list(SeqIO.parse(str(path),"fasta"))
    if not source:return {}
    seq={}
    for record in source:
        name=record.id
        if name.startswith("_R_"):name=name[3:]
        name=name.replace("_"," ")
        if name not in expected or name in seq or "\t" in name:
            raise ValueError("Unexpected or duplicate original published GenBank source taxon in MAFFT")
        sequence=str(record.seq).upper().replace("U","T")
        if any(z not in "ACGTN-?RYSWKMBDHV" for z in sequence):
            raise ValueError("Unrecognized MAFFT original nucleotide")
        seq[name]=sequence
    if len({len(z) for z in seq.values()})!=1:
        raise ValueError("MAFFT must produce one alignment length per original rDNA region")
    return seq

def trim_marker(aligned,*,occupancy,min_tips,min_information,kind):
    if kind not in ("ITS","ETS"):
        raise ValueError("Unrecognized linked nuclear rDNA region")
    if len(aligned)<min_tips:return {},{"status":"INSUFFICIENT_GENUINE_SOURCE_TIPS","source_tips":len(aligned)}
    names=sorted(aligned)
    length=len(aligned[names[0]])
    if not 180<=length<=MAX_NUCLEOTIDE_ALIGNMENT:
        raise ValueError("Original MSA region has implausible nucleotide alignment length")
    retained=[]
    variable=informative=0
    for i in range(length):
        calls=[aligned[n][i] for n in names if aligned[n][i] in "ACGT"]
        if len(calls)/len(names)<occupancy:continue
        counts=Counter(calls)
        variable+=len(counts)>=2
        informative+=sum(c>=2 for c in counts.values())>=2
        retained.append(i)
    if len(retained)<180 or informative<min_information:
        return {},{
            "status":"INSUFFICIENT_HOMOLOGOUS_REAL_RDNA_SITES",
            "original_msa_sites":length,"passed_sites":len(retained),
            "variable_sites":variable,"parsimony_informative_sites":informative,
        }
    sequences={name:"".join(aligned[name][i] if aligned[name][i] in "ACGT"
                            else "N" for i in retained) for name in names}
    return sequences,{
        "status":"GENBANK_SOURCE_RIBOSOMAL_REGION_MSA_QC_PASS",
        "region":kind,"original_msa_sites":length,
        "retained_sites":len(retained),"source_tips":len(names),
        "variable_sites":variable,
        "parsimony_informative_sites":informative,
        "original_site_occupancy_cutoff":occupancy,
        "ambiguous_real_genbank_bases_are_N_NOT_reference_fills":True,
    }

def original_states():
    known={}
    for evidence in rows(ORIGINAL_TRAITS_PATH):
        if evidence["source_posture"] in ("D","U","?","V"):
            known[evidence["taxon"]]={
                "posture":evidence["source_posture"],
                "source_url":evidence["source_url"],
                "source_stage":evidence["life_stage_scope"],
                "status":"ORIGINAL_SPECIES_AUTHORITY",
            }
    for evidence in rows(FOCAL_TRAIT_PATH):
        name=evidence["taxon"]
        posture=evidence["source_posture"]
        if name in known and known[name]["posture"]!=posture:
            raise ValueError("Conflicting original T3 published capitulum state")
        known[name]={"posture":posture,
                      "source_url":evidence["independent_flora_or_experiment_url"],
                      "source_stage":evidence["life_stage_scope"],
                      "status":"SUPPLEMENTAL_SPECIES_AUTHORITY"}
    return known

def prepare(mapping,its_file,ets_file,source_receipt,outdir):
    published=rows(SOURCE_TABLE_PATH)
    if len(published)!=ORIGINAL_RECORDS or len({r["original_published_scientific_name"] for r in published})!=ORIGINAL_NAMES:
        raise ValueError("Original Taiwania 163 GenBank accession/154 names lost")
    pairs=rows(mapping)
    receipt=json.loads(Path(source_receipt).read_text())
    if (receipt.get("original_163_source_rows")!=163
        or receipt.get("original_154_unique_published_taxa")!=154
        or receipt.get("source_rows_queried")!=163
        or receipt.get("source_failure_batches_not_ITS_absence")!=[]
        or not receipt["original_ITS_ETS_rDNA_spacers_physically_linked_not_independent_nuclear_orthologs"]):
        raise ValueError("Actual complete 2025 GenBank provenance not verified")
    names={r["published_source_species"] for r in pairs}
    if len(names)!=len(pairs):
        raise ValueError("Public source duplicated biological taxon in inferred species-tip registry")
    if not names<={r["original_published_scientific_name"] for r in published}:
        raise ValueError("Selected GenBank rDNA taxon not in original publication")
    ITS=load_fasta(its_file,names)
    ETS=load_fasta(ets_file,names)
    if set(ITS)!=names or not set(ETS)<=names:
        raise ValueError("Original exact-GenBank ITS vs ETS MAFFT taxa differ from registry")
    if len(ITS)!=receipt["n_unique_taxa_with_exact_published_organism_ITS"] or len(ETS)!=receipt["n_unique_taxa_with_exact_published_organism_both_ITS_and_ETS"]:
        raise ValueError("Unrelated sequence set substituted for actual original 163 accession source")
    a,r1=trim_marker(ITS,occupancy=OCCUPANCY_ITS,min_tips=MIN_ITS_TIPS,
                     min_information=MIN_ITS_INFORMATIVE,kind="ITS")
    b,r2=trim_marker(ETS,occupancy=OCCUPANCY_ETS,min_tips=MIN_ETS_TIPS,
                     min_information=MIN_ETS_INFORMATIVE,kind="ETS")
    d=Path(outdir);d.mkdir(parents=True,exist_ok=True)
    ready=(r1["status"].endswith("QC_PASS") and r2["status"].endswith("QC_PASS"))
    evidence=original_states()
    coverage=Counter(evidence.get(n,{}).get("posture","?") for n in names)
    if ready:
        dna={}
        for n in sorted(names):
            dna[n]=a[n]+b.get(n,"N"*len(next(iter(b.values()))))
        fasta=d/"T3_2025_AUTHENTIC_GENBANK_LINKED_ITS_ETS_CONCAT.fasta"
        with fasta.open("w",encoding="ascii") as f:
            for taxon,seq in dna.items():
                f.write(">"+taxon.replace(" ","_")+"\n")
                for i in range(0,len(seq),80):f.write(seq[i:i+80]+"\n")
        (d/"T3_2025_NUCLEAR_rDNA_SOURCE_PARTITIONS.txt").write_text(
            f"DNA, ITS = 1-{r1['retained_sites']}\n"
            f"DNA, ETS = {r1['retained_sites']+1}-{r1['retained_sites']+r2['retained_sites']}\n",
            encoding="utf-8")
        for marker,calls in (("ITS",a),("ETS",b)):
            with (d/f"T3_2025_{marker}_ACTUAL_GENBANK_SITE_QC_ALIGNED.fasta").open("w",encoding="ascii") as handle:
                for taxon,seq in calls.items():
                    handle.write(f">{taxon.replace(' ','_')}\n{seq}\n")
    report={
        "status":"T3_REAL_2025_LINKED_ITS_ETS_RIBOSOMAL_QC_READY" if ready
                 else "T3_REAL_2025_LINKED_ITS_ETS_RIBOSOMAL_QC_STOP",
        "source_original_163_rows":163,
        "source_original_154_distinct_scientific_names":154,
        "source_selected_authenticated_ITS_taxa":len(names),
        "source_selected_authenticated_ETS_taxa":len(ETS),
        "same_original_ITS_ETS_marker_record_ITS_ETS_genes_independent":False,
        "marker_site_occupancy": {"ITS":r1,"ETS":r2},
        "actual_2025_source_head_state_distribution":dict(coverage),
        "species_level_U_D_fixed_ancestral_transition_rate_computable":False,
        "source_original_GenBank_retrieval_digest":hashlib.sha256(Path(source_receipt).read_bytes()).hexdigest(),
        "tip_photo_camera_vertical_used":False,
        "tree_inference_origin_not_reused_2021_published_phylogeny":True,
        "rain_and_ancestral_gain_causality_not_calculated":True,
    }
    (d/"T3_2025_AUTHENTIC_ITS_ETS_MSA_QC_RECEIPT.json").write_text(
        json.dumps(report,ensure_ascii=False,sort_keys=True,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print("T3_GENBANK_REAL_ITS_ETS_ALIGNMENT_PREP="+json.dumps(report,sort_keys=True,ensure_ascii=False),flush=True)
    return report

def verify_tree(tree_path,mapping,prep,outdir):
    report=json.loads(Path(prep).read_text(encoding="utf-8"))
    if report["status"]!="T3_REAL_2025_LINKED_ITS_ETS_RIBOSOMAL_QC_READY":
        raise ValueError("Source-authentic ITS+ETS marker QC did not pass")
    taxa=rows(mapping)
    names={x["published_source_species"] for x in taxa}
    t=Phylo.read(str(tree_path),"newick")
    terminals=[x.name.replace("_"," ") for x in t.get_terminals()]
    if len(terminals)!=len(names) or len(set(terminals))!=len(names) or set(terminals)!=names:
        raise ValueError("Actual newly inferred IQTREE linked ribosomal tree does not contain original biology name tips")
    annotation=original_states()
    labels=[]
    for tip in t.get_terminals():
        name=tip.name.replace("_"," ")
        info=annotation.get(name,{})
        labels.append({
            "GenBank_original_published_species":name,
            "original_species_morphology_state":info.get("posture","?"),
            "species_flora_or_manipulation_source":info.get("source_url",""),
            "DNA_voucher_equals_flora_specimen":False,
            "linked_multicopy_rDNA_not_single_copy_independent_genes":True,
        })
    counts=Counter(x["original_species_morphology_state"] for x in labels)
    d={x.name for x in t.get_terminals()
       if annotation.get(x.name.replace("_"," "),{}).get("posture")=="D"}
    lca_clades=Counter(x.replace("_"," ").split()[0] for x in d)
    unique_genera=len(lca_clades)
    res={
        "status":"REAL_2025_PUBLISHED_GENBANK_ITS_ETS_IQTREE_OWN_NUCLEAR_RIBOSOMAL_GENEALOGY",
        "original_2025_genbank_linked_ITS_ETS_source_biological_tips":len(names),
        "posture_sources_tip_distribution":dict(counts),
        "phylogram_substitutions_per_site_not_dated":True,
        "phylogenetic_replicate_independent_gene_trees":False,
        "source_documented_nodding_tip_genera":dict(lca_clades),
        "n_source_different_genera_with_nodding_heads":unique_genera,
        "two_genera_are_not_proof_of_multiple_U_to_D_independent_origins":True,
        "no_rooted_true_Senecioneae_species_tree_or_ancestral_climate":True,
        "no_erect_morphology_control_on_source_genuine_rDNA_genealogy_proven":counts.get("U",0)==0,
        "new_original_IQTREE_tree_sha256":hashlib.sha256(Path(tree_path).read_bytes()).hexdigest(),
        "original_source_genbank_name_sequence_archive_sha256":report["source_original_GenBank_retrieval_digest"],
    }
    dpath=Path(outdir);dpath.mkdir(parents=True,exist_ok=True)
    with (dpath/"T3_2025_OWN_REAL_LINKED_ITS_ETS_GENBANK_TREE_TIP_EVIDENCE.csv").open(
        "w",newline="",encoding="utf-8") as h:
        w=csv.DictWriter(h,fieldnames=list(labels[0]))
        w.writeheader();w.writerows(labels)
    (dpath/"T3_2025_OWN_REAL_NUCLEAR_rDNA_ANCESTRY_AND_POSTURE.json").write_text(
        json.dumps(res,indent=2,ensure_ascii=False,sort_keys=True,allow_nan=False)+"\n",
        encoding="utf-8")
    print("T3_2025_LINKED_REAL_ITS_ETS_ML="+json.dumps(res,ensure_ascii=False,sort_keys=True),flush=True)
    return res

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--stage",required=True,choices=("prepare","tree"))
    p.add_argument("--mapping",type=Path,required=True)
    p.add_argument("--source-receipt",type=Path)
    p.add_argument("--ITS-msa",type=Path)
    p.add_argument("--ETS-msa",type=Path)
    p.add_argument("--prepared",type=Path)
    p.add_argument("--tree",type=Path)
    p.add_argument("--outdir",type=Path,required=True)
    a=p.parse_args()
    if a.stage=="prepare":
        if not a.source_receipt or not a.ITS_msa or not a.ETS_msa:
            raise ValueError("Source receipt and independently MAFFT aligned markers required")
        prepare(a.mapping,a.ITS_msa,a.ETS_msa,a.source_receipt,a.outdir)
    else:
        if not a.prepared or not a.tree:
            raise ValueError("True IQ-TREE and marker QC required")
        verify_tree(a.tree,a.mapping,a.prepared,a.outdir)

if __name__=="__main__":main()
