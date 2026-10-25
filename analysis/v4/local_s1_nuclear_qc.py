"""Strict, aggregate-only locus admission for an explicit local S1 panel.

Retains the existing 241-locus / 80% occupancy / zero-copy-warning / 100-locus
contract. Incomplete panels fail closed; completing a smaller panel is never
silently substituted for the requested panel. This is not tree acceptance.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

LOCUS_SHA = "d561c6e393b1964fdd4b3acf14fda8b10f2f43923b1074cd35f86bfed07ebf73"


def gate(loci, known, present, warned):
    if not known or len(loci) != 241 or len(set(loci)) != 241:
        raise ValueError("Explicit panel and frozen 241 loci required")
    if any((set(v) - known) for v in [*present.values(), *warned.values()]):
        raise ValueError("Unexpected sample in sequence or paralog data")
    minimum = math.ceil(.8 * len(known))
    rows = []
    for locus in sorted(loci):
        n = len(present.get(locus, set()))
        warning = len(warned.get(locus, set()))
        rows.append(dict(locus=locus, present=n, panel_size=len(known),
                         minimum_present=minimum, paralog_warning_samples=warning,
                         eligible=n >= minimum and warning == 0))
    admitted = [r['locus'] for r in rows if r['eligible']]
    return rows, admitted


def analyse(manifest, qc, locus_path):
    panel = list(csv.DictReader(manifest.open(encoding='utf-8-sig', newline='')))
    known = {r['run'] for r in panel}
    if not known or len(known) != len(panel):
        raise ValueError("Unique explicit run manifest required")
    names = set((qc/'sample_names.txt').read_text().splitlines())
    if names != known:
        raise ValueError("QC panel differs from requested manifest; do not prune missing samples")
    raw = locus_path.read_bytes()
    canonical = raw.replace(b'\r\n', b'\n')
    if hashlib.sha256(canonical).hexdigest() != LOCUS_SHA:
        raise ValueError("Frozen locus checksum mismatch")
    loci = raw.decode().splitlines()
    paralogs = list(csv.DictReader((qc/'paralog_report.tsv').open(), delimiter='\t'))
    if len(paralogs) != len(known) or {r['Species'] for r in paralogs} != known:
        raise ValueError("Paralog report panel mismatch")
    warned = {}
    for r in paralogs:
        for locus, value in r.items():
            if locus != 'Species' and int(value) > 1:
                warned.setdefault(locus, set()).add(r['Species'])
    present = {}
    for path in (qc/'retrieved_dna').glob('*.FNA'):
        names = [line[1:].split()[0] for line in path.read_text().splitlines()
                 if line.startswith('>')]
        if len(names) != len(set(names)):
            raise ValueError("Duplicate sample sequence at " + path.stem)
        present[path.stem] = set(names)
    rows, admitted = gate(loci, known, present, warned)
    summary = dict(status='LOCUS_QC_ONLY_NOT_TREE_ACCEPTANCE',
                   requested_samples=len(known), candidate_loci=241,
                   eligible_loci=len(admitted), minimum_loci=100,
                   enough_loci_for_tree_inputs=len(admitted) >= 100,
                   automatic_filter_relaxation=False,
                   locus_file_sha256=hashlib.sha256(raw).hexdigest(),
                   locus_LF_canonical_sha256=hashlib.sha256(canonical).hexdigest(),
                   tree_completed=False,
                   remaining_checks=['coding_sequence_warnings', 'root_references',
                                     'alignment_quality', 'tree_support_and_tip_identity'])
    return rows, admitted, summary


def main():
    p = argparse.ArgumentParser()
    for key in ('manifest', 'qc', 'loci', 'out'):
        p.add_argument('--'+key, type=Path, required=True)
    a = p.parse_args()
    rows, admitted, summary = analyse(a.manifest, a.qc, a.loci)
    a.out.mkdir(parents=True, exist_ok=True)
    with (a.out/'locus_qc.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    (a.out/'eligible_loci.txt').write_text(''.join(x+'\n' for x in admitted))
    (a.out/'summary.json').write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary))


if __name__ == '__main__':
    main()
