"""Local S1 concatenation: explicit tips, fixed loci, Carthamus-only root.

Generalises the checked EAzami c9870e9 20-tip concatenation contract to an
explicit S1 manifest. It never invents sequences for an entirely missing tip.
"""
import hashlib
from pathlib import Path

ROOT = 'OUTGROUP_saff'
REFERENCE_SHA = '77d510ef101d08a7a23a4df391d077d3b7f75482c66f7f4bea6d32cf290ced2c'


def fasta(path):
    rows = {}; name = None
    for raw in Path(path).read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith('>'):
            name = line[1:].split()[0]
            if name in rows:
                raise ValueError('Duplicate FASTA label: '+name)
            rows[name] = ''
        elif name is None:
            raise ValueError('Sequence before header')
        else:
            rows[name] += line.upper()
    if not rows or any(not seq or set(seq)-set('ACGTRYSWKMBDHVN-?') for seq in rows.values()):
        raise ValueError('Empty or invalid DNA FASTA')
    return rows


def write_fasta(path, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open('w') as f:
        for name, seq in rows.items():
            f.write('>'+name+'\n')
            for start in range(0, len(seq), 80):
                f.write(seq[start:start+80]+'\n')


def prepare_loci(loci, known, retrieved, reference, output):
    if hashlib.sha256(Path(reference).read_bytes()).hexdigest() != REFERENCE_SHA:
        raise ValueError('Reference digest mismatch')
    refs = fasta(reference)
    admitted = []; excluded = []
    for locus in loci:
        root_name = 'saff-'+locus
        if root_name not in refs:
            excluded.append(locus)
            continue
        records = fasta(Path(retrieved)/(locus+'.FNA'))
        if set(records)-known:
            raise ValueError('Unexpected recovered tip')
        records[ROOT] = refs[root_name]
        write_fasta(Path(output)/(locus+'.fasta'), records)
        admitted.append(locus)
    return admitted, excluded


def concatenate(loci, known, alignments):
    if len(loci) != len(set(loci)) or len(loci) < 100:
        raise ValueError('At least 100 unique admitted loci required')
    if not known or ROOT in known:
        raise ValueError('Explicit ingroup required')
    matrix = {tip: [] for tip in sorted(known)+[ROOT]}
    partitions = []; start = 1
    for locus in loci:
        rows = fasta(Path(alignments)/(locus+'.aln.fasta'))
        sizes = {len(seq) for seq in rows.values()}
        if len(sizes) != 1:
            raise ValueError('Ragged alignment: '+locus)
        if set(rows)-set(matrix) or ROOT not in rows:
            raise ValueError('Alignment tip or root mismatch')
        size = sizes.pop()
        for tip in matrix:
            matrix[tip].append(rows.get(tip, '-'*size))
        partitions.append(dict(locus=locus, start=start, end=start+size-1))
        start += size
    joined = {tip: ''.join(parts) for tip, parts in matrix.items()}
    if any(not (set(seq) & set('ACGT')) for seq in joined.values()):
        raise ValueError('An entirely missing tip cannot be placed into a nuclear tree')
    return joined, partitions
