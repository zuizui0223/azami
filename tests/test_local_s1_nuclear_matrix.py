import pytest
from analysis.v4.local_s1_nuclear_matrix import concatenate, write_fasta, ROOT
from analysis.v4.local_s1_nuclear_matrix import project_anchor_columns


def test_projection_keeps_exact_source_coordinates_and_never_adds_anchor_tips():
    rows = {'lett-L':'A-CNT', 'sunf-L':'AGC-T', 'saff-L':'ATCGT', 'A':'A-CGT'}
    projected, columns = project_anchor_columns(rows, 'L', {'A', 'missing'})
    assert columns == [0, 2, 4]
    assert projected == {ROOT:'ACT', 'A':'ACT'}
    assert 'missing' not in projected


def test_projection_rejects_unexpected_tip_and_absent_homology():
    rows = {'lett-L':'A-', 'sunf-L':'-C', 'saff-L':'AC', 'A':'AC'}
    with pytest.raises(ValueError, match='No shared'):
        project_anchor_columns(rows, 'L', {'A'})
    rows['foreign'] = 'AC'
    with pytest.raises(ValueError, match='identity mismatch'):
        project_anchor_columns(rows, 'L', {'A'})


def panel(tmp_path, last=None):
    loci = [f'L{i}' for i in range(100)]
    for locus in loci:
        write_fasta(tmp_path/(locus+'.aln.fasta'), {'A':'ACGT', ROOT:'AGGT'})
    if last:
        write_fasta(tmp_path/(loci[-1]+'.aln.fasta'), last)
    return loci


def test_never_adds_an_unsequenced_tip(tmp_path):
    loci = panel(tmp_path)
    with pytest.raises(ValueError, match='entirely missing'):
        concatenate(loci, {'A','B'}, tmp_path)


def test_missing_locus_is_gap_padded_but_partition_bounds_are_exact(tmp_path):
    loci = panel(tmp_path, {'A':'ACGT', 'B':'TGCA', ROOT:'AGGT'})
    matrix, partitions = concatenate(loci, {'A','B'}, tmp_path)
    assert matrix['B'] == '-'*396+'TGCA'
    assert partitions[-1] == {'locus':'L99', 'start':397, 'end':400}


def test_root_identity_and_alignment_shape_are_checked(tmp_path):
    loci = panel(tmp_path, {'A':'ACGT', ROOT:'AC'})
    with pytest.raises(ValueError, match='Ragged'):
        concatenate(loci, {'A'}, tmp_path)
    write_fasta(tmp_path/'L99.aln.fasta', {'A':'ACGT', 'foreign':'TGCA'})
    with pytest.raises(ValueError, match='root mismatch'):
        concatenate(loci, {'A'}, tmp_path)
