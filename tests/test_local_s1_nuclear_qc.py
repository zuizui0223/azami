import pytest
from analysis.v4.local_s1_nuclear_qc import gate, analyse


def test_occupancy_rounds_up_and_any_paralog_excludes_locus():
    loci = [f'L{i}' for i in range(241)]
    samples = set('ABCDEFGH')
    present = {'L0': set('ABCDEFG'), 'L1': set('ABCDEF'), 'L2': samples}
    rows, admitted = gate(loci, samples, present, {'L2': {'A'}})
    assert admitted == ['L0']
    assert all(row['minimum_present'] == 7 for row in rows)


def test_unknown_sample_cannot_inflate_occupancy():
    with pytest.raises(ValueError, match='Unexpected sample'):
        gate([f'L{i}' for i in range(241)], {'A', 'B'}, {'L0': {'A', 'alien'}}, {})


def test_unfinished_panel_is_not_silently_pruned(tmp_path):
    manifest = tmp_path/'manifest.csv'
    manifest.write_text('run\nA\nB\n')
    (tmp_path/'sample_names.txt').write_text('A\n')
    with pytest.raises(ValueError, match='do not prune'):
        analyse(manifest, tmp_path, tmp_path/'unused-loci.txt')
