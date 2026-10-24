from analysis.v4 import t3_blind_s1_orientation_qc_source_packet as m

def test_phenotype_blind_23_taxa_115_image_queue_protocol():
    assert m.N_FULL==131
    assert m.N_READY==23
    assert m.N_PER_SPECIES==5
    assert m.N_READY*m.N_PER_SPECIES==115
    assert m.SEED==20261024
    assert len(m.SHA_JSON)==64
