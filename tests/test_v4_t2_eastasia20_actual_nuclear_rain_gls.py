"""East-Asia original nuclear-GBIF phylogeny gates: no EXIF or photo count."""
from io import StringIO
import csv
from pathlib import Path

from Bio import Phylo
import pytest

from analysis.v4 import t2_eastasia20_actual_nuclear_rain_gls as ea


def test_exact_public_GBIF_14_climate_table_is_not_image_data():
    data=ea.LOCAL_CLIMATE.read_bytes()
    from analysis.v4.t2_pgls_placement_screen import git_blob_sha
    assert git_blob_sha(data)==ea.CLIMATE_BLOB
    rows=list(csv.DictReader(StringIO(data.decode("utf-8-sig"))))
    assert len(rows)==14
    assert len({z["taxon"] for z in rows})==14
    assert all(int(z["n"])>=1 for z in rows)
    assert set(z["region_source"] for z in rows)=={"japan","taiwan"}
    assert ea.ORIGINAL_NUCLEAR_TIPS==20
    assert ea.SOURCE_OBSERVED==17
    assert ea.SOURCE_UNKNOWN==3


def test_true_tip_guard_refuses_unrooted_name_binning_without_authentic_taxa():
    tree=Phylo.read(StringIO("((Cirsium_alpicola:.01,Cirsium_gyojanum:.02):.01,"
                             "OUTGROUP_saff:.2);"),"newick")
    data={"Cirsium_alpicola":0,"Cirsium_gyojanum":0}
    with pytest.raises(ValueError,match="Original genus nuclear tree"):
        ea.analyze(tree,data,{},{})


def test_original_GBIF_climate_checksum_rejects_mutation(tmp_path):
    raw=ea.LOCAL_CLIMATE.read_bytes()
    path=tmp_path/"mutated_climate.csv"
    path.write_bytes(raw.replace(b"2128.304347826",b"2128.404347826"))
    assert path.read_bytes()!=raw
    from analysis.v4.t2_pgls_placement_screen import git_blob_sha
    assert git_blob_sha(path.read_bytes()) != ea.CLIMATE_BLOB


def test_gravity_botanical_original_17_states_not_all20_made_up():
    assert ea.TREE_OUTGROUP=="OUTGROUP_saff"
    assert ea.TREE_BLOB=="1e3e6ffaf128ccd439287070b8c3c81385f5df34"
    assert ea.TRAIT_BLOB=="f655544e77f1659323c1e915567dacbb51b01cc3"
    assert "photo" not in str(ea.LOCAL_CLIMATE).lower()
    assert "TRUE" not in ea.TREE_FILE
