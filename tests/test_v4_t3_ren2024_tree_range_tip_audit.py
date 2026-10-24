"""Synthetic HTTP range and real-publication contract tests (not empirical tree data)."""
import io
import re
import zipfile

import pytest
from analysis.v4 import t3_ren2024_tree_range_tip_audit as audit


class FakeResponse:
    def __init__(self, status, data=b"", headers=None, payload=None):
        self.status_code, self.data = status, data
        self.headers, self.payload = headers or {}, payload
    def json(self): return self.payload
    def raise_for_status(self):
        if self.status_code >= 400: raise RuntimeError(str(self.status_code))
    def iter_content(self, chunk_size=262144):
        for i in range(0, len(self.data), chunk_size):
            yield self.data[i:i+chunk_size]
    def close(self): pass


class FakeRangeSession:
    def __init__(self, data, *, ignore_range=False, bad_length=False):
        self.data, self.ignore_range, self.bad_length = data, ignore_range, bad_length
        self.calls = []
    def get(self, url, headers=None, stream=False, timeout=None, allow_redirects=True):
        start, end = map(int, re.match(r"bytes=(\d+)-(\d+)", headers["Range"]).groups())
        self.calls.append((start, end))
        if self.ignore_range:
            return FakeResponse(200, b"NEVER CONSUME THIS FULL ZIP")
        size = len(self.data) + (1 if self.bad_length else 0)
        return FakeResponse(206, self.data[start:end+1],
                            {"Content-Range": f"bytes {start}-{end}/{size}"})


def small_zip():
    out = io.BytesIO()
    with zipfile.ZipFile(out, "w", compression=zipfile.ZIP_DEFLATED) as z:
        z.writestr("root/1_G727_hyb/2_astral_analysis/astral.tre",
                   "((Cremanthodium_arnicoides:0.1,Otherus_aa:0.2):0.1,"
                   "(Ligularia_lingiana:0.1,Xenophyllum_esquilachense:0.2):0.3);")
        z.writestr("root/notes/readme.md", "synthetic protocol test")
    return out.getvalue()


def test_http_206_exact_byte_range_and_unrooted_binomial_string_match():
    archive = small_zip()
    session = FakeRangeSession(archive)
    result = audit.zip_member_tree(session, "https://example.test/archive.zip",
                                   expected_bytes=len(archive))
    assert result["zip_internal_path"].endswith("2_astral_analysis/astral.tre")
    assert not result["downloaded_entire_original_zip"]
    assert session.calls
    source = [{"taxon": n, "source_posture": "D",
               "taxonomic_status": "no_conflict_identified"} for n in (
               "Cremanthodium arnicoides", "Ligularia lingiana",
               "Xenophyllum esquilachense", "Senecio humbertii")]
    matches = audit.tip_binomial_matches(result["newick"], source)
    assert matches["published_astral_tree_tip_count"] == 4
    assert matches["n_candidate_source_names_with_exact_binomial_string_match"] == 3
    assert not matches["sufficient_tree_and_traits_to_infer_independent_evolutionary_gains"]
    assert not matches["biological_rain_protection_causality_inferred"]


def test_http_200_full_archive_response_is_denied_before_reading():
    b = small_zip()
    with pytest.raises(audit.StrictRangeError, match="Range required"):
        audit.zip_member_tree(FakeRangeSession(b, ignore_range=True),
                              "https://example.test/archive.zip", expected_bytes=len(b))


def test_bad_content_range_is_denied():
    b = small_zip()
    with pytest.raises(audit.StrictRangeError, match="Content-Range"):
        audit.zip_member_tree(FakeRangeSession(b, bad_length=True),
                              "https://example.test/archive.zip", expected_bytes=len(b))


def test_missing_real_ast_member_stops_instead_of_inventing_tree():
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        z.writestr("nuclear_reads.fastq", "ACGT")
    b = buf.getvalue()
    with pytest.raises(ValueError, match="Missing/ambiguous"):
        audit.zip_member_tree(FakeRangeSession(b),
                              "https://example.test/archive.zip", expected_bytes=len(b))


def test_dryad_immutable_source_sha_file_id_and_transfer_caps():
    assert audit.ARCHIVE_ID == 3337673
    assert audit.VERSION_ID == 340361
    assert audit.ARCHIVE_SIZE == 1583779312
    assert audit.MAX_TOTAL_RANGE_BYTES < audit.ARCHIVE_SIZE // 20
    assert audit.ARCHIVE_SHA256 == (
        "962495cf8131fadd3a1cdf1ea88204432f2aaa81c9e6c0a97435756f0ad88139")
    class FakeMeta:
        def get(self, *args, **kwargs):
            return FakeResponse(200, payload={"_embedded": {"stash:files": [
                {"_links": {"self": {"href": "/api/v2/files/3337673"}},
                 "size": audit.ARCHIVE_SIZE,
                 "digestType": "sha-256", "digest": audit.ARCHIVE_SHA256}
            ]}})
    assert audit.check_official_manifest(FakeMeta())[
        "published_sha_verified_against_live_dryad_metadata"]
