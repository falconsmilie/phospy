from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from tests.support.fixture_byte_policy import (
    TEXT_FIXTURE_SUFFIXES,
    assert_lf_gitattributes_coverage,
    assert_text_fixture_matches_sha256,
    iter_importer_fixture_index_references,
    iter_manifest_governed_text_fixture_references,
)

ROOT = Path(__file__).resolve().parents[2]


def test_importer_fixture_index_references_canonical_raw_bytes() -> None:
    references = iter_importer_fixture_index_references(ROOT)

    assert {
        reference.path.relative_to(ROOT).as_posix() for reference in references
    } == {
        "tests/fixtures/fragpipe/ptmprophet_explicit_site_edge_cases.tsv",
        "tests/fixtures/fragpipe/ptmprophet_missing_required_start.tsv",
        "tests/fixtures/fragpipe/ptmprophet_peptide_position_edge_cases.tsv",
        "tests/fixtures/fragpipe/ptmprophet_sites.tsv",
        "tests/fixtures/maxquant/phospho_sty_sites_lfq_only.txt",
        "tests/fixtures/maxquant/phospho_sty_sites_multisite.txt",
        "tests/fixtures/maxquant/phospho_sty_sites_raw_and_lfq.txt",
        "tests/fixtures/maxquant/phospho_sty_sites_realistic_variants.txt",
        "tests/fixtures/maxquant/phospho_sty_sites_standard.txt",
    }
    for reference in references:
        assert reference.path.is_file()
        assert_text_fixture_matches_sha256(
            reference.path,
            expected_sha256=reference.expected_sha256,
            repo_root=ROOT,
        )


def test_manifest_governed_text_fixtures_have_lf_gitattributes_rules() -> None:
    references = iter_manifest_governed_text_fixture_references(ROOT)

    assert references, "expected manifest-governed text fixture references"
    assert_lf_gitattributes_coverage(
        ROOT,
        tuple(reference.path for reference in references),
    )


def test_fixture_lf_gitattributes_rules_remain_extension_specific() -> None:
    for line in (ROOT / ".gitattributes").read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        pattern, *attributes = stripped.split()
        if not pattern.startswith("tests/fixtures/"):
            continue
        if "text" not in attributes or "eol=lf" not in attributes:
            continue
        assert any(
            pattern.endswith(f"*{suffix}") for suffix in TEXT_FIXTURE_SUFFIXES
        ), (
            "tests/fixtures LF text rules must stay extension-specific so binary "
            f"fixtures are not promoted to text: {pattern}"
        )


@pytest.mark.parametrize(
    ("data", "violation"),
    (
        (b"header,value\r\nrow,1\r\n", "CRLF"),
        (b"header,value\rrow,1\r", "lone CR"),
        (b"header,value\nrow,1", "missing final newline"),
    ),
)
def test_text_fixture_byte_policy_reports_exact_newline_violation(
    tmp_path: Path,
    data: bytes,
    violation: str,
) -> None:
    fixture_path = tmp_path / "fixtures" / "bad_fixture.csv"
    fixture_path.parent.mkdir()
    fixture_path.write_bytes(data)

    with pytest.raises(AssertionError) as exc_info:
        assert_text_fixture_matches_sha256(
            fixture_path,
            expected_sha256=hashlib.sha256(data).hexdigest(),
            repo_root=tmp_path,
        )

    message = str(exc_info.value)
    assert "fixtures/bad_fixture.csv" in message
    assert violation in message


def test_text_fixture_byte_policy_reports_exact_digest_mismatch(
    tmp_path: Path,
) -> None:
    data = b"header,value\nrow,1\n"
    fixture_path = tmp_path / "fixtures" / "digest_mismatch.csv"
    fixture_path.parent.mkdir()
    fixture_path.write_bytes(data)
    expected_sha256 = "0" * 64

    with pytest.raises(AssertionError) as exc_info:
        assert_text_fixture_matches_sha256(
            fixture_path,
            expected_sha256=expected_sha256,
            repo_root=tmp_path,
        )

    message = str(exc_info.value)
    assert "fixtures/digest_mismatch.csv" in message
    assert "digest mismatch" in message
    assert f"expected sha256={expected_sha256}" in message
    assert f"actual sha256={hashlib.sha256(data).hexdigest()}" in message
