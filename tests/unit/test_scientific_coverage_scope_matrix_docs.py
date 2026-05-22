from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCIENTIFIC_COVERAGE_DOC = ROOT / "docs" / "scientific-coverage.md"


def _scientific_coverage_text() -> str:
    return SCIENTIFIC_COVERAGE_DOC.read_text(encoding="utf-8")


def test_scientific_coverage_doc_exists() -> None:
    assert SCIENTIFIC_COVERAGE_DOC.is_file()


def test_scope_categories_are_listed() -> None:
    text = _scientific_coverage_text()
    for label in (
        "parity-gated",
        "validated PhosPy implementation",
        "experimental",
        "open gap",
        "deliberate scope difference",
        "not planned",
    ):
        assert label in text


def test_scope_matrix_is_the_single_source_of_truth() -> None:
    text = _scientific_coverage_text()
    lowered = text.lower()
    assert "## Scientific Scope Matrix (Single Source Of Truth)" in text
    assert "full phosr package equivalence is not" in lowered
    assert "not claimed" in lowered


def test_scope_matrix_columns_and_required_rows_are_present() -> None:
    text = _scientific_coverage_text()
    lowered = text.lower()
    assert (
        "| Area | Scope category | Current executable support | "
        "Evidence and release checks | Limits and non-claims |"
    ) in text

    for row_name in (
        "differential analysis",
        "kinase scoring",
        "kinase activity scoring",
        "kinase prediction",
        "signalome analysis",
        "sequence context",
        "localisation handling",
        "missing values",
        "normalisation",
        "imputation",
        "batch correction / ruv",
        "enrichment",
        "visualisation",
        "supported bundled organisms and references",
        "full phosr package equivalence claim",
    ):
        assert f"| {row_name} |" in lowered
