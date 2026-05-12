from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCIENTIFIC_COVERAGE_DOC = ROOT / "docs" / "scientific-coverage.md"


def _scientific_coverage_text() -> str:
    return SCIENTIFIC_COVERAGE_DOC.read_text(encoding="utf-8")


def test_scientific_coverage_doc_exists() -> None:
    assert SCIENTIFIC_COVERAGE_DOC.is_file()


def test_scientific_confidence_labels_are_listed() -> None:
    text = _scientific_coverage_text()
    for label in (
        "PARITY_GATED_ACTIVE_SCIENCE",
        "PHOSPY_VALIDATED_SCIENCE",
        "SUPPORTED_CONTRACT_CHANGED",
        "OPEN_GAP",
    ):
        assert label in text


def test_parity_scope_terms_are_listed() -> None:
    lowered = _scientific_coverage_text().lower()
    for term in (
        "required parity",
        "deliberate scope difference",
        "useful future extension",
        "not planned",
    ):
        assert term in lowered


def test_scope_matrix_columns_and_required_rows_are_present() -> None:
    text = _scientific_coverage_text()
    lowered = text.lower()
    assert "| Area | Current confidence | Intended PhosR parity scope |" in text

    for row_name in (
        "input formats",
        "phosphosite representation",
        "site/flanking sequence",
        "localisation confidence",
        "replicate/condition modelling",
        "missing-value handling",
        "normalisation",
        "imputation",
        "batch correction",
        "differential phosphorylation",
        "kinase/substrate analysis",
        "motif/sequence-aware analysis",
        "enrichment analysis",
        "clustering/time-series",
        "visualisation",
        "reproducibility/reporting",
        "workflow composition/extensibility",
    ):
        assert f"| {row_name} |" in lowered
