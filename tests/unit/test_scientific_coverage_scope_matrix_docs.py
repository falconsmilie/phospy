from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCIENTIFIC_COVERAGE_DOC = ROOT / "docs" / "scientific-coverage.md"
API_GUIDE_DOC = ROOT / "docs" / "api" / "guide.md"


def _scientific_coverage_text() -> str:
    return SCIENTIFIC_COVERAGE_DOC.read_text(encoding="utf-8")


def _api_guide_text() -> str:
    return API_GUIDE_DOC.read_text(encoding="utf-8")


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


def test_fixed_effect_differential_limitations_are_documented() -> None:
    text = (_scientific_coverage_text() + "\n" + _api_guide_text()).lower()

    assert "fixed-effect covariates" in text
    assert "batch can be modelled" in text
    assert "not batch correction" in text
    assert "combat" in text
    assert "ruv" in text
    assert "removebatcheffect" in text
    assert "duplicatecorrelation" in text
    assert "mixed-effects" in text
    assert "fixed-block" in text
    assert "complete within-block contrast coverage" in text


def test_differential_provenance_docs_name_fixed_block_limitations() -> None:
    text = (_scientific_coverage_text() + "\n" + _api_guide_text()).lower()
    normalized = " ".join(text.split())

    assert 'paired_design_policy="fixed_block"' in normalized
    assert "fixed-block paired designs are supported only when" in normalized
    assert "block terms are ordinary fixed effects" in normalized
    assert "block terms are fixed effects" in normalized
    assert "not limma `duplicatecorrelation`" in normalized
    assert "not mixed-effects modelling" in normalized
    assert "no mixed effects are fitted" in normalized
    assert "not random subject modelling" in normalized
    assert "incomplete or partially covered blocks are rejected" in normalized
    assert "does not drop those blocks or samples" in normalized
    assert "simple unpaired workflows" in normalized
