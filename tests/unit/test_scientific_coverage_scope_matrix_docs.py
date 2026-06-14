from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCIENTIFIC_COVERAGE_DOC = ROOT / "docs" / "scientific-coverage.md"
API_GUIDE_DOC = ROOT / "docs" / "api" / "guide.md"
WORKFLOW_CONTRACTS_DOC = ROOT / "docs" / "workflow_contracts.md"
IMPORTERS_DOC = ROOT / "docs" / "importers.md"
ADR_0025_DOC = (
    ROOT
    / "docs"
    / "adr"
    / "adr_0025_competitive_phosphoproteomics_workflow_coverage.md"
)


def _scientific_coverage_text() -> str:
    return SCIENTIFIC_COVERAGE_DOC.read_text(encoding="utf-8")


def _api_guide_text() -> str:
    return API_GUIDE_DOC.read_text(encoding="utf-8")


def _workflow_contracts_text() -> str:
    return WORKFLOW_CONTRACTS_DOC.read_text(encoding="utf-8")


def _importers_text() -> str:
    return IMPORTERS_DOC.read_text(encoding="utf-8")


def _adr_0025_text() -> str:
    return ADR_0025_DOC.read_text(encoding="utf-8")


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
        "phosphosite importers",
        "missing values",
        "normalisation",
        "imputation",
        "batch correction: `linear_residualize_batch`",
        "ruv, combat, and `removebatcheffect` parity",
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


def test_batch_correction_scope_names_only_linear_residualize_batch_supported() -> None:
    text = (
        _scientific_coverage_text()
        + "\n"
        + _api_guide_text()
        + "\n"
        + _workflow_contracts_text()
    ).lower()
    normalized = " ".join(text.split())

    assert "| batch correction: `linear_residualize_batch` |" in normalized
    assert "| ruv, combat, and `removebatcheffect` parity | `open gap` |" in normalized
    assert "`linear_residualize_batch` fixed-effect residualisation" in normalized
    assert "preserves condition effects by design" in normalized
    assert "confounded batch/condition designs are rejected" in normalized
    assert "not combat" in normalized
    assert "not ruv" in normalized
    assert "not limma `removebatcheffect` parity" in normalized
    assert "not mixed-effects modelling" in normalized
    assert "does not solve all batch-effect problems" in normalized
    assert "do not interpret `ruv_readiness` as ruv support" in normalized


def test_protein_aware_preparation_scope_is_separate_from_modelling() -> None:
    text = (
        _scientific_coverage_text()
        + "\n"
        + _api_guide_text()
        + "\n"
        + _workflow_contracts_text()
    ).lower()
    normalized = " ".join(text.split())

    assert "| total-protein subtraction: `subtract_log_total` |" in normalized
    assert (
        "| protein-aware preparation | `validated phospy implementation` |"
        in normalized
    )
    assert (
        "| joint ptm/protein modelling and msstatsptm-style inference | `open gap` |"
    ) in normalized
    assert "`log2_phospho - log2_total`" in normalized
    assert "preparation-only" in normalized
    assert "does not modify phosphosite values" in normalized
    assert "does not subtract total protein" in normalized
    assert "does not run joint ptm/protein differential modelling" in normalized
    assert "does not adjust differential models" in normalized
    assert "does not claim msstatsptm-style inference" in normalized
    assert "current `differentialanalysisworkflow` does not consume" in normalized


def test_enrichment_scope_is_offline_ora_with_user_supplied_collections() -> None:
    text = (
        _scientific_coverage_text()
        + "\n"
        + _api_guide_text()
        + "\n"
        + _workflow_contracts_text()
    ).lower()
    normalized = " ".join(text.split())

    assert "| enrichment | `validated phospy implementation` |" in normalized
    assert "offline over-representation analysis" in normalized
    assert "caller-supplied `genesetcollection`, `ptmsetcollection`" in normalized
    assert "background universe is explicit and required" in normalized
    assert "go, kegg, reactome, ptm-sea" in normalized
    assert "are not bundled" in normalized
    assert "enrichr, gseapy, clusterprofiler" in normalized
    assert "not native core workflow behavior" in normalized
    assert "ora is not gsea, ssgsea, or ptm-sea support" in normalized
    assert "gene-level and site-level enrichment require explicit identifier" in (
        normalized
    )


def test_importer_scope_documents_targeted_maxquant_fragpipe_support() -> None:
    text = (
        _scientific_coverage_text()
        + "\n"
        + _api_guide_text()
        + "\n"
        + _importers_text()
        + "\n"
        + _adr_0025_text()
    ).lower()
    normalized = " ".join(text.split())

    assert "| phosphosite importers | `validated phospy implementation` |" in (
        normalized
    )
    assert "maxquantphosphositeimporter" in normalized
    assert "fragpipeptmprophetimporter" in normalized
    assert "mappedphosphositetableimporter" in normalized
    assert "phosphositeimportresult" in normalized
    assert "dataset-builder requests" in normalized
    assert "do not construct analysis-ready datasets" in normalized
    assert "infer sample groups" in normalized
    assert "infer contrasts" in normalized
    assert "infer batches or blocks" in normalized
    assert "infer differential design" in normalized
    assert "bypass builder validation" in normalized
    assert "not broad support for all vendor" in normalized
    assert "spectronaut/dia-nn support" in normalized
    assert "upstream statistical result import" in normalized
    assert "does not currently provide broad semantic importers" not in normalized


def test_adr_0025_current_state_names_implemented_support_without_parity_claims() -> (
    None
):
    normalized = " ".join(_adr_0025_text().lower().split())

    assert "current executable analysis/workflow lanes" in normalized
    assert "`enrichmentworkflow` for offline over-representation analysis" in (
        normalized
    )
    assert "maxquant phosphosite import" in normalized
    assert "fragpipe/philosopher/ptmprophet phosphosite import" in normalized
    assert '`paired_design_policy="fixed_block"`' in normalized
    assert "`linear_residualize_batch`" in normalized
    assert '`datasetproteinawarepreparationconfig(policy="prepare_model_inputs")`' in (
        normalized
    )
    assert "`kinaselibraryresource` / `kinaselibraryresourceloader`" in normalized
    assert "`ssgsea_substrate_enrichment_activity_v1`" in normalized
    assert "offline ora through `enrichmentworkflow`" in normalized

    assert "not combat" in normalized
    assert "not ruv" in normalized
    assert "does not bundle official kinase library data" in normalized
    assert "does not claim validated kinase library parity" in normalized
    assert "not ptm-sea parity" in normalized
    assert "does not imply gsea, ssgsea, or ptm-sea support" in normalized
    assert "does not claim msstatsptm-style" in normalized
    assert "duplicatecorrelation" in normalized
    assert "mixed-effects modelling" in normalized

    stale_fragments = (
        "batch-aware, block, paired, and repeated-measure modeling are not executable",
        "does not currently provide broad semantic importers",
    )
    for fragment in stale_fragments:
        assert fragment not in normalized
