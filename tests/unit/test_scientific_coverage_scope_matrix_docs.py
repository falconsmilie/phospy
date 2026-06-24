from __future__ import annotations

from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python < 3.11
    import tomli as tomllib  # type: ignore[no-redef]

ROOT = Path(__file__).resolve().parents[2]
README_DOC = ROOT / "README.md"
SCIENTIFIC_COVERAGE_DOC = ROOT / "docs" / "scientific-coverage.md"
API_GUIDE_DOC = ROOT / "docs" / "api" / "guide.md"
WORKFLOW_CONTRACTS_DOC = ROOT / "docs" / "workflow_contracts.md"
IMPORTERS_DOC = ROOT / "docs" / "importers.md"
PARITY_DOC = ROOT / "docs" / "parity.md"
MAINTENANCE_DOC = ROOT / "docs" / "maintenance.md"
CONTRIBUTING_DOC = ROOT / "docs" / "contributing.md"
PYPROJECT = ROOT / "pyproject.toml"
ADR_0025_DOC = (
    ROOT
    / "docs"
    / "adr"
    / "adr_0025_competitive_phosphoproteomics_workflow_coverage.md"
)
ADR_0027_DOC = (
    ROOT
    / "docs"
    / "adr"
    / "adr_0027_target_future_native_phosr_style_sps_ruv_iii_correction.md"
)


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _readme_text() -> str:
    return _read(README_DOC)


def _scientific_coverage_text() -> str:
    return _read(SCIENTIFIC_COVERAGE_DOC)


def _api_guide_text() -> str:
    return _read(API_GUIDE_DOC)


def _workflow_contracts_text() -> str:
    return _read(WORKFLOW_CONTRACTS_DOC)


def _importers_text() -> str:
    return _read(IMPORTERS_DOC)


def _parity_text() -> str:
    return _read(PARITY_DOC)


def _adr_0025_text() -> str:
    return _read(ADR_0025_DOC)


def _adr_0027_text() -> str:
    return _read(ADR_0027_DOC)


def _package_description() -> str:
    pyproject_data = tomllib.loads(_read(PYPROJECT))
    return pyproject_data["project"]["description"]


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


def test_package_metadata_does_not_claim_phosr_parity() -> None:
    description = _package_description()
    normalized = " ".join(description.lower().split())

    assert description == (
        "Python package for selected PhosR-inspired phosphoproteomics workflows"
    )
    assert "phosr-inspired" in normalized
    assert "parity" not in normalized
    assert "ruv" not in normalized
    assert "sps" not in normalized
    assert "ruv-iii" not in normalized
    assert "support" not in normalized


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
        "batch correction: `spsruvbatchcorrectionconfig`",
        "phosr-equivalent sps/ruv-iii, combat, and `removebatcheffect` parity",
        "enrichment",
        "visualisation",
        "supported bundled organisms and references",
        "full phosr package equivalence claim",
    ):
        assert f"| {row_name} |" in lowered


def test_docs_list_current_workflow_packages() -> None:
    text = (
        _readme_text()
        + "\n"
        + _scientific_coverage_text()
        + "\n"
        + _workflow_contracts_text()
        + "\n"
        + _parity_text()
        + "\n"
        + _adr_0025_text()
        + "\n"
        + _adr_0027_text()
    )
    normalized = " ".join(text.lower().split())

    for workflow in (
        "DifferentialAnalysisWorkflow",
        "EnrichmentWorkflow",
        "KinaseWorkflow",
        "SignalomeWorkflow",
    ):
        assert workflow.lower() in normalized

    readme = _readme_text().lower()
    assert "differential" in readme
    assert "enrichment" in readme
    assert "kinase" in readme
    assert "signalome" in readme
    assert "future enrichment" not in normalized
    assert "future-only enrichment" not in normalized
    assert "enrichment is future" not in normalized


def test_docs_importer_support_matches_public_importers() -> None:
    import phospy.io.readers as readers

    docs_text = (
        _scientific_coverage_text()
        + "\n"
        + _importers_text()
        + "\n"
        + _parity_text()
        + "\n"
        + _adr_0025_text()
    )
    normalized = " ".join(docs_text.lower().split())
    public_symbols = set(readers.__all__)

    supported_importers = {
        "MappedPhosphositeTableImporter",
        "MaxQuantPhosphositeImporter",
        "FragPipePTMProphetImporter",
    }
    unsupported_importers = {
        "SpectronautPhosphositeImporter",
        "DIANNPhosphositeImporter",
        "DiaNNPhosphositeImporter",
    }

    assert supported_importers <= public_symbols
    assert public_symbols.isdisjoint(unsupported_importers)
    for importer_name in supported_importers:
        assert importer_name.lower() in normalized

    assert "spectronaut and dia-nn phosphosite importers are not currently" in (
        normalized
    )
    assert "not spectronaut/dia-nn support" in normalized
    assert "not upstream statistical result import" in normalized


def test_docs_do_not_claim_ruv_support() -> None:
    text = (
        _readme_text()
        + "\n"
        + _scientific_coverage_text()
        + "\n"
        + _workflow_contracts_text()
        + "\n"
        + _parity_text()
        + "\n"
        + _adr_0025_text()
        + "\n"
        + _adr_0027_text()
    ).lower()
    normalized = " ".join(text.split())

    assert "| batch correction: `spsruvbatchcorrectionconfig` |" in normalized
    assert (
        "| phosr-equivalent sps/ruv-iii, combat, and `removebatcheffect` parity | `open gap` |"
        in normalized
    )
    assert "no phosr-equivalent sps/ruv-iii, combat, or limma" in normalized
    assert "native sps/ruv-style preprocessing correction" in normalized
    assert "do not interpret `ruv_readiness` as ruv support" in normalized
    assert "limited fixed-effect residualisation" in normalized
    assert "not combat" in normalized
    assert "not ruv" in normalized
    assert "not limma `removebatcheffect` parity" in normalized
    assert "not current phosr parity" in normalized

    for forbidden in (
        "supports ruv",
        "ruv is supported",
        "ruv correction is supported",
        "ruv-compatible correction is supported",
        "combat is supported",
        "removebatcheffect is supported",
    ):
        assert forbidden not in normalized


def test_docs_distinguish_linear_residualisation_from_ruv_sps() -> None:
    text = (
        _readme_text()
        + "\n"
        + _scientific_coverage_text()
        + "\n"
        + _parity_text()
        + "\n"
        + _adr_0027_text()
    ).lower()
    normalized = " ".join(text.split())

    assert "linear_residualize_batch" in normalized
    assert "fixed-effect residualisation" in normalized
    assert "not native sps/ruv-style correction" in normalized
    assert "not phosr-equivalent batch correction" in normalized
    assert "not equivalent to phosr-style ruv/sps correction" in normalized
    assert "no treatment of linear residualisation as equivalent" in normalized


def test_scientific_coverage_lists_sps_ruv_as_explicit_non_parity_support() -> None:
    normalized = " ".join(_scientific_coverage_text().lower().split())

    assert "### ruv/sps/ruv-iii batch-correction posture" in normalized
    assert "supported today:" in normalized
    assert "not supported today:" in normalized
    assert "`spsruvbatchcorrectionconfig`" in normalized
    assert "caller-supplied controls" in normalized
    assert "unwanted-factor count" in normalized
    assert "phosr-equivalent sps/ruv-iii batch correction" in normalized
    assert "not a claim of phosr-equivalent sps/ruv-iii correction" in normalized


def test_parity_docs_state_ruv_sps_is_future_work() -> None:
    normalized = " ".join(_parity_text().lower().split())

    assert "not currently parity-equivalent with phosr for sps/ruv-iii" in (normalized)
    assert "no sps control-selection fixtures" in normalized
    assert "no native ruv-iii correction-kernel parity fixtures" in normalized
    assert "no phosr `ruvphospho` corrected-output parity fixtures" in normalized
    assert "validated phospy implementation, not current phosr parity" in normalized


def test_adr_0027_records_future_preprocessing_direction_and_constraints() -> None:
    normalized = " ".join(_adr_0027_text().lower().split())

    assert "adr-0027" in normalized
    assert "status:** accepted" in normalized
    assert (
        "phospy targets native sps/ruv-style correction as a preprocessing/normalisation"
        in normalized
    )
    assert "not phosr-equivalent sps/ruv-iii parity" in normalized
    assert "preprocessing/normalisation" in normalized
    assert "must not be owned by differential analysis" in normalized
    assert "kinase analysis, enrichment, or signalome execution" in normalized
    assert "`spsruvbatchcorrectionconfig`" in normalized
    assert "`ruv_readiness` and similarly named diagnostics are report-only" in (
        normalized
    )
    assert "sps/control phosphosite selection" in normalized
    assert "unwanted-factor count" in normalized
    assert "temporary imputation rules" in normalized
    assert "restore missingness after correction" in normalized
    assert "carry an observation/imputation mask forward" in normalized
    assert "flag imputed positions in downstream outputs" in normalized
    assert "withhold features from downstream statistical testing" in normalized
    assert "temporary imputation for correction mechanics is not equivalent" in (
        normalized
    )


def test_readme_does_not_imply_phosr_batch_correction_parity() -> None:
    normalized = " ".join(_readme_text().lower().split())

    assert "not native ruv/sps/ruv-iii correction" in normalized
    assert "not phosr-equivalent batch correction" in normalized
    assert "readiness signals and do not apply correction" in normalized
    assert "future architecture commitment, not a current feature claim" in normalized


def test_docs_scientific_coverage_mentions_enrichment_scope() -> None:
    text = (
        _readme_text()
        + "\n"
        + _scientific_coverage_text()
        + "\n"
        + _parity_text()
        + "\n"
        + _adr_0025_text()
    ).lower()
    normalized = " ".join(text.split())

    assert "| enrichment | `validated phospy implementation` |" in normalized
    assert "`enrichmentworkflow` runs offline over-representation analysis" in (
        normalized
    )
    assert "caller-supplied `genesetcollection`, `ptmsetcollection`" in normalized
    assert "background universe is explicit and required" in normalized
    assert "not a phosr enrichment parity lane" in normalized
    assert "ora is not gsea, ssgsea, or ptm-sea support" in normalized
    assert "future enrichment" not in normalized


def test_pyright_docs_match_current_include_and_strict_scope() -> None:
    pyproject_data = tomllib.loads(_read(PYPROJECT))
    pyright_config = pyproject_data["tool"]["pyright"]
    include_paths = tuple(pyright_config["include"])
    strict_paths = tuple(pyright_config["strict"])
    docs_text = (_read(MAINTENANCE_DOC) + "\n" + _read(CONTRIBUTING_DOC)).lower()

    expected_include_paths = (
        "src/phospy/api",
        "src/phospy/errors",
        "src/phospy/frames",
        "src/phospy/io",
        "src/phospy/policies",
        "src/phospy/provenance",
        "src/phospy/science",
        "src/phospy/tables",
        "src/phospy/validation",
        "src/phospy/workflows",
    )
    strict_dataset_model = "src/phospy/science/datasets/models.py"

    assert include_paths == expected_include_paths
    assert strict_dataset_model in strict_paths
    for include_path in expected_include_paths:
        assert f"`{include_path}`" in docs_text
    assert f"`{strict_dataset_model}`" in docs_text
    assert "already strict-checked" in docs_text


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


def test_multiple_testing_finite_denominator_and_differential_validation_are_documented() -> (
    None
):
    normalized = " ".join(_scientific_coverage_text().lower().split())

    assert "shared multiple-testing helper adjusts only finite p-values" in (normalized)
    assert "denominator is the number of finite p-values passed to the helper" in (
        normalized
    )
    assert "differential execution validates generated `p.value` values are finite" in (
        normalized
    )
    assert "withheld rows receive missing" in normalized
    assert "are excluded from the benjamini-hochberg denominator" in normalized


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


def test_batch_correction_scope_names_supported_preprocessing_methods() -> None:
    text = (
        _scientific_coverage_text()
        + "\n"
        + _api_guide_text()
        + "\n"
        + _workflow_contracts_text()
    ).lower()
    normalized = " ".join(text.split())

    assert "| batch correction: `linear_residualize_batch` |" in normalized
    assert "| batch correction: `spsruvbatchcorrectionconfig` |" in normalized
    assert (
        "| phosr-equivalent sps/ruv-iii, combat, and `removebatcheffect` parity | `open gap` |"
        in normalized
    )
    assert "`linear_residualize_batch` fixed-effect residualisation" in normalized
    assert "preserves condition effects by design" in normalized
    assert "confounded batch/condition designs are rejected" in normalized
    assert "not combat" in normalized
    assert "not ruv" in normalized
    assert "not limma `removebatcheffect` parity" in normalized
    assert "not mixed-effects modelling" in normalized
    assert "does not solve all batch-effect problems" in normalized
    assert "caller-supplied controls" in normalized
    assert "correction remains in dataset preprocessing" in normalized
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
