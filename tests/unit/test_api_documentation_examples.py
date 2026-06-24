from __future__ import annotations

import ast
import re
from pathlib import Path

import pandas as pd

from phospy import AnalysisReadyDatasetBuilder
from phospy.api import (
    Contrast,
    ControlSiteSet,
    CorrectionMissingnessPolicy,
    DatasetBatchCorrectionConfig,
    DatasetBuildRequest,
    DatasetIntensityTransformConfig,
    DatasetLocalisationConfig,
    DatasetMissingDataConfig,
    DatasetNormalisationConfig,
    DatasetPreprocessingConfig,
    DifferentialAnalysisRequest,
    EnrichmentConfig,
    EnrichmentWorkflow,
    EnrichmentWorkflowRequest,
    ExperimentalDesign,
    GeneSetCollection,
    IntensityScaleKind,
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflowRequest,
    KinaseWorkflowResult,
    ObservationMask,
    Organism,
    OriginallyMissingCellTracking,
    ReferenceBundle,
    ReferencePreset,
    SampleDesignRecord,
    SignalomeClusteringConfig,
    SignalomeConfig,
    SignalomeOutputConfig,
    SignalomePerformanceConfig,
    SignalomeScientificConfig,
    SignalomeValidationConfig,
    SignalomeWorkflowRequest,
    SpsRuvBatchCorrectionConfig,
    TemporaryImputationMethod,
    TemporaryImputationPolicy,
)
from phospy.api.results import KinasePredictionResult, KinaseScoringResult

ROOT = Path(__file__).resolve().parents[2]
README = ROOT / "README.md"
API_DOCS_DIR = ROOT / "docs" / "api"
WORKFLOW_DOCS_DIR = API_DOCS_DIR
DATASET_BUILD_DOC = API_DOCS_DIR / "dataset-build-workflow.md"
DIFFERENTIAL_DOC = WORKFLOW_DOCS_DIR / "differential-analysis.md"
ENRICHMENT_DOC = WORKFLOW_DOCS_DIR / "enrichment.md"
KINASE_DOC = WORKFLOW_DOCS_DIR / "kinase.md"
SIGNALOME_DOC = WORKFLOW_DOCS_DIR / "signalome.md"


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _normalise_whitespace(text: str) -> str:
    return " ".join(text.split())


def _iter_python_code_blocks(source: str) -> tuple[str, ...]:
    return tuple(
        match.group("code").strip()
        for match in re.finditer(
            r"```python\s*\n(?P<code>.*?)\n```",
            source,
            flags=re.DOTALL,
        )
    )


def _parse_python_code_blocks(source: str) -> tuple[ast.Module, ...]:
    parsed_blocks: list[ast.Module] = []
    for block in _iter_python_code_blocks(source):
        try:
            parsed_blocks.append(ast.parse(block))
        except SyntaxError as exc:  # pragma: no cover - assertion context only
            raise AssertionError(f"invalid documented Python example: {exc}") from exc
    return tuple(parsed_blocks)


def _imported_names(source: str, module_name: str) -> set[str]:
    names: set[str] = set()
    for tree in _parse_python_code_blocks(source):
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.module == module_name:
                names.update(alias.name for alias in node.names)
    return names


def _assert_python_imports(
    source: str,
    module_name: str,
    names: tuple[str, ...],
    *,
    context: str,
) -> None:
    imported = _imported_names(source, module_name)
    missing = sorted(set(names) - imported)
    assert not missing, f"{context} missing imports from {module_name}: {missing}"


def _assert_python_imports_absent(
    source: str,
    module_name: str,
    names: tuple[str, ...],
    *,
    context: str,
) -> None:
    imported = _imported_names(source, module_name)
    present = sorted(set(names) & imported)
    assert not present, f"{context} documents unsupported imports: {present}"


def _call_name(call: ast.Call) -> str | None:
    if isinstance(call.func, ast.Name):
        return call.func.id
    return None


def _is_workflow_run_call(call: ast.Call, workflow_name: str) -> bool:
    return (
        isinstance(call.func, ast.Attribute)
        and call.func.attr == "run"
        and isinstance(call.func.value, ast.Call)
        and isinstance(call.func.value.func, ast.Name)
        and call.func.value.func.id == workflow_name
    )


def _assert_python_call(source: str, call_name: str, *, context: str) -> None:
    assert any(
        _call_name(call) == call_name
        for tree in _parse_python_code_blocks(source)
        for call in ast.walk(tree)
        if isinstance(call, ast.Call)
    ), f"{context} must include a documented {call_name}(...) example"


def _assert_python_run_call(
    source: str,
    workflow_name: str,
    *,
    context: str,
) -> None:
    assert any(
        _is_workflow_run_call(call, workflow_name)
        for tree in _parse_python_code_blocks(source)
        for call in ast.walk(tree)
        if isinstance(call, ast.Call)
    ), f"{context} must include {workflow_name}().run(...)"


_MISSING = object()


def _literal_value(node: ast.AST) -> object:
    if isinstance(node, ast.Constant):
        return node.value
    if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name):
        return f"{node.value.id}.{node.attr}"
    return _MISSING


def _assert_python_call_keyword(
    source: str,
    call_name: str,
    keyword_name: str,
    expected_value: object = _MISSING,
    *,
    context: str,
) -> None:
    for tree in _parse_python_code_blocks(source):
        for call in ast.walk(tree):
            if not isinstance(call, ast.Call) or _call_name(call) != call_name:
                continue
            for keyword in call.keywords:
                if keyword.arg != keyword_name:
                    continue
                if expected_value is _MISSING:
                    return
                if _literal_value(keyword.value) == expected_value:
                    return

    expected = (
        keyword_name
        if expected_value is _MISSING
        else f"{keyword_name}={expected_value!r}"
    )
    raise AssertionError(f"{context} must document {call_name}({expected}, ...)")


def _assert_python_constant(source: str, expected: object, *, context: str) -> None:
    assert any(
        isinstance(node, ast.Constant) and node.value == expected
        for tree in _parse_python_code_blocks(source)
        for node in ast.walk(tree)
    ), f"{context} must include literal {expected!r} in a Python example"


def _attribute_chain(node: ast.AST) -> tuple[str, ...] | None:
    if isinstance(node, ast.Name):
        return (node.id,)
    if isinstance(node, ast.Attribute):
        parent = _attribute_chain(node.value)
        if parent is not None:
            return (*parent, node.attr)
    return None


def _assert_python_attribute_chain(
    source: str,
    chain: tuple[str, ...],
    *,
    context: str,
) -> None:
    assert any(
        _attribute_chain(node) == chain
        for tree in _parse_python_code_blocks(source)
        for node in ast.walk(tree)
        if isinstance(node, ast.Attribute)
    ), f"{context} must document {'.'.join(chain)}"


def _assert_documented_terms(
    source: str,
    expected: tuple[str, ...],
    *,
    context: str,
) -> None:
    missing = [term for term in expected if term not in source]
    assert not missing, f"{context} missing public documentation terms: {missing}"


def _iter_documentation_statements(source: str) -> tuple[str, ...]:
    statements: list[str] = []
    for block in re.split(r"\n\s*\n", source):
        stripped_block = block.strip()
        if not stripped_block:
            continue
        if stripped_block.startswith("|"):
            chunks = stripped_block.splitlines()
        else:
            chunks = re.split(r"\n(?=\s*[-*]\s+)", stripped_block)
        for chunk in chunks:
            normalised = _normalise_whitespace(chunk)
            if not normalised or normalised.startswith("```"):
                continue
            statements.extend(
                sentence
                for sentence in re.split(r"(?<=[.!?])\s+", normalised)
                if sentence
            )
    return tuple(statements)


_NEGATED_SCOPE_PATTERN = re.compile(
    r"\b(no|not|without|unsupported|report-only)\b|"
    r"\b(does|do|is|are)\s+not\b",
    flags=re.IGNORECASE,
)
_FAILURE_PATTERN = re.compile(r"\bfail(?:s|ed|ing)?\b", flags=re.IGNORECASE)


def _statement_contains_term(statement: str, term: str) -> bool:
    return bool(
        re.search(
            rf"(?<![a-z0-9_]){re.escape(term.lower())}(?![a-z0-9_])",
            statement.lower(),
        )
    )


def _assert_negated_scope_statement(
    source: str,
    unsupported_terms: tuple[str, ...],
    *,
    context: str,
) -> None:
    statements = _iter_documentation_statements(source)
    missing = []
    for term in unsupported_terms:
        if not any(
            _statement_contains_term(statement, term)
            and _NEGATED_SCOPE_PATTERN.search(statement)
            for statement in statements
        ):
            missing.append(term)

    assert not missing, (
        f"{context} must state unsupported scope with negation for: {missing}"
    )


def _assert_statement_contains_all(
    source: str,
    required_terms: tuple[str, ...],
    *,
    context: str,
) -> None:
    assert any(
        all(_statement_contains_term(statement, term) for term in required_terms)
        for statement in _iter_documentation_statements(source)
    ), f"{context} must include a statement containing {required_terms}"


def _assert_failure_boundary_statement(
    source: str,
    required_terms: tuple[str, ...],
    *,
    context: str,
) -> None:
    assert any(
        _FAILURE_PATTERN.search(statement)
        and all(_statement_contains_term(statement, term) for term in required_terms)
        for statement in _iter_documentation_statements(source)
    ), f"{context} must include a failure statement containing {required_terms}"


def _tiny_inputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    phospho = pd.DataFrame(
        {
            "A_1": [100.0, 70.0],
            "A_2": [101.0, 72.0],
            "B_1": [120.0, 80.0],
            "B_2": [118.0, 82.0],
        },
        index=["TSC2;S939;", "GSK3B;S9;"],
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["TSC2", "GSK3B"],
            "site": ["S939", "S9"],
            "site_sequence": [
                "FDDTPEKDSFRARSTSLNERPKSLRIARAPK",
                "ATMSGRPRTTSFAESSKPVQQPSAFGQAAAL",
            ],
            "protein_id": ["TSC2", "GSK3B"],
            "localisation_confidence": [0.95, 0.92],
        },
        index=phospho.index.copy(),
    )
    return phospho, site_metadata


def _build_dataset():
    phospho, site_metadata = _tiny_inputs()
    return AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            preprocessing_config=DatasetPreprocessingConfig(
                intensity_transform=DatasetIntensityTransformConfig(
                    policy="log2",
                    pseudocount=1.0,
                ),
                localisation=DatasetLocalisationConfig(
                    mode="require_threshold",
                    confidence_column="localisation_confidence",
                    min_confidence=0.75,
                ),
            ),
        )
    )


def _minimal_kinase_result() -> KinaseWorkflowResult:
    dataset = _build_dataset()
    site_keys = dataset.phospho.index.astype(str).tolist()
    display_ids = dataset.site_metadata.loc[:, "display_id"].astype(str).tolist()
    index = pd.Index(site_keys, name="site_key")
    kinases = pd.Index(["K1"], name="kinase")
    score_matrix = pd.DataFrame([[0.9], [0.8]], index=index, columns=kinases)
    prediction_matrix = pd.DataFrame([[0.7], [0.6]], index=index, columns=kinases)
    references = ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["K1", "K1"],
                "substrate_site": display_ids,
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": ["A" * 31, "C" * 31]},
            index=pd.Index(display_ids, name="site_id"),
        ),
    )
    return KinaseWorkflowResult(
        dataset=dataset,
        references=references,
        scoring_result=KinaseScoringResult(profile_scores=score_matrix),
        prediction_result=KinasePredictionResult(pred_mat=prediction_matrix),
        activity_result=None,
    )


def test_api_docs_public_imports_are_valid() -> None:
    namespace: dict[str, object] = {}

    exec(
        """from phospy import (
    AnalysisReadyDatasetBuilder,
    AnalysisReadyPhosphoDataset,
    DifferentialAnalysisWorkflow,
    KinaseWorkflow,
    SignalomeWorkflow,
)""",
        namespace,
    )
    exec(
        """from phospy.api import (
    DatasetBuildRequest,
    ExperimentalDesign,
    Contrast,
    SampleDesignRecord,
    DatasetPreprocessingConfig,
    DifferentialAnalysisRequest,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
    ReferencePreset,
    SignalomeConfig,
    SignalomeWorkflowRequest,
    UnsupportedInputFormatError,
    WorkflowValidationError,
)""",
        namespace,
    )
    exec(
        """from phospy.api import (
    EnrichmentConfig,
    EnrichmentWorkflow,
    EnrichmentWorkflowRequest,
    GeneSetCollection,
)""",
        namespace,
    )

    assert "AnalysisReadyDatasetBuilder" in namespace
    assert "DifferentialAnalysisWorkflow" in namespace
    assert "EnrichmentWorkflow" in namespace
    assert "DatasetBuildRequest" in namespace
    assert "DifferentialAnalysisRequest" in namespace
    assert "EnrichmentWorkflowRequest" in namespace
    assert "SignalomeWorkflowRequest" in namespace


def test_each_public_workflow_has_dedicated_api_page_with_contract_classes() -> None:
    expected = {
        DIFFERENTIAL_DOC: (
            "DifferentialAnalysisWorkflow",
            "DifferentialAnalysisRequest",
            "DifferentialAnalysisConfig",
            "DifferentialAnalysisResult",
        ),
        ENRICHMENT_DOC: (
            "EnrichmentWorkflow",
            "EnrichmentWorkflowRequest",
            "EnrichmentConfig",
            "EnrichmentWorkflowResult",
        ),
        KINASE_DOC: (
            "KinaseWorkflow",
            "KinaseWorkflowRequest",
            "KinaseScoringConfig",
            "KinasePredictionConfig",
            "KinaseActivityConfig",
            "KinaseWorkflowResult",
        ),
        SIGNALOME_DOC: (
            "SignalomeWorkflow",
            "SignalomeWorkflowRequest",
            "SignalomeConfig",
            "SignalomeWorkflowResult",
        ),
    }

    for path, class_names in expected.items():
        assert path.exists(), f"missing workflow API page: {path}"
        source = _read(path)
        for class_name in class_names:
            assert class_name in source, f"{class_name} missing from {path}"


def test_api_docs_enrichment_example_uses_public_api_and_runs_offline() -> None:
    source = _read(ENRICHMENT_DOC)

    _assert_python_imports(
        source,
        "phospy.api",
        (
            "EnrichmentConfig",
            "EnrichmentWorkflow",
            "EnrichmentWorkflowRequest",
            "GeneSetCollection",
        ),
        context="enrichment public API example",
    )
    _assert_python_call(
        source,
        "GeneSetCollection",
        context="enrichment public API example",
    )
    _assert_python_call_keyword(
        source,
        "EnrichmentWorkflowRequest",
        "background_universe",
        context="enrichment public API example",
    )
    _assert_python_run_call(
        source,
        "EnrichmentWorkflow",
        context="enrichment public API example",
    )
    _assert_statement_contains_all(
        source,
        ("provenance", "offline", "online-resource"),
        context="enrichment provenance scope",
    )
    _assert_negated_scope_statement(
        source,
        ("GO", "KEGG", "Reactome", "PTM-SEA", "Enrichr", "gseapy"),
        context="enrichment bundled/online resource scope",
    )
    _assert_negated_scope_statement(
        source,
        ("GSEA", "ssGSEA", "PTM-SEA"),
        context="enrichment method scope",
    )

    collection = GeneSetCollection(
        sets={
            "kinase_response": ("AKT1", "MAPK1", "MTOR"),
            "cell_cycle": ("CDK1", "CDK2", "MAPK1"),
        },
        identifier_kind="gene_symbol",
        term_names={
            "kinase_response": "Kinase response",
            "cell_cycle": "Cell cycle",
        },
        source_name="example in-memory gene sets",
        source_version="2026-06",
    )
    request = EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind="gene_symbol",
        set_collection=collection,
        selected_identifiers=("AKT1", "MAPK1"),
        background_universe=("AKT1", "MAPK1", "MTOR", "CDK1", "CDK2"),
        config=EnrichmentConfig(
            method="over_representation",
            multiple_testing_correction="benjamini_hochberg",
        ),
    )

    result = EnrichmentWorkflow().run(request)

    assert tuple(result.table["term_id"]) == ("kinase_response", "cell_cycle")
    assert result.provenance is not None
    assert result.provenance.workflow_parameters["background_universe_size"] == 5


def test_api_docs_dataset_build_request_example_is_constructible() -> None:
    phospho, site_metadata = _tiny_inputs()
    preprocessing = DatasetPreprocessingConfig(
        intensity_transform=DatasetIntensityTransformConfig(
            policy="log2",
            pseudocount=1.0,
        ),
        normalisation=DatasetNormalisationConfig(policy="median_center"),
        missing_data=DatasetMissingDataConfig(policy="forbid"),
    )

    request = DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        preprocessing_config=preprocessing,
    )

    assert request.organism is Organism.RAT
    assert request.preprocessing_config.normalisation.policy == "median_center"


def test_api_docs_batch_correction_example_is_constructible() -> None:
    source = _read(DATASET_BUILD_DOC)

    _assert_python_call_keyword(
        source,
        "DatasetBatchCorrectionConfig",
        "method",
        "linear_residualize_batch",
        context="batch-correction configuration example",
    )
    _assert_documented_terms(
        source,
        ("dataset.preprocessing_report.batch_correction",),
        context="batch-correction report contract",
    )
    _assert_statement_contains_all(
        source,
        ("confounding", "status"),
        context="batch-correction report fields",
    )
    _assert_negated_scope_statement(
        source,
        ("ComBat", "RUV", "removeBatchEffect", "mixed-effects modelling"),
        context="batch-correction scientific scope",
    )

    phospho = pd.DataFrame(
        {
            "sample_1": [10.0, 2.0],
            "sample_2": [15.0, 7.0],
            "sample_3": [14.0, 1.0],
            "sample_4": [19.0, 6.0],
        },
        index=["MAPK14;Y182;", "AKT1;T308;"],
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": [
                "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA",
                "AAAAAAAAAAAAAAATAAAAAAAAAAAAAAA",
            ],
            "display_id": ["MAPK14;Y182;", "AKT1;T308;"],
            "organism": ["rat", "rat"],
            "protein_namespace": ["protein_id", "protein_id"],
            "protein_identifier": ["MAPK14", "AKT1"],
            "protein_id": ["MAPK14", "AKT1"],
            "localisation_confidence": [0.95, 0.92],
        },
        index=phospho.index.copy(),
    )
    sample_metadata = pd.DataFrame(
        {
            "batch": ["run_1", "run_2", "run_1", "run_2"],
            "condition": ["control", "control", "treated", "treated"],
        },
        index=phospho.columns.copy(),
    )
    preprocessing = DatasetPreprocessingConfig(
        batch_correction=DatasetBatchCorrectionConfig(
            method="linear_residualize_batch",
            batch_column="batch",
            condition_column="condition",
            preserve_condition_effects=True,
        )
    )

    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=sample_metadata,
            organism=Organism.RAT,
            input_intensity_scale=IntensityScaleKind.LOG2,
            preprocessing_config=preprocessing,
        )
    )

    assert dataset.preprocessing_report is not None
    report = dataset.preprocessing_report.batch_correction
    assert report is not None
    assert report.status == "applied"
    assert report.method == "linear_residualize_batch"
    assert report.confounding_check_status == "passed"
    assert report.batch_levels == ("run_1", "run_2")
    assert report.condition_levels == ("control", "treated")


def test_api_docs_sps_ruv_batch_correction_example_is_explicit() -> None:
    source = _read(DATASET_BUILD_DOC)

    _assert_python_call_keyword(
        source,
        "SpsRuvBatchCorrectionConfig",
        "n_unwanted_factors",
        1,
        context="SPS/RUV batch-correction configuration example",
    )
    _assert_python_call_keyword(
        source,
        "SpsRuvBatchCorrectionConfig",
        "provenance_enabled",
        True,
        context="SPS/RUV batch-correction provenance setting",
    )
    _assert_documented_terms(
        source,
        (
            "ControlSiteSet.from_site_keys",
            "CorrectionMissingnessPolicy",
            "ObservationMask",
            "caller-supplied",
            "Minimal valid native-correction example",
            "Rejected unsafe example",
            "TemporaryImputationMethod.ROW_MEDIAN_TEMPORARY",
            "provenance",
        ),
        context="SPS/RUV explicit controls and policy",
    )

    config = SpsRuvBatchCorrectionConfig(
        control_site_set=ControlSiteSet.from_site_keys(("site_a", "site_c")),
        batch_column="batch",
        condition_columns=("condition",),
        replicate_column="replicate",
        missingness_policy=CorrectionMissingnessPolicy(),
        n_unwanted_factors=1,
        diagnostics_enabled=True,
        provenance_enabled=True,
    )

    preprocessing = DatasetPreprocessingConfig(batch_correction=config)
    assert preprocessing.batch_correction is config

    mask = ObservationMask(
        feature_ids=("site_a", "site_b", "site_c"),
        sample_ids=("sample_1", "sample_2", "sample_3", "sample_4"),
        originally_missing_cells=(("site_b", "sample_2"),),
    )
    missingness_policy = CorrectionMissingnessPolicy(
        temporary_imputation=TemporaryImputationPolicy(
            allowed=True,
            method=TemporaryImputationMethod.ROW_MEDIAN_TEMPORARY,
            method_parameters={"min_observed_values": 2},
        ),
        originally_missing_cells_tracked_by=(
            OriginallyMissingCellTracking.OBSERVATION_MASK
        ),
        observation_mask=mask,
    )
    assert missingness_policy.observation_mask is mask


def test_api_docs_protein_aware_preparation_boundary_is_documented() -> None:
    source = _read(DATASET_BUILD_DOC)

    _assert_python_call_keyword(
        source,
        "DatasetProteinAwarePreparationConfig",
        "policy",
        "prepare_model_inputs",
        context="protein-aware preparation configuration example",
    )
    _assert_python_attribute_chain(
        source,
        ("dataset", "protein_aware_preparation"),
        context="protein-aware preparation result contract",
    )
    _assert_python_attribute_chain(
        source,
        ("report", "site_eligibility_dataframe"),
        context="protein-aware preparation report contract",
    )
    _assert_negated_scope_statement(
        source,
        (
            "phosphosite matrix",
            "subtract total protein",
            "normalise intensities",
            "differential analysis",
            "MSstatsPTM-style inference",
        ),
        context="protein-aware preparation boundary",
    )


def test_api_docs_differential_request_example_is_constructible() -> None:
    dataset = _build_dataset()
    assert dataset.intensity_scale_state.kind.value == "log2"
    assert dataset.intensity_scale_state.is_established

    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_1",
                condition="control",
                biological_replicate_id="control_r1",
            ),
            SampleDesignRecord(
                sample_id="A_2",
                condition="control",
                biological_replicate_id="control_r2",
            ),
            SampleDesignRecord(
                sample_id="B_1",
                condition="treatment",
                biological_replicate_id="treatment_r1",
            ),
            SampleDesignRecord(
                sample_id="B_2",
                condition="treatment",
                biological_replicate_id="treatment_r2",
            ),
        )
    )
    contrasts = (
        Contrast(
            name="treatment_vs_control",
            numerator_condition="treatment",
            denominator_condition="control",
        ),
    )
    request = DifferentialAnalysisRequest(
        dataset=dataset,
        design=design,
        contrasts=contrasts,
    )

    assert request.design.samples[0].sample_id == "A_1"
    assert request.contrasts[0].name == "treatment_vs_control"
    assert request.config.minimum_condition_replicates == 2


def test_api_docs_kinase_request_example_is_constructible() -> None:
    dataset = _build_dataset()
    request = KinaseWorkflowRequest(
        dataset=dataset,
        references=ReferencePreset.AUTO,
        scoring_config=KinaseScoringConfig(
            min_substrates=2,
            include_diagnostic_scoring_tables=False,
            profile_missing_value_strategy="strict",
        ),
        prediction_config=KinasePredictionConfig(
            mode="deterministic_ranking",
            top_k=30,
            deterministic_max_selected_kinases=10,
        ),
        activity_config=KinaseActivityConfig(
            enabled=True,
            method="simplified_weighted_substrate_activity",
            threshold=0.6,
            min_substrates=3,
            top_n_substrates=20,
        ),
        site_sequence_conflict_policy="prefer_reference",
    )

    assert request.references is ReferencePreset.AUTO
    assert request.prediction_config.mode == "deterministic_ranking"
    assert request.site_sequence_conflict_policy == "prefer_reference"


def test_api_docs_signalome_request_example_is_constructible() -> None:
    kinase_result = _minimal_kinase_result()
    config = SignalomeConfig(
        scientific=SignalomeScientificConfig(
            substrate_support_cutoff=0.5,
            assignment_policy="cutoff_binary",
        ),
        clustering=SignalomeClusteringConfig(
            module_count=None,
            module_selection_primary_correlation_threshold=0.5,
            module_selection_fallback_correlation_threshold=0.1,
            module_selection_max_clusters=10,
            candidate_scoring_policy="full",
            clustering_engine="scipy_hierarchical",
        ),
        validation=SignalomeValidationConfig(
            score_preconditioning_policy="error_on_drop",
            allow_mixed_total_protein_quantitative_meaning=False,
        ),
        output=SignalomeOutputConfig(
            network_correlation_threshold=0.5,
            network_policy="signed",
        ),
        performance=SignalomePerformanceConfig(
            max_exact_tree_sites=2000,
            max_full_candidate_scoring_sites=2000,
        ),
    )
    request = SignalomeWorkflowRequest(
        kinase_result=kinase_result,
        config=config,
    )

    assert request.kinase_result.dataset.organism is Organism.RAT
    assert request.config.output.network_policy == "signed"


def test_api_docs_differential_import_route_uses_supported_public_path() -> None:
    source = _read(DIFFERENTIAL_DOC)

    _assert_python_imports(
        source,
        "phospy",
        ("DifferentialAnalysisWorkflow",),
        context="differential workflow import route",
    )
    _assert_python_imports(
        source,
        "phospy.api",
        ("DifferentialAnalysisRequest",),
        context="differential workflow import route",
    )
    _assert_python_run_call(
        source,
        "DifferentialAnalysisWorkflow",
        context="differential workflow import route",
    )
    _assert_python_imports_absent(
        source,
        "phospy",
        ("DifferentialAnalysis",),
        context="differential workflow import route",
    )
    _assert_python_imports_absent(
        source,
        "phospy.api",
        ("DifferentialAnalysis",),
        context="differential workflow import route",
    )
    _assert_negated_scope_statement(
        source,
        (
            "from phospy import DifferentialAnalysis",
            "from phospy.api import DifferentialAnalysis",
        ),
        context="unsupported differential import routes",
    )


def test_readme_and_differential_docs_keep_scientific_scope_contracts() -> None:
    readme_source = _read(README)
    differential_source = _read(DIFFERENTIAL_DOC)

    assert "minimum_condition_replicates=1" not in readme_source
    assert "minimum_condition_replicates=1" not in differential_source
    _assert_python_run_call(
        readme_source,
        "KinaseWorkflow",
        context="README scientific-scope workflow example",
    )
    _assert_python_call(
        readme_source,
        "KinaseWorkflowRequest",
        context="README scientific-scope workflow example",
    )
    _assert_python_constant(
        readme_source,
        "site_sequence",
        context="README scientific-scope workflow example",
    )
    _assert_python_call_keyword(
        readme_source,
        "KinaseWorkflowRequest",
        "references",
        "ReferencePreset.AUTO",
        context="README scientific-scope workflow example",
    )
    _assert_documented_terms(
        readme_source,
        ("`linear_residualize_batch`",),
        context="README batch-correction method contract",
    )
    _assert_statement_contains_all(
        readme_source,
        ("confounded", "batch/condition"),
        context="README batch confounding warning",
    )
    _assert_negated_scope_statement(
        readme_source,
        ("ComBat", "RUV", "removeBatchEffect", "mixed-effects"),
        context="README batch-correction scope",
    )
    for sample_id in (
        "control_rep1",
        "control_rep2",
        "treatment_rep1",
        "treatment_rep2",
    ):
        _assert_python_constant(
            differential_source,
            sample_id,
            context="differential minimal example sample design",
        )
    _assert_statement_contains_all(
        differential_source,
        ("DatasetIntensityTransformConfig", "log2"),
        context="differential log2 dataset preparation guidance",
    )


def test_public_workflow_docs_make_localisation_policy_explicit() -> None:
    readme_source = _read(README)
    dataset_source = _read(DATASET_BUILD_DOC)
    differential_source = _read(DIFFERENTIAL_DOC)
    kinase_source = _read(KINASE_DOC)
    signalome_source = _read(SIGNALOME_DOC)

    _assert_python_call_keyword(
        readme_source,
        "DatasetLocalisationConfig",
        "confidence_column",
        "localisation_confidence",
        context="README localisation policy example",
    )
    _assert_python_call_keyword(
        readme_source,
        "DatasetLocalisationConfig",
        "min_confidence",
        0.75,
        context="README localisation policy example",
    )

    _assert_python_call_keyword(
        dataset_source,
        "DatasetLocalisationConfig",
        "mode",
        "require_threshold",
        context="dataset localisation policy example",
    )
    _assert_python_call_keyword(
        dataset_source,
        "DatasetLocalisationConfig",
        "confidence_column",
        "localisation_confidence",
        context="dataset localisation policy example",
    )
    _assert_python_call_keyword(
        dataset_source,
        "DatasetLocalisationConfig",
        "min_confidence",
        0.75,
        context="dataset localisation policy example",
    )
    _assert_failure_boundary_statement(
        dataset_source,
        ("missing",),
        context="dataset localisation missing-value boundary",
    )
    _assert_failure_boundary_statement(
        dataset_source,
        ("invalid",),
        context="dataset localisation invalid-value boundary",
    )
    _assert_failure_boundary_statement(
        dataset_source,
        ("min_confidence",),
        context="dataset localisation threshold boundary",
    )

    _assert_python_call_keyword(
        differential_source,
        "DatasetLocalisationConfig",
        "confidence_column",
        "localisation_confidence",
        context="differential localisation prerequisite",
    )
    _assert_failure_boundary_statement(
        differential_source,
        ("low-confidence", "phosphosite"),
        context="differential localisation prerequisite",
    )

    _assert_python_call(
        kinase_source,
        "DatasetLocalisationConfig",
        context="kinase localisation prerequisite",
    )
    _assert_failure_boundary_statement(
        kinase_source,
        ("localisation", "missing"),
        context="kinase localisation prerequisite",
    )

    _assert_python_call(
        signalome_source,
        "DatasetLocalisationConfig",
        context="signalome localisation prerequisite",
    )
    _assert_failure_boundary_statement(
        signalome_source,
        ("localisation", "missing"),
        context="signalome localisation prerequisite",
    )


def test_kinase_docs_explain_reference_display_ambiguity_policy() -> None:
    source = _read(KINASE_DOC)

    _assert_documented_terms(
        source,
        (
            "reference_display_ambiguity_policy",
            '"error"',
            '"allow_with_diagnostics"',
        ),
        context="reference display ambiguity policy contract",
    )
    _assert_statement_contains_all(
        source,
        ("diagnostics", "matched", "site_key"),
        context="reference display ambiguity diagnostics",
    )
    _assert_negated_scope_statement(
        source,
        ("collapse", "duplicate display labels"),
        context="reference display ambiguity non-collapse policy",
    )
