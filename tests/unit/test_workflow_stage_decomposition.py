from __future__ import annotations

import ast
import inspect
from pathlib import Path

import pandas as pd
import pytest

from phospy.api import (
    Contrast,
    DifferentialAnalysisConfig,
    EnrichmentConfig,
    ExperimentalDesign,
    GeneSetCollection,
    SampleDesignRecord,
)
from phospy.api.configs import ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL
from phospy.errors import WorkflowBoundaryError
from phospy.science.design.models import PAIRED_DESIGN_POLICY_REJECT
from phospy.science.differential.executor import (
    DifferentialAnalysisExecutor as DifferentialComputationExecutor,
)
from phospy.science.differential.linear_model import decompose_differential_design
from phospy.science.differential.models import (
    DIFFERENTIAL_RESULT_STATUS_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_TESTED,
)
from phospy.science.differential.models import (
    DifferentialAnalysisRequest as DifferentialComputationRequest,
)
from phospy.science.enrichment.ora import OraConfig, OraEngine
from phospy.workflows.differential.design_assembly import (
    DifferentialExecutionDesignAssembler,
)
from phospy.workflows.differential.eligibility import (
    DifferentialComputationEligibilityResolver,
    DifferentialExecutionEligibilityResolution,
)
from phospy.workflows.differential.fitting import DifferentialModelFitter
from phospy.workflows.differential.models import (
    DifferentialFeatureEligibilityInputs,
    InterpretedDifferentialAnalysisRequest,
    ResolvedDifferentialExecutionConfig,
)
from phospy.workflows.differential.provenance import (
    DifferentialWorkflowProvenanceAssembler,
)
from phospy.workflows.differential.result_assembly import DifferentialResultAssembler
from phospy.workflows.enrichment.models import (
    EnrichmentIdentifierSemantics,
    InterpretedEnrichmentWorkflowRequest,
)
from phospy.workflows.enrichment.ora_execution import EnrichmentOraRunner
from phospy.workflows.enrichment.provenance import EnrichmentRunProvenanceAssembler
from phospy.workflows.enrichment.result_assembly import EnrichmentResultAssembler
from phospy.workflows.enrichment.set_filtering import EnrichmentSetExecutionPreparer
from tests.support.site_keys import protein_site_key_index, site_key_context_columns

ROOT = Path(__file__).resolve().parents[2]


def _design_frame() -> pd.DataFrame:
    frame = pd.DataFrame(
        {
            "A": [1.0, 1.0, 0.0, 0.0],
            "B": [0.0, 0.0, 1.0, 1.0],
        },
        index=pd.Index(["A_1", "A_2", "B_1", "B_2"], name="sample"),
    )
    frame.columns = pd.Index(frame.columns, name="coefficient")
    return frame


def _contrast_frame() -> pd.DataFrame:
    frame = pd.DataFrame(
        {"B_vs_A": [-1.0, 1.0]},
        index=pd.Index(["A", "B"], name="coefficient"),
    )
    frame.columns = pd.Index(frame.columns, name="contrast")
    return frame


def _differential_matrix() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "A_1": [1.0, 2.0],
            "A_2": [1.2, 2.2],
            "B_1": [2.0, 1.8],
            "B_2": [2.2, 2.0],
        },
        index=protein_site_key_index(
            protein_identifiers=["MAPK14", "GSK3B"],
            sites=["Y182", "S9"],
        ),
    )


def _computation_request(
    matrix: pd.DataFrame | None = None,
) -> DifferentialComputationRequest:
    return DifferentialComputationRequest(
        matrix=_differential_matrix() if matrix is None else matrix,
        design=_design_frame(),
        contrasts=_contrast_frame(),
    )


def _identity_metadata(index: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "site_key": index.astype(str).tolist(),
            "display_id": ["MAPK14;Y182;", "GSK3B;S9;"][: len(index)],
            **site_key_context_columns(index),
            "gene_symbol": ["MAPK14", "GSK3B"][: len(index)],
            "site": ["Y182", "S9"][: len(index)],
        },
        index=index.copy(),
    )


def _interpreted_differential_request(
    *,
    computation_request: DifferentialComputationRequest | None = None,
    feature_eligibility_inputs: DifferentialFeatureEligibilityInputs | None = None,
) -> InterpretedDifferentialAnalysisRequest:
    request = (
        _computation_request() if computation_request is None else computation_request
    )
    decomposition = request.design_decomposition
    config = DifferentialAnalysisConfig()
    return InterpretedDifferentialAnalysisRequest(
        computation_request=request,
        result_identity_metadata=_identity_metadata(request.matrix.index),
        config=config,
        execution_config=ResolvedDifferentialExecutionConfig(
            technical_replicate_policy=config.technical_replicate_policy,
            paired_design_policy=config.paired_design_policy,
            imputed_value_policy=config.imputed_value_policy,
            imputed_value_max_fraction=config.imputed_value_max_fraction,
            allow_design_subset=config.allow_design_subset,
            allow_suspicious_declared_input_scale=(
                config.allow_suspicious_declared_input_scale
            ),
            reliability_profile=config.reliability_profile,
            minimum_condition_replicates=config.minimum_condition_replicates,
            empirical_bayes=config.empirical_bayes,
            multiple_testing_method=config.multiple_testing.method,
        ),
        design_rank=int(decomposition.rank),
        residual_degrees_of_freedom=float(decomposition.residual_degrees_of_freedom),
        design_decomposition=decomposition,
        workflow_provenance={"input_intensity_scale": "log2"},
        feature_eligibility_inputs=feature_eligibility_inputs,
        normalisation_state="not_recorded",
    )


def _feature_eligibility(
    *,
    index: pd.Index,
    testable_feature_ids: tuple[str, ...],
) -> DifferentialFeatureEligibilityInputs:
    metadata = pd.DataFrame(
        {
            "site_key": index.astype(str).tolist(),
            DIFFERENTIAL_RESULT_STATUS_COLUMN: [
                DIFFERENTIAL_RESULT_STATUS_TESTED
                if str(label) in set(testable_feature_ids)
                else "withheld_other"
                for label in index
            ],
            DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN: [
                "tested"
                if str(label) in set(testable_feature_ids)
                else "withheld for unit test"
                for label in index
            ],
        },
        index=index.copy(),
    )
    return DifferentialFeatureEligibilityInputs(
        feature_metadata=metadata,
        result_status=metadata[DIFFERENTIAL_RESULT_STATUS_COLUMN],
        testable_feature_ids=testable_feature_ids,
        attach_to_result_tables=True,
    )


def _enrichment_request(
    *,
    config: EnrichmentConfig | None = None,
) -> InterpretedEnrichmentWorkflowRequest:
    set_collection = GeneSetCollection(
        sets={
            "PASS": ("AKT1", "MAPK1"),
            "SMALL": ("MTOR",),
        },
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        source_name="unit_test",
    )
    resolved_config = EnrichmentConfig() if config is None else config
    method_config = OraConfig(
        selected_outside_background_policy=(
            resolved_config.selected_outside_background_policy
        ),
        set_outside_background_policy=(
            resolved_config.set_member_outside_background_policy
        ),
        multiple_testing_correction=resolved_config.multiple_testing_correction,
    )
    return InterpretedEnrichmentWorkflowRequest(
        selected_identifiers=("AKT1", "MAPK1"),
        background_universe=("AKT1", "MAPK1", "MTOR"),
        set_collection=set_collection,
        method_config=method_config,
        identifier_semantics=EnrichmentIdentifierSemantics(
            identifier_column="gene_symbol",
            identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
            collection_kind=set_collection.collection_kind,
            analysis_level="gene",
        ),
        config=resolved_config,
        selected_identifier_source="selected_identifiers",
        method_metadata={"method": resolved_config.method},
        background_summary={
            "source": "explicit",
            "provided_identifier_count": 3,
            "universe_size": 3,
            "selected_identifier_count": 2,
            "selected_identifier_input_count": 2,
            "selected_identifier_source": "selected_identifiers",
        },
        set_collection_summary={
            "collection_kind": set_collection.collection_kind,
            "identifier_kind": set_collection.identifier_kind,
            "set_count": len(set_collection.enrichment_sets),
            "member_count": 3,
            "distinct_member_count": 3,
            "source_name": "unit_test",
            "source_version": None,
        },
        diagnostics={"interpreter": {"selected_identifiers_prepared": 2}},
        selected_identifier_input_count=2,
        background_identifier_input_count=3,
        selected_identifier_provenance=None,
        background_identifier_provenance=None,
    )


def test_differential_design_assembler_owns_execution_design_metadata() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="A_1", condition="A"),
            SampleDesignRecord(sample_id="A_2", condition="A"),
            SampleDesignRecord(sample_id="B_1", condition="B"),
            SampleDesignRecord(sample_id="B_2", condition="B"),
        )
    )
    contrasts = (
        Contrast(name="B_vs_A", numerator_condition="B", denominator_condition="A"),
    )
    design_frame = _design_frame()
    result = DifferentialExecutionDesignAssembler().run(
        design=design,
        contrasts=contrasts,
        design_aligned=design_frame,
        contrasts_aligned=_contrast_frame(),
        design_build_result=None,
        paired_design_policy=PAIRED_DESIGN_POLICY_REJECT,
        design_decomposition=decompose_differential_design(
            design_frame.to_numpy(dtype=float)
        ),
    )

    assert result.description == "condition-only fixed-effect design"
    assert result.sample_order == ("A_1", "A_2", "B_1", "B_2")
    assert result.condition_contrast_vectors[0].coefficients == (
        ("A", -1.0),
        ("B", 1.0),
    )


def test_differential_eligibility_resolver_filters_execution_matrix() -> None:
    request = _computation_request()
    decomposition = request.design_decomposition
    first_feature = str(request.matrix.index[0])
    eligibility = _feature_eligibility(
        index=request.matrix.index,
        testable_feature_ids=(first_feature,),
    )

    result = DifferentialComputationEligibilityResolver().run(
        _interpreted_differential_request(
            computation_request=request,
            feature_eligibility_inputs=eligibility,
        )
    )

    assert result.computation_request.matrix.index.tolist() == [first_feature]
    assert result.computation_request.design_decomposition is decomposition
    assert result.input_feature_ids == tuple(request.matrix.index.astype(str))
    assert result.multiple_testing_feature_ids == (first_feature,)


def test_differential_eligibility_failure_reports_executor_stage() -> None:
    request = _computation_request()
    eligibility = _feature_eligibility(
        index=request.matrix.index,
        testable_feature_ids=(),
    )

    with pytest.raises(
        WorkflowBoundaryError,
        match="differential.executor.feature_eligibility_testable_features",
    ):
        DifferentialComputationEligibilityResolver().run(
            _interpreted_differential_request(
                computation_request=request,
                feature_eligibility_inputs=eligibility,
            )
        )


def test_differential_model_fitter_delegates_to_science_executor() -> None:
    calls: list[object] = []
    expected = object()

    class Executor:
        def run(self, request: object) -> object:
            calls.append(request)
            return expected

    request = _computation_request()
    result = DifferentialModelFitter(
        computation_executor=Executor(),  # type: ignore[arg-type]
    ).run(request)

    assert result is expected
    assert calls == [request]


def test_differential_decomposition_identity_is_preserved_across_fit_and_diagnostics() -> (
    None
):
    interpreted = _interpreted_differential_request()
    eligibility = DifferentialComputationEligibilityResolver().run(interpreted)
    computation_result = DifferentialComputationExecutor().run(
        eligibility.computation_request
    )

    assert interpreted.computation_request.design_decomposition is (
        interpreted.design_decomposition
    )
    assert eligibility.computation_request.design_decomposition is (
        interpreted.design_decomposition
    )
    assert computation_result.design_decomposition is interpreted.design_decomposition

    public_result = DifferentialResultAssembler().run(
        request=interpreted,
        computation_result=computation_result,
        eligibility=eligibility,
        workflow_provenance={"row_attrition_metrics": {}},
    )

    assert public_result.diagnostics.rank == interpreted.design_decomposition.rank
    assert public_result.diagnostics.singular_values == (
        interpreted.design_decomposition.singular_values
    )
    assert public_result.diagnostics.condition_number == pytest.approx(
        interpreted.design_decomposition.condition_number
    )


def test_differential_provenance_assembler_records_row_attrition() -> None:
    provenance = DifferentialWorkflowProvenanceAssembler().run(
        workflow_provenance={"input_intensity_scale": "log2"},
        input_feature_ids=("site1", "site2", "site3"),
        model_fit_feature_ids=("site1", "site2"),
        failed_model_fit_feature_ids=("site2",),
        multiple_testing_feature_ids=("site1",),
    )

    assert provenance["row_attrition_metrics"] == {
        "input_sites": 3,
        "sites_retained_for_model_fitting": 2,
        "sites_excluded_before_testing": 1,
        "sites_with_failed_model_fit": 1,
        "sites_included_in_multiple_testing_family": 1,
    }
    records = provenance["row_attrition"]["records"]  # type: ignore[index]
    assert [record["stage"] for record in records] == [
        "differential_feature_eligibility",
        "differential_model_fit",
    ]


def test_differential_result_assembler_owns_public_result_shape() -> None:
    request = _interpreted_differential_request()
    computation_result = DifferentialComputationExecutor().run(
        request.computation_request
    )
    eligibility = DifferentialExecutionEligibilityResolution(
        computation_request=request.computation_request,
        feature_eligibility_inputs=None,
        input_feature_ids=tuple(request.computation_request.matrix.index.astype(str)),
        model_fit_feature_ids=tuple(
            request.computation_request.matrix.index.astype(str)
        ),
        failed_model_fit_feature_ids=(),
        multiple_testing_feature_ids=tuple(
            request.computation_request.matrix.index.astype(str)
        ),
    )

    result = DifferentialResultAssembler().run(
        request=request,
        computation_result=computation_result,
        eligibility=eligibility,
        workflow_provenance={"row_attrition_metrics": {}},
    )

    table = result.table_for("B_vs_A")
    assert "display_id" in table.columns
    assert "P.Value" in table.columns
    assert result.diagnostics.model_type == "moderated_ols_fixed_effect"


def test_enrichment_set_preparer_owns_min_max_set_filtering() -> None:
    request = _enrichment_request(config=EnrichmentConfig(min_set_size=2))

    result = EnrichmentSetExecutionPreparer().run(request)

    assert result.tested_set_count == 1
    assert result.dropped_sets[0].set_id == "SMALL"
    assert result.dropped_sets[0].reason == "below_min_set_size"


def test_enrichment_ora_runner_handles_empty_filtered_collection() -> None:
    request = _enrichment_request(config=EnrichmentConfig(min_set_size=10))
    set_filter = EnrichmentSetExecutionPreparer().run(request)

    result = EnrichmentOraRunner().run(
        request=request,
        set_size_filter_result=set_filter,
    )

    assert result.records == ()
    assert result.selected_identifiers == ("AKT1", "MAPK1")


def test_enrichment_result_assembler_owns_public_records_and_diagnostics() -> None:
    request = _enrichment_request()
    set_filter = EnrichmentSetExecutionPreparer().run(request)
    ora_result = OraEngine().run(
        selected_identifiers=request.selected_identifiers,
        background_universe=request.background_universe,
        enrichment_sets=request.set_collection,
        config=request.method_config,
    )

    result = EnrichmentResultAssembler().run(
        request=request,
        ora_result=ora_result,
        set_size_filter_result=set_filter,
    )

    assert result.records[0].term_id == "PASS"
    assert result.diagnostics["ora"]["record_count"] == 2
    assert result.background_summary["source"] == "explicit"


def test_enrichment_provenance_assembler_owns_run_fingerprints() -> None:
    request = _enrichment_request(config=EnrichmentConfig(min_set_size=2))
    set_filter = EnrichmentSetExecutionPreparer().run(request)
    ora_result = EnrichmentOraRunner().run(
        request=request,
        set_size_filter_result=set_filter,
    )
    public_result = EnrichmentResultAssembler().run(
        request=request,
        ora_result=ora_result,
        set_size_filter_result=set_filter,
    )

    provenance = EnrichmentRunProvenanceAssembler().run(
        request=request,
        ora_result=ora_result,
        result_table=public_result.table,
        set_size_filter_result=set_filter,
    )

    assert provenance.workflow_name == "enrichment"
    assert provenance.workflow_parameters["number_of_tests"] == 1
    assert provenance.workflow_parameters["set_size_filter"] == {
        "min_set_size": 2,
        "max_set_size": None,
        "applied_after_background_intersection": True,
    }


def test_ora_kernel_has_no_workflow_or_output_serialization_imports() -> None:
    import phospy.science.enrichment.ora as ora

    tree = ast.parse(inspect.getsource(ora))
    imported_modules: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported_modules.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module is not None:
            imported_modules.add(node.module)

    assert "pandas" not in imported_modules
    assert "phospy.provenance" not in imported_modules
    assert "phospy.contracts.results" not in imported_modules
    assert not any(module.startswith("phospy.workflows") for module in imported_modules)


def test_workflow_executors_no_longer_build_public_structures_directly() -> None:
    from phospy.workflows.differential.executor import DifferentialAnalysisExecutor
    from phospy.workflows.enrichment.executor import EnrichmentWorkflowExecutor

    differential_source = inspect.getsource(DifferentialAnalysisExecutor)
    enrichment_source = inspect.getsource(EnrichmentWorkflowExecutor)

    assert "DifferentialModelDiagnostics(" not in differential_source
    assert "DifferentialAnalysisResult._from_owned" not in differential_source
    assert "EnrichmentResultRecord(" not in enrichment_source
    assert "RunProvenance(" not in enrichment_source
    assert "fingerprint_table(" not in enrichment_source


def test_differential_execution_stages_do_not_rebuild_design_decomposition() -> None:
    from phospy.workflows.differential import eligibility, fitting, result_assembly
    from phospy.workflows.differential.executor import DifferentialAnalysisExecutor

    sources = (
        inspect.getsource(DifferentialComputationExecutor),
        inspect.getsource(DifferentialAnalysisExecutor),
        inspect.getsource(eligibility),
        inspect.getsource(fitting),
        inspect.getsource(result_assembly),
    )

    for source in sources:
        assert "decompose_differential_design" not in source


def test_decomposed_workflow_modules_do_not_create_import_cycles() -> None:
    package_dirs = (
        ROOT / "src" / "phospy" / "workflows" / "differential",
        ROOT / "src" / "phospy" / "workflows" / "enrichment",
    )
    module_paths = [
        path for package_dir in package_dirs for path in package_dir.glob("*.py")
    ]
    modules = {
        "phospy." + ".".join(path.relative_to(ROOT / "src").with_suffix("").parts): path
        for path in module_paths
    }
    graph = {module: set[str]() for module in modules}
    for module, path in modules.items():
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for imported_module in _internal_imports(tree):
            if imported_module in modules:
                graph[module].add(imported_module)
                continue
            prefix = _nearest_known_module(imported_module, modules)
            if prefix is not None:
                graph[module].add(prefix)

    cycle = _find_cycle(graph)

    assert cycle is None


def _internal_imports(tree: ast.AST) -> set[str]:
    imports: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module is not None:
            imports.add(node.module)
            imports.update(f"{node.module}.{alias.name}" for alias in node.names)
    return {
        imported
        for imported in imports
        if imported.startswith("phospy.workflows.differential")
        or imported.startswith("phospy.workflows.enrichment")
    }


def _nearest_known_module(
    imported_module: str,
    modules: dict[str, Path],
) -> str | None:
    parts = imported_module.split(".")
    for end in range(len(parts), 0, -1):
        candidate = ".".join(parts[:end])
        if candidate in modules:
            return candidate
    return None


def _find_cycle(graph: dict[str, set[str]]) -> list[str] | None:
    visiting: set[str] = set()
    visited: set[str] = set()
    stack: list[str] = []

    def visit(module: str) -> list[str] | None:
        if module in visited:
            return None
        if module in visiting:
            start = stack.index(module)
            return stack[start:] + [module]
        visiting.add(module)
        stack.append(module)
        for dependency in graph[module]:
            cycle = visit(dependency)
            if cycle is not None:
                return cycle
        stack.pop()
        visiting.remove(module)
        visited.add(module)
        return None

    for module in graph:
        cycle = visit(module)
        if cycle is not None:
            return cycle
    return None
