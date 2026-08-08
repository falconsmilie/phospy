from __future__ import annotations

import ast
import inspect
from typing import Any, cast

import pandas as pd
import pytest

import phospy.science.datasets.preprocessing.pipeline as pipeline_module
import phospy.science.datasets.preprocessing.stage_registry as stage_registry_module
from phospy.advanced.configs import (
    DatasetComparisonBuildingConfig,
    DatasetIntensityTransformConfig,
    DatasetNormalisationConfig,
    DatasetSiteMatrixConfig,
)
from phospy.api.configs import (
    DatasetLocalisationConfig,
    DatasetPreprocessingConfig,
)
from phospy.contracts.configs import DATASET_LOCALISATION_MODE_ALLOW_MISSING_WITH_WAIVER
from phospy.errors.build import DatasetBuildError
from phospy.provenance.models import (
    PREPROCESSING_EXTERNAL_NONDETERMINISM_CAVEAT_CODE,
)
from phospy.science.datasets.builders.preprocessing import DatasetPreprocessor
from phospy.science.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION,
    DATASET_PREPROCESSING_STAGE_COMPARISONS,
    DATASET_PREPROCESSING_STAGE_GROUP_COVERAGE_FILTER,
    DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
    DATASET_PREPROCESSING_STAGE_LOCALISATION,
    DATASET_PREPROCESSING_STAGE_MISSING_DATA,
    DATASET_PREPROCESSING_STAGE_NORMALISATION,
    DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
    DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
    DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
    PreprocessingPlan,
    PreprocessingStageResult,
    PreprocessingState,
)
from phospy.science.datasets.preprocessing.pipeline import PreprocessingPipeline
from phospy.science.datasets.preprocessing.stage_contract import (
    DeterminismKind,
    PreprocessingStageFactoryContext,
    PreprocessingStagePlanInterpreter,
    PreprocessingStageRegistrationMetadata,
)
from phospy.science.datasets.preprocessing.stage_registry import (
    PreprocessingStageMetadata,
    build_registered_preprocessing_stage_instances,
    get_preprocessing_stage_metadata,
    list_registered_preprocessing_stages,
    resolve_builder_provenance_stage_order,
    resolve_registered_preprocessing_stages,
)
from phospy.science.datasets.preprocessing.stages.batch_correction import (
    BATCH_CORRECTION_STAGE_CONTRACT,
)
from phospy.science.transformations.quantitative_contracts import (
    QuantitativeOperationContract,
    preserve_quantitative_contract,
)
from tests.support.site_keys import site_key_index_from_display_ids

_DISPLAY_IDS = ["MAPK14;Y182;", "AKT1;T308;"]


def _preserve_quantitative_contract() -> QuantitativeOperationContract:
    return preserve_quantitative_contract()


def _phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "sample_b": [2.0, 4.0],
        },
        index=pd.Index(_DISPLAY_IDS, name="site_id"),
    )


def _site_metadata(index: pd.Index) -> pd.DataFrame:
    site_keys = site_key_index_from_display_ids(
        _DISPLAY_IDS,
        protein_namespace="gene_symbol",
    )
    return pd.DataFrame(
        {
            "site_key": site_keys.astype(str).tolist(),
            "display_id": _DISPLAY_IDS,
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": ["SEQ_A", "SEQ_R"],
            "localisation_confidence": [0.95, 0.91],
        },
        index=index.copy(),
    )


def _sample_metadata(columns: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {"comparison_group": ["group_a", "group_b"]},
        index=columns.copy(),
    )


class _NoOpPreQuantitativeContractStage:
    def validate_before_quantitative_contract(
        self,
        state: PreprocessingState,
    ) -> None:
        del state
        return None


def _plan_with_multiple_stages() -> PreprocessingPlan:
    return PreprocessingPlan.from_config(
        DatasetPreprocessingConfig(
            intensity_transform=DatasetIntensityTransformConfig(
                policy="log2",
                pseudocount=1.0,
            ),
            normalisation=DatasetNormalisationConfig(policy="median_center"),
            site_matrix=DatasetSiteMatrixConfig(policy="build_from_metadata"),
            comparisons=DatasetComparisonBuildingConfig(policy="sample_metadata_pairs"),
            localisation=DatasetLocalisationConfig(
                mode=DATASET_LOCALISATION_MODE_ALLOW_MISSING_WITH_WAIVER,
                waiver_reason="test waiver",
            ),
        )
    )


def test_every_preprocessing_stage_has_registry_metadata() -> None:
    expected = {
        DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
        DATASET_PREPROCESSING_STAGE_LOCALISATION,
        DATASET_PREPROCESSING_STAGE_GROUP_COVERAGE_FILTER,
        DATASET_PREPROCESSING_STAGE_MISSING_DATA,
        DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
        DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION,
        DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
        DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
        DATASET_PREPROCESSING_STAGE_NORMALISATION,
        DATASET_PREPROCESSING_STAGE_COMPARISONS,
    }
    observed = {item.stage_key for item in list_registered_preprocessing_stages()}
    assert observed == expected


def test_registered_stage_keys_are_unique() -> None:
    registered = list_registered_preprocessing_stages()
    stage_keys = [metadata.stage_key for metadata in registered]
    assert len(stage_keys) == len(set(stage_keys))


def test_registered_stage_order_is_stable() -> None:
    observed = tuple(
        metadata.stage_key for metadata in list_registered_preprocessing_stages()
    )
    assert observed == (
        DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
        DATASET_PREPROCESSING_STAGE_LOCALISATION,
        DATASET_PREPROCESSING_STAGE_GROUP_COVERAGE_FILTER,
        DATASET_PREPROCESSING_STAGE_MISSING_DATA,
        DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
        DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION,
        DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
        DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
        DATASET_PREPROCESSING_STAGE_NORMALISATION,
        DATASET_PREPROCESSING_STAGE_COMPARISONS,
    )


def test_every_registered_stage_has_required_metadata_contract_fields() -> None:
    plan = _plan_with_multiple_stages()
    for metadata in list_registered_preprocessing_stages():
        assert metadata.stage_key
        assert metadata.display_label
        assert metadata.provenance_stage_key
        assert callable(metadata.operation_name)
        assert callable(metadata.serialize_parameters)
        assert isinstance(metadata.serialize_parameters(plan), dict)
        assert isinstance(metadata.resolve_determinism_kind(plan), DeterminismKind)
        assert isinstance(metadata.diagnostics_metadata, dict)
        assert metadata.diagnostics_metadata
        contract = metadata.resolve_quantitative_contract(plan)
        assert isinstance(contract, QuantitativeOperationContract)
        assert contract.accepted_input_scale_kinds
        assert contract.accepted_quantitative_meanings
        assert contract.output_scale_transition.output_scale_by_input
        assert contract.output_meaning_transition.output_meaning_by_input
        assert contract.required_evidence


def test_stage_contract_separates_static_metadata_from_plan_interpretation() -> None:
    contract = get_preprocessing_stage_metadata(
        DATASET_PREPROCESSING_STAGE_MISSING_DATA
    )

    static_metadata = contract.registration_metadata
    plan_interpreter = contract.plan_interpreter

    assert isinstance(static_metadata, PreprocessingStageRegistrationMetadata)
    assert static_metadata.stage_key == DATASET_PREPROCESSING_STAGE_MISSING_DATA
    assert static_metadata.consumed_input_tables == contract.consumed_input_tables
    assert isinstance(plan_interpreter, PreprocessingStagePlanInterpreter)
    assert plan_interpreter.operation_name is contract.operation_name


def test_duplicate_override_stage_keys_fail_registry_resolution() -> None:
    duplicate_key = DATASET_PREPROCESSING_STAGE_MISSING_DATA
    duplicate_entries = (
        PreprocessingStageMetadata(
            stage_key=duplicate_key,
            display_label="missing_data_override_a",
            operation_name=lambda _plan: "forbid",
            serialize_parameters=lambda _plan: {},
            consumed_input_tables=("dataset.phospho",),
            produced_output_tables=("dataset.phospho",),
            quantitative_contract=_preserve_quantitative_contract(),
            diagnostics_metadata={"known_diagnostics_fields": ("policy",)},
        ),
        PreprocessingStageMetadata(
            stage_key=duplicate_key,
            display_label="missing_data_override_b",
            operation_name=lambda _plan: "forbid",
            serialize_parameters=lambda _plan: {},
            consumed_input_tables=("dataset.phospho",),
            produced_output_tables=("dataset.phospho",),
            quantitative_contract=_preserve_quantitative_contract(),
            diagnostics_metadata={"known_diagnostics_fields": ("policy",)},
        ),
    )

    with pytest.raises(DatasetBuildError, match="duplicate stage key"):
        resolve_registered_preprocessing_stages(duplicate_entries)


def test_registry_rejects_quantitative_stage_without_semantic_contract() -> None:
    missing_contract = PreprocessingStageMetadata(
        stage_key="custom_quantitative_stage",
        display_label="custom_quantitative_stage",
        operation_name=lambda _plan: "custom_quantitative_stage",
        serialize_parameters=lambda _plan: {},
        consumed_input_tables=("dataset.phospho",),
        produced_output_tables=("dataset.phospho",),
        diagnostics_metadata={"known_diagnostics_fields": ("policy",)},
    )

    with pytest.raises(
        DatasetBuildError, match="missing quantitative semantic contract"
    ):
        resolve_registered_preprocessing_stages((missing_contract,))


def test_stage_metadata_rejects_unknown_consumed_table_key() -> None:
    with pytest.raises(DatasetBuildError, match="unknown table key"):
        PreprocessingStageMetadata(
            stage_key="custom_stage",
            display_label="custom_stage",
            operation_name=lambda _plan: "custom_stage",
            serialize_parameters=lambda _plan: {},
            consumed_input_tables=("dataset.sampl_metadata",),
            produced_output_tables=("dataset.phospho",),
            quantitative_contract=_preserve_quantitative_contract(),
            diagnostics_metadata={"known_diagnostics_fields": ("policy",)},
        )


def test_stage_metadata_rejects_unknown_produced_table_key() -> None:
    with pytest.raises(DatasetBuildError, match="unknown table key"):
        PreprocessingStageMetadata(
            stage_key="custom_stage",
            display_label="custom_stage",
            operation_name=lambda _plan: "custom_stage",
            serialize_parameters=lambda _plan: {},
            consumed_input_tables=("dataset.phospho",),
            produced_output_tables=("report.row_audt",),
            quantitative_contract=_preserve_quantitative_contract(),
            diagnostics_metadata={"known_diagnostics_fields": ("policy",)},
        )


def test_known_optional_missing_table_is_skipped_in_trace_fingerprints() -> None:
    class _OptionalSampleMetadataStage(_NoOpPreQuantitativeContractStage):
        stage_key = "optional_sample_metadata_stage"

        def run(self, state: PreprocessingState) -> PreprocessingStageResult:
            return PreprocessingStageResult(state=state)

    phospho = _phospho()
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=_site_metadata(phospho.index),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(stage_order=("optional_sample_metadata_stage",)),
    )
    metadata = PreprocessingStageMetadata(
        stage_key="optional_sample_metadata_stage",
        display_label="optional_sample_metadata_stage",
        operation_name=lambda _plan: "optional_sample_metadata_stage",
        serialize_parameters=lambda _plan: {},
        consumed_input_tables=("dataset.sample_metadata",),
        produced_output_tables=("dataset.phospho",),
        quantitative_contract=_preserve_quantitative_contract(),
        diagnostics_metadata={"known_diagnostics_fields": ("policy",)},
    )

    _, trace = PreprocessingPipeline(
        stage_registry=(_OptionalSampleMetadataStage(),),
        stage_contract_registry=(metadata,),
    ).run_with_trace(state)

    assert len(trace) == 1
    assert trace[0].consumed_input_tables == ()
    assert tuple(item.name for item in trace[0].produced_output_tables) == (
        "dataset.phospho",
    )


def test_pipeline_trace_parameters_and_operation_come_from_registry() -> None:
    phospho = _phospho()
    plan = _plan_with_multiple_stages()
    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=_site_metadata(phospho.index),
        sample_metadata=_sample_metadata(phospho.columns),
        total=None,
        plan=plan,
    )
    trace = preprocessed.preprocessing_trace or ()
    assert trace

    for entry in trace:
        metadata = get_preprocessing_stage_metadata(entry.stage)
        assert entry.operation == metadata.operation_name(plan)
        assert dict(entry.parameters) == metadata.serialize_parameters(plan)


def test_pipeline_and_builder_use_same_stage_labels_and_operations() -> None:
    phospho = _phospho()
    plan = _plan_with_multiple_stages()
    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=_site_metadata(phospho.index),
        sample_metadata=_sample_metadata(phospho.columns),
        total=None,
        plan=plan,
    )
    trace = preprocessed.preprocessing_trace or ()
    operations = preprocessed.preprocessing_operations
    assert operations is not None

    for entry in trace:
        metadata = get_preprocessing_stage_metadata(entry.stage)
        matching = operations.loc[operations.loc[:, "stage"] == metadata.display_label]
        assert matching.shape[0] == 1
        assert str(matching.iloc[0]["operation"]) == entry.operation


def test_provenance_operations_are_resolved_from_registry_metadata() -> None:
    phospho = _phospho()
    plan = _plan_with_multiple_stages()
    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=_site_metadata(phospho.index),
        sample_metadata=_sample_metadata(phospho.columns),
        total=None,
        plan=plan,
    )
    operations = preprocessed.preprocessing_operations
    assert operations is not None

    expected_metadata = {
        metadata.display_label: metadata
        for metadata in resolve_builder_provenance_stage_order(plan)
    }
    for record in operations.to_dict(orient="records"):
        label = str(record["stage"])
        metadata = expected_metadata.get(label)
        assert metadata is not None
        assert str(record["operation"]) == metadata.operation_name(plan)


def test_unknown_stage_metadata_fails_with_clear_error() -> None:
    class _UnknownStage(_NoOpPreQuantitativeContractStage):
        stage_key = "unknown_stage"

        def run(self, state: PreprocessingState) -> PreprocessingStageResult:
            return PreprocessingStageResult(state=state)

    phospho = _phospho()
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=_site_metadata(phospho.index),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(stage_order=("unknown_stage",)),
    )

    with pytest.raises(
        DatasetBuildError,
        match="without registered stage contract",
    ):
        PreprocessingPipeline(stage_registry=(_UnknownStage(),)).run_with_trace(state)


def test_custom_stage_valid_only_under_non_default_plan() -> None:
    observed_plans: list[tuple[str, ...]] = []

    class NonDefaultOnlyStage(_NoOpPreQuantitativeContractStage):
        stage_key = "non_default_only_stage"

        def run(self, state: PreprocessingState) -> PreprocessingStageResult:
            return PreprocessingStageResult(
                state=state,
                diagnostics={"diagnostics": {"policy": "non_default_only"}},
            )

    def _require_actual_custom_plan(plan: PreprocessingPlan) -> None:
        observed_plans.append(tuple(plan.stage_order))
        if tuple(plan.stage_order) != ("non_default_only_stage",):
            raise AssertionError("custom stage was interpreted with an unrelated plan")

    def _operation(plan: PreprocessingPlan) -> str:
        _require_actual_custom_plan(plan)
        return "non_default_only"

    def _parameters(plan: PreprocessingPlan) -> dict[str, object]:
        _require_actual_custom_plan(plan)
        return {"stage_order": list(plan.stage_order)}

    def _contract(plan: PreprocessingPlan) -> QuantitativeOperationContract:
        _require_actual_custom_plan(plan)
        return _preserve_quantitative_contract()

    def _build_stage(
        _context: PreprocessingStageFactoryContext,
    ) -> NonDefaultOnlyStage:
        return NonDefaultOnlyStage()

    contract = PreprocessingStageMetadata(
        stage_key="non_default_only_stage",
        display_label="non_default_only_stage",
        operation_name=_operation,
        serialize_parameters=_parameters,
        consumed_input_tables=("dataset.phospho",),
        produced_output_tables=("dataset.phospho",),
        stage_factory=_build_stage,
        quantitative_contract=_contract,
        diagnostics_metadata={"known_diagnostics_fields": ("policy",)},
    )

    _, trace = PreprocessingPipeline(
        stage_contract_registry=(contract,)
    ).run_with_trace(_fake_stage_state("non_default_only_stage"))

    assert trace[0].operation == "non_default_only"
    assert observed_plans
    assert set(observed_plans) == {("non_default_only_stage",)}


def test_stage_pre_quantitative_contract_validator_is_invoked_explicitly() -> None:
    calls: list[str] = []

    class ExplicitValidatorStage:
        stage_key = "explicit_validator_stage"

        def validate_before_quantitative_contract(
            self,
            state: PreprocessingState,
        ) -> None:
            del state
            calls.append("validate_before_quantitative_contract")

        def run(self, state: PreprocessingState) -> PreprocessingStageResult:
            calls.append("run")
            return PreprocessingStageResult(state=state)

    def _build_stage(
        _context: PreprocessingStageFactoryContext,
    ) -> ExplicitValidatorStage:
        return ExplicitValidatorStage()

    contract = PreprocessingStageMetadata(
        stage_key="explicit_validator_stage",
        display_label="explicit_validator_stage",
        operation_name=lambda _plan: "explicit_validator_stage",
        serialize_parameters=lambda _plan: {},
        consumed_input_tables=("dataset.phospho",),
        produced_output_tables=("dataset.phospho",),
        stage_factory=_build_stage,
        quantitative_contract=_preserve_quantitative_contract(),
        diagnostics_metadata={"known_diagnostics_fields": ("policy",)},
    )

    PreprocessingPipeline(stage_contract_registry=(contract,)).run_with_trace(
        _fake_stage_state("explicit_validator_stage")
    )

    assert calls == ["validate_before_quantitative_contract", "run"]


def test_missing_required_lifecycle_method_fails_at_composition() -> None:
    class MissingPreContractValidatorStage:
        stage_key = "missing_pre_contract_validator_stage"

        def run(self, state: PreprocessingState) -> PreprocessingStageResult:
            return PreprocessingStageResult(state=state)

    def _build_stage(
        _context: PreprocessingStageFactoryContext,
    ) -> MissingPreContractValidatorStage:
        return MissingPreContractValidatorStage()

    contract = PreprocessingStageMetadata(
        stage_key="missing_pre_contract_validator_stage",
        display_label="missing_pre_contract_validator_stage",
        operation_name=lambda _plan: "missing_pre_contract_validator_stage",
        serialize_parameters=lambda _plan: {},
        consumed_input_tables=("dataset.phospho",),
        produced_output_tables=("dataset.phospho",),
        stage_factory=_build_stage,
        quantitative_contract=_preserve_quantitative_contract(),
        diagnostics_metadata={"known_diagnostics_fields": ("policy",)},
    )

    with pytest.raises(
        DatasetBuildError,
        match="validate_before_quantitative_contract.*run",
    ):
        PreprocessingPipeline(stage_contract_registry=(contract,))


def test_plan_table_dependency_validation_uses_actual_stage_order() -> None:
    calls: list[str] = []

    class RowAuditConsumerStage(_NoOpPreQuantitativeContractStage):
        stage_key = "row_audit_consumer_stage"

        def run(self, state: PreprocessingState) -> PreprocessingStageResult:
            calls.append("run")
            return PreprocessingStageResult(state=state)

    def _build_stage(
        _context: PreprocessingStageFactoryContext,
    ) -> RowAuditConsumerStage:
        return RowAuditConsumerStage()

    contract = PreprocessingStageMetadata(
        stage_key="row_audit_consumer_stage",
        display_label="row_audit_consumer_stage",
        operation_name=lambda _plan: "row_audit_consumer_stage",
        serialize_parameters=lambda _plan: {},
        consumed_input_tables=("report.row_audit",),
        produced_output_tables=("dataset.phospho",),
        stage_factory=_build_stage,
        quantitative_contract=_preserve_quantitative_contract(),
        diagnostics_metadata={"known_diagnostics_fields": ("policy",)},
    )

    with pytest.raises(DatasetBuildError, match="invalid stage table dependencies"):
        PreprocessingPipeline(stage_contract_registry=(contract,)).run_with_trace(
            _fake_stage_state("row_audit_consumer_stage")
        )

    assert calls == []


def test_pipeline_uses_explicit_pre_contract_hook_without_hidden_getattr() -> None:
    source = inspect.getsource(pipeline_module)

    assert "getattr(" not in source
    assert "validate_before_quantitative_contract" in source


def test_registered_stage_factories_expose_run_method() -> None:
    context = PreprocessingStageFactoryContext()
    for metadata in list_registered_preprocessing_stages():
        if metadata.stage_factory is None:
            continue
        stage = metadata.stage_factory(context)
        assert callable(stage.validate_before_quantitative_contract)
        assert callable(stage.run)


def test_dependency_free_stage_factory_receives_uniform_context() -> None:
    observed_contexts: list[PreprocessingStageFactoryContext] = []

    class ContextFreeStage(_NoOpPreQuantitativeContractStage):
        stage_key = "context_free_stage"

        def run(self, state: PreprocessingState) -> PreprocessingStageResult:
            return PreprocessingStageResult(state=state)

    def _build_context_free_stage(
        context: PreprocessingStageFactoryContext,
    ) -> ContextFreeStage:
        observed_contexts.append(context)
        return ContextFreeStage()

    contract = PreprocessingStageMetadata(
        stage_key="context_free_stage",
        display_label="context_free_stage",
        operation_name=lambda _plan: "context_free",
        serialize_parameters=lambda _plan: {},
        consumed_input_tables=("dataset.phospho",),
        produced_output_tables=("dataset.phospho",),
        stage_factory=_build_context_free_stage,
        quantitative_contract=_preserve_quantitative_contract(),
        diagnostics_metadata={"known_diagnostics_fields": ("policy",)},
    )
    context = PreprocessingStageFactoryContext()

    instances = build_registered_preprocessing_stage_instances(
        (contract,),
        context=context,
    )

    assert len(instances) == 1
    assert observed_contexts == [context]
    assert instances[0].stage_key == "context_free_stage"


def test_dependency_bearing_stage_receives_collaborator_from_uniform_context() -> None:
    class FakeCollaborator:
        marker = "from_context"

    class CollaboratorStage(_NoOpPreQuantitativeContractStage):
        stage_key = "collaborator_stage"

        def __init__(self, collaborator: object) -> None:
            self.collaborator = collaborator

        def run(self, state: PreprocessingState) -> PreprocessingStageResult:
            return PreprocessingStageResult(state=state)

    def _build_collaborator_stage(
        context: PreprocessingStageFactoryContext,
    ) -> CollaboratorStage:
        return CollaboratorStage(context.batch_correction_runner)

    contract = PreprocessingStageMetadata(
        stage_key="collaborator_stage",
        display_label="collaborator_stage",
        operation_name=lambda _plan: "collaborator_stage",
        serialize_parameters=lambda _plan: {},
        consumed_input_tables=("dataset.phospho",),
        produced_output_tables=("dataset.phospho",),
        stage_factory=_build_collaborator_stage,
        quantitative_contract=_preserve_quantitative_contract(),
        diagnostics_metadata={"known_diagnostics_fields": ("policy",)},
    )
    collaborator = FakeCollaborator()
    context = PreprocessingStageFactoryContext(
        batch_correction_runner=cast(Any, collaborator),
    )

    (stage,) = build_registered_preprocessing_stage_instances(
        (contract,),
        context=context,
    )

    assert cast(CollaboratorStage, stage).collaborator is collaborator


def test_pipeline_builds_dependency_bearing_custom_stage_without_registry_branch() -> (
    None
):
    class FakeCollaborator:
        marker = "future_collaborator"

    class CollaboratorStage(_NoOpPreQuantitativeContractStage):
        stage_key = "future_dependency_stage"

        def __init__(self, collaborator: FakeCollaborator) -> None:
            self._collaborator = collaborator

        def run(self, state: PreprocessingState) -> PreprocessingStageResult:
            return PreprocessingStageResult(
                state=state,
                diagnostics={
                    "diagnostics": {"collaborator": self._collaborator.marker}
                },
            )

    def _build_collaborator_stage(
        context: PreprocessingStageFactoryContext,
    ) -> CollaboratorStage:
        return CollaboratorStage(
            cast(FakeCollaborator, context.batch_correction_runner)
        )

    contract = PreprocessingStageMetadata(
        stage_key="future_dependency_stage",
        display_label="future_dependency_stage",
        operation_name=lambda _plan: "future_dependency_stage",
        serialize_parameters=lambda _plan: {},
        consumed_input_tables=("dataset.phospho",),
        produced_output_tables=("dataset.phospho",),
        stage_factory=_build_collaborator_stage,
        quantitative_contract=_preserve_quantitative_contract(),
        diagnostics_metadata={"known_diagnostics_fields": ("collaborator",)},
    )
    phospho = _phospho()
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=_site_metadata(phospho.index),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(stage_order=("future_dependency_stage",)),
    )

    _, trace = PreprocessingPipeline(
        stage_contract_registry=(contract,),
        batch_correction_runner=cast(Any, FakeCollaborator()),
    ).run_with_trace(state)

    assert trace[0].diagnostics["collaborator"] == "future_collaborator"


def test_batch_correction_constructs_and_runs_through_uniform_context_path() -> None:
    phospho = _phospho()
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=_site_metadata(phospho.index),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(
            stage_order=(DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION,),
        ),
    )

    _, trace = PreprocessingPipeline(
        stage_contract_registry=(BATCH_CORRECTION_STAGE_CONTRACT,),
    ).run_with_trace(state)

    assert len(trace) == 1
    assert trace[0].stage == DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION
    assert trace[0].operation == "none"
    assert trace[0].diagnostics["status"] == "disabled"


def test_registry_instance_builder_contains_no_concrete_batch_special_case() -> None:
    source = inspect.getsource(
        stage_registry_module.build_registered_preprocessing_stage_instances
    )
    tree = ast.parse(source)

    assert "BatchCorrectionStage" not in source
    assert "BATCH_CORRECTION_STAGE_CONTRACT" not in source
    assert "batch_correction" not in source
    assert "factory is" not in source
    assert "stage_key ==" not in source
    for node in ast.walk(tree):
        if isinstance(node, ast.If):
            assert "stage_key" not in ast.unparse(node.test)


def test_pipeline_rejects_both_registry_aliases() -> None:
    with pytest.raises(
        DatasetBuildError,
        match="stage_contract_registry.*stage_metadata_registry.*aliases.*only one",
    ):
        PreprocessingPipeline(
            stage_contract_registry=(),
            stage_metadata_registry=(),
        )


def test_legacy_stage_metadata_registry_alias_warns_and_still_works() -> None:
    class FakeStage(_NoOpPreQuantitativeContractStage):
        stage_key = "legacy_alias_stage"

        def run(self, state: PreprocessingState) -> PreprocessingStageResult:
            return PreprocessingStageResult(state=state)

    contract = PreprocessingStageMetadata(
        stage_key="legacy_alias_stage",
        display_label="legacy_alias_stage",
        operation_name=lambda _plan: "legacy_alias_stage",
        serialize_parameters=lambda _plan: {},
        consumed_input_tables=("dataset.phospho",),
        produced_output_tables=("dataset.phospho",),
        quantitative_contract=_preserve_quantitative_contract(),
        diagnostics_metadata={"known_diagnostics_fields": ("policy",)},
    )
    state = _fake_stage_state("legacy_alias_stage")

    with pytest.warns(DeprecationWarning, match="stage_metadata_registry"):
        pipeline = PreprocessingPipeline(
            stage_registry=(FakeStage(),),
            stage_metadata_registry=(contract,),
        )

    _, trace = pipeline.run_with_trace(state)

    assert trace[0].stage == "legacy_alias_stage"


def test_custom_stage_registration_is_stage_owned() -> None:
    class FakeStage(_NoOpPreQuantitativeContractStage):
        stage_key = "fake_stage"

        def run(self, state: PreprocessingState) -> PreprocessingStageResult:
            return PreprocessingStageResult(
                state=state,
                diagnostics={
                    "notes": "stage executed",
                    "diagnostics": {"policy": "fake"},
                },
            )

    def _build_fake_stage(_context: PreprocessingStageFactoryContext) -> FakeStage:
        return FakeStage()

    fake_contract = PreprocessingStageMetadata(
        stage_key="fake_stage",
        display_label="fake_stage",
        provenance_stage="fake_stage",
        operation_name=lambda _plan: "fake",
        serialize_parameters=lambda _plan: {"mode": "test"},
        consumed_input_tables=("dataset.phospho",),
        produced_output_tables=("dataset.phospho",),
        stage_factory=_build_fake_stage,
        quantitative_contract=_preserve_quantitative_contract(),
        diagnostics_metadata={"known_diagnostics_fields": ("policy",)},
    )
    state = PreprocessingState(
        phospho=_phospho(),
        site_metadata=_site_metadata(_phospho().index),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(stage_order=("fake_stage",)),
    )

    _, trace = PreprocessingPipeline(
        stage_contract_registry=(fake_contract,),
    ).run_with_trace(state)

    assert len(trace) == 1
    assert trace[0].stage == "fake_stage"
    assert trace[0].operation == "fake"
    assert trace[0].parameters == {"mode": "test"}
    assert trace[0].diagnostics["policy"] == "fake"


def _fake_stage_state(stage_key: str) -> PreprocessingState:
    phospho = _phospho()
    return PreprocessingState(
        phospho=phospho,
        site_metadata=_site_metadata(phospho.index),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(stage_order=(stage_key,)),
    )


def _fake_stage_contract(
    *,
    stage_key: str,
    determinism_kind: DeterminismKind,
) -> PreprocessingStageMetadata:
    return PreprocessingStageMetadata(
        stage_key=stage_key,
        display_label=stage_key,
        provenance_stage=stage_key,
        operation_name=lambda _plan: "fake_operation",
        serialize_parameters=lambda _plan: {"mode": "test"},
        consumed_input_tables=("dataset.phospho",),
        produced_output_tables=("dataset.phospho",),
        determinism_kind=determinism_kind,
        quantitative_contract=_preserve_quantitative_contract(),
        diagnostics_metadata={"known_diagnostics_fields": ("random_seed",)},
    )


def test_deterministic_stage_contract_runs_without_seed() -> None:
    class FakeDeterministicStage(_NoOpPreQuantitativeContractStage):
        stage_key = "fake_deterministic"

        def run(self, state: PreprocessingState) -> PreprocessingStageResult:
            return PreprocessingStageResult(state=state)

    _, trace = PreprocessingPipeline(
        stage_registry=(FakeDeterministicStage(),),
        stage_contract_registry=(
            _fake_stage_contract(
                stage_key="fake_deterministic",
                determinism_kind=DeterminismKind.DETERMINISTIC,
            ),
        ),
    ).run_with_trace(_fake_stage_state("fake_deterministic"))

    assert trace[0].determinism is DeterminismKind.DETERMINISTIC
    assert trace[0].random_seed is None
    assert trace[0].reproducibility_caveats == ()


def test_seeded_stochastic_stage_contract_runs_with_explicit_seed() -> None:
    class FakeSeededStage(_NoOpPreQuantitativeContractStage):
        stage_key = "fake_seeded"

        def run(self, state: PreprocessingState) -> PreprocessingStageResult:
            return PreprocessingStageResult(
                state=state,
                diagnostics={"diagnostics": {"random_seed": 123}},
            )

    _, trace = PreprocessingPipeline(
        stage_registry=(FakeSeededStage(),),
        stage_contract_registry=(
            _fake_stage_contract(
                stage_key="fake_seeded",
                determinism_kind=DeterminismKind.SEEDED_STOCHASTIC,
            ),
        ),
    ).run_with_trace(_fake_stage_state("fake_seeded"))

    assert trace[0].determinism is DeterminismKind.SEEDED_STOCHASTIC
    assert trace[0].random_seed == 123
    assert trace[0].reproducibility_caveats == ()


def test_seeded_stochastic_stage_contract_without_seed_fails() -> None:
    class FakeSeededStage(_NoOpPreQuantitativeContractStage):
        stage_key = "fake_seeded_without_seed"

        def run(self, state: PreprocessingState) -> PreprocessingStageResult:
            return PreprocessingStageResult(state=state)

    with pytest.raises(
        DatasetBuildError, match="did not record an explicit random seed"
    ):
        PreprocessingPipeline(
            stage_registry=(FakeSeededStage(),),
            stage_contract_registry=(
                _fake_stage_contract(
                    stage_key="fake_seeded_without_seed",
                    determinism_kind=DeterminismKind.SEEDED_STOCHASTIC,
                ),
            ),
        ).run_with_trace(_fake_stage_state("fake_seeded_without_seed"))


def test_externally_nondeterministic_stage_records_reproducibility_caveat() -> None:
    class FakeExternalStage(_NoOpPreQuantitativeContractStage):
        stage_key = "fake_external"

        def run(self, state: PreprocessingState) -> PreprocessingStageResult:
            return PreprocessingStageResult(state=state)

    _, trace = PreprocessingPipeline(
        stage_registry=(FakeExternalStage(),),
        stage_contract_registry=(
            _fake_stage_contract(
                stage_key="fake_external",
                determinism_kind=DeterminismKind.EXTERNALLY_NONDETERMINISTIC,
            ),
        ),
    ).run_with_trace(_fake_stage_state("fake_external"))

    record = trace[0]
    assert record.determinism is DeterminismKind.EXTERNALLY_NONDETERMINISTIC
    assert record.is_deterministic is False
    assert len(record.reproducibility_caveats) == 1
    caveat = record.reproducibility_caveats[0]
    assert caveat.code == PREPROCESSING_EXTERNAL_NONDETERMINISM_CAVEAT_CODE
    assert caveat.severity == "warning"
    assert caveat.details["determinism_kind"] == "externally_nondeterministic"
