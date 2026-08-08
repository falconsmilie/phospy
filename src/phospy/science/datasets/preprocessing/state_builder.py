"""Processing-state construction collaborator for dataset preprocessing builder."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from phospy.errors.build import DatasetBuildError
from phospy.errors.input import PhosPyInputError
from phospy.errors.transformations import InvalidTransformationStateError
from phospy.provenance.models import TableFingerprint
from phospy.science.datasets.preprocessing.diagnostics import ProcessingTraceDiagnostics
from phospy.science.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_MISSING_DATA,
    DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
    DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
    PreprocessingPlan,
    PreprocessingStageExecution,
)
from phospy.science.datasets.preprocessing.policy_models import (
    ComparisonBuildingPolicy,
    SiteMatrixPolicy,
    TotalProteinCorrectionPolicy,
)
from phospy.science.datasets.preprocessing.quantitative_evidence import (
    QuantitativeOperationEvidenceContext,
    QuantitativeOperationEvidenceValidator,
)
from phospy.science.datasets.preprocessing.stage_registry import (
    get_preprocessing_stage_metadata,
)
from phospy.science.datasets.processing_state import (
    TOTAL_PROTEIN_CORRECTION_DIAGNOSTICS_SCHEMA_VERSION_V1,
    ComparisonState,
    DatasetProcessingState,
    MissingDataDiagnostics,
    MissingDataState,
    NormalisationState,
    RuvReadinessState,
    SiteMatrixState,
    SiteSequenceResolutionState,
    TotalProteinCorrectionDiagnostics,
    TotalProteinCorrectionState,
)
from phospy.science.transformations._authority import (
    dataset_quantitative_meaning_transition_authority,
)
from phospy.science.transformations.models import (
    QUANTITATIVE_MEANING_OPERATION_CALLER_DECLARATION,
    QUANTITATIVE_MEANING_OPERATION_SCALE_CONTRACT_INFERENCE,
    QUANTITATIVE_MEANING_USER_DECLARED_CAVEAT_CODE,
    IntensityScaleState,
    QuantitativeMeaning,
    QuantitativeMeaningEvidenceMode,
    QuantitativeMeaningTransitionProvenance,
    caller_declarable_quantitative_meaning_values,
    default_quantitative_meaning_for_scale_kind,
    is_caller_declarable_quantitative_meaning,
)
from phospy.science.transformations.quantitative_contracts import (
    QuantitativeContractState,
    QuantitativeOperationContract,
)

_DATASET_BUILDER_QUANTITATIVE_MEANING_PRODUCER = (
    "phospy.science.datasets.preprocessing.state_builder"
)


@dataclass(frozen=True, slots=True)
class _PreprocessingTraceFacts:
    comparison_pairs: tuple[tuple[str, str], ...] | None
    resolved_total_policy: TotalProteinCorrectionPolicy
    total_correction_applied: bool
    parsed: ProcessingTraceDiagnostics
    correction_diagnostics: dict[str, object] | None
    missing_data_diagnostics: dict[str, object] | None
    site_sequence_resolution_diagnostics: dict[str, object] | None


@dataclass(frozen=True, slots=True)
class _SiteSequenceResolutionDefaults:
    configured: bool
    mode: str | None
    flank_size: int | None
    conflict_policy: str | None


@dataclass(frozen=True, slots=True)
class _QuantitativeMeaningResolution:
    intensity_scale_state: IntensityScaleState
    default_quantitative_meaning: str


@dataclass(frozen=True, slots=True)
class _MatrixCompletenessResolution:
    correction_diagnostics: dict[str, object]
    typed_correction_diagnostics: TotalProteinCorrectionDiagnostics
    typed_missing_data_diagnostics: MissingDataDiagnostics | None
    output_missing_cell_count: int
    imputed_cell_count: int
    default_imputation_input_scale: str | None
    default_imputation_operation_order: str | None

    @property
    def complete_matrix(self) -> bool:
        return self.output_missing_cell_count == 0

    @property
    def imputed(self) -> bool:
        return self.imputed_cell_count > 0


@dataclass(frozen=True, slots=True)
class _BatchCorrectionAndRuvReadinessResolution:
    ruv_readiness: RuvReadinessState


@dataclass(frozen=True, slots=True)
class _ProcessingStateComponents:
    intensity_scale: IntensityScaleState
    site_sequence_resolution: SiteSequenceResolutionState
    missing_data: MissingDataState
    normalisation: NormalisationState
    total_protein_correction: TotalProteinCorrectionState
    site_matrix: SiteMatrixState
    comparisons: ComparisonState
    ruv_readiness: RuvReadinessState


class DatasetProcessingStateBuilder:
    """Build `DatasetProcessingState` from a preprocessing plan and execution trace."""

    def build(
        self,
        *,
        plan: PreprocessingPlan,
        intensity_scale_state: IntensityScaleState,
        explicit_quantitative_meaning: QuantitativeMeaning | None = None,
        preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None = None,
        final_phospho: pd.DataFrame | None = None,
        final_site_metadata: pd.DataFrame | None = None,
        final_sample_metadata: pd.DataFrame | None = None,
    ) -> DatasetProcessingState:
        trace_facts = _collect_preprocessing_trace_facts(
            plan=plan,
            preprocessing_trace=preprocessing_trace,
        )
        site_sequence_defaults = _resolve_site_sequence_resolution_defaults(
            plan=plan,
            trace_facts=trace_facts,
        )
        quantitative_state = _resolve_quantitative_meaning(
            plan=plan,
            intensity_scale_state=intensity_scale_state,
            trace_facts=trace_facts,
            explicit_quantitative_meaning=explicit_quantitative_meaning,
            preprocessing_trace=preprocessing_trace,
        )
        matrix_state = _resolve_intensity_scale_and_matrix_completeness_state(
            plan=plan,
            trace_facts=trace_facts,
            quantitative_state=quantitative_state,
        )
        batch_and_ruv_state = _resolve_batch_correction_and_ruv_readiness(
            plan=plan,
            trace_facts=trace_facts,
            matrix_state=matrix_state,
            final_phospho=final_phospho,
            final_site_metadata=final_site_metadata,
            final_sample_metadata=final_sample_metadata,
        )
        _validate_processing_state_cross_invariants(
            quantitative_state=quantitative_state,
            matrix_state=matrix_state,
            batch_and_ruv_state=batch_and_ruv_state,
        )
        components = _assemble_processing_state_components(
            plan=plan,
            trace_facts=trace_facts,
            site_sequence_defaults=site_sequence_defaults,
            quantitative_state=quantitative_state,
            matrix_state=matrix_state,
            batch_and_ruv_state=batch_and_ruv_state,
        )
        return _assemble_final_processing_state(components)


def _collect_preprocessing_trace_facts(
    *,
    plan: PreprocessingPlan,
    preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None,
) -> _PreprocessingTraceFacts:
    comparison_pairs = (
        None
        if plan.comparison_pairs is None
        else tuple((str(left), str(right)) for left, right in plan.comparison_pairs)
    )
    resolved_total_policy = plan.total_protein_correction_policy
    parsed = ProcessingTraceDiagnostics.from_trace(preprocessing_trace)
    return _PreprocessingTraceFacts(
        comparison_pairs=comparison_pairs,
        resolved_total_policy=resolved_total_policy,
        total_correction_applied=(
            resolved_total_policy is not TotalProteinCorrectionPolicy.NONE
        ),
        parsed=parsed,
        correction_diagnostics=parsed.total_protein_correction,
        missing_data_diagnostics=parsed.missing_data,
        site_sequence_resolution_diagnostics=parsed.site_sequence_resolution,
    )


def _resolve_site_sequence_resolution_defaults(
    *,
    plan: PreprocessingPlan,
    trace_facts: _PreprocessingTraceFacts,
) -> _SiteSequenceResolutionDefaults:
    configured = bool(
        trace_facts.parsed.resolve_optional_bool(
            trace_facts.site_sequence_resolution_diagnostics,
            stage=DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
            key="configured",
            default=bool(plan.site_sequence_resolution_enabled),
        )
    )
    return _SiteSequenceResolutionDefaults(
        configured=configured,
        mode=plan.site_sequence_resolution_mode.value if configured else None,
        flank_size=(
            int(plan.site_sequence_resolution_flank_size) if configured else None
        ),
        conflict_policy=(
            plan.site_sequence_resolution_conflict_policy.value if configured else None
        ),
    )


def _resolve_quantitative_meaning(
    *,
    plan: PreprocessingPlan,
    intensity_scale_state: IntensityScaleState,
    trace_facts: _PreprocessingTraceFacts,
    explicit_quantitative_meaning: QuantitativeMeaning | None,
    preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None,
) -> _QuantitativeMeaningResolution:
    resolved_state = _resolve_quantitative_meaning_state(
        plan=plan,
        intensity_scale_state=intensity_scale_state,
        total_correction_policy=trace_facts.resolved_total_policy,
        correction_diagnostics=trace_facts.correction_diagnostics,
        explicit_quantitative_meaning=explicit_quantitative_meaning,
        preprocessing_trace=preprocessing_trace,
    )
    quantitative_meaning = resolved_state.quantity
    if quantitative_meaning is None:
        raise DatasetBuildError("intensity-scale state is missing quantitative meaning")
    return _QuantitativeMeaningResolution(
        intensity_scale_state=resolved_state,
        default_quantitative_meaning=quantitative_meaning.value,
    )


def _resolve_intensity_scale_and_matrix_completeness_state(
    *,
    plan: PreprocessingPlan,
    trace_facts: _PreprocessingTraceFacts,
    quantitative_state: _QuantitativeMeaningResolution,
) -> _MatrixCompletenessResolution:
    correction_diagnostics = _correction_diagnostics_with_defaults(
        trace_facts=trace_facts,
        quantitative_state=quantitative_state,
    )
    trace_facts.parsed.validate_site_sequence_resolution_payload(
        trace_facts.site_sequence_resolution_diagnostics
    )
    typed_correction_diagnostics = TotalProteinCorrectionDiagnostics.from_payload(
        correction_diagnostics,
        field_name="dataset processing state total_protein_correction.diagnostics",
    )
    typed_missing_data_diagnostics = (
        None
        if trace_facts.missing_data_diagnostics is None
        else MissingDataDiagnostics.from_payload(
            trace_facts.missing_data_diagnostics,
            field_name="dataset processing state missing_data.diagnostics",
        )
    )
    output_missing_cell_count = trace_facts.parsed.resolve_optional_int(
        trace_facts.missing_data_diagnostics,
        stage=DATASET_PREPROCESSING_STAGE_MISSING_DATA,
        key="output_missing_cell_count",
        default=0,
    )
    imputed_cell_count = trace_facts.parsed.resolve_optional_int(
        trace_facts.missing_data_diagnostics,
        stage=DATASET_PREPROCESSING_STAGE_MISSING_DATA,
        key="imputed_cell_count",
        default=0,
    )
    return _MatrixCompletenessResolution(
        correction_diagnostics=correction_diagnostics,
        typed_correction_diagnostics=typed_correction_diagnostics,
        typed_missing_data_diagnostics=typed_missing_data_diagnostics,
        output_missing_cell_count=output_missing_cell_count,
        imputed_cell_count=imputed_cell_count,
        default_imputation_input_scale=(
            None
            if plan.missing_data_input_scale is None
            else plan.missing_data_input_scale.value
        ),
        default_imputation_operation_order=(
            plan.missing_data_imputation_operation_order
        ),
    )


def _correction_diagnostics_with_defaults(
    *,
    trace_facts: _PreprocessingTraceFacts,
    quantitative_state: _QuantitativeMeaningResolution,
) -> dict[str, object]:
    if trace_facts.correction_diagnostics is not None:
        return trace_facts.correction_diagnostics
    return {
        "diagnostics_schema_version": (
            TOTAL_PROTEIN_CORRECTION_DIAGNOSTICS_SCHEMA_VERSION_V1
        ),
        "policy": trace_facts.resolved_total_policy.value,
        "requested_policy": trace_facts.resolved_total_policy.value,
        "resolved_policy": trace_facts.resolved_total_policy.value,
        "quantitative_meaning": quantitative_state.default_quantitative_meaning,
    }


def _resolve_batch_correction_and_ruv_readiness(
    *,
    plan: PreprocessingPlan,
    trace_facts: _PreprocessingTraceFacts,
    matrix_state: _MatrixCompletenessResolution,
    final_phospho: pd.DataFrame | None,
    final_site_metadata: pd.DataFrame | None,
    final_sample_metadata: pd.DataFrame | None,
) -> _BatchCorrectionAndRuvReadinessResolution:
    return _BatchCorrectionAndRuvReadinessResolution(
        ruv_readiness=_resolve_ruv_readiness_state(
            plan=plan,
            final_phospho=final_phospho,
            final_site_metadata=final_site_metadata,
            final_sample_metadata=final_sample_metadata,
            missing_data_diagnostics=trace_facts.missing_data_diagnostics,
            default_matrix_complete=matrix_state.complete_matrix,
        )
    )


def _validate_processing_state_cross_invariants(
    *,
    quantitative_state: _QuantitativeMeaningResolution,
    matrix_state: _MatrixCompletenessResolution,
    batch_and_ruv_state: _BatchCorrectionAndRuvReadinessResolution,
) -> None:
    if quantitative_state.intensity_scale_state.quantity is None:
        raise DatasetBuildError("intensity-scale state is missing quantitative meaning")
    if matrix_state.complete_matrix != (matrix_state.output_missing_cell_count == 0):
        raise DatasetBuildError(
            "dataset processing state matrix completeness disagrees with missing-data diagnostics"
        )
    ruv_readiness = batch_and_ruv_state.ruv_readiness
    if ruv_readiness.ready and ruv_readiness.requires_complete_matrix:
        if not ruv_readiness.matrix_complete:
            raise DatasetBuildError(
                "dataset processing state RUV readiness cannot be ready when the "
                "matrix is incomplete"
            )


def _assemble_processing_state_components(
    *,
    plan: PreprocessingPlan,
    trace_facts: _PreprocessingTraceFacts,
    site_sequence_defaults: _SiteSequenceResolutionDefaults,
    quantitative_state: _QuantitativeMeaningResolution,
    matrix_state: _MatrixCompletenessResolution,
    batch_and_ruv_state: _BatchCorrectionAndRuvReadinessResolution,
) -> _ProcessingStateComponents:
    return _ProcessingStateComponents(
        intensity_scale=quantitative_state.intensity_scale_state,
        site_sequence_resolution=_resolve_site_sequence_resolution_state(
            trace_facts=trace_facts,
            defaults=site_sequence_defaults,
        ),
        missing_data=_resolve_missing_data_state(
            plan=plan,
            trace_facts=trace_facts,
            matrix_state=matrix_state,
        ),
        normalisation=NormalisationState(policy=plan.normalisation_policy.value),
        total_protein_correction=_resolve_total_protein_correction_state(
            trace_facts=trace_facts,
            quantitative_state=quantitative_state,
            matrix_state=matrix_state,
        ),
        site_matrix=SiteMatrixState(
            policy=plan.site_matrix_policy.value,
            constructed=plan.site_matrix_policy is SiteMatrixPolicy.BUILD_FROM_METADATA,
            missing_data_policy=plan.site_matrix_missing_data_policy.value,
            minimum_observed_values=plan.site_matrix_minimum_observed_values,
            duplicate_site_policy=plan.site_matrix_duplicate_site_policy.value,
        ),
        comparisons=ComparisonState(
            policy=plan.comparison_building_policy.value,
            sample_group_column=plan.comparison_sample_group_column,
            pairs=(
                None
                if plan.comparison_building_policy is ComparisonBuildingPolicy.NONE
                else trace_facts.comparison_pairs
            ),
        ),
        ruv_readiness=batch_and_ruv_state.ruv_readiness,
    )


def _resolve_site_sequence_resolution_state(
    *,
    trace_facts: _PreprocessingTraceFacts,
    defaults: _SiteSequenceResolutionDefaults,
) -> SiteSequenceResolutionState:
    diagnostics = trace_facts.site_sequence_resolution_diagnostics
    parsed = trace_facts.parsed
    return SiteSequenceResolutionState(
        configured=defaults.configured,
        mode=parsed.resolve_optional_string(
            diagnostics,
            stage=DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
            key="mode",
            default=defaults.mode,
        ),
        flank_size=parsed.resolve_optional_nullable_int(
            diagnostics,
            stage=DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
            key="flank_size",
            default=defaults.flank_size,
        ),
        fasta_source_path=parsed.resolve_optional_string(
            diagnostics,
            stage=DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
            key="fasta_source_path",
            default=None,
        ),
        fasta_source_label=parsed.resolve_optional_string(
            diagnostics,
            stage=DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
            key="fasta_source_label",
            default=None,
        ),
        fasta_sha256=parsed.resolve_optional_string(
            diagnostics,
            stage=DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
            key="fasta_sha256",
            default=None,
        ),
        resolver_version=parsed.resolve_optional_string(
            diagnostics,
            stage=DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
            key="resolver_version",
            default=None,
        ),
        resolved_site_count=parsed.resolve_optional_int(
            diagnostics,
            stage=DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
            key="resolved_site_count",
            default=0,
        ),
        unresolved_site_count=parsed.resolve_optional_int(
            diagnostics,
            stage=DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
            key="unresolved_site_count",
            default=0,
        ),
        unresolved_counts_by_reason=parsed.resolve_optional_mapping_int(
            diagnostics,
            stage=DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
            key="unresolved_counts_by_reason",
        ),
        filled_missing_count=parsed.resolve_optional_int(
            diagnostics,
            stage=DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
            key="filled_missing_count",
            default=0,
        ),
        replaced_existing_count=parsed.resolve_optional_int(
            diagnostics,
            stage=DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
            key="replaced_existing_count",
            default=0,
        ),
        preserved_existing_count=parsed.resolve_optional_int(
            diagnostics,
            stage=DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
            key="preserved_existing_count",
            default=0,
        ),
        existing_sequence_conflict_count=parsed.resolve_optional_int(
            diagnostics,
            stage=DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
            key="existing_sequence_conflict_count",
            default=0,
        ),
        conflict_policy=parsed.resolve_optional_string(
            diagnostics,
            stage=DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
            key="conflict_policy",
            default=defaults.conflict_policy,
        ),
        row_diagnostics=parsed.resolve_site_sequence_row_diagnostics(diagnostics),
    )


def _resolve_missing_data_state(
    *,
    plan: PreprocessingPlan,
    trace_facts: _PreprocessingTraceFacts,
    matrix_state: _MatrixCompletenessResolution,
) -> MissingDataState:
    return MissingDataState(
        policy=plan.missing_data_policy,
        min_observed_values=plan.missing_data_min_observed_values,
        complete_matrix=matrix_state.complete_matrix,
        imputed=matrix_state.imputed,
        diagnostics=matrix_state.typed_missing_data_diagnostics,
        has_missing_values=(matrix_state.output_missing_cell_count > 0),
        missing_value_count=matrix_state.output_missing_cell_count,
        imputation_input_scale=trace_facts.parsed.resolve_optional_string(
            trace_facts.missing_data_diagnostics,
            stage=DATASET_PREPROCESSING_STAGE_MISSING_DATA,
            key="imputation_input_scale",
            default=matrix_state.default_imputation_input_scale,
        ),
        imputation_operation_order=trace_facts.parsed.resolve_optional_string(
            trace_facts.missing_data_diagnostics,
            stage=DATASET_PREPROCESSING_STAGE_MISSING_DATA,
            key="imputation_operation_order",
            default=matrix_state.default_imputation_operation_order,
        ),
    )


def _resolve_total_protein_correction_state(
    *,
    trace_facts: _PreprocessingTraceFacts,
    quantitative_state: _QuantitativeMeaningResolution,
    matrix_state: _MatrixCompletenessResolution,
) -> TotalProteinCorrectionState:
    default_formula = (
        "log2_phospho - log2_total"
        if trace_facts.resolved_total_policy
        is TotalProteinCorrectionPolicy.SUBTRACT_LOG_TOTAL
        else None
    )
    default_input_scale = (
        "log2"
        if trace_facts.resolved_total_policy
        is TotalProteinCorrectionPolicy.SUBTRACT_LOG_TOTAL
        else None
    )
    default_output_scale = (
        "log2_ratio"
        if trace_facts.resolved_total_policy
        is TotalProteinCorrectionPolicy.SUBTRACT_LOG_TOTAL
        else None
    )
    return TotalProteinCorrectionState(
        policy=trace_facts.resolved_total_policy,
        applied=trace_facts.total_correction_applied,
        formula=trace_facts.parsed.resolve_optional_string(
            matrix_state.correction_diagnostics,
            stage=DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
            key="formula",
            default=default_formula,
        ),
        requires_log_scale=trace_facts.parsed.resolve_optional_bool(
            matrix_state.correction_diagnostics,
            stage=DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
            key="requires_log_scale",
            default=bool(trace_facts.total_correction_applied),
        ),
        input_scale=trace_facts.parsed.resolve_optional_string(
            matrix_state.correction_diagnostics,
            stage=DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
            key="input_scale",
            default=default_input_scale,
        ),
        output_scale=trace_facts.parsed.resolve_optional_string(
            matrix_state.correction_diagnostics,
            stage=DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
            key="output_scale",
            default=default_output_scale,
        ),
        quantitative_meaning=_resolve_total_correction_quantitative_meaning(
            trace_facts.parsed.resolve_optional_string(
                matrix_state.correction_diagnostics,
                stage=DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
                key="quantitative_meaning",
                default=quantitative_state.default_quantitative_meaning,
            )
        ),
        diagnostics=matrix_state.typed_correction_diagnostics,
    )


def _assemble_final_processing_state(
    components: _ProcessingStateComponents,
) -> DatasetProcessingState:
    return DatasetProcessingState(
        intensity_scale=components.intensity_scale,
        site_sequence_resolution=components.site_sequence_resolution,
        missing_data=components.missing_data,
        normalisation=components.normalisation,
        total_protein_correction=components.total_protein_correction,
        site_matrix=components.site_matrix,
        comparisons=components.comparisons,
        ruv_readiness=components.ruv_readiness,
    )


def _resolve_quantitative_meaning_state(
    *,
    plan: PreprocessingPlan,
    intensity_scale_state: IntensityScaleState,
    total_correction_policy: TotalProteinCorrectionPolicy,
    correction_diagnostics: dict[str, object] | None,
    explicit_quantitative_meaning: QuantitativeMeaning | None = None,
    preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None = None,
) -> IntensityScaleState:
    resolved_state = intensity_scale_state
    if explicit_quantitative_meaning is not None:
        resolved_state = _apply_caller_declared_quantitative_meaning(
            intensity_scale_state=resolved_state,
            target=explicit_quantitative_meaning,
        )
    elif resolved_state.quantity is None:
        resolved_state = _infer_quantitative_meaning_from_scale_contract(
            intensity_scale_state=resolved_state
        )
    resolved_state = _apply_quantitative_operation_contract_transitions(
        plan=plan,
        plan_total_correction_policy=total_correction_policy,
        intensity_scale_state=resolved_state,
        preprocessing_trace=preprocessing_trace,
    )
    _reject_diagnostic_quantitative_meaning_mismatch(
        correction_diagnostics=correction_diagnostics,
        resolved_state=resolved_state,
    )
    return resolved_state


def _apply_caller_declared_quantitative_meaning(
    *,
    intensity_scale_state: IntensityScaleState,
    target: QuantitativeMeaning,
) -> IntensityScaleState:
    if not is_caller_declarable_quantitative_meaning(target):
        allowed = caller_declarable_quantitative_meaning_values()
        raise DatasetBuildError(
            "dataset build request quantitative_meaning may only declare direct "
            "input meanings: " + ", ".join(allowed) + f"; got {target.value!r}"
        )
    provenance = QuantitativeMeaningTransitionProvenance(
        source_quantity=intensity_scale_state.quantity,
        target_quantity=target,
        operation_id=QUANTITATIVE_MEANING_OPERATION_CALLER_DECLARATION,
        producer_id=_DATASET_BUILDER_QUANTITATIVE_MEANING_PRODUCER,
        evidence_mode=QuantitativeMeaningEvidenceMode.DECLARED_BY_CALLER,
        parameters={
            "request_field": "DatasetBuildRequest.quantitative_meaning",
            "declared_quantitative_meaning": target.value,
        },
        diagnostic_caveat_codes=(QUANTITATIVE_MEANING_USER_DECLARED_CAVEAT_CODE,),
    )
    try:
        return intensity_scale_state.transition_quantitative_meaning(
            target_quantity=target,
            provenance=provenance,
            authority=dataset_quantitative_meaning_transition_authority(),
        )
    except InvalidTransformationStateError as exc:
        raise PhosPyInputError(
            "dataset build request quantitative_meaning is incompatible with the "
            f"established intensity scale: {exc}"
        ) from exc


def _infer_quantitative_meaning_from_scale_contract(
    *,
    intensity_scale_state: IntensityScaleState,
) -> IntensityScaleState:
    target = default_quantitative_meaning_for_scale_kind(intensity_scale_state.kind)
    provenance = QuantitativeMeaningTransitionProvenance(
        source_quantity=intensity_scale_state.quantity,
        target_quantity=target,
        operation_id=QUANTITATIVE_MEANING_OPERATION_SCALE_CONTRACT_INFERENCE,
        producer_id=_DATASET_BUILDER_QUANTITATIVE_MEANING_PRODUCER,
        evidence_mode=QuantitativeMeaningEvidenceMode.INFERRED_FROM_SCALE_CONTRACT,
        parameters={
            "scale_kind": intensity_scale_state.kind.value,
            "scale_label": intensity_scale_state.label,
        },
    )
    return intensity_scale_state.transition_quantitative_meaning(
        target_quantity=target,
        provenance=provenance,
        authority=dataset_quantitative_meaning_transition_authority(),
    )


def _apply_quantitative_operation_contract_transitions(
    *,
    plan: PreprocessingPlan,
    intensity_scale_state: IntensityScaleState,
    plan_total_correction_policy: TotalProteinCorrectionPolicy,
    preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None,
) -> IntensityScaleState:
    _require_trace_for_contract_emitting_plan_stages(
        plan=plan,
        preprocessing_trace=preprocessing_trace,
    )
    resolved_state = intensity_scale_state
    evidence_validator = QuantitativeOperationEvidenceValidator()
    for stage in preprocessing_trace or ():
        contract = _resolve_trace_quantitative_contract(stage=stage, plan=plan)
        current_quantity = resolved_state.quantity
        if current_quantity is None:
            raise DatasetBuildError(
                "quantitative operation contract transition requires an established "
                f"input quantitative meaning: stage={stage.stage!r}"
            )
        contract_state = QuantitativeContractState(
            scale_kind=resolved_state.kind,
            meaning=current_quantity,
        )
        evidence_validator.validate(
            QuantitativeOperationEvidenceContext(
                stage=stage.stage,
                operation=stage.operation,
                quantitative_contract=contract,
                trace_record=stage,
                input_quantitative_state=contract_state,
                input_intensity_scale_state=resolved_state,
                input_quantitative_meaning_evidence_mode=(
                    _resolve_quantitative_meaning_evidence_mode(resolved_state)
                ),
            )
        )
        if not contract.emits_quantitative_meaning_state_event:
            continue
        target_state = contract.validate_and_transition(
            contract_state,
            stage=stage.stage,
            operation=stage.operation,
            evidence=stage.quantitative_transition_evidence,
        )
        if target_state.meaning is current_quantity:
            continue
        provenance = _build_contract_quantitative_meaning_provenance(
            stage=stage,
            contract=contract,
            source=current_quantity,
            target=target_state.meaning,
        )
        resolved_state = resolved_state.transition_quantitative_meaning(
            target_quantity=target_state.meaning,
            provenance=provenance,
            authority=dataset_quantitative_meaning_transition_authority(),
        )
    if (
        plan_total_correction_policy is TotalProteinCorrectionPolicy.SUBTRACT_LOG_TOTAL
        and resolved_state.quantity
        not in {
            QuantitativeMeaning.PHOSPHO_TOTAL_LOG_RATIO,
            QuantitativeMeaning.MIXED_PHOSPHO_TOTAL_LOG_RATIO_AND_PHOSPHOSITE_LOG_ABUNDANCE,
        }
    ):
        raise DatasetBuildError(
            "dataset preprocessing plan requested total-protein correction but no "
            "quantitative operation contract transition changed the output meaning"
        )
    return resolved_state


def _resolve_quantitative_meaning_evidence_mode(
    state: IntensityScaleState,
) -> QuantitativeMeaningEvidenceMode | None:
    provenance = state.quantitative_meaning_provenance
    if provenance is None:
        return None
    evidence_mode = provenance.evidence_mode
    if isinstance(evidence_mode, QuantitativeMeaningEvidenceMode):
        return evidence_mode
    return QuantitativeMeaningEvidenceMode(str(evidence_mode))


def _require_trace_for_contract_emitting_plan_stages(
    *,
    plan: PreprocessingPlan,
    preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None,
) -> None:
    observed_stages = {stage.stage for stage in preprocessing_trace or ()}
    stage_keys = list(plan.stage_order)
    if (
        plan.total_protein_correction_policy
        is TotalProteinCorrectionPolicy.SUBTRACT_LOG_TOTAL
        and DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION not in stage_keys
    ):
        stage_keys.append(DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION)
    for stage_key in stage_keys:
        if stage_key in observed_stages:
            continue
        metadata = get_preprocessing_stage_metadata(stage_key)
        interpreted = metadata.interpret(plan)
        if (
            interpreted.quantitative_contract.emits_quantitative_meaning_state_event
            and interpreted.stage not in observed_stages
        ):
            raise DatasetBuildError(
                "quantitative operation contract transition requires a "
                f"{interpreted.stage} preprocessing trace stage with fingerprints"
            )


def _resolve_trace_quantitative_contract(
    *,
    stage: PreprocessingStageExecution,
    plan: PreprocessingPlan,
) -> QuantitativeOperationContract:
    if stage.quantitative_contract is not None:
        return stage.quantitative_contract
    metadata = get_preprocessing_stage_metadata(stage.stage)
    interpreted = metadata.interpret(plan)
    return interpreted.quantitative_contract


def _build_contract_quantitative_meaning_provenance(
    *,
    stage: PreprocessingStageExecution,
    contract: QuantitativeOperationContract,
    source: QuantitativeMeaning,
    target: QuantitativeMeaning,
) -> QuantitativeMeaningTransitionProvenance:
    operation_id = contract.operation_id
    producer_id = contract.producer_id
    output_table = contract.provenance_output_table
    if operation_id is None or producer_id is None or output_table is None:
        raise DatasetBuildError(
            "quantitative operation contract is missing semantic provenance "
            f"metadata: stage={stage.stage!r}, operation={stage.operation!r}"
        )
    _require_stage_input_fingerprint_names(
        stage,
        required_names=contract.provenance_input_tables,
    )
    output_fingerprint = _require_stage_output_fingerprint(
        stage,
        table_name=output_table,
    )
    parameters = _contract_transition_parameters(
        stage=stage,
        contract=contract,
    )
    return QuantitativeMeaningTransitionProvenance(
        source_quantity=source,
        target_quantity=target,
        operation_id=operation_id,
        producer_id=producer_id,
        evidence_mode=contract.evidence_mode,
        parameters=parameters,
        input_table_fingerprints=stage.consumed_input_tables,
        output_table_fingerprint=output_fingerprint,
        trace_id=f"{stage.stage}:{stage.operation}:{stage.output_hash}",
        diagnostic_caveat_codes=contract.caveat_codes(
            target=target,
            evidence=stage.quantitative_transition_evidence,
        ),
    )


def _contract_transition_parameters(
    *,
    stage: PreprocessingStageExecution,
    contract: QuantitativeOperationContract,
) -> dict[str, object]:
    parameters = dict(stage.parameters)
    diagnostics = stage.diagnostics
    for key in contract.provenance_diagnostic_fields:
        if key in diagnostics:
            parameters[f"diagnostics.{key}"] = diagnostics[key]
    evidence = stage.quantitative_transition_evidence
    if evidence is not None:
        for key, value in evidence.to_payload().items():
            parameters[f"evidence.{key}"] = value
    parameters["semantic_contract"] = {
        "preserves_abundance": bool(contract.preserves_abundance),
        "negative_domain_policy": contract.negative_domain_policy.value,
        "reversibility": contract.reversibility.value,
        "information_loss": contract.information_loss.value,
        "emits_state_transition_event": bool(contract.emits_state_transition_event),
        "required_evidence": sorted(item.value for item in contract.required_evidence),
    }
    return parameters


def _require_stage_output_fingerprint(
    stage: PreprocessingStageExecution,
    *,
    table_name: str,
) -> TableFingerprint:
    for fingerprint in stage.produced_output_tables:
        if fingerprint.name == table_name:
            return fingerprint
    raise DatasetBuildError(
        "quantitative operation contract transition requires output "
        f"fingerprint for {table_name!r}"
    )


def _require_stage_input_fingerprint_names(
    stage: PreprocessingStageExecution,
    *,
    required_names: tuple[str, ...],
) -> None:
    observed = {fingerprint.name for fingerprint in stage.consumed_input_tables}
    missing = tuple(name for name in required_names if name not in observed)
    if missing:
        raise DatasetBuildError(
            "quantitative operation contract transition requires input "
            "fingerprints for " + ", ".join(repr(name) for name in missing)
        )


def _reject_diagnostic_quantitative_meaning_mismatch(
    *,
    correction_diagnostics: dict[str, object] | None,
    resolved_state: IntensityScaleState,
) -> None:
    diagnostic_meaning = ProcessingTraceDiagnostics.resolve_optional_string(
        correction_diagnostics,
        stage=DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
        key="quantitative_meaning",
        default=None,
    )
    if diagnostic_meaning is None:
        return
    try:
        diagnostic_quantity = QuantitativeMeaning(diagnostic_meaning)
    except ValueError as exc:
        supported = ", ".join(member.value for member in QuantitativeMeaning)
        raise DatasetBuildError(
            "dataset preprocessing total_protein_correction diagnostics "
            "quantitative_meaning must be one of: "
            f"{supported}; got {diagnostic_meaning!r}"
        ) from exc
    resolved_quantity = resolved_state.quantity
    if resolved_quantity is None:
        raise DatasetBuildError("intensity-scale state is missing quantitative meaning")
    if diagnostic_quantity is resolved_quantity:
        return
    raise DatasetBuildError(
        "dataset preprocessing total_protein_correction diagnostics "
        "quantitative_meaning does not match the typed quantitative contract "
        "transition output: "
        f"diagnostics={diagnostic_meaning!r}, contract={resolved_quantity.value!r}"
    )


def _resolve_ruv_readiness_state(
    *,
    plan: PreprocessingPlan,
    final_phospho: pd.DataFrame | None,
    final_site_metadata: pd.DataFrame | None,
    final_sample_metadata: pd.DataFrame | None,
    missing_data_diagnostics: dict[str, object] | None,
    default_matrix_complete: bool,
) -> RuvReadinessState:
    enabled = bool(plan.ruv_readiness_enabled)
    matrix_complete = _resolve_matrix_completeness(
        final_phospho=final_phospho,
        default_matrix_complete=default_matrix_complete,
    )
    control_feature_column = str(plan.ruv_readiness_control_feature_column).strip()
    replicate_group_column = str(plan.ruv_readiness_replicate_group_column).strip()
    batch_column = plan.ruv_readiness_batch_column

    control_feature_count = _count_control_features(
        site_metadata=final_site_metadata,
        control_feature_column=control_feature_column,
    )
    replicate_group_count = _count_distinct_non_missing(
        sample_metadata=final_sample_metadata,
        column=replicate_group_column,
    )
    batch_count = (
        None
        if batch_column is None
        else _count_distinct_non_missing(
            sample_metadata=final_sample_metadata,
            column=batch_column,
        )
    )

    imputation_method_id = ProcessingTraceDiagnostics.resolve_optional_string(
        missing_data_diagnostics,
        stage=DATASET_PREPROCESSING_STAGE_MISSING_DATA,
        key="imputation_method_id",
        default=None,
    )
    missingness_mask_hash = ProcessingTraceDiagnostics.resolve_optional_string(
        missing_data_diagnostics,
        stage=DATASET_PREPROCESSING_STAGE_MISSING_DATA,
        key="missingness_mask_hash",
        default=None,
    )
    missingness_mask_preserved = missingness_mask_hash is not None
    reasons: list[str] = []
    if not enabled:
        reasons.append("not configured")
    else:
        if not matrix_complete:
            reasons.append("matrix contains missing values")
        if missing_data_diagnostics is None:
            reasons.append("missing-data diagnostics unavailable")
        if missingness_mask_hash is None:
            reasons.append("missingness_mask_hash unavailable")
        if final_site_metadata is None:
            reasons.append("site metadata unavailable")
        else:
            if control_feature_column not in final_site_metadata.columns:
                reasons.append("control feature column missing")
            if control_feature_count < 1:
                reasons.append("no control features present")
        if final_sample_metadata is None:
            reasons.append("sample metadata unavailable")
        else:
            if replicate_group_column not in final_sample_metadata.columns:
                reasons.append("replicate group column missing")
            if replicate_group_count < 2:
                reasons.append("insufficient replicate groups")
            if batch_column is not None:
                if batch_column not in final_sample_metadata.columns:
                    reasons.append("batch column missing")
                elif (batch_count or 0) < 1:
                    reasons.append("no batch values present")

    return RuvReadinessState(
        enabled=enabled,
        ready=enabled and len(reasons) == 0,
        reasons=tuple(reasons),
        control_feature_column=control_feature_column,
        replicate_group_column=replicate_group_column,
        batch_column=batch_column,
        control_feature_count=control_feature_count,
        replicate_group_count=replicate_group_count,
        batch_count=batch_count,
        requires_complete_matrix=True,
        matrix_complete=matrix_complete,
        imputation_method_id=imputation_method_id,
        missingness_mask_preserved=missingness_mask_preserved,
    )


def _resolve_matrix_completeness(
    *,
    final_phospho: pd.DataFrame | None,
    default_matrix_complete: bool,
) -> bool:
    if final_phospho is None:
        return bool(default_matrix_complete)
    return int(final_phospho.isna().to_numpy().sum()) == 0


def _count_control_features(
    *,
    site_metadata: pd.DataFrame | None,
    control_feature_column: str,
) -> int:
    if site_metadata is None or control_feature_column not in site_metadata.columns:
        return 0
    series = site_metadata.loc[:, control_feature_column]
    return int(series.map(_is_truthy_control_feature).sum())


def _is_truthy_control_feature(value: object) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return bool(value)
    if isinstance(value, str):
        normalized = value.strip().lower()
        return normalized in {"1", "true", "t", "yes", "y"}
    return False


def _count_distinct_non_missing(
    *,
    sample_metadata: pd.DataFrame | None,
    column: str,
) -> int:
    if sample_metadata is None or column not in sample_metadata.columns:
        return 0
    values = (
        sample_metadata.loc[:, column]
        .astype("string")
        .str.strip()
        .replace("", pd.NA)
        .dropna()
    )
    return int(values.nunique())


def _resolve_total_correction_quantitative_meaning(
    quantitative_meaning: str | None,
) -> QuantitativeMeaning | None:
    if quantitative_meaning is None:
        return None
    try:
        return QuantitativeMeaning(quantitative_meaning)
    except ValueError as exc:
        supported = ", ".join(member.value for member in QuantitativeMeaning)
        raise DatasetBuildError(
            "dataset preprocessing total_protein_correction quantitative_meaning "
            f"must be one of: {supported}; got {quantitative_meaning!r}"
        ) from exc
