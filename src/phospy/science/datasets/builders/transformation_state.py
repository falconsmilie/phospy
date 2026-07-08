"""Transformation-state resolution collaborators for dataset builder execution."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import pandas as pd

from phospy.errors.build import DatasetBuildError
from phospy.errors.transformations import TransformationStateEstablishmentError
from phospy.science.datasets.builders.contracts import (
    DatasetIntensityScaleResolverContract,
    InterpretedDatasetBuildRequest,
    PreprocessedDatasetBuildTables,
)
from phospy.science.datasets.builders.preprocessing import (
    build_dataset_processing_state,
)
from phospy.science.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
    PreprocessingPlan,
    PreprocessingStageExecution,
)
from phospy.science.datasets.preprocessing.policy_models import IntensityTransformPolicy
from phospy.science.datasets.processing_state import DatasetProcessingState
from phospy.science.transformations.models import (
    DeclaredIntensityScaleDiagnosticPolicy,
    IntensityScaleEstablishmentMode,
    IntensityScaleEvidenceLevel,
    IntensityScaleKind,
    IntensityScaleState,
    IntensityTransformationEvent,
    MatrixIntensityScaleState,
)


@dataclass(frozen=True, slots=True)
class ResolvedDatasetTransformationState:
    """Post-preprocessing transformation state and resolved output matrices."""

    phospho: pd.DataFrame
    total: pd.DataFrame | None
    intensity_scale_state: IntensityScaleState
    processing_state: DatasetProcessingState
    quantitative_meaning: str
    intensity_scale_establishment: Mapping[str, object]


@dataclass(frozen=True, slots=True)
class _DeclaredInputIntensityScaleResolution:
    state: IntensityScaleState
    establishment_mode: IntensityScaleEstablishmentMode
    input_declaration_source: str | None
    establishment_parameters: Mapping[str, object]
    establishment_transformer_name: str | None


class DatasetTransformationStateResolver:
    """Resolve final dataset transformation state from preprocessing outputs."""

    def __init__(
        self,
        *,
        intensity_scale_resolver: DatasetIntensityScaleResolverContract,
    ) -> None:
        self._intensity_scale_resolver = intensity_scale_resolver

    def run(
        self,
        *,
        request: InterpretedDatasetBuildRequest,
        preprocessed: PreprocessedDatasetBuildTables,
        validated_site_metadata: pd.DataFrame,
    ) -> ResolvedDatasetTransformationState:
        intensity_transformation_event = (
            preprocessed.intensity_transformation_event
            if preprocessed.intensity_transformation_event is not None
            else _resolve_intensity_transformation_event(
                preprocessed.preprocessing_trace
            )
        )
        declared_input_scale_resolution = _resolve_declared_input_intensity_scale_resolution(
            preprocessing_plan=request.preprocessing_plan,
            preprocessing_trace=preprocessed.preprocessing_trace,
            intensity_transformation_event=intensity_transformation_event,
            declared_input_scale_kind=request.declared_input_intensity_scale_kind,
            declared_input_scale_source=request.declared_input_intensity_scale_source,
            has_total_matrix=preprocessed.total is not None,
        )
        declared_input_scale_state = (
            None
            if declared_input_scale_resolution is None
            else declared_input_scale_resolution.state
        )
        expected_scale_kind = _resolve_expected_intensity_scale_kind(
            request.preprocessing_plan,
            declared_input_scale_kind=request.declared_input_intensity_scale_kind,
        )
        resolved = self._intensity_scale_resolver.run(
            phospho=preprocessed.phospho,
            total=preprocessed.total,
            expected_scale_kind=expected_scale_kind,
            declared_input_scale_state=declared_input_scale_state,
            declared_input_establishment_mode=(
                None
                if declared_input_scale_resolution is None
                else declared_input_scale_resolution.establishment_mode
            ),
            input_declaration_source=(
                None
                if declared_input_scale_resolution is None
                else declared_input_scale_resolution.input_declaration_source
            ),
            scale_establishment_parameters=(
                None
                if declared_input_scale_resolution is None
                else dict(declared_input_scale_resolution.establishment_parameters)
            ),
            establishment_transformer_name=(
                None
                if declared_input_scale_resolution is None
                else declared_input_scale_resolution.establishment_transformer_name
            ),
            establishment_trace_id=_resolve_scale_establishment_trace_id(
                preprocessed.preprocessing_trace
            ),
            declared_scale_diagnostic_policy=_resolve_declared_scale_diagnostic_policy(
                request
            ),
        )
        if not resolved.intensity_scale_state.is_established:
            raise TransformationStateEstablishmentError(
                "intensity-scale resolver returned a non-established intensity scale "
                "state; this violates the dataset boundary contract"
            )
        processing_state = build_dataset_processing_state(
            plan=request.preprocessing_plan,
            intensity_scale_state=resolved.intensity_scale_state,
            explicit_quantitative_meaning=request.quantitative_meaning,
            preprocessing_trace=preprocessed.preprocessing_trace,
            final_phospho=resolved.phospho,
            final_site_metadata=validated_site_metadata,
            final_sample_metadata=preprocessed.sample_metadata,
        )
        intensity_scale_state = processing_state.intensity_scale
        quantitative_meaning = intensity_scale_state.quantity
        if quantitative_meaning is None:
            raise DatasetBuildError(
                "intensity-scale state is missing quantitative meaning"
            )
        establishment_provenance = intensity_scale_state.establishment_provenance
        if establishment_provenance is None:
            raise TransformationStateEstablishmentError(
                "intensity-scale state is established but missing establishment "
                "provenance"
            )
        return ResolvedDatasetTransformationState(
            phospho=resolved.phospho,
            total=resolved.total,
            intensity_scale_state=intensity_scale_state,
            processing_state=processing_state,
            quantitative_meaning=quantitative_meaning.value,
            intensity_scale_establishment=establishment_provenance.to_payload(),
        )


def _resolve_expected_intensity_scale_kind(
    preprocessing_plan: PreprocessingPlan,
    *,
    declared_input_scale_kind: IntensityScaleKind | None = None,
) -> IntensityScaleKind | None:
    if preprocessing_plan.intensity_transform_policy is IntensityTransformPolicy.LOG2:
        return IntensityScaleKind.LOG2
    if declared_input_scale_kind is not None:
        return declared_input_scale_kind
    return None


def _resolve_declared_scale_diagnostic_policy(
    request: InterpretedDatasetBuildRequest,
) -> DeclaredIntensityScaleDiagnosticPolicy:
    if request.allow_suspicious_declared_input_intensity_scale:
        return DeclaredIntensityScaleDiagnosticPolicy.WARN
    if request.declared_input_intensity_scale_kind is IntensityScaleKind.LOG2:
        return DeclaredIntensityScaleDiagnosticPolicy.ERROR
    return DeclaredIntensityScaleDiagnosticPolicy.WARN


def _resolve_declared_input_intensity_scale_resolution(
    *,
    preprocessing_plan: PreprocessingPlan,
    preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None,
    intensity_transformation_event: IntensityTransformationEvent | None,
    declared_input_scale_kind: IntensityScaleKind | None,
    declared_input_scale_source: str | None,
    has_total_matrix: bool,
) -> _DeclaredInputIntensityScaleResolution | None:
    intensity_transform_stage = _resolve_intensity_transform_stage(preprocessing_trace)
    if intensity_transformation_event is not None:
        return _DeclaredInputIntensityScaleResolution(
            state=_build_intensity_scale_state_from_event(
                intensity_transformation_event=intensity_transformation_event,
                has_total_matrix=has_total_matrix,
            ),
            establishment_mode=_resolve_event_establishment_mode(
                intensity_transformation_event
            ),
            input_declaration_source=None,
            establishment_parameters=_resolve_transformed_establishment_parameters(
                intensity_transform_stage=intensity_transform_stage,
                intensity_transformation_event=intensity_transformation_event,
            ),
            establishment_transformer_name=intensity_transformation_event.transformer_name,
        )
    if declared_input_scale_kind is not None:
        return _DeclaredInputIntensityScaleResolution(
            state=_build_declared_intensity_scale_state(
                kind=declared_input_scale_kind,
                has_total_matrix=has_total_matrix,
                established_by=(
                    "phospy.science.datasets.builders.executor.input_intensity_scale"
                ),
            ),
            establishment_mode=IntensityScaleEstablishmentMode.DECLARED,
            input_declaration_source=(
                declared_input_scale_source
                or "dataset_build_request.input_intensity_scale"
            ),
            establishment_parameters={
                "declared_scale_kind": declared_input_scale_kind.value,
            },
            establishment_transformer_name=None,
        )
    if (
        intensity_transform_stage is None
        and preprocessing_plan.intensity_transform_policy
        is IntensityTransformPolicy.LOG2
    ):
        return None
    return None


def _resolve_intensity_transform_stage(
    preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None,
) -> PreprocessingStageExecution | None:
    if preprocessing_trace is None:
        return None
    for stage in preprocessing_trace:
        if stage.stage != DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM:
            continue
        return stage
    return None


def _resolve_intensity_transformation_event(
    preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None,
) -> IntensityTransformationEvent | None:
    if preprocessing_trace is None:
        return None
    event: IntensityTransformationEvent | None = None
    for stage in preprocessing_trace:
        if stage.intensity_transformation_event is None:
            continue
        event = stage.intensity_transformation_event
    return event


def _build_intensity_scale_state_from_event(
    *,
    intensity_transformation_event: IntensityTransformationEvent,
    has_total_matrix: bool,
) -> IntensityScaleState:
    output_scale = intensity_transformation_event.output_scale
    total_state = (
        None
        if not has_total_matrix
        else MatrixIntensityScaleState(
            kind=output_scale.kind,
            transformed=output_scale.transformed,
            established_by=output_scale.established_by,
        )
    )
    return IntensityScaleState(
        phospho=output_scale,
        total=total_state,
    )


def _resolve_event_establishment_mode(
    event: IntensityTransformationEvent,
) -> IntensityScaleEstablishmentMode:
    if event.evidence_level is IntensityScaleEvidenceLevel.DECLARED_BY_USER:
        return IntensityScaleEstablishmentMode.DECLARED
    normalized_kind = event.transformation_kind.lower().replace("-", "_")
    if normalized_kind in {"identity", "passthrough", "pass_through"}:
        return IntensityScaleEstablishmentMode.IDENTITY
    return IntensityScaleEstablishmentMode.TRANSFORMED


def _resolve_transformed_establishment_parameters(
    *,
    intensity_transform_stage: PreprocessingStageExecution | None,
    intensity_transformation_event: IntensityTransformationEvent,
) -> Mapping[str, object]:
    operation = (
        intensity_transformation_event.transformation_kind
        if intensity_transform_stage is None
        else str(intensity_transform_stage.operation)
    )
    parameters: dict[str, object] = {
        "operation": operation,
        "transformation_kind": intensity_transformation_event.transformation_kind,
        "input_scale": intensity_transformation_event.input_scale.kind.value,
        "output_scale": intensity_transformation_event.output_scale.kind.value,
    }
    if intensity_transform_stage is not None:
        parameters.update(dict(intensity_transform_stage.parameters))
    if intensity_transformation_event.pseudocount is not None:
        parameters["pseudocount"] = intensity_transformation_event.pseudocount
    if intensity_transformation_event.input_fingerprint is not None:
        parameters["input_fingerprint"] = (
            intensity_transformation_event.input_fingerprint
        )
    if intensity_transformation_event.output_fingerprint is not None:
        parameters["output_fingerprint"] = (
            intensity_transformation_event.output_fingerprint
        )
    diagnostics: Mapping[str, object] = {}
    if intensity_transform_stage is not None and isinstance(
        intensity_transform_stage.diagnostics,
        Mapping,
    ):
        diagnostics = intensity_transform_stage.diagnostics
    for optional_key in ("affected_matrices",):
        value = diagnostics.get(optional_key)
        if value is not None:
            parameters[optional_key] = value
    return parameters


def _build_declared_intensity_scale_state(
    *,
    kind: IntensityScaleKind,
    has_total_matrix: bool,
    established_by: str,
) -> IntensityScaleState:
    if kind is IntensityScaleKind.LOG2:
        phospho_state = MatrixIntensityScaleState.log2(established_by=established_by)
        if has_total_matrix:
            return IntensityScaleState(
                phospho=phospho_state,
                total=MatrixIntensityScaleState.log2(established_by=established_by),
            )
        return IntensityScaleState(phospho=phospho_state, total=None)
    phospho_state = MatrixIntensityScaleState.linear(established_by=established_by)
    if has_total_matrix:
        return IntensityScaleState(
            phospho=phospho_state,
            total=MatrixIntensityScaleState.linear(established_by=established_by),
        )
    return IntensityScaleState(phospho=phospho_state, total=None)


def _resolve_scale_establishment_trace_id(
    preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None,
) -> str | None:
    if preprocessing_trace is None or not preprocessing_trace:
        return None
    final_stage = preprocessing_trace[-1]
    return f"{final_stage.stage}:{final_stage.operation}:{final_stage.output_hash}"
