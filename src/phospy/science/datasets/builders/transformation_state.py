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
    quantitative_meaning_provenance: Mapping[str, object]


@dataclass(frozen=True, slots=True)
class _DeclaredInputIntensityScaleResolution:
    state: IntensityScaleState
    establishment_mode: IntensityScaleEstablishmentMode
    evidence_level: IntensityScaleEvidenceLevel
    input_declaration_source: str | None
    establishment_parameters: Mapping[str, object]
    establishment_transformer_name: str | None


@dataclass(frozen=True, slots=True)
class _SourcedIntensityTransformationEvent:
    source: str
    event: IntensityTransformationEvent


_INTENSITY_TRANSFORMATION_EVENT_DEDUP_FIELDS = (
    "transformer_name",
    "input_scale",
    "output_scale",
    "evidence_level",
    "transformation_kind",
    "pseudocount",
    "input_fingerprint",
    "output_fingerprint",
)


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
        intensity_transformation_events = _resolve_intensity_transformation_events(
            preprocessed=preprocessed
        )
        observed_intensity_transformation_event = (
            _select_intensity_transformation_event(
                intensity_transformation_events,
                evidence_level=IntensityScaleEvidenceLevel.OBSERVED_TRANSFORMATION,
            )
        )
        inferred_intensity_transformation_event = (
            _select_intensity_transformation_event(
                intensity_transformation_events,
                evidence_level=IntensityScaleEvidenceLevel.INFERRED_FROM_METADATA,
            )
        )
        _reject_missing_observed_intensity_transformation_event(
            preprocessing_plan=request.preprocessing_plan,
            preprocessing_trace=preprocessed.preprocessing_trace,
            observed_intensity_transformation_event=(
                observed_intensity_transformation_event
            ),
        )
        declared_input_scale_resolution = _resolve_declared_input_intensity_scale_resolution(
            preprocessing_plan=request.preprocessing_plan,
            preprocessing_trace=preprocessed.preprocessing_trace,
            observed_intensity_transformation_event=(
                observed_intensity_transformation_event
            ),
            inferred_intensity_transformation_event=(
                inferred_intensity_transformation_event
            ),
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
            scale_establishment_evidence_level=(
                None
                if declared_input_scale_resolution is None
                else declared_input_scale_resolution.evidence_level
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
        quantitative_meaning_provenance = (
            intensity_scale_state.quantitative_meaning_provenance
        )
        if quantitative_meaning_provenance is None:
            raise TransformationStateEstablishmentError(
                "intensity-scale state is missing quantitative meaning provenance"
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
            quantitative_meaning_provenance=(
                quantitative_meaning_provenance.to_payload()
            ),
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
    observed_intensity_transformation_event: IntensityTransformationEvent | None,
    inferred_intensity_transformation_event: IntensityTransformationEvent | None,
    declared_input_scale_kind: IntensityScaleKind | None,
    declared_input_scale_source: str | None,
    has_total_matrix: bool,
) -> _DeclaredInputIntensityScaleResolution | None:
    intensity_transform_stage = _resolve_intensity_transform_stage(preprocessing_trace)
    if observed_intensity_transformation_event is not None:
        return _DeclaredInputIntensityScaleResolution(
            state=_build_intensity_scale_state_from_event(
                intensity_transformation_event=observed_intensity_transformation_event,
                has_total_matrix=has_total_matrix,
            ),
            establishment_mode=_resolve_event_establishment_mode(
                observed_intensity_transformation_event
            ),
            evidence_level=observed_intensity_transformation_event.evidence_level,
            input_declaration_source=None,
            establishment_parameters=_resolve_transformed_establishment_parameters(
                intensity_transform_stage=intensity_transform_stage,
                intensity_transformation_event=observed_intensity_transformation_event,
                has_total_matrix=has_total_matrix,
            ),
            establishment_transformer_name=(
                observed_intensity_transformation_event.transformer_name
            ),
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
            evidence_level=IntensityScaleEvidenceLevel.DECLARED_BY_USER,
            input_declaration_source=(
                declared_input_scale_source
                or "dataset_build_request.input_intensity_scale"
            ),
            establishment_parameters={
                "declared_scale_kind": declared_input_scale_kind.value,
            },
            establishment_transformer_name=None,
        )
    if inferred_intensity_transformation_event is not None:
        return _DeclaredInputIntensityScaleResolution(
            state=_build_intensity_scale_state_from_event(
                intensity_transformation_event=inferred_intensity_transformation_event,
                has_total_matrix=has_total_matrix,
            ),
            establishment_mode=_resolve_event_establishment_mode(
                inferred_intensity_transformation_event
            ),
            evidence_level=inferred_intensity_transformation_event.evidence_level,
            input_declaration_source=None,
            establishment_parameters=_resolve_transformed_establishment_parameters(
                intensity_transform_stage=intensity_transform_stage,
                intensity_transformation_event=inferred_intensity_transformation_event,
                has_total_matrix=has_total_matrix,
            ),
            establishment_transformer_name=(
                inferred_intensity_transformation_event.transformer_name
            ),
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


def _resolve_intensity_transformation_events(
    *,
    preprocessed: PreprocessedDatasetBuildTables,
) -> tuple[IntensityTransformationEvent, ...]:
    events: list[_SourcedIntensityTransformationEvent] = []
    for index, stage in enumerate(preprocessed.preprocessing_trace or ()):
        if stage.intensity_transformation_event is None:
            continue
        source = (
            f"preprocessed.preprocessing_trace[{index}].intensity_transformation_event"
        )
        events.append(
            _SourcedIntensityTransformationEvent(
                source=source,
                event=_require_intensity_transformation_event(
                    stage.intensity_transformation_event,
                    source=source,
                ),
            )
        )
    if preprocessed.intensity_transformation_event is not None:
        events.append(
            _SourcedIntensityTransformationEvent(
                source="preprocessed.intensity_transformation_event",
                event=_require_intensity_transformation_event(
                    preprocessed.intensity_transformation_event,
                    source="preprocessed.intensity_transformation_event",
                ),
            )
        )
    return _deduplicate_or_reject_intensity_transformation_events(tuple(events))


def _deduplicate_or_reject_intensity_transformation_events(
    events: tuple[_SourcedIntensityTransformationEvent, ...],
) -> tuple[IntensityTransformationEvent, ...]:
    selected_by_evidence: dict[
        IntensityScaleEvidenceLevel, _SourcedIntensityTransformationEvent
    ] = {}
    unique_events: list[IntensityTransformationEvent] = []
    for sourced in events:
        evidence_level = sourced.event.evidence_level
        existing = selected_by_evidence.get(evidence_level)
        if existing is None:
            selected_by_evidence[evidence_level] = sourced
            unique_events.append(sourced.event)
            continue
        differing_fields = _differing_intensity_transformation_event_fields(
            existing.event,
            sourced.event,
        )
        if differing_fields:
            _raise_intensity_transformation_event_conflict(
                evidence_level=evidence_level,
                first_source=existing.source,
                conflicting_source=sourced.source,
                differing_fields=differing_fields,
                first_event=existing.event,
                conflicting_event=sourced.event,
            )
        # Duplicate event from trace/top-level provenance surfaces.
    return tuple(unique_events)


def _intensity_transformation_event_comparison_payload(
    event: IntensityTransformationEvent,
) -> dict[str, object]:
    payload = event.to_payload()
    return {
        field_name: payload[field_name]
        for field_name in _INTENSITY_TRANSFORMATION_EVENT_DEDUP_FIELDS
    }


def _differing_intensity_transformation_event_fields(
    left: IntensityTransformationEvent,
    right: IntensityTransformationEvent,
) -> tuple[str, ...]:
    left_payload = _intensity_transformation_event_comparison_payload(left)
    right_payload = _intensity_transformation_event_comparison_payload(right)
    return tuple(
        field_name
        for field_name in _INTENSITY_TRANSFORMATION_EVENT_DEDUP_FIELDS
        if left_payload[field_name] != right_payload[field_name]
    )


def _raise_intensity_transformation_event_conflict(
    *,
    evidence_level: IntensityScaleEvidenceLevel,
    first_source: str,
    conflicting_source: str,
    differing_fields: tuple[str, ...],
    first_event: IntensityTransformationEvent,
    conflicting_event: IntensityTransformationEvent,
) -> None:
    field_label = "/".join(differing_fields)
    raise DatasetBuildError(
        "dataset preprocessing intensity transformation event conflict: "
        f"evidence_level={evidence_level.value!r} appears more than once with "
        f"different {field_label}; sources: {first_source} and "
        f"{conflicting_source}; transformer_names: "
        f"first={first_event.transformer_name!r}, "
        f"conflicting={conflicting_event.transformer_name!r}"
    )


def _require_intensity_transformation_event(
    value: object,
    *,
    source: str,
) -> IntensityTransformationEvent:
    if isinstance(value, IntensityTransformationEvent):
        return value
    raise DatasetBuildError(
        "dataset preprocessing intensity transformation event parse error: "
        f"{source} must be IntensityTransformationEvent or None, got "
        f"{value!r} ({type(value).__name__})"
    )


def _select_intensity_transformation_event(
    events: tuple[IntensityTransformationEvent, ...],
    *,
    evidence_level: IntensityScaleEvidenceLevel,
) -> IntensityTransformationEvent | None:
    matching_events: list[IntensityTransformationEvent] = []
    for event in events:
        if event.evidence_level is evidence_level:
            matching_events.append(event)
    if not matching_events:
        return None
    first_event = matching_events[0]
    for index, event in enumerate(matching_events[1:], start=1):
        differing_fields = _differing_intensity_transformation_event_fields(
            first_event,
            event,
        )
        if differing_fields:
            _raise_intensity_transformation_event_conflict(
                evidence_level=evidence_level,
                first_source="events[0]",
                conflicting_source=f"events[{index}]",
                differing_fields=differing_fields,
                first_event=first_event,
                conflicting_event=event,
            )
    return first_event


def _reject_missing_observed_intensity_transformation_event(
    *,
    preprocessing_plan: PreprocessingPlan,
    preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None,
    observed_intensity_transformation_event: IntensityTransformationEvent | None,
) -> None:
    if observed_intensity_transformation_event is not None:
        return
    intensity_transform_stage = _resolve_intensity_transform_stage(preprocessing_trace)
    if intensity_transform_stage is None:
        return
    operation = str(intensity_transform_stage.operation).strip()
    if (
        preprocessing_plan.intensity_transform_policy
        is not IntensityTransformPolicy.LOG2
        and operation != IntensityTransformPolicy.LOG2.value
    ):
        return
    raise TransformationStateEstablishmentError(
        "preprocessing intensity_transform policy='log2' did not emit a typed "
        "observed IntensityTransformationEvent; dataset intensity scale state "
        "cannot be established from diagnostics-only transformation metadata"
    )


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
    if event.evidence_level is IntensityScaleEvidenceLevel.INFERRED_FROM_METADATA:
        return IntensityScaleEstablishmentMode.DERIVED
    normalized_kind = event.transformation_kind.lower().replace("-", "_")
    if normalized_kind in {"identity", "passthrough", "pass_through"}:
        return IntensityScaleEstablishmentMode.IDENTITY
    return IntensityScaleEstablishmentMode.TRANSFORMED


def _resolve_transformed_establishment_parameters(
    *,
    intensity_transform_stage: PreprocessingStageExecution | None,
    intensity_transformation_event: IntensityTransformationEvent,
    has_total_matrix: bool,
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
        "affected_matrices": ["phospho", "total"] if has_total_matrix else ["phospho"],
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
