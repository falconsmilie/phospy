"""Transformation-state resolution collaborators for dataset builder execution."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import TypedDict

import pandas as pd

from phospy.errors.build import DatasetBuildError
from phospy.errors.transformations import TransformationStateEstablishmentError
from phospy.science.datasets.builders.contracts import (
    InterpretedDatasetBuildRequest,
    PreprocessedDatasetBuildTables,
)
from phospy.science.datasets.builders.preprocessing import (
    build_dataset_processing_state,
)
from phospy.science.datasets.builders.transformation_resolver import (
    DatasetIntensityScaleResolver,
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
    IntensityScaleKind,
    IntensityScaleState,
    MatrixIntensityScaleState,
    QuantitativeMeaning,
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


class _DeclaredScaleResolverKwargs(TypedDict, total=False):
    declared_input_establishment_mode: IntensityScaleEstablishmentMode
    declared_scale_diagnostic_policy: DeclaredIntensityScaleDiagnosticPolicy
    input_declaration_source: str | None
    scale_establishment_parameters: Mapping[str, object]
    establishment_transformer_name: str | None
    establishment_trace_id: str | None


class DatasetTransformationStateResolver:
    """Resolve final dataset transformation state from preprocessing outputs."""

    def __init__(
        self,
        *,
        intensity_scale_resolver: DatasetIntensityScaleResolver,
    ) -> None:
        self._intensity_scale_resolver = intensity_scale_resolver

    def run(
        self,
        *,
        request: InterpretedDatasetBuildRequest,
        preprocessed: PreprocessedDatasetBuildTables,
        validated_site_metadata: pd.DataFrame,
    ) -> ResolvedDatasetTransformationState:
        declared_input_scale_resolution = _resolve_declared_input_intensity_scale_resolution(
            preprocessing_plan=request.preprocessing_plan,
            preprocessing_trace=preprocessed.preprocessing_trace,
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
        declared_scale_kwargs = _resolve_declared_scale_resolver_kwargs(
            resolver=self._intensity_scale_resolver,
            resolution=declared_input_scale_resolution,
            preprocessing_trace=preprocessed.preprocessing_trace,
            declared_scale_diagnostic_policy=(
                _resolve_declared_scale_diagnostic_policy(request)
            ),
        )
        if _resolver_supports_declared_input_scale_state(
            self._intensity_scale_resolver
        ):
            resolved = self._intensity_scale_resolver.run(
                phospho=preprocessed.phospho,
                total=preprocessed.total,
                expected_scale_kind=expected_scale_kind,
                declared_input_scale_state=declared_input_scale_state,
                **declared_scale_kwargs,
            )
        else:
            resolved = self._intensity_scale_resolver.run(
                phospho=preprocessed.phospho,
                total=preprocessed.total,
                expected_scale_kind=expected_scale_kind,
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
    declared_input_scale_kind: IntensityScaleKind | None,
    declared_input_scale_source: str | None,
    has_total_matrix: bool,
) -> _DeclaredInputIntensityScaleResolution | None:
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
    intensity_transform_stage = _resolve_intensity_transform_stage(preprocessing_trace)
    transformed_state = _resolve_transformed_state_from_intensity_stage(
        intensity_transform_stage=intensity_transform_stage,
        has_total_matrix=has_total_matrix,
    )
    if transformed_state is not None:
        return _DeclaredInputIntensityScaleResolution(
            state=transformed_state,
            establishment_mode=IntensityScaleEstablishmentMode.TRANSFORMED,
            input_declaration_source=None,
            establishment_parameters=_resolve_transformed_establishment_parameters(
                intensity_transform_stage=intensity_transform_stage
            ),
            establishment_transformer_name=_resolve_transformed_transformer_name(
                intensity_transform_stage=intensity_transform_stage
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


def _resolve_transformed_state_from_intensity_stage(
    *,
    intensity_transform_stage: PreprocessingStageExecution | None,
    has_total_matrix: bool,
) -> IntensityScaleState | None:
    if intensity_transform_stage is None:
        return None
    diagnostics = intensity_transform_stage.diagnostics
    if not isinstance(diagnostics, Mapping):
        return None
    transformer_state = diagnostics.get("transformer_state")
    if not isinstance(transformer_state, Mapping):
        return None
    phospho_payload = transformer_state.get("phospho")
    if not isinstance(phospho_payload, Mapping):
        raise DatasetBuildError(
            "intensity-transform stage diagnostics is missing transformer_state.phospho"
        )
    phospho_state = _matrix_state_from_payload(
        phospho_payload,
        field_name="transformer_state.phospho",
    )
    total_payload = transformer_state.get("total")
    if total_payload is None:
        if has_total_matrix:
            raise DatasetBuildError(
                "intensity-transform stage diagnostics omitted transformer_state.total "
                "for a dataset with total input matrix"
            )
        total_state = None
    else:
        if not has_total_matrix:
            raise DatasetBuildError(
                "intensity-transform stage diagnostics reported transformer_state.total "
                "for a dataset without total input matrix"
            )
        if not isinstance(total_payload, Mapping):
            raise DatasetBuildError(
                "intensity-transform stage diagnostics contains invalid "
                "transformer_state.total payload"
            )
        total_state = _matrix_state_from_payload(
            total_payload,
            field_name="transformer_state.total",
        )
    quantity_payload = transformer_state.get("quantity")
    quantity = _resolve_transformer_state_quantity(
        quantity_payload,
        field_name="transformer_state.quantity",
    )
    return IntensityScaleState(
        phospho=phospho_state,
        total=total_state,
        quantity=quantity,
    )


def _matrix_state_from_payload(
    payload: Mapping[str, object],
    *,
    field_name: str,
) -> MatrixIntensityScaleState:
    kind_payload = payload.get("kind")
    transformed_payload = payload.get("transformed")
    established_by_payload = payload.get("established_by")
    if not isinstance(kind_payload, str) or not kind_payload.strip():
        raise DatasetBuildError(
            f"intensity-transform stage diagnostics {field_name}.kind must be a "
            "non-empty string"
        )
    if not isinstance(transformed_payload, bool):
        raise DatasetBuildError(
            f"intensity-transform stage diagnostics {field_name}.transformed must be a "
            "boolean"
        )
    if (
        not isinstance(established_by_payload, str)
        or not established_by_payload.strip()
    ):
        raise DatasetBuildError(
            f"intensity-transform stage diagnostics {field_name}.established_by "
            "must be a non-empty string"
        )
    try:
        kind = IntensityScaleKind(str(kind_payload))
    except ValueError as exc:
        supported = ", ".join(item.value for item in IntensityScaleKind)
        raise DatasetBuildError(
            f"intensity-transform stage diagnostics {field_name}.kind must be one of: "
            f"{supported}"
        ) from exc
    return MatrixIntensityScaleState(
        kind=kind,
        transformed=transformed_payload,
        established_by=established_by_payload,
    )


def _resolve_transformer_state_quantity(
    value: object,
    *,
    field_name: str,
) -> QuantitativeMeaning | None:
    if value is None:
        return None
    normalized = str(value).strip()
    if not normalized:
        return None
    try:
        return QuantitativeMeaning(normalized)
    except ValueError as exc:
        supported = ", ".join(item.value for item in QuantitativeMeaning)
        raise DatasetBuildError(
            f"intensity-transform stage diagnostics {field_name} must be one of: "
            f"{supported}"
        ) from exc


def _resolve_transformed_establishment_parameters(
    *,
    intensity_transform_stage: PreprocessingStageExecution | None,
) -> Mapping[str, object]:
    if intensity_transform_stage is None:
        return {}
    diagnostics = (
        {}
        if not isinstance(intensity_transform_stage.diagnostics, Mapping)
        else intensity_transform_stage.diagnostics
    )
    parameters: dict[str, object] = {
        "operation": str(intensity_transform_stage.operation),
        **dict(intensity_transform_stage.parameters),
    }
    for optional_key in ("pseudocount", "affected_matrices"):
        value = diagnostics.get(optional_key)
        if value is not None:
            parameters[optional_key] = value
    return parameters


def _resolve_transformed_transformer_name(
    *,
    intensity_transform_stage: PreprocessingStageExecution | None,
) -> str | None:
    if intensity_transform_stage is None:
        return None
    diagnostics = intensity_transform_stage.diagnostics
    if not isinstance(diagnostics, Mapping):
        return None
    value = diagnostics.get("transformer_name")
    if value is None:
        return None
    normalized = str(value).strip()
    if not normalized:
        return None
    return normalized


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


def _resolver_supports_declared_input_scale_state(resolver: object) -> bool:
    run_method = getattr(resolver, "run", None)
    if run_method is None:
        return False
    code_object = getattr(run_method, "__code__", None)
    if code_object is None:
        return False
    return "declared_input_scale_state" in code_object.co_varnames


def _resolve_declared_scale_resolver_kwargs(
    *,
    resolver: object,
    resolution: _DeclaredInputIntensityScaleResolution | None,
    preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None,
    declared_scale_diagnostic_policy: DeclaredIntensityScaleDiagnosticPolicy,
) -> _DeclaredScaleResolverKwargs:
    run_method = getattr(resolver, "run", None)
    if run_method is None:
        return _DeclaredScaleResolverKwargs()
    code_object = getattr(run_method, "__code__", None)
    if code_object is None:
        return _DeclaredScaleResolverKwargs()
    supported_parameters = set(code_object.co_varnames)
    kwargs = _DeclaredScaleResolverKwargs()
    if "declared_scale_diagnostic_policy" in supported_parameters:
        kwargs["declared_scale_diagnostic_policy"] = declared_scale_diagnostic_policy
    if resolution is not None:
        if "declared_input_establishment_mode" in supported_parameters:
            kwargs["declared_input_establishment_mode"] = resolution.establishment_mode
        if "input_declaration_source" in supported_parameters:
            kwargs["input_declaration_source"] = resolution.input_declaration_source
        if "scale_establishment_parameters" in supported_parameters:
            kwargs["scale_establishment_parameters"] = dict(
                resolution.establishment_parameters
            )
        if "establishment_transformer_name" in supported_parameters:
            kwargs["establishment_transformer_name"] = (
                resolution.establishment_transformer_name
            )
    if "establishment_trace_id" in supported_parameters:
        kwargs["establishment_trace_id"] = _resolve_scale_establishment_trace_id(
            preprocessing_trace
        )
    return kwargs


def _resolve_scale_establishment_trace_id(
    preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None,
) -> str | None:
    if preprocessing_trace is None or not preprocessing_trace:
        return None
    final_stage = preprocessing_trace[-1]
    return f"{final_stage.stage}:{final_stage.operation}:{final_stage.output_hash}"
