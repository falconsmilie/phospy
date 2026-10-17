"""Capability-owned method selection and preprocessing-stage execution."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import cast

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.preparation_records.batch_correction import (
    BATCH_CORRECTION_CONFOUNDING_NOT_APPLICABLE,
    BATCH_CORRECTION_DESIGN_PRESERVATION_PRESERVE_CONDITION_EFFECTS,
    BATCH_CORRECTION_STATUS_DISABLED,
    BatchCorrectionDiagnostics,
    BatchCorrectionPolicy,
    BatchCorrectionReport,
)
from phospy.preprocessing.batch_correction.contracts.configuration import (
    DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH,
    DATASET_BATCH_CORRECTION_METHOD_NONE,
    SPS_RUV_BATCH_CORRECTION_METHODS,
    DatasetBatchCorrectionConfig,
)
from phospy.preprocessing.batch_correction.contracts.output import (
    CorrectedPreprocessingOutput,
)
from phospy.preprocessing.batch_correction.contracts.stage import (
    BatchCorrectionAdequacyValidatorProtocol,
    BatchDesignMetadataValidatorProtocol,
    SpsRuvStyleBatchCorrectionRunner,
)
from phospy.preprocessing.batch_correction.execution.linear import (
    LinearResidualizeBatchCorrectionEngine,
)
from phospy.preprocessing.batch_correction.metadata import (
    BatchCorrectionMetadataResolver,
    ResolvedBatchCorrectionMetadata,
)
from phospy.preprocessing.batch_correction.provenance.assembly import (
    build_native_batch_correction_provenance,
)
from phospy.preprocessing.control_sites.contracts import ControlSiteSet
from phospy.provenance.models._shared import JsonValue
from phospy.provenance.models.batch_correction import BatchCorrectionProvenance


@dataclass(frozen=True, slots=True)
class BatchCorrectionStageExecutionRequest:
    """Passive input assembled by the core pipeline stage adapter."""

    phospho: pd.DataFrame
    sample_metadata: pd.DataFrame | None
    site_metadata: pd.DataFrame
    batch_correction_method: str
    batch_correction_batch_column: str
    batch_correction_condition_column: str
    batch_correction_condition_columns: tuple[str, ...]
    batch_correction_preserve_condition_effects: bool
    batch_correction_internal_request: object | None
    batch_correction_control_site_set: object | None
    batch_correction_missingness_policy: object | None
    imputation_observation_mask: pd.DataFrame | None
    stage_order: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class BatchCorrectionStageExecutionResult:
    """Passive capability result returned to core pipeline composition."""

    corrected_matrix: pd.DataFrame
    output_observation_mask: pd.DataFrame | None
    metadata: ResolvedBatchCorrectionMetadata | None
    report: BatchCorrectionReport
    diagnostics: Mapping[str, JsonValue]
    provenance: BatchCorrectionProvenance | None


class _MissingBatchDesignMetadataValidator:
    def run(self, **_: object) -> object:
        raise PhosPyInputError(
            "BatchCorrectionStage requires an injected BatchDesignMetadataValidator"
        )


class _MissingBatchCorrectionAdequacyValidator:
    def run(self, **_: object) -> object:
        raise PhosPyInputError(
            "BatchCorrectionStage requires an injected BatchCorrectionAdequacyValidator"
        )


class BatchCorrectionStageExecutor:
    """Select and execute one supported correction method using passive inputs."""

    def __init__(
        self,
        *,
        metadata_resolver: BatchCorrectionMetadataResolver | None = None,
        metadata_validator: BatchDesignMetadataValidatorProtocol | None = None,
        adequacy_validator: BatchCorrectionAdequacyValidatorProtocol | None = None,
        engine: LinearResidualizeBatchCorrectionEngine | None = None,
        sps_ruv_runner: SpsRuvStyleBatchCorrectionRunner | None = None,
    ) -> None:
        self._metadata_resolver = metadata_resolver or BatchCorrectionMetadataResolver()
        self._metadata_validator = (
            metadata_validator or _MissingBatchDesignMetadataValidator()
        )
        self._adequacy_validator = (
            adequacy_validator or _MissingBatchCorrectionAdequacyValidator()
        )
        self._engine = engine or LinearResidualizeBatchCorrectionEngine()
        self._sps_ruv_runner = sps_ruv_runner

    def run(
        self,
        request: BatchCorrectionStageExecutionRequest,
    ) -> BatchCorrectionStageExecutionResult:
        method = resolve_batch_correction_method(request.batch_correction_method)
        if method == DATASET_BATCH_CORRECTION_METHOD_NONE:
            report = _build_disabled_report(request)
            return BatchCorrectionStageExecutionResult(
                corrected_matrix=request.phospho,
                output_observation_mask=request.imputation_observation_mask,
                metadata=None,
                report=report,
                diagnostics={
                    "method": method,
                    "status": BATCH_CORRECTION_STATUS_DISABLED,
                    "matrix_shape_before": list(report.matrix_shape_before or ()),
                    "matrix_shape_after": list(report.matrix_shape_after or ()),
                },
                provenance=None,
            )
        if method in SPS_RUV_BATCH_CORRECTION_METHODS:
            return self._run_sps_ruv_style(request, method=method)
        if method != DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH:
            raise PhosPyInputError(
                "dataset preprocessing plan contains unsupported "
                f"batch_correction_method={method!r}"
            )

        try:
            metadata = self._resolve_validated_linear_metadata(request)
            result = self._engine.run(
                phospho=request.phospho,
                batch_labels=metadata.batch_labels,
                condition_labels=metadata.condition_labels,
                config=DatasetBatchCorrectionConfig(
                    method=method,
                    batch_column=request.batch_correction_batch_column,
                    condition_column=request.batch_correction_condition_column,
                    preserve_condition_effects=True,
                ),
            )
        except PhosPyInputError as exc:
            raise PhosPyInputError(
                _validation_failure_message(request, method=method, reason=str(exc))
            ) from exc
        provenance = build_native_batch_correction_provenance(
            input_matrix=request.phospho,
            output_matrix=result.corrected_matrix,
            plan=request,
            report=result.report,
            metadata=metadata,
            diagnostics=result.diagnostics,
            warnings=result.report.warnings,
            observation_mask=request.imputation_observation_mask,
            source="native_preprocessing_stage",
        )
        return BatchCorrectionStageExecutionResult(
            corrected_matrix=result.corrected_matrix,
            output_observation_mask=request.imputation_observation_mask,
            metadata=metadata,
            report=result.report,
            diagnostics=cast(dict[str, JsonValue], dict(result.diagnostics)),
            provenance=provenance,
        )

    def validate_linear(
        self,
        request: BatchCorrectionStageExecutionRequest,
    ) -> None:
        method = resolve_batch_correction_method(request.batch_correction_method)
        if method != DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH:
            return
        try:
            self._resolve_validated_linear_metadata(request)
        except PhosPyInputError as exc:
            raise PhosPyInputError(
                _validation_failure_message(request, method=method, reason=str(exc))
            ) from exc

    def _resolve_validated_linear_metadata(
        self,
        request: BatchCorrectionStageExecutionRequest,
    ) -> ResolvedBatchCorrectionMetadata:
        self._metadata_validator.run(
            phospho=request.phospho,
            sample_metadata=request.sample_metadata,
            batch_column=request.batch_correction_batch_column,
            condition_columns=(request.batch_correction_condition_column,),
            context="dataset build request preprocessing_config.batch_correction",
        )
        metadata = self._metadata_resolver.run(
            phospho=request.phospho,
            sample_metadata=request.sample_metadata,
            batch_column=request.batch_correction_batch_column,
            condition_column=request.batch_correction_condition_column,
        )
        self._adequacy_validator.run(
            batch_by_sample=metadata.batch_by_sample,
            condition_by_sample=metadata.condition_by_sample,
            sample_order=metadata.sample_order,
            preserve_condition_effects=(
                request.batch_correction_preserve_condition_effects
            ),
        )
        return metadata

    def _run_sps_ruv_style(
        self,
        request: BatchCorrectionStageExecutionRequest,
        *,
        method: str,
    ) -> BatchCorrectionStageExecutionResult:
        internal_request = request.batch_correction_internal_request
        control_site_set = request.batch_correction_control_site_set
        missingness_policy = request.batch_correction_missingness_policy
        if (
            internal_request is None
            or control_site_set is None
            or missingness_policy is None
        ):
            raise PhosPyInputError(
                "dataset preprocessing plan contains incomplete SPS/RUV-style "
                "batch-correction configuration"
            )
        if not isinstance(control_site_set, ControlSiteSet):
            raise PhosPyInputError(
                "dataset preprocessing plan SPS/RUV-style batch correction requires "
                "a ControlSiteSet"
            )
        if self._sps_ruv_runner is None:
            raise PhosPyInputError(
                "dataset preprocessing plan SPS/RUV-style batch correction requires "
                "an injected SpsRuvStyleBatchCorrectionRunner"
            )
        try:
            result = self._sps_ruv_runner.run(
                phospho=request.phospho,
                config=internal_request,
                sample_metadata=request.sample_metadata,
                control_site_set=control_site_set,
                missingness_policy=missingness_policy,
                upstream_observation_mask=request.imputation_observation_mask,
                site_metadata=request.site_metadata,
            )
        except PhosPyInputError as exc:
            raise PhosPyInputError(
                _validation_failure_message(request, method=method, reason=str(exc))
            ) from exc

        corrected_output = result.corrected_preprocessing_output
        if not isinstance(corrected_output, CorrectedPreprocessingOutput):
            raise PhosPyInputError(
                "SPS/RUV-style batch correction did not produce complete "
                "analysis-ready preprocessing output; check missingness policy and "
                "observation-mask diagnostics"
            )
        return BatchCorrectionStageExecutionResult(
            corrected_matrix=corrected_output.corrected_matrix,
            output_observation_mask=corrected_output.output_observation_mask,
            metadata=None,
            report=corrected_output.batch_correction_report,
            diagnostics=dict(result.diagnostics),
            provenance=result.provenance,
        )


def resolve_batch_correction_method(value: object) -> str:
    """Resolve the configured method while preserving the historical no-op default."""

    method = str(value).strip()
    if not method:
        return DATASET_BATCH_CORRECTION_METHOD_NONE
    return method


def _build_disabled_report(
    request: BatchCorrectionStageExecutionRequest,
) -> BatchCorrectionReport:
    shape = (int(request.phospho.shape[0]), int(request.phospho.shape[1]))
    return BatchCorrectionReport(
        status=BATCH_CORRECTION_STATUS_DISABLED,
        policy=BatchCorrectionPolicy(
            method=DATASET_BATCH_CORRECTION_METHOD_NONE,
            batch_column=request.batch_correction_batch_column,
            condition_column=request.batch_correction_condition_column,
            condition_columns=request.batch_correction_condition_columns,
            design_preservation_policy=(
                BATCH_CORRECTION_DESIGN_PRESERVATION_PRESERVE_CONDITION_EFFECTS
            ),
            preserve_condition_effects=(
                request.batch_correction_preserve_condition_effects
            ),
        ),
        diagnostics=BatchCorrectionDiagnostics(
            confounding_check_status=BATCH_CORRECTION_CONFOUNDING_NOT_APPLICABLE,
            matrix_shape_before=shape,
            matrix_shape_after=shape,
            limitations=("batch correction disabled by preprocessing configuration",),
        ),
    )


def _validation_failure_message(
    request: BatchCorrectionStageExecutionRequest,
    *,
    method: str,
    reason: str,
) -> str:
    return (
        "batch correction validation failed before correction execution; "
        f"method={method!r}; "
        f"input_shape={[int(request.phospho.shape[0]), int(request.phospho.shape[1])]}; "
        f"batch_column={request.batch_correction_batch_column!r}; "
        f"condition_column={request.batch_correction_condition_column!r}; "
        f"stage_order={list(request.stage_order)!r}; "
        f"reason: {reason}"
    )


__all__ = [
    "BatchCorrectionStageExecutionRequest",
    "BatchCorrectionStageExecutionResult",
    "BatchCorrectionStageExecutor",
    "resolve_batch_correction_method",
]
