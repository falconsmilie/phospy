"""Batch-correction preprocessing intent configuration."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from phospy.contracts.configs.preprocessing.correction_missingness import (
    CorrectionMissingnessPolicy,
    TemporaryImputationMethod,
)
from phospy.contracts.configs.preprocessing.internal_batch_correction import (
    InternalBatchCorrectionControlSiteMode,
    InternalBatchCorrectionControlSiteSource,
    InternalBatchCorrectionImputationPolicy,
    InternalBatchCorrectionMethod,
    InternalBatchCorrectionMissingValuePolicy,
    InternalBatchCorrectionRequest,
    InternalBatchCorrectionStageOrder,
)
from phospy.errors.input import PhosPyInputError
from phospy.validation.common.config_values import require_non_empty_string
from phospy.validation.common.numbers import require_int_at_least
from phospy.validation.configs.preprocessing import (
    reject_unsupported_ruv_iii_style_method,
    validate_batch_correction_config,
)

DATASET_BATCH_CORRECTION_METHOD_NONE = "none"
DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH = "linear_residualize_batch"
DatasetBatchCorrectionMethod = Literal["none", "linear_residualize_batch"]
DATASET_BATCH_CORRECTION_METHODS = frozenset(
    {
        DATASET_BATCH_CORRECTION_METHOD_NONE,
        DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH,
    }
)
DATASET_BATCH_CORRECTION_METHOD_SPS_RUV_STYLE = "sps_ruv_style"
DATASET_BATCH_CORRECTION_METHOD_CONTROL_SITE_RUV_STYLE = "control_site_ruv_style"
DATASET_BATCH_CORRECTION_METHOD_RUV_III_STYLE = "ruv_iii_style"
SpsRuvBatchCorrectionMethod = Literal[
    "sps_ruv_style",
    "control_site_ruv_style",
    "ruv_iii_style",
]
SPS_RUV_BATCH_CORRECTION_METHODS = frozenset(
    {
        DATASET_BATCH_CORRECTION_METHOD_SPS_RUV_STYLE,
        DATASET_BATCH_CORRECTION_METHOD_CONTROL_SITE_RUV_STYLE,
        DATASET_BATCH_CORRECTION_METHOD_RUV_III_STYLE,
    }
)


@dataclass(frozen=True, slots=True)
class DatasetBatchCorrectionConfig:
    """Public batch-correction intent for dataset preprocessing.

    - `"none"`: do not request batch correction.
    - `"linear_residualize_batch"`: run fixed-effect residualisation of batch
      terms while preserving condition effects by design during dataset
      preprocessing.

    `"linear_residualize_batch"` is not ComBat, not RUV, and not limma
    `removeBatchEffect` parity. Dataset-build execution resolves sample
    metadata, validates design adequacy, applies correction, and records a
    typed preprocessing report.
    """

    method: DatasetBatchCorrectionMethod = DATASET_BATCH_CORRECTION_METHOD_NONE
    batch_column: str = "batch"
    condition_column: str = "condition"
    preserve_condition_effects: Literal[True] = True

    def __post_init__(self) -> None:
        validate_batch_correction_config(
            method=self.method,
            batch_column=self.batch_column,
            condition_column=self.condition_column,
            preserve_condition_effects=self.preserve_condition_effects,
            supported_methods=DATASET_BATCH_CORRECTION_METHODS,
        )


@dataclass(frozen=True, slots=True)
class SpsRuvBatchCorrectionConfig:
    """Explicit public config for native SPS/RUV-style preprocessing correction.

    The caller must supply controls and correction policy metadata. This
    configuration never selects controls, fetches online resources, or permits
    correction without provenance. `replicate_column`, when provided for the
    native lane, is recorded for provenance and diagnostics only; it does not
    enable replicate-aware RUV-III correction semantics. The `ruv_iii_style`
    label is retained for compatibility and roadmap clarity, but it is not
    executable until replicate-aware RUV-III semantics are implemented.
    """

    control_site_set: object
    batch_column: str
    condition_columns: tuple[str, ...]
    missingness_policy: CorrectionMissingnessPolicy
    n_unwanted_factors: int
    method: SpsRuvBatchCorrectionMethod = DATASET_BATCH_CORRECTION_METHOD_SPS_RUV_STYLE
    replicate_column: str | None = None
    stage_order: InternalBatchCorrectionStageOrder = (
        InternalBatchCorrectionStageOrder.AFTER_MISSING_DATA_BEFORE_DOWNSTREAM
    )
    diagnostics_enabled: bool = True
    provenance_enabled: Literal[True] = True

    def __post_init__(self) -> None:
        method = InternalBatchCorrectionMethod.parse(
            self.method,
            field_name="dataset build request preprocessing_config.batch_correction.method",
        )
        if method.value not in SPS_RUV_BATCH_CORRECTION_METHODS:
            raise PhosPyInputError(
                "dataset build request preprocessing_config.batch_correction."
                "method must be one of: "
                + ", ".join(sorted(SPS_RUV_BATCH_CORRECTION_METHODS))
            )
        reject_unsupported_ruv_iii_style_method(
            method,
            field_name=(
                "dataset build request preprocessing_config.batch_correction.method"
            ),
        )
        _require_control_site_set(self.control_site_set)
        batch_column = require_non_empty_string(
            self.batch_column,
            field_name=(
                "dataset build request preprocessing_config.batch_correction."
                "batch_column"
            ),
            error_type=PhosPyInputError,
        )
        condition_columns = _require_condition_columns(self.condition_columns)
        replicate_column = (
            None
            if self.replicate_column is None
            else require_non_empty_string(
                self.replicate_column,
                field_name=(
                    "dataset build request preprocessing_config.batch_correction."
                    "replicate_column"
                ),
                error_type=PhosPyInputError,
                when_provided=True,
            )
        )
        if not isinstance(self.missingness_policy, CorrectionMissingnessPolicy):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.batch_correction."
                "missingness_policy must be a CorrectionMissingnessPolicy"
            )
        n_unwanted_factors = require_int_at_least(
            self.n_unwanted_factors,
            field_name=(
                "dataset build request preprocessing_config.batch_correction."
                "n_unwanted_factors"
            ),
            minimum=1,
            error_type=PhosPyInputError,
        )
        stage_order = InternalBatchCorrectionStageOrder.parse(
            self.stage_order,
            field_name=(
                "dataset build request preprocessing_config.batch_correction."
                "stage_order"
            ),
        )
        if not isinstance(self.diagnostics_enabled, bool):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.batch_correction."
                "diagnostics_enabled must be a bool"
            )
        if self.provenance_enabled is not True:
            raise PhosPyInputError(
                "dataset build request preprocessing_config.batch_correction."
                "provenance_enabled must be True; native correction cannot run "
                "without provenance"
            )

        object.__setattr__(self, "method", method.value)
        object.__setattr__(self, "batch_column", batch_column)
        object.__setattr__(self, "condition_columns", condition_columns)
        object.__setattr__(self, "replicate_column", replicate_column)
        object.__setattr__(self, "n_unwanted_factors", n_unwanted_factors)
        object.__setattr__(self, "stage_order", stage_order)

    def to_internal_request(self) -> InternalBatchCorrectionRequest:
        """Return the protected internal workflow request contract."""

        missing_value_policy, imputation_policy = _resolve_internal_missingness(
            self.missingness_policy
        )
        return InternalBatchCorrectionRequest(
            method=InternalBatchCorrectionMethod(str(self.method)),
            batch_column=self.batch_column,
            condition_columns=self.condition_columns,
            replicate_column=self.replicate_column,
            control_site_source=InternalBatchCorrectionControlSiteSource.CALLER_SUPPLIED,
            control_site_mode=_resolve_control_site_mode(self.control_site_set),
            missing_value_policy=missing_value_policy,
            imputation_policy=imputation_policy,
            n_unwanted_factors=self.n_unwanted_factors,
            stage_order=self.stage_order,
            diagnostics_enabled=self.diagnostics_enabled,
        )


DatasetPreprocessingBatchCorrectionConfig = (
    DatasetBatchCorrectionConfig | SpsRuvBatchCorrectionConfig
)


def _require_control_site_set(value: object) -> None:
    from phospy.science.datasets.preprocessing.control_sites import ControlSiteSet

    if not isinstance(value, ControlSiteSet):
        raise PhosPyInputError(
            "dataset build request preprocessing_config.batch_correction."
            "control_site_set must be a ControlSiteSet"
        )
    if not value.annotations:
        raise PhosPyInputError(
            "dataset build request preprocessing_config.batch_correction."
            "control_site_set must contain caller-supplied control annotations"
        )


def _require_condition_columns(value: object) -> tuple[str, ...]:
    if isinstance(value, str) or not isinstance(value, tuple | list):
        raise PhosPyInputError(
            "dataset build request preprocessing_config.batch_correction."
            "condition_columns must be a non-empty sequence of non-empty strings"
        )
    resolved = tuple(
        require_non_empty_string(
            item,
            field_name=(
                "dataset build request preprocessing_config.batch_correction."
                "condition_columns[]"
            ),
            error_type=PhosPyInputError,
        )
        for item in value
    )
    if not resolved:
        raise PhosPyInputError(
            "dataset build request preprocessing_config.batch_correction."
            "condition_columns must be a non-empty sequence of non-empty strings"
        )
    if len(set(resolved)) != len(resolved):
        raise PhosPyInputError(
            "dataset build request preprocessing_config.batch_correction."
            "condition_columns must not contain duplicates"
        )
    return resolved


def _resolve_internal_missingness(
    policy: CorrectionMissingnessPolicy,
) -> tuple[
    InternalBatchCorrectionMissingValuePolicy,
    InternalBatchCorrectionImputationPolicy,
]:
    method = policy.temporary_imputation.method
    if policy.temporary_imputation.allowed:
        return (
            InternalBatchCorrectionMissingValuePolicy.ALLOW_TEMPORARY_IMPUTATION,
            _internal_imputation_method(method),
        )
    return (
        InternalBatchCorrectionMissingValuePolicy.REJECT_MISSING,
        InternalBatchCorrectionImputationPolicy.NONE,
    )


def _internal_imputation_method(
    method: TemporaryImputationMethod,
) -> InternalBatchCorrectionImputationPolicy:
    if method is TemporaryImputationMethod.ROW_MEDIAN_TEMPORARY:
        return InternalBatchCorrectionImputationPolicy.ROW_MEDIAN_TEMPORARY
    if method is TemporaryImputationMethod.MINPROB_TEMPORARY:
        return InternalBatchCorrectionImputationPolicy.MINPROB_TEMPORARY
    if method is TemporaryImputationMethod.KNN_TEMPORARY:
        return InternalBatchCorrectionImputationPolicy.KNN_TEMPORARY
    raise PhosPyInputError(
        "dataset build request preprocessing_config.batch_correction."
        "missingness_policy temporary imputation method is unsupported for "
        "native SPS/RUV-style correction"
    )


def _resolve_control_site_mode(
    control_site_set: object,
) -> InternalBatchCorrectionControlSiteMode:
    annotations = getattr(control_site_set, "annotations", ())
    if any(
        getattr(annotation, "weight", None) is not None for annotation in annotations
    ):
        return InternalBatchCorrectionControlSiteMode.WEIGHT_COLUMN
    return InternalBatchCorrectionControlSiteMode.SITE_KEY_LIST


__all__ = [
    "DATASET_BATCH_CORRECTION_METHOD_CONTROL_SITE_RUV_STYLE",
    "DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH",
    "DATASET_BATCH_CORRECTION_METHOD_NONE",
    "DATASET_BATCH_CORRECTION_METHOD_RUV_III_STYLE",
    "DATASET_BATCH_CORRECTION_METHOD_SPS_RUV_STYLE",
    "DATASET_BATCH_CORRECTION_METHODS",
    "SPS_RUV_BATCH_CORRECTION_METHODS",
    "DatasetBatchCorrectionConfig",
    "DatasetBatchCorrectionMethod",
    "DatasetPreprocessingBatchCorrectionConfig",
    "SpsRuvBatchCorrectionConfig",
    "SpsRuvBatchCorrectionMethod",
]
