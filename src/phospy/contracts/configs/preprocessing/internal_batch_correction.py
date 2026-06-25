"""Internal future SPS/RUV-style batch-correction request contracts.

These contracts describe preprocessing intent only. They do not select
controls, estimate unwanted factors, modify matrices, fetch resources, or claim
equivalence to another implementation.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

from phospy.policies import PolicyEnum
from phospy.validation.configs.preprocessing import (
    validate_internal_batch_correction_request,
)


class InternalBatchCorrectionMethod(PolicyEnum):
    """Known correction intent labels; not every label is executable."""

    SPS_RUV_STYLE = "sps_ruv_style"
    CONTROL_SITE_RUV_STYLE = "control_site_ruv_style"
    RUV_III_STYLE = "ruv_iii_style"


class InternalBatchCorrectionControlSiteSource(PolicyEnum):
    """Where a future correction request expects control-site labels to come from."""

    CALLER_SUPPLIED = "caller_supplied"
    DATASET_METADATA = "dataset_metadata"
    PRECOMPUTED_SELECTION = "precomputed_selection"
    PACKAGED_REFERENCE = "packaged_reference"


class InternalBatchCorrectionControlSiteMode(PolicyEnum):
    """Shape of the requested control-site representation."""

    SITE_KEY_LIST = "site_key_list"
    BOOLEAN_COLUMN = "boolean_column"
    WEIGHT_COLUMN = "weight_column"
    SELECTION_RESULT = "selection_result"


class InternalBatchCorrectionMissingValuePolicy(PolicyEnum):
    """Requested missing-value handling boundary for a future correction."""

    REJECT_MISSING = "reject_missing"
    ALLOW_TEMPORARY_IMPUTATION = "allow_temporary_imputation"
    USE_EXISTING_IMPUTATION_PROVENANCE = "use_existing_imputation_provenance"


class InternalBatchCorrectionImputationPolicy(PolicyEnum):
    """Requested temporary imputation policy label, without implementation."""

    NONE = "none"
    ROW_MEDIAN_TEMPORARY = "row_median_temporary"
    MINPROB_TEMPORARY = "minprob_temporary"
    KNN_TEMPORARY = "knn_temporary"


class InternalBatchCorrectionStageOrder(PolicyEnum):
    """Requested preprocessing stage placement label."""

    AFTER_INTENSITY_TRANSFORM_BEFORE_MISSING_DATA = (
        "after_intensity_transform_before_missing_data"
    )
    AFTER_MISSING_DATA_BEFORE_DOWNSTREAM = "after_missing_data_before_downstream"
    AFTER_TOTAL_PROTEIN_CORRECTION_BEFORE_DOWNSTREAM = (
        "after_total_protein_correction_before_downstream"
    )


INTERNAL_BATCH_CORRECTION_METHODS = frozenset(InternalBatchCorrectionMethod)
INTERNAL_BATCH_CORRECTION_CONTROL_SITE_SOURCES = frozenset(
    InternalBatchCorrectionControlSiteSource
)
INTERNAL_BATCH_CORRECTION_CONTROL_SITE_MODES = frozenset(
    InternalBatchCorrectionControlSiteMode
)
INTERNAL_BATCH_CORRECTION_MISSING_VALUE_POLICIES = frozenset(
    InternalBatchCorrectionMissingValuePolicy
)
INTERNAL_BATCH_CORRECTION_IMPUTATION_POLICIES = frozenset(
    InternalBatchCorrectionImputationPolicy
)
INTERNAL_BATCH_CORRECTION_STAGE_ORDERS = frozenset(InternalBatchCorrectionStageOrder)
SUPPORTED_INTERNAL_BATCH_CORRECTION_STAGE_ORDER = (
    InternalBatchCorrectionStageOrder.AFTER_MISSING_DATA_BEFORE_DOWNSTREAM
)
SUPPORTED_INTERNAL_BATCH_CORRECTION_EXECUTED_STAGE_ORDER = (
    "missing_data",
    "batch_correction",
    "downstream_workflows",
)


@dataclass(frozen=True, slots=True)
class InternalBatchCorrectionRequest:
    """Internal request contract for future native SPS/RUV-style correction.

    Construction validates local scalar values and stores typed enum members.
    It intentionally does not validate scientific eligibility, select or
    resolve controls, estimate unwanted factors, impute values, or modify any
    phosphosite matrix.
    """

    method: InternalBatchCorrectionMethod
    batch_column: str
    condition_columns: tuple[str, ...]
    replicate_column: str | None
    control_site_source: InternalBatchCorrectionControlSiteSource
    control_site_mode: InternalBatchCorrectionControlSiteMode
    missing_value_policy: InternalBatchCorrectionMissingValuePolicy
    imputation_policy: InternalBatchCorrectionImputationPolicy
    n_unwanted_factors: int | None
    stage_order: InternalBatchCorrectionStageOrder
    diagnostics_enabled: bool

    def __post_init__(self) -> None:
        missing_policy = InternalBatchCorrectionMissingValuePolicy
        existing_imputation = missing_policy.USE_EXISTING_IMPUTATION_PROVENANCE.value
        no_imputation = InternalBatchCorrectionImputationPolicy.NONE.value
        validated = validate_internal_batch_correction_request(
            method=self.method,
            batch_column=self.batch_column,
            condition_columns=self.condition_columns,
            replicate_column=self.replicate_column,
            control_site_source=self.control_site_source,
            control_site_mode=self.control_site_mode,
            missing_value_policy=self.missing_value_policy,
            imputation_policy=self.imputation_policy,
            n_unwanted_factors=self.n_unwanted_factors,
            stage_order=self.stage_order,
            diagnostics_enabled=self.diagnostics_enabled,
            method_type=InternalBatchCorrectionMethod,
            control_site_source_type=InternalBatchCorrectionControlSiteSource,
            control_site_mode_type=InternalBatchCorrectionControlSiteMode,
            missing_value_policy_type=InternalBatchCorrectionMissingValuePolicy,
            imputation_policy_type=InternalBatchCorrectionImputationPolicy,
            stage_order_type=InternalBatchCorrectionStageOrder,
            missing_value_policy_reject=missing_policy.REJECT_MISSING.value,
            missing_value_policy_allow_temporary_imputation=(
                missing_policy.ALLOW_TEMPORARY_IMPUTATION.value
            ),
            missing_value_policy_use_existing_imputation_provenance=(
                existing_imputation
            ),
            imputation_policy_none=no_imputation,
        )
        if validated["method"] is InternalBatchCorrectionMethod.CONTROL_SITE_RUV_STYLE:
            warnings.warn(
                "control_site_ruv_style is a deprecated internal compatibility "
                "alias; use sps_ruv_style for native PhosPy SPS/RUV-style "
                "correction.",
                DeprecationWarning,
                stacklevel=2,
            )
            validated["method"] = InternalBatchCorrectionMethod.SPS_RUV_STYLE
        for field_name, value in validated.items():
            object.__setattr__(self, field_name, value)


__all__ = [
    "INTERNAL_BATCH_CORRECTION_CONTROL_SITE_MODES",
    "INTERNAL_BATCH_CORRECTION_CONTROL_SITE_SOURCES",
    "INTERNAL_BATCH_CORRECTION_IMPUTATION_POLICIES",
    "INTERNAL_BATCH_CORRECTION_METHODS",
    "INTERNAL_BATCH_CORRECTION_MISSING_VALUE_POLICIES",
    "INTERNAL_BATCH_CORRECTION_STAGE_ORDERS",
    "SUPPORTED_INTERNAL_BATCH_CORRECTION_EXECUTED_STAGE_ORDER",
    "SUPPORTED_INTERNAL_BATCH_CORRECTION_STAGE_ORDER",
    "InternalBatchCorrectionControlSiteMode",
    "InternalBatchCorrectionControlSiteSource",
    "InternalBatchCorrectionImputationPolicy",
    "InternalBatchCorrectionMethod",
    "InternalBatchCorrectionMissingValuePolicy",
    "InternalBatchCorrectionRequest",
    "InternalBatchCorrectionStageOrder",
]
