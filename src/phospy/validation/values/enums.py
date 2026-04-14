from __future__ import annotations

from ...errors import PhospyValidationError
from ...internal.types import (
    DuplicateSiteStrategy,
    KinaseProfileMissingValueStrategy,
    PredictionSvmMode,
    PredictionTraceFormat,
    PredictionTraceLevel,
    SignalomeKinaseNetworkPolicy,
    SignalomeModuleSelectionStrategy,
)


def validate_duplicate_site_strategy(
    value: DuplicateSiteStrategy,
) -> DuplicateSiteStrategy:
    """Validate the configured duplicate-site handling strategy."""

    if value not in {
        "max_mean_signal",
        "first",
        "aggregate_mean",
        "aggregate_median",
        "error",
    }:
        msg = (
            "duplicate_site_strategy must be one of: 'max_mean_signal', 'first', "
            "'aggregate_mean', 'aggregate_median', 'error'"
        )
        raise PhospyValidationError(msg)
    return value


def validate_missing_value_strategy(
    value: KinaseProfileMissingValueStrategy,
) -> KinaseProfileMissingValueStrategy:
    """Validate the configured kinase-profile missing-value strategy."""

    if value not in {"propagate_any_missing", "median_skipna"}:
        msg = (
            "missing_value_strategy must be one of: 'propagate_any_missing', "
            "'median_skipna'"
        )
        raise PhospyValidationError(msg)
    return value


def validate_module_selection_strategy(
    value: SignalomeModuleSelectionStrategy,
) -> SignalomeModuleSelectionStrategy:
    """Validate the configured signalome module-selection strategy."""

    if value not in {"correlation_thresholds", "single_module"}:
        msg = (
            "module_selection_strategy must be one of: 'correlation_thresholds', "
            "'single_module'"
        )
        raise PhospyValidationError(msg)
    return value


def validate_kinase_network_policy(
    value: SignalomeKinaseNetworkPolicy,
) -> SignalomeKinaseNetworkPolicy:
    """Validate the configured signalome kinase-network edge policy."""

    if value not in {"positive_only", "absolute_threshold", "signed"}:
        msg = (
            "kinase_network_policy must be one of: 'positive_only', "
            "'absolute_threshold', 'signed'"
        )
        raise PhospyValidationError(msg)
    return value


def validate_svm_mode(value: PredictionSvmMode) -> PredictionSvmMode:
    """Validate the configured prediction SVM mode."""

    if value not in {"default", "r_parity"}:
        msg = "svm_mode must be one of: 'default', 'r_parity'"
        raise PhospyValidationError(msg)
    return value


def validate_trace_level(value: PredictionTraceLevel) -> PredictionTraceLevel:
    """Validate the configured prediction trace level."""

    if value not in {"none", "summary", "full"}:
        msg = "trace_level must be one of: 'none', 'summary', 'full'"
        raise PhospyValidationError(msg)
    return value


def validate_trace_format(value: PredictionTraceFormat) -> PredictionTraceFormat:
    """Validate the configured prediction trace sink format."""

    if value not in {"csv", "parquet"}:
        msg = "trace_sink_format must be one of: 'csv', 'parquet'"
        raise PhospyValidationError(msg)
    return value


__all__ = [
    "validate_duplicate_site_strategy",
    "validate_kinase_network_policy",
    "validate_missing_value_strategy",
    "validate_module_selection_strategy",
    "validate_svm_mode",
    "validate_trace_format",
    "validate_trace_level",
]
