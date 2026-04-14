from __future__ import annotations

from ...errors import PhospyValidationError
from ...internal.types import (
    DUPLICATE_SITE_STRATEGIES,
    KINASE_PROFILE_MISSING_VALUE_STRATEGIES,
    PREDICTION_SVM_MODES,
    PREDICTION_TRACE_FORMATS,
    PREDICTION_TRACE_LEVELS,
    SIGNALOME_ASSIGNMENT_POLICIES,
    SIGNALOME_KINASE_NETWORK_POLICIES,
    SIGNALOME_MODULE_SELECTION_STRATEGIES,
    SITE_MATRIX_MISSING_DATA_POLICIES,
    DuplicateSiteStrategy,
    KinaseProfileMissingValueStrategy,
    PredictionSvmMode,
    PredictionTraceFormat,
    PredictionTraceLevel,
    SignalomeAssignmentPolicy,
    SignalomeKinaseNetworkPolicy,
    SignalomeModuleSelectionStrategy,
    SiteMatrixMissingDataPolicy,
)


def validate_duplicate_site_strategy(
    value: DuplicateSiteStrategy,
) -> DuplicateSiteStrategy:
    """Validate the configured duplicate-site handling strategy."""

    if value not in DUPLICATE_SITE_STRATEGIES:
        msg = "duplicate_site_strategy must be one of: " + ", ".join(
            f"'{token}'" for token in DUPLICATE_SITE_STRATEGIES
        )
        raise PhospyValidationError(msg)
    return value


def validate_site_matrix_missing_data_policy(
    value: SiteMatrixMissingDataPolicy,
) -> SiteMatrixMissingDataPolicy:
    """Validate the configured site-matrix missing-data policy."""

    if value not in SITE_MATRIX_MISSING_DATA_POLICIES:
        msg = "missing_data_policy must be one of: " + ", ".join(
            f"'{token}'" for token in SITE_MATRIX_MISSING_DATA_POLICIES
        )
        raise PhospyValidationError(msg)
    return value


def validate_missing_value_strategy(
    value: KinaseProfileMissingValueStrategy,
) -> KinaseProfileMissingValueStrategy:
    """Validate the configured kinase-profile missing-value strategy."""

    if value not in KINASE_PROFILE_MISSING_VALUE_STRATEGIES:
        msg = "missing_value_strategy must be one of: " + ", ".join(
            f"'{token}'" for token in KINASE_PROFILE_MISSING_VALUE_STRATEGIES
        )
        raise PhospyValidationError(msg)
    return value


def validate_module_selection_strategy(
    value: SignalomeModuleSelectionStrategy,
) -> SignalomeModuleSelectionStrategy:
    """Validate the configured signalome module-selection strategy."""

    if value not in SIGNALOME_MODULE_SELECTION_STRATEGIES:
        msg = "module_selection_strategy must be one of: " + ", ".join(
            f"'{token}'" for token in SIGNALOME_MODULE_SELECTION_STRATEGIES
        )
        raise PhospyValidationError(msg)
    return value


def validate_kinase_network_policy(
    value: SignalomeKinaseNetworkPolicy,
) -> SignalomeKinaseNetworkPolicy:
    """Validate the configured signalome kinase-network edge policy."""

    if value not in SIGNALOME_KINASE_NETWORK_POLICIES:
        msg = "kinase_network_policy must be one of: " + ", ".join(
            f"'{token}'" for token in SIGNALOME_KINASE_NETWORK_POLICIES
        )
        raise PhospyValidationError(msg)
    return value


def validate_signalome_assignment_policy(
    value: SignalomeAssignmentPolicy,
) -> SignalomeAssignmentPolicy:
    """Validate signalome assignment propagation policy."""

    if value not in SIGNALOME_ASSIGNMENT_POLICIES:
        msg = "assignment_policy must be one of: " + ", ".join(
            f"'{token}'" for token in SIGNALOME_ASSIGNMENT_POLICIES
        )
        raise PhospyValidationError(msg)
    return value


def validate_svm_mode(value: PredictionSvmMode) -> PredictionSvmMode:
    """Validate the configured prediction SVM mode."""

    if value not in PREDICTION_SVM_MODES:
        msg = "svm_mode must be one of: " + ", ".join(
            f"'{token}'" for token in PREDICTION_SVM_MODES
        )
        raise PhospyValidationError(msg)
    return value


def validate_trace_level(value: PredictionTraceLevel) -> PredictionTraceLevel:
    """Validate the configured prediction trace level."""

    if value not in PREDICTION_TRACE_LEVELS:
        msg = "trace_level must be one of: " + ", ".join(
            f"'{token}'" for token in PREDICTION_TRACE_LEVELS
        )
        raise PhospyValidationError(msg)
    return value


def validate_trace_format(value: PredictionTraceFormat) -> PredictionTraceFormat:
    """Validate the configured prediction trace sink format."""

    if value not in PREDICTION_TRACE_FORMATS:
        msg = "trace_sink_format must be one of: " + ", ".join(
            f"'{token}'" for token in PREDICTION_TRACE_FORMATS
        )
        raise PhospyValidationError(msg)
    return value


__all__ = [
    "validate_duplicate_site_strategy",
    "validate_kinase_network_policy",
    "validate_missing_value_strategy",
    "validate_module_selection_strategy",
    "validate_site_matrix_missing_data_policy",
    "validate_signalome_assignment_policy",
    "validate_svm_mode",
    "validate_trace_format",
    "validate_trace_level",
]
