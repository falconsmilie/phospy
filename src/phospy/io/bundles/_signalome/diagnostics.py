"""Signalome diagnostics payload conversion helpers."""

from __future__ import annotations

from phospy.contracts.configs import SIGNALOME_SCORE_PRECONDITIONING_POLICIES
from phospy.errors.input import PhosPyInputError
from phospy.io.bundles._shared.primitives import (
    require_bool,
    require_float,
    require_mapping,
    require_str,
)
from phospy.io.bundles._signalome.primitives import (
    _parse_optional_int,
    _reject_unsupported_fields,
    _require_fields,
    _require_int,
)
from phospy.provenance.serialization import (
    table_fingerprint_from_payload,
    table_fingerprint_to_payload,
)
from phospy.science.signalomes.models import (
    SIGNALOME_MODULE_SELECTION_STABILITY_STATUS_NOT_COMPUTABLE,
    SIGNALOME_MODULE_SELECTION_STABILITY_STATUS_STABLE,
    SIGNALOME_MODULE_SELECTION_STABILITY_STATUS_UNSTABLE,
    SIGNALOME_MODULE_SELECTION_STRATEGY_CORRELATION_THRESHOLDS,
    SIGNALOME_MODULE_SELECTION_STRATEGY_EXPLICIT_MODULE_COUNT,
    SignalomeAlignmentDiagnostics,
    SignalomeAlignmentInputDiagnostics,
    SignalomeClusterCandidateScore,
    SignalomeClusteringPreparationDiagnostics,
    SignalomeModuleSelectionAssignmentSimilaritySummary,
    SignalomeModuleSelectionDiagnostics,
    SignalomeModuleSelectionStabilityReport,
    SignalomeModuleSelectionThresholdSensitivity,
    SignalomeModuleSelectionThresholdSensitivityRecord,
    SignalomeNetworkCorrelationDiagnostics,
    SignalomeScorePreconditioningDiagnostics,
    default_signalome_clustering_preparation_diagnostics,
    default_signalome_module_selection_stability_report,
)


def signalome_module_selection_diagnostics_to_payload(
    diagnostics: SignalomeModuleSelectionDiagnostics,
) -> dict[str, object]:
    return {
        "strategy": str(diagnostics.strategy),
        "selected_module_count": int(diagnostics.selected_module_count),
        "requested_module_count": (
            None
            if diagnostics.requested_module_count is None
            else int(diagnostics.requested_module_count)
        ),
        "threshold_used": (
            None
            if diagnostics.threshold_used is None
            else float(diagnostics.threshold_used)
        ),
        "max_clusters_evaluated": int(diagnostics.max_clusters_evaluated),
        "candidate_scores": {
            str(cluster_count): {
                "min_median_correlation": float(score.min_median_correlation),
                "mean_median_correlation": float(score.mean_median_correlation),
            }
            for cluster_count, score in diagnostics.candidate_scores.items()
        },
        "reason": str(diagnostics.reason),
        "zero_variance_profile_count": int(diagnostics.zero_variance_profile_count),
        "near_constant_profile_count": int(diagnostics.near_constant_profile_count),
        "excluded_from_correlation_count": int(
            diagnostics.excluded_from_correlation_count
        ),
        "stability_report": signalome_module_selection_stability_report_to_payload(
            diagnostics.stability_report
        ),
    }


def signalome_module_selection_diagnostics_from_payload(
    payload: object,
    *,
    scope: str,
) -> SignalomeModuleSelectionDiagnostics:
    diagnostics_payload = require_mapping(
        payload,
        field_name=f"{scope}.module_selection_diagnostics",
    )
    diagnostics_field_name = f"{scope}.module_selection_diagnostics"
    allowed_fields = frozenset(
        {
            "strategy",
            "selected_module_count",
            "requested_module_count",
            "threshold_used",
            "max_clusters_evaluated",
            "candidate_scores",
            "reason",
            "zero_variance_profile_count",
            "near_constant_profile_count",
            "excluded_from_correlation_count",
            "stability_report",
        }
    )
    _reject_unsupported_fields(
        diagnostics_payload,
        field_name=diagnostics_field_name,
        allowed_fields=allowed_fields,
    )
    required_fields = allowed_fields - {"stability_report"}
    _require_fields(
        diagnostics_payload,
        field_name=diagnostics_field_name,
        required_fields=required_fields,
    )
    strategy = require_str(
        diagnostics_payload.get("strategy"),
        field_name=f"{scope}.module_selection_diagnostics.strategy",
    )
    if strategy not in {
        SIGNALOME_MODULE_SELECTION_STRATEGY_CORRELATION_THRESHOLDS,
        SIGNALOME_MODULE_SELECTION_STRATEGY_EXPLICIT_MODULE_COUNT,
    }:
        allowed = ", ".join(
            (
                SIGNALOME_MODULE_SELECTION_STRATEGY_CORRELATION_THRESHOLDS,
                SIGNALOME_MODULE_SELECTION_STRATEGY_EXPLICIT_MODULE_COUNT,
            )
        )
        raise PhosPyInputError(
            f"{scope}.module_selection_diagnostics.strategy must be one of: {allowed}"
        )
    candidate_scores_payload = require_mapping(
        diagnostics_payload.get("candidate_scores"),
        field_name=f"{scope}.module_selection_diagnostics.candidate_scores",
    )
    candidate_scores: dict[int, SignalomeClusterCandidateScore] = {}
    for cluster_count_raw, score_payload in candidate_scores_payload.items():
        score_mapping = require_mapping(
            score_payload,
            field_name=(
                f"{scope}.module_selection_diagnostics.candidate_scores."
                f"{cluster_count_raw}"
            ),
        )
        candidate_scores[int(cluster_count_raw)] = SignalomeClusterCandidateScore(
            min_median_correlation=require_float(
                score_mapping.get("min_median_correlation"),
                field_name=(
                    f"{scope}.module_selection_diagnostics.candidate_scores."
                    f"{cluster_count_raw}.min_median_correlation"
                ),
            ),
            mean_median_correlation=require_float(
                score_mapping.get("mean_median_correlation"),
                field_name=(
                    f"{scope}.module_selection_diagnostics.candidate_scores."
                    f"{cluster_count_raw}.mean_median_correlation"
                ),
            ),
        )
    requested_module_count = _parse_optional_int(
        diagnostics_payload.get("requested_module_count"),
        field_name=f"{scope}.module_selection_diagnostics.requested_module_count",
    )
    threshold_used_raw = diagnostics_payload.get("threshold_used")
    threshold_used = (
        None
        if threshold_used_raw is None
        else require_float(
            threshold_used_raw,
            field_name=f"{scope}.module_selection_diagnostics.threshold_used",
        )
    )
    return SignalomeModuleSelectionDiagnostics(
        strategy=strategy,  # type: ignore[arg-type]
        selected_module_count=_require_int(
            diagnostics_payload.get("selected_module_count"),
            field_name=f"{scope}.module_selection_diagnostics.selected_module_count",
        ),
        requested_module_count=requested_module_count,
        threshold_used=threshold_used,
        max_clusters_evaluated=_require_int(
            diagnostics_payload.get("max_clusters_evaluated"),
            field_name=f"{scope}.module_selection_diagnostics.max_clusters_evaluated",
        ),
        candidate_scores=candidate_scores,
        reason=require_str(
            diagnostics_payload.get("reason"),
            field_name=f"{scope}.module_selection_diagnostics.reason",
        ),
        zero_variance_profile_count=_require_int(
            diagnostics_payload.get("zero_variance_profile_count"),
            field_name=(
                f"{scope}.module_selection_diagnostics.zero_variance_profile_count"
            ),
        ),
        near_constant_profile_count=_require_int(
            diagnostics_payload.get("near_constant_profile_count"),
            field_name=(
                f"{scope}.module_selection_diagnostics.near_constant_profile_count"
            ),
        ),
        excluded_from_correlation_count=_require_int(
            diagnostics_payload.get("excluded_from_correlation_count"),
            field_name=(
                f"{scope}.module_selection_diagnostics.excluded_from_correlation_count"
            ),
        ),
        stability_report=signalome_module_selection_stability_report_from_payload(
            diagnostics_payload.get("stability_report"),
            scope=scope,
        ),
    )


def signalome_module_selection_stability_report_to_payload(
    report: SignalomeModuleSelectionStabilityReport,
) -> dict[str, object]:
    return {
        "evaluation_method": str(report.evaluation_method),
        "evaluation_version": str(report.evaluation_version),
        "seed_policy": str(report.seed_policy),
        "random_seed": None if report.random_seed is None else int(report.random_seed),
        "perturbation_count": int(report.perturbation_count),
        "selected_count_frequency": {
            str(count): int(frequency)
            for count, frequency in report.selected_count_frequency.items()
        },
        "assignment_similarity_metric": str(report.assignment_similarity_metric),
        "assignment_similarity": _assignment_similarity_to_payload(
            report.assignment_similarity
        ),
        "threshold_sensitivity": _threshold_sensitivity_to_payload(
            report.threshold_sensitivity
        ),
        "status": str(report.status),
        "limitations": list(report.limitations),
        "not_computable_reason": report.not_computable_reason,
        "base_selected_module_count": int(report.base_selected_module_count),
        "input_site_count": int(report.input_site_count),
        "input_dimension_count": int(report.input_dimension_count),
    }


def signalome_module_selection_stability_report_from_payload(
    payload: object,
    *,
    scope: str,
) -> SignalomeModuleSelectionStabilityReport:
    if payload is None:
        return default_signalome_module_selection_stability_report()
    field_name = f"{scope}.module_selection_diagnostics.stability_report"
    report_payload = require_mapping(payload, field_name=field_name)
    allowed_fields = frozenset(
        {
            "evaluation_method",
            "evaluation_version",
            "seed_policy",
            "random_seed",
            "perturbation_count",
            "selected_count_frequency",
            "assignment_similarity_metric",
            "assignment_similarity",
            "threshold_sensitivity",
            "status",
            "limitations",
            "not_computable_reason",
            "base_selected_module_count",
            "input_site_count",
            "input_dimension_count",
        }
    )
    _reject_unsupported_fields(
        report_payload,
        field_name=field_name,
        allowed_fields=allowed_fields,
    )
    _require_fields(
        report_payload,
        field_name=field_name,
        required_fields=allowed_fields,
    )
    status = require_str(
        report_payload.get("status"), field_name=f"{field_name}.status"
    )
    if status not in {
        SIGNALOME_MODULE_SELECTION_STABILITY_STATUS_STABLE,
        SIGNALOME_MODULE_SELECTION_STABILITY_STATUS_UNSTABLE,
        SIGNALOME_MODULE_SELECTION_STABILITY_STATUS_NOT_COMPUTABLE,
    }:
        raise PhosPyInputError(
            f"{field_name}.status must be one of: stable, unstable, not_computable"
        )
    return SignalomeModuleSelectionStabilityReport(
        evaluation_method=require_str(
            report_payload.get("evaluation_method"),
            field_name=f"{field_name}.evaluation_method",
        ),
        evaluation_version=require_str(
            report_payload.get("evaluation_version"),
            field_name=f"{field_name}.evaluation_version",
        ),
        seed_policy=require_str(
            report_payload.get("seed_policy"),
            field_name=f"{field_name}.seed_policy",
        ),
        random_seed=_parse_optional_int(
            report_payload.get("random_seed"),
            field_name=f"{field_name}.random_seed",
        ),
        perturbation_count=_require_int(
            report_payload.get("perturbation_count"),
            field_name=f"{field_name}.perturbation_count",
        ),
        selected_count_frequency=_frequency_from_payload(
            report_payload.get("selected_count_frequency"),
            field_name=f"{field_name}.selected_count_frequency",
        ),
        assignment_similarity_metric=require_str(
            report_payload.get("assignment_similarity_metric"),
            field_name=f"{field_name}.assignment_similarity_metric",
        ),
        assignment_similarity=_assignment_similarity_from_payload(
            report_payload.get("assignment_similarity"),
            field_name=f"{field_name}.assignment_similarity",
        ),
        threshold_sensitivity=_threshold_sensitivity_from_payload(
            report_payload.get("threshold_sensitivity"),
            field_name=f"{field_name}.threshold_sensitivity",
        ),
        status=status,  # type: ignore[arg-type]
        limitations=_string_tuple_from_payload(
            report_payload.get("limitations"),
            field_name=f"{field_name}.limitations",
        ),
        not_computable_reason=_optional_str(
            report_payload.get("not_computable_reason"),
            field_name=f"{field_name}.not_computable_reason",
        ),
        base_selected_module_count=_require_int(
            report_payload.get("base_selected_module_count"),
            field_name=f"{field_name}.base_selected_module_count",
        ),
        input_site_count=_require_int(
            report_payload.get("input_site_count"),
            field_name=f"{field_name}.input_site_count",
        ),
        input_dimension_count=_require_int(
            report_payload.get("input_dimension_count"),
            field_name=f"{field_name}.input_dimension_count",
        ),
    )


def _assignment_similarity_to_payload(
    summary: SignalomeModuleSelectionAssignmentSimilaritySummary,
) -> dict[str, object]:
    return {
        "metric": str(summary.metric),
        "evaluated_perturbations": int(summary.evaluated_perturbations),
        "minimum": None if summary.minimum is None else float(summary.minimum),
        "median": None if summary.median is None else float(summary.median),
        "mean": None if summary.mean is None else float(summary.mean),
        "maximum": None if summary.maximum is None else float(summary.maximum),
    }


def _assignment_similarity_from_payload(
    payload: object,
    *,
    field_name: str,
) -> SignalomeModuleSelectionAssignmentSimilaritySummary:
    mapping = require_mapping(payload, field_name=field_name)
    allowed_fields = frozenset(
        {"metric", "evaluated_perturbations", "minimum", "median", "mean", "maximum"}
    )
    _reject_unsupported_fields(
        mapping, field_name=field_name, allowed_fields=allowed_fields
    )
    _require_fields(mapping, field_name=field_name, required_fields=allowed_fields)
    return SignalomeModuleSelectionAssignmentSimilaritySummary(
        metric=require_str(mapping.get("metric"), field_name=f"{field_name}.metric"),
        evaluated_perturbations=_require_int(
            mapping.get("evaluated_perturbations"),
            field_name=f"{field_name}.evaluated_perturbations",
        ),
        minimum=_optional_float(
            mapping.get("minimum"), field_name=f"{field_name}.minimum"
        ),
        median=_optional_float(
            mapping.get("median"), field_name=f"{field_name}.median"
        ),
        mean=_optional_float(mapping.get("mean"), field_name=f"{field_name}.mean"),
        maximum=_optional_float(
            mapping.get("maximum"), field_name=f"{field_name}.maximum"
        ),
    )


def _threshold_sensitivity_to_payload(
    sensitivity: SignalomeModuleSelectionThresholdSensitivity,
) -> dict[str, object]:
    return {
        "method": str(sensitivity.method),
        "version": str(sensitivity.version),
        "records": [
            {
                "primary_threshold": float(record.primary_threshold),
                "fallback_threshold": float(record.fallback_threshold),
                "selected_module_count": int(record.selected_module_count),
                "threshold_used": (
                    None
                    if record.threshold_used is None
                    else float(record.threshold_used)
                ),
            }
            for record in sensitivity.records
        ],
        "selected_count_frequency": {
            str(count): int(frequency)
            for count, frequency in sensitivity.selected_count_frequency.items()
        },
        "disagrees_with_selected_count": bool(
            sensitivity.disagrees_with_selected_count
        ),
    }


def _threshold_sensitivity_from_payload(
    payload: object,
    *,
    field_name: str,
) -> SignalomeModuleSelectionThresholdSensitivity:
    mapping = require_mapping(payload, field_name=field_name)
    allowed_fields = frozenset(
        {
            "method",
            "version",
            "records",
            "selected_count_frequency",
            "disagrees_with_selected_count",
        }
    )
    _reject_unsupported_fields(
        mapping, field_name=field_name, allowed_fields=allowed_fields
    )
    _require_fields(mapping, field_name=field_name, required_fields=allowed_fields)
    records_raw = mapping.get("records")
    if not isinstance(records_raw, list | tuple):
        raise PhosPyInputError(f"{field_name}.records must be an array")
    return SignalomeModuleSelectionThresholdSensitivity(
        method=require_str(mapping.get("method"), field_name=f"{field_name}.method"),
        version=require_str(
            mapping.get("version"),
            field_name=f"{field_name}.version",
        ),
        records=tuple(
            _threshold_sensitivity_record_from_payload(
                item,
                field_name=f"{field_name}.records[{position}]",
            )
            for position, item in enumerate(records_raw)
        ),
        selected_count_frequency=_frequency_from_payload(
            mapping.get("selected_count_frequency"),
            field_name=f"{field_name}.selected_count_frequency",
        ),
        disagrees_with_selected_count=require_bool(
            mapping.get("disagrees_with_selected_count"),
            field_name=f"{field_name}.disagrees_with_selected_count",
        ),
    )


def _threshold_sensitivity_record_from_payload(
    payload: object,
    *,
    field_name: str,
) -> SignalomeModuleSelectionThresholdSensitivityRecord:
    mapping = require_mapping(payload, field_name=field_name)
    allowed_fields = frozenset(
        {
            "primary_threshold",
            "fallback_threshold",
            "selected_module_count",
            "threshold_used",
        }
    )
    _reject_unsupported_fields(
        mapping, field_name=field_name, allowed_fields=allowed_fields
    )
    _require_fields(mapping, field_name=field_name, required_fields=allowed_fields)
    return SignalomeModuleSelectionThresholdSensitivityRecord(
        primary_threshold=require_float(
            mapping.get("primary_threshold"),
            field_name=f"{field_name}.primary_threshold",
        ),
        fallback_threshold=require_float(
            mapping.get("fallback_threshold"),
            field_name=f"{field_name}.fallback_threshold",
        ),
        selected_module_count=_require_int(
            mapping.get("selected_module_count"),
            field_name=f"{field_name}.selected_module_count",
        ),
        threshold_used=_optional_float(
            mapping.get("threshold_used"),
            field_name=f"{field_name}.threshold_used",
        ),
    )


def _frequency_from_payload(payload: object, *, field_name: str) -> dict[int, int]:
    mapping = require_mapping(payload, field_name=field_name)
    return {
        int(key): _require_int(value, field_name=f"{field_name}.{key}")
        for key, value in mapping.items()
    }


def _optional_float(value: object, *, field_name: str) -> float | None:
    if value is None:
        return None
    return require_float(value, field_name=field_name)


def _optional_str(value: object, *, field_name: str) -> str | None:
    if value is None:
        return None
    return require_str(value, field_name=field_name)


def _string_tuple_from_payload(value: object, *, field_name: str) -> tuple[str, ...]:
    if not isinstance(value, list | tuple):
        raise PhosPyInputError(f"{field_name} must be an array")
    return tuple(
        require_str(item, field_name=f"{field_name}[{position}]")
        for position, item in enumerate(value)
    )


def signalome_clustering_preparation_diagnostics_to_payload(
    diagnostics: SignalomeClusteringPreparationDiagnostics,
) -> dict[str, object]:
    return {
        "preparation_policy_id": str(diagnostics.preparation_policy_id),
        "input_dimension_count": int(diagnostics.input_dimension_count),
        "retained_dimension_count": int(diagnostics.retained_dimension_count),
        "retained_dimension_labels": list(diagnostics.retained_dimension_labels),
        "dropped_fully_missing_dimension_count": int(
            diagnostics.dropped_fully_missing_dimension_count
        ),
        "dropped_fully_missing_dimension_labels": list(
            diagnostics.dropped_fully_missing_dimension_labels
        ),
        "dropped_fully_missing_dimension_preview": list(
            diagnostics.dropped_fully_missing_dimension_preview
        ),
        "dropped_fully_missing_value_count": int(
            diagnostics.dropped_fully_missing_value_count
        ),
        "non_finite_input_value_count": int(diagnostics.non_finite_input_value_count),
        "missing_after_non_finite_normalization_count": int(
            diagnostics.missing_after_non_finite_normalization_count
        ),
        "imputed_value_count": int(diagnostics.imputed_value_count),
        "imputed_value_counts_by_dimension": {
            str(key): int(value)
            for key, value in diagnostics.imputed_value_counts_by_dimension.items()
        },
        "prepared_matrix_fingerprint": (
            None
            if diagnostics.prepared_matrix_fingerprint is None
            else table_fingerprint_to_payload(diagnostics.prepared_matrix_fingerprint)
        ),
    }


def signalome_clustering_preparation_diagnostics_from_payload(
    payload: object,
    *,
    scope: str,
) -> SignalomeClusteringPreparationDiagnostics:
    if payload is None:
        return default_signalome_clustering_preparation_diagnostics()
    diagnostics_payload = require_mapping(
        payload,
        field_name=f"{scope}.clustering_preparation_diagnostics",
    )
    diagnostics_field_name = f"{scope}.clustering_preparation_diagnostics"
    allowed_fields = frozenset(
        {
            "preparation_policy_id",
            "input_dimension_count",
            "retained_dimension_count",
            "retained_dimension_labels",
            "dropped_fully_missing_dimension_count",
            "dropped_fully_missing_dimension_labels",
            "dropped_fully_missing_dimension_preview",
            "dropped_fully_missing_value_count",
            "non_finite_input_value_count",
            "missing_after_non_finite_normalization_count",
            "imputed_value_count",
            "imputed_value_counts_by_dimension",
            "prepared_matrix_fingerprint",
        }
    )
    _reject_unsupported_fields(
        diagnostics_payload,
        field_name=diagnostics_field_name,
        allowed_fields=allowed_fields,
    )
    _require_fields(
        diagnostics_payload,
        field_name=diagnostics_field_name,
        required_fields=allowed_fields,
    )
    imputation_counts_payload = require_mapping(
        diagnostics_payload.get("imputed_value_counts_by_dimension"),
        field_name=(
            f"{scope}.clustering_preparation_diagnostics."
            "imputed_value_counts_by_dimension"
        ),
    )
    fingerprint_payload = diagnostics_payload.get("prepared_matrix_fingerprint")
    return SignalomeClusteringPreparationDiagnostics(
        preparation_policy_id=require_str(
            diagnostics_payload.get("preparation_policy_id"),
            field_name=(
                f"{scope}.clustering_preparation_diagnostics.preparation_policy_id"
            ),
        ),
        input_dimension_count=_require_int(
            diagnostics_payload.get("input_dimension_count"),
            field_name=(
                f"{scope}.clustering_preparation_diagnostics.input_dimension_count"
            ),
        ),
        retained_dimension_count=_require_int(
            diagnostics_payload.get("retained_dimension_count"),
            field_name=(
                f"{scope}.clustering_preparation_diagnostics.retained_dimension_count"
            ),
        ),
        retained_dimension_labels=_require_string_tuple(
            diagnostics_payload.get("retained_dimension_labels"),
            field_name=(
                f"{scope}.clustering_preparation_diagnostics.retained_dimension_labels"
            ),
        ),
        dropped_fully_missing_dimension_count=_require_int(
            diagnostics_payload.get("dropped_fully_missing_dimension_count"),
            field_name=(
                f"{scope}.clustering_preparation_diagnostics."
                "dropped_fully_missing_dimension_count"
            ),
        ),
        dropped_fully_missing_dimension_labels=_require_string_tuple(
            diagnostics_payload.get("dropped_fully_missing_dimension_labels"),
            field_name=(
                f"{scope}.clustering_preparation_diagnostics."
                "dropped_fully_missing_dimension_labels"
            ),
        ),
        dropped_fully_missing_dimension_preview=_require_string_tuple(
            diagnostics_payload.get("dropped_fully_missing_dimension_preview"),
            field_name=(
                f"{scope}.clustering_preparation_diagnostics."
                "dropped_fully_missing_dimension_preview"
            ),
        ),
        dropped_fully_missing_value_count=_require_int(
            diagnostics_payload.get("dropped_fully_missing_value_count"),
            field_name=(
                f"{scope}.clustering_preparation_diagnostics."
                "dropped_fully_missing_value_count"
            ),
        ),
        non_finite_input_value_count=_require_int(
            diagnostics_payload.get("non_finite_input_value_count"),
            field_name=(
                f"{scope}.clustering_preparation_diagnostics."
                "non_finite_input_value_count"
            ),
        ),
        missing_after_non_finite_normalization_count=_require_int(
            diagnostics_payload.get("missing_after_non_finite_normalization_count"),
            field_name=(
                f"{scope}.clustering_preparation_diagnostics."
                "missing_after_non_finite_normalization_count"
            ),
        ),
        imputed_value_count=_require_int(
            diagnostics_payload.get("imputed_value_count"),
            field_name=(
                f"{scope}.clustering_preparation_diagnostics.imputed_value_count"
            ),
        ),
        imputed_value_counts_by_dimension={
            str(key): _require_int(
                value,
                field_name=(
                    f"{scope}.clustering_preparation_diagnostics."
                    f"imputed_value_counts_by_dimension.{key}"
                ),
            )
            for key, value in imputation_counts_payload.items()
        },
        prepared_matrix_fingerprint=(
            None
            if fingerprint_payload is None
            else table_fingerprint_from_payload(
                require_mapping(
                    fingerprint_payload,
                    field_name=(
                        f"{scope}.clustering_preparation_diagnostics."
                        "prepared_matrix_fingerprint"
                    ),
                )
            )
        ),
    )


def signalome_score_preconditioning_diagnostics_to_payload(
    diagnostics: SignalomeScorePreconditioningDiagnostics,
) -> dict[str, object]:
    return {
        "policy": str(diagnostics.policy),
        "input_row_count": int(diagnostics.input_row_count),
        "dropped_all_missing_row_count": int(diagnostics.dropped_all_missing_row_count),
        "retained_row_count": int(diagnostics.retained_row_count),
    }


def signalome_score_preconditioning_diagnostics_from_payload(
    payload: object,
    *,
    scope: str,
) -> SignalomeScorePreconditioningDiagnostics:
    diagnostics_payload = require_mapping(
        payload,
        field_name=f"{scope}.score_preconditioning_diagnostics",
    )
    diagnostics_field_name = f"{scope}.score_preconditioning_diagnostics"
    allowed_fields = frozenset(
        {
            "policy",
            "input_row_count",
            "dropped_all_missing_row_count",
            "retained_row_count",
        }
    )
    _reject_unsupported_fields(
        diagnostics_payload,
        field_name=diagnostics_field_name,
        allowed_fields=allowed_fields,
    )
    _require_fields(
        diagnostics_payload,
        field_name=diagnostics_field_name,
        required_fields=allowed_fields,
    )
    policy = require_str(
        diagnostics_payload.get("policy"),
        field_name=f"{scope}.score_preconditioning_diagnostics.policy",
    )
    if policy not in SIGNALOME_SCORE_PRECONDITIONING_POLICIES:
        allowed = ", ".join(sorted(SIGNALOME_SCORE_PRECONDITIONING_POLICIES))
        raise PhosPyInputError(
            f"{scope}.score_preconditioning_diagnostics.policy must be one of: "
            f"{allowed}"
        )
    input_row_count = _require_int(
        diagnostics_payload.get("input_row_count"),
        field_name=f"{scope}.score_preconditioning_diagnostics.input_row_count",
    )
    dropped_all_missing_row_count = _require_int(
        diagnostics_payload.get("dropped_all_missing_row_count"),
        field_name=(
            f"{scope}.score_preconditioning_diagnostics.dropped_all_missing_row_count"
        ),
    )
    retained_row_count = _require_int(
        diagnostics_payload.get("retained_row_count"),
        field_name=f"{scope}.score_preconditioning_diagnostics.retained_row_count",
    )
    if (
        dropped_all_missing_row_count < 0
        or retained_row_count < 0
        or input_row_count < 0
        or dropped_all_missing_row_count + retained_row_count != input_row_count
    ):
        raise PhosPyInputError(
            f"{scope}.score_preconditioning_diagnostics counts must be non-negative "
            "and satisfy dropped_all_missing_row_count + retained_row_count = input_row_count"
        )
    return SignalomeScorePreconditioningDiagnostics(
        policy=policy,  # type: ignore[arg-type]
        input_row_count=input_row_count,
        dropped_all_missing_row_count=dropped_all_missing_row_count,
        retained_row_count=retained_row_count,
    )


def signalome_network_correlation_diagnostics_to_payload(
    diagnostics: SignalomeNetworkCorrelationDiagnostics,
) -> dict[str, object]:
    return {
        "total_candidate_correlations": int(diagnostics.total_candidate_correlations),
        "finite_correlations": int(diagnostics.finite_correlations),
        "undefined_correlations": int(diagnostics.undefined_correlations),
        "constant_profile_correlations": int(diagnostics.constant_profile_correlations),
        "insufficient_observation_correlations": int(
            diagnostics.insufficient_observation_correlations
        ),
        "missing_value_correlations": int(diagnostics.missing_value_correlations),
        "non_finite_value_correlations": int(diagnostics.non_finite_value_correlations),
        "edges_created": int(diagnostics.edges_created),
        "edges_skipped_non_finite_correlation": int(
            diagnostics.edges_skipped_non_finite_correlation
        ),
        "edges_skipped_below_threshold": int(diagnostics.edges_skipped_below_threshold),
        "edges_skipped_insufficient_paired_observations": int(
            diagnostics.edges_skipped_insufficient_paired_observations
        ),
        "edges_skipped_constant_profile": int(
            diagnostics.edges_skipped_constant_profile
        ),
        "edges_skipped_missing_score": int(diagnostics.edges_skipped_missing_score),
        "edges_skipped_non_finite_score": int(
            diagnostics.edges_skipped_non_finite_score
        ),
        "edges_skipped_undefined_correlation": int(
            diagnostics.edges_skipped_undefined_correlation
        ),
    }


def signalome_network_correlation_diagnostics_from_payload(
    payload: object,
    *,
    scope: str,
) -> SignalomeNetworkCorrelationDiagnostics:
    diagnostics_payload = require_mapping(
        payload,
        field_name=f"{scope}.network_correlation_diagnostics",
    )
    diagnostics_field_name = f"{scope}.network_correlation_diagnostics"
    allowed_fields = frozenset(
        {
            "total_candidate_correlations",
            "finite_correlations",
            "undefined_correlations",
            "constant_profile_correlations",
            "insufficient_observation_correlations",
            "missing_value_correlations",
            "non_finite_value_correlations",
            "edges_created",
            "edges_skipped_non_finite_correlation",
            "edges_skipped_below_threshold",
            "edges_skipped_insufficient_paired_observations",
            "edges_skipped_constant_profile",
            "edges_skipped_missing_score",
            "edges_skipped_non_finite_score",
            "edges_skipped_undefined_correlation",
        }
    )
    required_fields = frozenset(
        {
            "total_candidate_correlations",
            "finite_correlations",
            "undefined_correlations",
            "constant_profile_correlations",
            "insufficient_observation_correlations",
            "missing_value_correlations",
            "non_finite_value_correlations",
            "edges_created",
            "edges_skipped_non_finite_correlation",
        }
    )
    _reject_unsupported_fields(
        diagnostics_payload,
        field_name=diagnostics_field_name,
        allowed_fields=allowed_fields,
    )
    _require_fields(
        diagnostics_payload,
        field_name=diagnostics_field_name,
        required_fields=required_fields,
    )
    return SignalomeNetworkCorrelationDiagnostics(
        total_candidate_correlations=_require_int(
            diagnostics_payload.get("total_candidate_correlations"),
            field_name=(
                f"{scope}.network_correlation_diagnostics.total_candidate_correlations"
            ),
        ),
        finite_correlations=_require_int(
            diagnostics_payload.get("finite_correlations"),
            field_name=f"{scope}.network_correlation_diagnostics.finite_correlations",
        ),
        undefined_correlations=_require_int(
            diagnostics_payload.get("undefined_correlations"),
            field_name=(
                f"{scope}.network_correlation_diagnostics.undefined_correlations"
            ),
        ),
        constant_profile_correlations=_require_int(
            diagnostics_payload.get("constant_profile_correlations"),
            field_name=(
                f"{scope}.network_correlation_diagnostics.constant_profile_correlations"
            ),
        ),
        insufficient_observation_correlations=_require_int(
            diagnostics_payload.get("insufficient_observation_correlations"),
            field_name=(
                f"{scope}.network_correlation_diagnostics."
                "insufficient_observation_correlations"
            ),
        ),
        missing_value_correlations=_require_int(
            diagnostics_payload.get("missing_value_correlations"),
            field_name=(
                f"{scope}.network_correlation_diagnostics.missing_value_correlations"
            ),
        ),
        non_finite_value_correlations=_require_int(
            diagnostics_payload.get("non_finite_value_correlations"),
            field_name=(
                f"{scope}.network_correlation_diagnostics.non_finite_value_correlations"
            ),
        ),
        edges_created=_require_int(
            diagnostics_payload.get("edges_created"),
            field_name=f"{scope}.network_correlation_diagnostics.edges_created",
        ),
        edges_skipped_non_finite_correlation=_require_int(
            diagnostics_payload.get("edges_skipped_non_finite_correlation"),
            field_name=(
                f"{scope}.network_correlation_diagnostics."
                "edges_skipped_non_finite_correlation"
            ),
        ),
        edges_skipped_below_threshold=_optional_diagnostics_int(
            diagnostics_payload,
            key="edges_skipped_below_threshold",
            field_name=(
                f"{scope}.network_correlation_diagnostics.edges_skipped_below_threshold"
            ),
        ),
        edges_skipped_insufficient_paired_observations=_optional_diagnostics_int(
            diagnostics_payload,
            key="edges_skipped_insufficient_paired_observations",
            field_name=(
                f"{scope}.network_correlation_diagnostics."
                "edges_skipped_insufficient_paired_observations"
            ),
        ),
        edges_skipped_constant_profile=_optional_diagnostics_int(
            diagnostics_payload,
            key="edges_skipped_constant_profile",
            field_name=(
                f"{scope}.network_correlation_diagnostics."
                "edges_skipped_constant_profile"
            ),
        ),
        edges_skipped_missing_score=_optional_diagnostics_int(
            diagnostics_payload,
            key="edges_skipped_missing_score",
            field_name=(
                f"{scope}.network_correlation_diagnostics.edges_skipped_missing_score"
            ),
        ),
        edges_skipped_non_finite_score=_optional_diagnostics_int(
            diagnostics_payload,
            key="edges_skipped_non_finite_score",
            field_name=(
                f"{scope}.network_correlation_diagnostics."
                "edges_skipped_non_finite_score"
            ),
        ),
        edges_skipped_undefined_correlation=_optional_diagnostics_int(
            diagnostics_payload,
            key="edges_skipped_undefined_correlation",
            field_name=(
                f"{scope}.network_correlation_diagnostics."
                "edges_skipped_undefined_correlation"
            ),
        ),
    )


def _optional_diagnostics_int(
    payload: object,
    *,
    key: str,
    field_name: str,
) -> int:
    mapping = require_mapping(payload, field_name=field_name.rsplit(".", 1)[0])
    if key not in mapping:
        return 0
    return _require_int(mapping.get(key), field_name=field_name)


def signalome_alignment_diagnostics_to_payload(
    diagnostics: SignalomeAlignmentDiagnostics,
) -> dict[str, object]:
    return {
        "dataset_sites": _alignment_input_to_payload(diagnostics.dataset_sites),
        "prediction_score_sites": _alignment_input_to_payload(
            diagnostics.prediction_score_sites
        ),
        "downstream_score_sites": _alignment_input_to_payload(
            diagnostics.downstream_score_sites
        ),
        "kinases": _alignment_input_to_payload(diagnostics.kinases),
        "protein_group_ids": _alignment_input_to_payload(diagnostics.protein_group_ids),
    }


def signalome_alignment_diagnostics_from_payload(
    payload: object,
    *,
    scope: str,
) -> SignalomeAlignmentDiagnostics:
    diagnostics_payload = require_mapping(
        payload,
        field_name=f"{scope}.alignment_diagnostics",
    )
    diagnostics_field_name = f"{scope}.alignment_diagnostics"
    required_fields = frozenset(
        {"dataset_sites", "prediction_score_sites", "downstream_score_sites", "kinases"}
    )
    allowed_fields = required_fields | frozenset(
        {"protein_group_ids", "protein_identifiers"}
    )
    _reject_unsupported_fields(
        diagnostics_payload,
        field_name=diagnostics_field_name,
        allowed_fields=allowed_fields,
    )
    _require_fields(
        diagnostics_payload,
        field_name=diagnostics_field_name,
        required_fields=required_fields,
    )
    protein_group_payload = diagnostics_payload.get("protein_group_ids")
    legacy_protein_identifier_payload = diagnostics_payload.get("protein_identifiers")
    if protein_group_payload is None and legacy_protein_identifier_payload is None:
        raise PhosPyInputError(
            f"{diagnostics_field_name} is missing required field(s): "
            "protein_group_ids (legacy alias: protein_identifiers)"
        )
    protein_group_ids = _alignment_input_from_payload(
        (
            protein_group_payload
            if protein_group_payload is not None
            else legacy_protein_identifier_payload
        ),
        field_name=f"{scope}.alignment_diagnostics.protein_group_ids",
    )
    if (
        protein_group_payload is not None
        and legacy_protein_identifier_payload is not None
    ):
        legacy_protein_identifiers = _alignment_input_from_payload(
            legacy_protein_identifier_payload,
            field_name=f"{scope}.alignment_diagnostics.protein_identifiers",
        )
        if legacy_protein_identifiers != protein_group_ids:
            raise PhosPyInputError(
                f"{diagnostics_field_name} has conflicting protein_group_ids "
                "and legacy protein_identifiers diagnostics"
            )
    return SignalomeAlignmentDiagnostics(
        dataset_sites=_alignment_input_from_payload(
            diagnostics_payload.get("dataset_sites"),
            field_name=f"{scope}.alignment_diagnostics.dataset_sites",
        ),
        prediction_score_sites=_alignment_input_from_payload(
            diagnostics_payload.get("prediction_score_sites"),
            field_name=f"{scope}.alignment_diagnostics.prediction_score_sites",
        ),
        downstream_score_sites=_alignment_input_from_payload(
            diagnostics_payload.get("downstream_score_sites"),
            field_name=f"{scope}.alignment_diagnostics.downstream_score_sites",
        ),
        kinases=_alignment_input_from_payload(
            diagnostics_payload.get("kinases"),
            field_name=f"{scope}.alignment_diagnostics.kinases",
        ),
        protein_group_ids=protein_group_ids,
    )


def _alignment_input_to_payload(
    diagnostics: SignalomeAlignmentInputDiagnostics,
) -> dict[str, object]:
    return {
        "provided_count": int(diagnostics.provided_count),
        "retained_count": int(diagnostics.retained_count),
        "dropped_count": int(diagnostics.dropped_count),
        "dropped_reasons": {
            str(reason): int(count)
            for reason, count in diagnostics.dropped_reasons.items()
        },
    }


def _alignment_input_from_payload(
    payload: object,
    *,
    field_name: str,
) -> SignalomeAlignmentInputDiagnostics:
    input_payload = require_mapping(payload, field_name=field_name)
    allowed_fields = frozenset(
        {
            "provided_count",
            "retained_count",
            "dropped_count",
            "dropped_reasons",
        }
    )
    _reject_unsupported_fields(
        input_payload,
        field_name=field_name,
        allowed_fields=allowed_fields,
    )
    _require_fields(
        input_payload,
        field_name=field_name,
        required_fields=allowed_fields,
    )
    provided_count = _require_int(
        input_payload.get("provided_count"),
        field_name=f"{field_name}.provided_count",
    )
    retained_count = _require_int(
        input_payload.get("retained_count"),
        field_name=f"{field_name}.retained_count",
    )
    dropped_count = _require_int(
        input_payload.get("dropped_count"),
        field_name=f"{field_name}.dropped_count",
    )
    if provided_count < 0 or retained_count < 0 or dropped_count < 0:
        raise PhosPyInputError(f"{field_name} counts must be non-negative integers")
    if provided_count != retained_count + dropped_count:
        raise PhosPyInputError(
            f"{field_name} must satisfy provided_count = retained_count + dropped_count"
        )
    dropped_reasons_payload = require_mapping(
        input_payload.get("dropped_reasons"),
        field_name=f"{field_name}.dropped_reasons",
    )
    dropped_reasons = {
        str(reason): _require_int(
            count,
            field_name=f"{field_name}.dropped_reasons.{reason}",
        )
        for reason, count in dropped_reasons_payload.items()
    }
    negative_reasons = sorted(
        reason for reason, count in dropped_reasons.items() if count < 0
    )
    if negative_reasons:
        joined = ", ".join(negative_reasons)
        raise PhosPyInputError(
            f"{field_name}.dropped_reasons contains negative count(s): {joined}"
        )
    return SignalomeAlignmentInputDiagnostics(
        provided_count=provided_count,
        retained_count=retained_count,
        dropped_count=dropped_count,
        dropped_reasons=dropped_reasons,
    )


def _require_string_tuple(value: object, *, field_name: str) -> tuple[str, ...]:
    if not isinstance(value, list | tuple):
        raise PhosPyInputError(f"{field_name} must be an array")
    return tuple(require_str(item, field_name=f"{field_name}[]") for item in value)
