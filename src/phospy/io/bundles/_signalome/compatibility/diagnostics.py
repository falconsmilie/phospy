"""Signalome diagnostics payload conversion helpers."""

from __future__ import annotations

from phospy.api.configs import SIGNALOME_SCORE_PRECONDITIONING_POLICIES
from phospy.errors.input import PhosPyInputError
from phospy.io.bundles._shared.primitives import (
    require_float,
    require_mapping,
    require_str,
)
from phospy.io.bundles._signalome.compatibility.primitives import (
    _parse_optional_int,
    _reject_unsupported_fields,
    _require_fields,
    _require_int,
)
from phospy.science.signalomes.models import (
    SIGNALOME_MODULE_SELECTION_STRATEGY_CORRELATION_THRESHOLDS,
    SIGNALOME_MODULE_SELECTION_STRATEGY_EXPLICIT_MODULE_COUNT,
    SignalomeAlignmentDiagnostics,
    SignalomeAlignmentInputDiagnostics,
    SignalomeClusterCandidateScore,
    SignalomeModuleSelectionDiagnostics,
    SignalomeNetworkCorrelationDiagnostics,
    SignalomeScorePreconditioningDiagnostics,
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
    )


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
        "protein_identifiers": _alignment_input_to_payload(
            diagnostics.protein_identifiers
        ),
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
    allowed_fields = frozenset(
        {
            "dataset_sites",
            "prediction_score_sites",
            "downstream_score_sites",
            "kinases",
            "protein_identifiers",
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
        protein_identifiers=_alignment_input_from_payload(
            diagnostics_payload.get("protein_identifiers"),
            field_name=f"{scope}.alignment_diagnostics.protein_identifiers",
        ),
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
