"""Isolated migration/compatibility helpers for signalome bundle loading."""

from __future__ import annotations

import ast
from collections.abc import Mapping

from phospy.api.configs import (
    SIGNALOME_ASSIGNMENT_POLICIES,
    SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY,
    SIGNALOME_KINASE_NETWORK_POLICIES,
    SIGNALOME_KINASE_NETWORK_POLICY_SIGNED,
    SIGNALOME_MODULE_SELECTION_FALLBACK_THRESHOLD_DEFAULT,
    SIGNALOME_MODULE_SELECTION_MAX_CLUSTERS_DEFAULT,
    SIGNALOME_MODULE_SELECTION_PRIMARY_THRESHOLD_DEFAULT,
    SIGNALOME_SCORE_PRECONDITIONING_POLICIES,
    SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT,
    SignalomeConfig,
)
from phospy.errors.input import PhosPyInputError
from phospy.io.bundles._shared.primitives import (
    require_float,
    require_mapping,
    require_str,
)
from phospy.signalomes.models import (
    SIGNALOME_MODULE_SELECTION_STRATEGY_CORRELATION_THRESHOLDS,
    SIGNALOME_MODULE_SELECTION_STRATEGY_EXPLICIT_MODULE_COUNT,
    SignalomeClusterCandidateScore,
    SignalomeModuleSelectionDiagnostics,
    SignalomeNetworkCorrelationDiagnostics,
    SignalomeScorePreconditioningDiagnostics,
    default_signalome_module_selection_diagnostics,
    default_signalome_network_correlation_diagnostics,
    default_signalome_score_preconditioning_diagnostics,
)


def signalome_config_from_payload_with_compatibility_support(
    payload: Mapping[str, object],
    *,
    scope: str,
) -> SignalomeConfig:
    """Parse signalome config payload with explicit compatibility cutoff fallback."""

    compatibility_cutoff = payload.get("signalome_cutoff")
    substrate_support_cutoff = payload.get("substrate_support_cutoff")
    network_correlation_threshold = payload.get("network_correlation_threshold")
    if substrate_support_cutoff is None and network_correlation_threshold is None:
        substrate_support_cutoff = compatibility_cutoff
        network_correlation_threshold = compatibility_cutoff
    module_selection_max_clusters_raw = payload.get("module_selection_max_clusters")
    module_selection_max_clusters = (
        SIGNALOME_MODULE_SELECTION_MAX_CLUSTERS_DEFAULT
        if module_selection_max_clusters_raw is None
        else _parse_optional_int(
            module_selection_max_clusters_raw,
            field_name=f"{scope}.signalome_config.module_selection_max_clusters",
        )
    )
    module_count = _parse_optional_int(
        payload.get("module_count"),
        field_name=f"{scope}.signalome_config.module_count",
    )
    assignment_policy_raw = payload.get("assignment_policy")
    assignment_policy = (
        SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY
        if assignment_policy_raw is None
        else require_str(
            assignment_policy_raw,
            field_name=f"{scope}.signalome_config.assignment_policy",
        )
    )
    if assignment_policy not in SIGNALOME_ASSIGNMENT_POLICIES:
        allowed = ", ".join(sorted(SIGNALOME_ASSIGNMENT_POLICIES))
        raise PhosPyInputError(
            f"{scope}.signalome_config.assignment_policy must be one of: {allowed}"
        )
    network_policy_raw = payload.get("network_policy")
    network_policy_field = "network_policy"
    if network_policy_raw is None:
        network_policy_raw = payload.get("kinase_network_policy")
        network_policy_field = "kinase_network_policy"
    network_policy = (
        SIGNALOME_KINASE_NETWORK_POLICY_SIGNED
        if network_policy_raw is None
        else require_str(
            network_policy_raw,
            field_name=f"{scope}.signalome_config.{network_policy_field}",
        )
    )
    if network_policy not in SIGNALOME_KINASE_NETWORK_POLICIES:
        allowed = ", ".join(sorted(SIGNALOME_KINASE_NETWORK_POLICIES))
        raise PhosPyInputError(
            f"{scope}.signalome_config.network_policy must be one of: {allowed}"
        )
    score_preconditioning_policy_raw = payload.get("score_preconditioning_policy")
    score_preconditioning_policy = (
        SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT
        if score_preconditioning_policy_raw is None
        else require_str(
            score_preconditioning_policy_raw,
            field_name=f"{scope}.signalome_config.score_preconditioning_policy",
        )
    )
    if score_preconditioning_policy not in SIGNALOME_SCORE_PRECONDITIONING_POLICIES:
        allowed = ", ".join(sorted(SIGNALOME_SCORE_PRECONDITIONING_POLICIES))
        raise PhosPyInputError(
            f"{scope}.signalome_config.score_preconditioning_policy must be one of: {allowed}"
        )
    return SignalomeConfig(
        substrate_support_cutoff=require_float(
            substrate_support_cutoff,
            field_name=f"{scope}.signalome_config.substrate_support_cutoff",
        ),
        network_correlation_threshold=require_float(
            network_correlation_threshold,
            field_name=f"{scope}.signalome_config.network_correlation_threshold",
        ),
        network_policy=network_policy,  # type: ignore[arg-type]
        assignment_policy=assignment_policy,  # type: ignore[arg-type]
        score_preconditioning_policy=score_preconditioning_policy,  # type: ignore[arg-type]
        module_count=module_count,
        module_selection_primary_correlation_threshold=require_float(
            payload.get("module_selection_primary_correlation_threshold")
            if payload.get("module_selection_primary_correlation_threshold") is not None
            else SIGNALOME_MODULE_SELECTION_PRIMARY_THRESHOLD_DEFAULT,
            field_name=(
                f"{scope}.signalome_config."
                "module_selection_primary_correlation_threshold"
            ),
        ),
        module_selection_fallback_correlation_threshold=require_float(
            payload.get("module_selection_fallback_correlation_threshold")
            if payload.get("module_selection_fallback_correlation_threshold")
            is not None
            else SIGNALOME_MODULE_SELECTION_FALLBACK_THRESHOLD_DEFAULT,
            field_name=(
                f"{scope}.signalome_config."
                "module_selection_fallback_correlation_threshold"
            ),
        ),
        module_selection_max_clusters=module_selection_max_clusters,
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


def signalome_module_selection_diagnostics_from_payload_with_compatibility_support(
    payload: object,
    *,
    scope: str,
) -> SignalomeModuleSelectionDiagnostics:
    if payload is None:
        return default_signalome_module_selection_diagnostics()
    diagnostics_payload = require_mapping(
        payload,
        field_name=f"{scope}.module_selection_diagnostics",
    )
    strategy = require_str(
        diagnostics_payload.get("strategy"),
        field_name=f"{scope}.module_selection_diagnostics.strategy",
    )
    if strategy not in {
        SIGNALOME_MODULE_SELECTION_STRATEGY_CORRELATION_THRESHOLDS,
        SIGNALOME_MODULE_SELECTION_STRATEGY_EXPLICIT_MODULE_COUNT,
    }:
        strategy = SIGNALOME_MODULE_SELECTION_STRATEGY_CORRELATION_THRESHOLDS
    candidate_scores_payload = require_mapping(
        diagnostics_payload.get("candidate_scores", {}),
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
            diagnostics_payload.get("zero_variance_profile_count", 0),
            field_name=(
                f"{scope}.module_selection_diagnostics.zero_variance_profile_count"
            ),
        ),
        near_constant_profile_count=_require_int(
            diagnostics_payload.get("near_constant_profile_count", 0),
            field_name=(
                f"{scope}.module_selection_diagnostics.near_constant_profile_count"
            ),
        ),
        excluded_from_correlation_count=_require_int(
            diagnostics_payload.get("excluded_from_correlation_count", 0),
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


def signalome_score_preconditioning_diagnostics_from_payload_with_compatibility_support(
    payload: object,
    *,
    scope: str,
) -> SignalomeScorePreconditioningDiagnostics:
    if payload is None:
        return default_signalome_score_preconditioning_diagnostics()
    diagnostics_payload = require_mapping(
        payload,
        field_name=f"{scope}.score_preconditioning_diagnostics",
    )
    policy = require_str(
        diagnostics_payload.get(
            "policy", SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT
        ),
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


def signalome_network_correlation_diagnostics_from_payload_with_compatibility_support(
    payload: object,
    *,
    scope: str,
) -> SignalomeNetworkCorrelationDiagnostics:
    if payload is None:
        return default_signalome_network_correlation_diagnostics()
    diagnostics_payload = require_mapping(
        payload,
        field_name=f"{scope}.network_correlation_diagnostics",
    )
    return SignalomeNetworkCorrelationDiagnostics(
        total_candidate_correlations=_require_int(
            diagnostics_payload.get("total_candidate_correlations", 0),
            field_name=(
                f"{scope}.network_correlation_diagnostics.total_candidate_correlations"
            ),
        ),
        finite_correlations=_require_int(
            diagnostics_payload.get("finite_correlations", 0),
            field_name=f"{scope}.network_correlation_diagnostics.finite_correlations",
        ),
        undefined_correlations=_require_int(
            diagnostics_payload.get("undefined_correlations", 0),
            field_name=(
                f"{scope}.network_correlation_diagnostics.undefined_correlations"
            ),
        ),
        constant_profile_correlations=_require_int(
            diagnostics_payload.get("constant_profile_correlations", 0),
            field_name=(
                f"{scope}.network_correlation_diagnostics.constant_profile_correlations"
            ),
        ),
        insufficient_observation_correlations=_require_int(
            diagnostics_payload.get("insufficient_observation_correlations", 0),
            field_name=(
                f"{scope}.network_correlation_diagnostics."
                "insufficient_observation_correlations"
            ),
        ),
        missing_value_correlations=_require_int(
            diagnostics_payload.get("missing_value_correlations", 0),
            field_name=(
                f"{scope}.network_correlation_diagnostics.missing_value_correlations"
            ),
        ),
        non_finite_value_correlations=_require_int(
            diagnostics_payload.get("non_finite_value_correlations", 0),
            field_name=(
                f"{scope}.network_correlation_diagnostics.non_finite_value_correlations"
            ),
        ),
        edges_created=_require_int(
            diagnostics_payload.get("edges_created", 0),
            field_name=f"{scope}.network_correlation_diagnostics.edges_created",
        ),
        edges_skipped_non_finite_correlation=_require_int(
            diagnostics_payload.get("edges_skipped_non_finite_correlation", 0),
            field_name=(
                f"{scope}.network_correlation_diagnostics."
                "edges_skipped_non_finite_correlation"
            ),
        ),
    )


def normalize_module_assignments_table(table):
    """Normalize compatibility tuple/list/dict-serialized signalome assignment fields."""

    normalized = table.copy(deep=True)
    candidate_columns = [
        str(column)
        for column in normalized.columns
        if str(column).endswith("_candidates")
    ]
    for candidates_column in candidate_columns:
        candidates_index = normalized.columns.get_loc(candidates_column)
        candidates = (
            normalized.loc[:, candidates_column]
            .map(_parse_kinase_candidates)
            .astype(object)
        )
        normalized = normalized.drop(columns=[candidates_column])
        normalized.insert(candidates_index, candidates_column, candidates)
    weight_columns = [
        str(column) for column in normalized.columns if str(column).endswith("_weights")
    ]
    for weight_column in weight_columns:
        weight_index = normalized.columns.get_loc(weight_column)
        weights = (
            normalized.loc[:, weight_column].map(_parse_kinase_weights).astype(object)
        )
        normalized = normalized.drop(columns=[weight_column])
        normalized.insert(weight_index, weight_column, weights)
    return normalized


def _parse_kinase_candidates(value: object) -> tuple[str, ...]:
    if isinstance(value, tuple):
        return tuple(str(item) for item in value)
    if isinstance(value, list):
        return tuple(str(item) for item in value)
    raw = str(value).strip()
    if raw == "" or raw.lower() == "nan":
        return ()
    try:
        parsed = ast.literal_eval(raw)
    except (ValueError, SyntaxError):
        return (raw,)
    if isinstance(parsed, tuple):
        return tuple(str(item) for item in parsed)
    if isinstance(parsed, list):
        return tuple(str(item) for item in parsed)
    return (str(parsed),)


def _parse_kinase_weights(value: object) -> tuple[tuple[str, float], ...]:
    if isinstance(value, dict):
        return tuple((str(key), float(weight)) for key, weight in value.items())
    if isinstance(value, (tuple, list)):
        return _normalize_kinase_weight_pairs(value)
    raw = str(value).strip()
    if raw == "" or raw.lower() == "nan":
        return ()
    try:
        parsed = ast.literal_eval(raw)
    except (ValueError, SyntaxError):
        return ()
    if isinstance(parsed, dict):
        return tuple((str(key), float(weight)) for key, weight in parsed.items())
    if isinstance(parsed, (tuple, list)):
        return _normalize_kinase_weight_pairs(parsed)
    return ()


def _normalize_kinase_weight_pairs(
    values: tuple[object, ...] | list[object],
) -> tuple[tuple[str, float], ...]:
    normalized_pairs: list[tuple[str, float]] = []
    for value in values:
        if not isinstance(value, (tuple, list)) or len(value) != 2:
            continue
        kinase, weight = value
        try:
            normalized_pairs.append((str(kinase), float(weight)))
        except (TypeError, ValueError):
            continue
    return tuple(normalized_pairs)


def _parse_optional_int(value: object, *, field_name: str) -> int | None:
    if value is None:
        return None
    if isinstance(value, bool):
        raise PhosPyInputError(f"{field_name} must be an int")
    if isinstance(value, int):
        return int(value)
    if isinstance(value, float) and float(value).is_integer():
        return int(value)
    raise PhosPyInputError(f"{field_name} must be an int")


def _require_int(value: object, *, field_name: str) -> int:
    if isinstance(value, bool):
        raise PhosPyInputError(f"{field_name} must be an int")
    if isinstance(value, int):
        return int(value)
    if isinstance(value, float) and float(value).is_integer():
        return int(value)
    raise PhosPyInputError(f"{field_name} must be an int")
