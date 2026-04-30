"""Signalome bundle payload parsing and normalization helpers."""

from __future__ import annotations

import ast
from collections.abc import Mapping

from phospy.api.configs import (
    SIGNALOME_ASSIGNMENT_POLICIES,
    SIGNALOME_CANDIDATE_SCORING_POLICIES,
    SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON,
    SIGNALOME_CLUSTERING_ENGINES,
    SIGNALOME_KINASE_NETWORK_POLICIES,
    SIGNALOME_SCORE_PRECONDITIONING_POLICIES,
    SIGNALOME_TREE_ENGINES,
    SignalomeClusteringConfig,
    SignalomeConfig,
    SignalomeOutputConfig,
    SignalomePerformanceConfig,
    SignalomeScientificConfig,
    SignalomeValidationConfig,
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
    SignalomeAlignmentDiagnostics,
    SignalomeAlignmentInputDiagnostics,
    SignalomeClusterCandidateScore,
    SignalomeModuleSelectionDiagnostics,
    SignalomeNetworkCorrelationDiagnostics,
    SignalomeScorePreconditioningDiagnostics,
)

_SIGNALOME_CONFIG_ALLOWED_FIELDS = frozenset(
    {"scientific", "clustering", "validation", "output", "performance"}
)
_SIGNALOME_CONFIG_REQUIRED_FIELDS = _SIGNALOME_CONFIG_ALLOWED_FIELDS

_SCIENTIFIC_ALLOWED_FIELDS = frozenset(
    {"substrate_support_cutoff", "assignment_policy"}
)
_SCIENTIFIC_REQUIRED_FIELDS = _SCIENTIFIC_ALLOWED_FIELDS

_CLUSTERING_ALLOWED_FIELDS = frozenset(
    {
        "module_count",
        "module_selection_primary_correlation_threshold",
        "module_selection_fallback_correlation_threshold",
        "module_selection_max_clusters",
        "tree_engine",
        "candidate_scoring_policy",
        "clustering_engine",
    }
)
_CLUSTERING_REQUIRED_FIELDS = frozenset(
    field for field in _CLUSTERING_ALLOWED_FIELDS if field not in {"clustering_engine"}
)

_VALIDATION_ALLOWED_FIELDS = frozenset({"score_preconditioning_policy"})
_VALIDATION_REQUIRED_FIELDS = _VALIDATION_ALLOWED_FIELDS

_OUTPUT_ALLOWED_FIELDS = frozenset({"network_correlation_threshold", "network_policy"})
_OUTPUT_REQUIRED_FIELDS = _OUTPUT_ALLOWED_FIELDS

_PERFORMANCE_ALLOWED_FIELDS = frozenset(
    {"max_exact_tree_sites", "max_full_candidate_scoring_sites"}
)
_PERFORMANCE_REQUIRED_FIELDS = _PERFORMANCE_ALLOWED_FIELDS


def signalome_config_from_payload(
    payload: Mapping[str, object],
    *,
    scope: str,
) -> SignalomeConfig:
    """Parse signalome config payload."""

    config_field_name = f"{scope}.signalome_config"
    _reject_unsupported_fields(
        payload,
        field_name=config_field_name,
        allowed_fields=_SIGNALOME_CONFIG_ALLOWED_FIELDS,
    )
    _require_fields(
        payload,
        field_name=config_field_name,
        required_fields=_SIGNALOME_CONFIG_REQUIRED_FIELDS,
    )

    scientific_payload = require_mapping(
        payload.get("scientific"),
        field_name=f"{scope}.signalome_config.scientific",
    )
    _reject_unsupported_fields(
        scientific_payload,
        field_name=f"{scope}.signalome_config.scientific",
        allowed_fields=_SCIENTIFIC_ALLOWED_FIELDS,
    )
    _require_fields(
        scientific_payload,
        field_name=f"{scope}.signalome_config.scientific",
        required_fields=_SCIENTIFIC_REQUIRED_FIELDS,
    )
    assignment_policy = require_str(
        scientific_payload.get("assignment_policy"),
        field_name=f"{scope}.signalome_config.scientific.assignment_policy",
    )
    if assignment_policy not in SIGNALOME_ASSIGNMENT_POLICIES:
        allowed = ", ".join(sorted(SIGNALOME_ASSIGNMENT_POLICIES))
        raise PhosPyInputError(
            f"{scope}.signalome_config.scientific.assignment_policy must be one of: "
            f"{allowed}"
        )

    clustering_payload = require_mapping(
        payload.get("clustering"),
        field_name=f"{scope}.signalome_config.clustering",
    )
    _reject_unsupported_fields(
        clustering_payload,
        field_name=f"{scope}.signalome_config.clustering",
        allowed_fields=_CLUSTERING_ALLOWED_FIELDS,
    )
    _require_fields(
        clustering_payload,
        field_name=f"{scope}.signalome_config.clustering",
        required_fields=_CLUSTERING_REQUIRED_FIELDS,
    )
    tree_engine = require_str(
        clustering_payload.get("tree_engine"),
        field_name=f"{scope}.signalome_config.clustering.tree_engine",
    )
    if tree_engine not in SIGNALOME_TREE_ENGINES:
        allowed = ", ".join(sorted(SIGNALOME_TREE_ENGINES))
        raise PhosPyInputError(
            f"{scope}.signalome_config.clustering.tree_engine must be one of: {allowed}"
        )
    candidate_scoring_policy = require_str(
        clustering_payload.get("candidate_scoring_policy"),
        field_name=f"{scope}.signalome_config.clustering.candidate_scoring_policy",
    )
    if candidate_scoring_policy not in SIGNALOME_CANDIDATE_SCORING_POLICIES:
        allowed = ", ".join(sorted(SIGNALOME_CANDIDATE_SCORING_POLICIES))
        raise PhosPyInputError(
            f"{scope}.signalome_config.clustering.candidate_scoring_policy "
            "must be one of: "
            f"{allowed}"
        )
    clustering_engine = clustering_payload.get("clustering_engine")
    if clustering_engine is None:
        # Keep historical bundle replay deterministic: payloads created before
        # `clustering.clustering_engine` was serialized are interpreted as the
        # legacy exact backend rather than today's public default.
        clustering_engine = SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON
    else:
        clustering_engine = require_str(
            clustering_engine,
            field_name=f"{scope}.signalome_config.clustering.clustering_engine",
        )
    if clustering_engine not in SIGNALOME_CLUSTERING_ENGINES:
        allowed = ", ".join(sorted(SIGNALOME_CLUSTERING_ENGINES))
        raise PhosPyInputError(
            f"{scope}.signalome_config.clustering.clustering_engine must be one "
            f"of: {allowed}"
        )
    module_selection_max_clusters = _require_int(
        clustering_payload.get("module_selection_max_clusters"),
        field_name=f"{scope}.signalome_config.clustering.module_selection_max_clusters",
    )
    module_count = _parse_optional_int(
        clustering_payload.get("module_count"),
        field_name=f"{scope}.signalome_config.clustering.module_count",
    )

    validation_payload = require_mapping(
        payload.get("validation"),
        field_name=f"{scope}.signalome_config.validation",
    )
    _reject_unsupported_fields(
        validation_payload,
        field_name=f"{scope}.signalome_config.validation",
        allowed_fields=_VALIDATION_ALLOWED_FIELDS,
    )
    _require_fields(
        validation_payload,
        field_name=f"{scope}.signalome_config.validation",
        required_fields=_VALIDATION_REQUIRED_FIELDS,
    )
    score_preconditioning_policy = require_str(
        validation_payload.get("score_preconditioning_policy"),
        field_name=(
            f"{scope}.signalome_config.validation.score_preconditioning_policy"
        ),
    )
    if score_preconditioning_policy not in SIGNALOME_SCORE_PRECONDITIONING_POLICIES:
        allowed = ", ".join(sorted(SIGNALOME_SCORE_PRECONDITIONING_POLICIES))
        raise PhosPyInputError(
            f"{scope}.signalome_config.validation.score_preconditioning_policy "
            f"must be one of: {allowed}"
        )

    output_payload = require_mapping(
        payload.get("output"),
        field_name=f"{scope}.signalome_config.output",
    )
    _reject_unsupported_fields(
        output_payload,
        field_name=f"{scope}.signalome_config.output",
        allowed_fields=_OUTPUT_ALLOWED_FIELDS,
    )
    _require_fields(
        output_payload,
        field_name=f"{scope}.signalome_config.output",
        required_fields=_OUTPUT_REQUIRED_FIELDS,
    )
    network_policy = require_str(
        output_payload.get("network_policy"),
        field_name=f"{scope}.signalome_config.output.network_policy",
    )
    if network_policy not in SIGNALOME_KINASE_NETWORK_POLICIES:
        allowed = ", ".join(sorted(SIGNALOME_KINASE_NETWORK_POLICIES))
        raise PhosPyInputError(
            f"{scope}.signalome_config.output.network_policy must be one of: {allowed}"
        )

    performance_payload = require_mapping(
        payload.get("performance"),
        field_name=f"{scope}.signalome_config.performance",
    )
    _reject_unsupported_fields(
        performance_payload,
        field_name=f"{scope}.signalome_config.performance",
        allowed_fields=_PERFORMANCE_ALLOWED_FIELDS,
    )
    _require_fields(
        performance_payload,
        field_name=f"{scope}.signalome_config.performance",
        required_fields=_PERFORMANCE_REQUIRED_FIELDS,
    )
    max_exact_tree_sites = _require_int(
        performance_payload.get("max_exact_tree_sites"),
        field_name=f"{scope}.signalome_config.performance.max_exact_tree_sites",
    )
    max_full_candidate_scoring_sites = _require_int(
        performance_payload.get("max_full_candidate_scoring_sites"),
        field_name=(
            f"{scope}.signalome_config.performance.max_full_candidate_scoring_sites"
        ),
    )

    return SignalomeConfig(
        scientific=SignalomeScientificConfig(
            substrate_support_cutoff=require_float(
                scientific_payload.get("substrate_support_cutoff"),
                field_name=(
                    f"{scope}.signalome_config.scientific.substrate_support_cutoff"
                ),
            ),
            assignment_policy=assignment_policy,  # type: ignore[arg-type]
        ),
        clustering=SignalomeClusteringConfig(
            module_count=module_count,
            module_selection_primary_correlation_threshold=require_float(
                clustering_payload.get(
                    "module_selection_primary_correlation_threshold"
                ),
                field_name=(
                    f"{scope}.signalome_config.clustering."
                    "module_selection_primary_correlation_threshold"
                ),
            ),
            module_selection_fallback_correlation_threshold=require_float(
                clustering_payload.get(
                    "module_selection_fallback_correlation_threshold"
                ),
                field_name=(
                    f"{scope}.signalome_config.clustering."
                    "module_selection_fallback_correlation_threshold"
                ),
            ),
            module_selection_max_clusters=module_selection_max_clusters,
            tree_engine=tree_engine,  # type: ignore[arg-type]
            candidate_scoring_policy=candidate_scoring_policy,  # type: ignore[arg-type]
            clustering_engine=clustering_engine,  # type: ignore[arg-type]
        ),
        validation=SignalomeValidationConfig(
            score_preconditioning_policy=score_preconditioning_policy,  # type: ignore[arg-type]
        ),
        output=SignalomeOutputConfig(
            network_correlation_threshold=require_float(
                output_payload.get("network_correlation_threshold"),
                field_name=(
                    f"{scope}.signalome_config.output.network_correlation_threshold"
                ),
            ),
            network_policy=network_policy,  # type: ignore[arg-type]
        ),
        performance=SignalomePerformanceConfig(
            max_exact_tree_sites=max_exact_tree_sites,
            max_full_candidate_scoring_sites=max_full_candidate_scoring_sites,
        ),
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


def normalize_module_assignments_table(table):
    """Normalize tuple/list/dict-serialized signalome assignment fields."""

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
    raise PhosPyInputError(f"{field_name} must be an int")


def _resolve_optional_alias_string(
    *,
    payload: Mapping[str, object],
    canonical_key: str,
    alias_key: str,
    scope: str,
) -> str | None:
    canonical_value = payload.get(canonical_key)
    alias_value = payload.get(alias_key)
    if canonical_value is None and alias_value is None:
        return None
    canonical = (
        None
        if canonical_value is None
        else require_str(
            canonical_value,
            field_name=f"{scope}.signalome_config.{canonical_key}",
        )
    )
    alias = (
        None
        if alias_value is None
        else require_str(
            alias_value,
            field_name=f"{scope}.signalome_config.{alias_key}",
        )
    )
    if canonical is not None and alias is not None and canonical != alias:
        raise PhosPyInputError(
            f"{scope}.signalome_config.{canonical_key} conflicts with "
            f"{scope}.signalome_config.{alias_key}; provide matching values."
        )
    return canonical if canonical is not None else alias


def _require_int(value: object, *, field_name: str) -> int:
    if isinstance(value, bool):
        raise PhosPyInputError(f"{field_name} must be an int")
    if isinstance(value, int):
        return int(value)
    if isinstance(value, float) and float(value).is_integer():
        return int(value)
    raise PhosPyInputError(f"{field_name} must be an int")


def _reject_unsupported_fields(
    payload: Mapping[str, object],
    *,
    field_name: str,
    allowed_fields: frozenset[str],
) -> None:
    unknown_fields = sorted(
        str(key) for key in payload.keys() if str(key) not in allowed_fields
    )
    if unknown_fields:
        unknown = ", ".join(unknown_fields)
        raise PhosPyInputError(f"{field_name} contains unsupported field(s): {unknown}")


def _require_fields(
    payload: Mapping[str, object],
    *,
    field_name: str,
    required_fields: frozenset[str],
) -> None:
    missing_fields = sorted(
        str(key) for key in required_fields if str(key) not in payload
    )
    if not missing_fields:
        return
    missing = ", ".join(missing_fields)
    raise PhosPyInputError(f"{field_name} is missing required field(s): {missing}")
