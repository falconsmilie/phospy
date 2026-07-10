"""Signalome bundle config payload parsing helpers."""

from __future__ import annotations

from collections.abc import Mapping

from phospy.api.configs import (
    SIGNALOME_ASSIGNMENT_POLICIES,
    SIGNALOME_CANDIDATE_SCORING_POLICIES,
    SIGNALOME_CLUSTERING_ENGINES,
    SIGNALOME_KINASE_NETWORK_POLICIES,
    SIGNALOME_SCORE_PRECONDITIONING_POLICIES,
    ReferenceContextCompatibilityPolicy,
    SignalomeClusteringConfig,
    SignalomeConfig,
    SignalomeOutputConfig,
    SignalomePerformanceConfig,
    SignalomeScientificConfig,
    SignalomeValidationConfig,
)
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
        "candidate_scoring_policy",
        "clustering_engine",
    }
)
_CLUSTERING_REQUIRED_FIELDS = _CLUSTERING_ALLOWED_FIELDS

_VALIDATION_ALLOWED_FIELDS = frozenset(
    {
        "score_preconditioning_policy",
        "allow_mixed_total_protein_quantitative_meaning",
        "reference_context_compatibility_policy",
    }
)
_VALIDATION_REQUIRED_FIELDS = _VALIDATION_ALLOWED_FIELDS

_OUTPUT_ALLOWED_FIELDS = frozenset(
    {
        "network_correlation_threshold",
        "network_policy",
        "network_min_paired_finite_observations",
    }
)
_OUTPUT_REQUIRED_FIELDS = frozenset({"network_correlation_threshold", "network_policy"})

_PERFORMANCE_ALLOWED_FIELDS = frozenset(
    {"max_exact_tree_sites", "max_full_candidate_scoring_sites"}
)
_PERFORMANCE_REQUIRED_FIELDS = _PERFORMANCE_ALLOWED_FIELDS


def _parse_reference_context_compatibility_policy(
    value: str, *, field_name: str
) -> ReferenceContextCompatibilityPolicy:
    try:
        return ReferenceContextCompatibilityPolicy(value)
    except ValueError as exc:
        allowed = ", ".join(
            sorted(str(policy) for policy in ReferenceContextCompatibilityPolicy)
        )
        raise PhosPyInputError(f"{field_name} must be one of: {allowed}") from exc


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
    clustering_engine = require_str(
        clustering_payload.get("clustering_engine"),
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
    allow_mixed_total_protein_quantitative_meaning = require_bool(
        validation_payload.get("allow_mixed_total_protein_quantitative_meaning"),
        field_name=(
            f"{scope}.signalome_config.validation."
            "allow_mixed_total_protein_quantitative_meaning"
        ),
    )
    reference_context_compatibility_policy = (
        _parse_reference_context_compatibility_policy(
            require_str(
                validation_payload.get("reference_context_compatibility_policy"),
                field_name=(
                    f"{scope}.signalome_config.validation."
                    "reference_context_compatibility_policy"
                ),
            ),
            field_name=(
                f"{scope}.signalome_config.validation."
                "reference_context_compatibility_policy"
            ),
        )
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
    network_min_paired_finite_observations = _parse_optional_int(
        output_payload.get("network_min_paired_finite_observations"),
        field_name=(
            f"{scope}.signalome_config.output.network_min_paired_finite_observations"
        ),
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
            candidate_scoring_policy=candidate_scoring_policy,  # type: ignore[arg-type]
            clustering_engine=clustering_engine,  # type: ignore[arg-type]
        ),
        validation=SignalomeValidationConfig(
            score_preconditioning_policy=score_preconditioning_policy,  # type: ignore[arg-type]
            allow_mixed_total_protein_quantitative_meaning=(
                allow_mixed_total_protein_quantitative_meaning
            ),
            reference_context_compatibility_policy=(
                reference_context_compatibility_policy
            ),
        ),
        output=SignalomeOutputConfig(
            network_correlation_threshold=require_float(
                output_payload.get("network_correlation_threshold"),
                field_name=(
                    f"{scope}.signalome_config.output.network_correlation_threshold"
                ),
            ),
            network_policy=network_policy,  # type: ignore[arg-type]
            network_min_paired_finite_observations=(
                network_min_paired_finite_observations
            ),
        ),
        performance=SignalomePerformanceConfig(
            max_exact_tree_sites=max_exact_tree_sites,
            max_full_candidate_scoring_sites=max_full_candidate_scoring_sites,
        ),
    )
