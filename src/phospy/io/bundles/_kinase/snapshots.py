"""Config snapshot contract for kinase bundles."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, cast

from phospy.contracts.configs import (
    KINASE_ACTIVITY_KSEA_P_VALUE_METHODS,
    KINASE_ACTIVITY_METHODS,
    KINASE_ACTIVITY_SSGSEA_DEFAULT_ADJUST_P_VALUES,
    KINASE_ACTIVITY_SSGSEA_DEFAULT_MIN_SUBSTRATES,
    KINASE_ACTIVITY_SSGSEA_DEFAULT_PERMUTATIONS,
    KINASE_ACTIVITY_SSGSEA_DEFAULT_RANDOM_SEED,
    KINASE_ACTIVITY_SSGSEA_RANKING_DIRECTION_DESCENDING,
    KINASE_ACTIVITY_SSGSEA_RANKING_DIRECTIONS,
    KINASE_ADAPTIVE_POLICIES,
    KINASE_ATTRITION_POLICY_ON_VIOLATION_MODES,
    KINASE_PREDICTION_MODES,
    KINASE_PROFILE_MISSING_VALUE_STRATEGIES,
    KINASE_PROFILE_SELF_INCLUSION_POLICIES,
    KINASE_PROFILE_SELF_INCLUSION_POLICY_ALLOW,
    KINASE_RELIABILITY_PROFILES,
    KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED,
    KINASE_SCORING_MODES,
    KinaseActivityConfig,
    KinaseActivityMethod,
    KinaseActivityPValueMethod,
    KinaseActivitySsgseaRankingDirection,
    KinaseAdaptivePolicy,
    KinaseAttritionPolicy,
    KinaseAttritionViolationMode,
    KinasePredictionConfig,
    KinasePredictionMode,
    KinaseProfileMissingValueStrategy,
    KinaseReliabilityProfile,
    KinaseScoringConfig,
    KinaseScoringMode,
    LocalisationRequirement,
    ProfileSelfInclusionPolicy,
    ReferenceContextCompatibilityPolicy,
)
from phospy.errors.input import PhosPyInputError
from phospy.io.bundles._shared.primitives import (
    require_bool,
    require_float,
    require_int,
    require_mapping,
    require_str,
)

if TYPE_CHECKING:
    from phospy.contracts.requests import KinaseWorkflowRequest

_SCORING_CONFIG_ALLOWED_FIELDS = frozenset(
    {
        "min_substrates",
        "scoring_mode",
        "include_diagnostic_scoring_tables",
        "include_substrate_contributions",
        "profile_missing_value_strategy",
        "profile_self_inclusion_policy",
        "reliability_profile",
        "requested_reliability_profile",
        "reference_context_compatibility_policy",
        "attrition_policy",
        "localisation_requirement",
        "allow_mixed_total_protein_quantitative_meaning",
    }
)
_SCORING_CONFIG_REQUIRED_FIELDS = frozenset(
    {
        "min_substrates",
        "include_diagnostic_scoring_tables",
        "profile_missing_value_strategy",
    }
)
_ATTRITION_POLICY_ALLOWED_FIELDS = frozenset(
    {
        "minimum_reference_overlap_fraction",
        "minimum_sequence_supported_fraction",
        "minimum_scored_fraction",
        "on_violation",
    }
)
_ATTRITION_POLICY_REQUIRED_FIELDS = _ATTRITION_POLICY_ALLOWED_FIELDS
_LOCALISATION_REQUIREMENT_ALLOWED_FIELDS = frozenset(
    {
        "policy",
        "require_present",
        "minimum_probability",
    }
)
_LOCALISATION_REQUIREMENT_REQUIRED_FIELDS = frozenset(
    {
        "require_present",
        "minimum_probability",
    }
)
_PREDICTION_CONFIG_ALLOWED_FIELDS = frozenset(
    {
        "top_k",
        "deterministic_max_selected_kinases",
        "adaptive_ensemble_runs",
        "mode",
        "adaptive_policy",
        "n_iterations",
        "random_state",
    }
)
_ACTIVITY_CONFIG_ALLOWED_FIELDS = frozenset(
    {
        "enabled",
        "method",
        "threshold",
        "min_substrates",
        "top_n_substrates",
        "ksea_min_substrates",
        "ksea_evidence_threshold",
        "ksea_p_value_method",
        "ksea_adjust_p_values",
        "ssgsea_min_substrates",
        "ssgsea_ranking_direction",
        "ssgsea_permutations",
        "ssgsea_random_seed",
        "ssgsea_adjust_p_values",
    }
)
_ACTIVITY_CONFIG_REQUIRED_FIELDS = frozenset(
    {
        "enabled",
        "method",
        "threshold",
        "min_substrates",
        "top_n_substrates",
        "ksea_min_substrates",
        "ksea_evidence_threshold",
        "ksea_p_value_method",
        "ksea_adjust_p_values",
    }
)
_CONFIG_SNAPSHOT_ALLOWED_FIELDS = frozenset(
    {"scoring_config", "prediction_config", "activity_config"}
)


@dataclass(frozen=True, slots=True)
class KinaseWorkflowConfigSnapshot:
    """Serializable snapshot of the kinase workflow configuration."""

    scoring_config: KinaseScoringConfig
    prediction_config: KinasePredictionConfig
    activity_config: KinaseActivityConfig | None
    _include_scoring_mode: bool = field(default=True, repr=False, compare=False)
    _include_attrition_policy: bool = field(default=False, repr=False, compare=False)

    @classmethod
    def from_request(
        cls, request: KinaseWorkflowRequest
    ) -> KinaseWorkflowConfigSnapshot:
        """Create a config snapshot from a workflow request."""

        from phospy.contracts.requests import KinaseWorkflowRequest

        if not isinstance(request, KinaseWorkflowRequest):
            raise PhosPyInputError(
                "config snapshot request must be a KinaseWorkflowRequest"
            )
        return cls(
            scoring_config=request.scoring_config,
            prediction_config=request.prediction_config,
            activity_config=request.activity_config,
        )

    def to_payload(self) -> dict[str, object]:
        """Return a manifest-safe JSON payload for this config snapshot."""

        activity_payload: dict[str, object] | None
        if self.activity_config is None:
            activity_payload = None
        else:
            activity_payload = {
                "enabled": self.activity_config.enabled,
                "method": self.activity_config.method,
                "threshold": float(self.activity_config.threshold),
                "min_substrates": self.activity_config.min_substrates,
                "top_n_substrates": self.activity_config.top_n_substrates,
                "ksea_min_substrates": self.activity_config.ksea_min_substrates,
                "ksea_evidence_threshold": self.activity_config.ksea_evidence_threshold,
                "ksea_p_value_method": self.activity_config.ksea_p_value_method,
                "ksea_adjust_p_values": self.activity_config.ksea_adjust_p_values,
                "ssgsea_min_substrates": self.activity_config.ssgsea_min_substrates,
                "ssgsea_ranking_direction": (
                    self.activity_config.ssgsea_ranking_direction
                ),
                "ssgsea_permutations": self.activity_config.ssgsea_permutations,
                "ssgsea_random_seed": self.activity_config.ssgsea_random_seed,
                "ssgsea_adjust_p_values": self.activity_config.ssgsea_adjust_p_values,
            }
        scoring_payload: dict[str, object] = {
            "reliability_profile": str(
                self.scoring_config.effective_reliability_profile
            ),
            "requested_reliability_profile": (
                None
                if self.scoring_config.requested_reliability_profile is None
                else str(self.scoring_config.requested_reliability_profile)
            ),
            "min_substrates": self.scoring_config.min_substrates,
            "include_diagnostic_scoring_tables": (
                self.scoring_config.include_diagnostic_scoring_tables
            ),
            "profile_missing_value_strategy": (
                self.scoring_config.profile_missing_value_strategy
            ),
            "profile_self_inclusion_policy": str(
                self.scoring_config.profile_self_inclusion_policy
            ),
            "allow_mixed_total_protein_quantitative_meaning": (
                self.scoring_config.allow_mixed_total_protein_quantitative_meaning
            ),
            "reference_context_compatibility_policy": str(
                self.scoring_config.reference_context_compatibility_policy
            ),
            "localisation_requirement": _localisation_requirement_to_payload(
                self.scoring_config.localisation_requirement
            ),
            "attrition_policy": _attrition_policy_to_payload(
                self.scoring_config.attrition_policy
            ),
        }
        if self.scoring_config.include_substrate_contributions:
            scoring_payload["include_substrate_contributions"] = True
        if self._include_scoring_mode:
            scoring_payload["scoring_mode"] = self.scoring_config.scoring_mode
        return {
            "scoring_config": scoring_payload,
            "prediction_config": {
                "top_k": self.prediction_config.top_k,
                "deterministic_max_selected_kinases": (
                    self.prediction_config.deterministic_max_selected_kinases
                ),
                "adaptive_ensemble_runs": self.prediction_config.adaptive_ensemble_runs,
                "mode": self.prediction_config.mode,
                "adaptive_policy": self.prediction_config.adaptive_policy,
                "n_iterations": self.prediction_config.n_iterations,
                "random_state": self.prediction_config.random_state,
            },
            "activity_config": activity_payload,
        }

    @classmethod
    def from_payload(
        cls, payload: Mapping[str, object]
    ) -> KinaseWorkflowConfigSnapshot:
        """Create a config snapshot from a decoded JSON payload."""

        scope = "config snapshot"
        _reject_unsupported_fields(
            payload,
            field_name=scope,
            allowed_fields=_CONFIG_SNAPSHOT_ALLOWED_FIELDS,
        )
        _require_fields(
            payload,
            field_name=scope,
            required_fields=_CONFIG_SNAPSHOT_ALLOWED_FIELDS,
        )
        scoring_payload = require_mapping(
            payload.get("scoring_config"),
            field_name=f"{scope}.scoring_config",
        )
        _reject_unsupported_fields(
            scoring_payload,
            field_name=f"{scope}.scoring_config",
            allowed_fields=_SCORING_CONFIG_ALLOWED_FIELDS,
        )
        _require_fields(
            scoring_payload,
            field_name=f"{scope}.scoring_config",
            required_fields=_SCORING_CONFIG_REQUIRED_FIELDS,
        )
        prediction_payload = require_mapping(
            payload.get("prediction_config"),
            field_name=f"{scope}.prediction_config",
        )
        _reject_unsupported_fields(
            prediction_payload,
            field_name=f"{scope}.prediction_config",
            allowed_fields=_PREDICTION_CONFIG_ALLOWED_FIELDS,
        )
        _require_fields(
            prediction_payload,
            field_name=f"{scope}.prediction_config",
            required_fields=_PREDICTION_CONFIG_ALLOWED_FIELDS,
        )
        activity_raw = payload.get("activity_config")
        if activity_raw is None:
            activity_config = None
        else:
            activity_payload = require_mapping(
                activity_raw,
                field_name=f"{scope}.activity_config",
            )
            _reject_unsupported_fields(
                activity_payload,
                field_name=f"{scope}.activity_config",
                allowed_fields=_ACTIVITY_CONFIG_ALLOWED_FIELDS,
            )
            _require_fields(
                activity_payload,
                field_name=f"{scope}.activity_config",
                required_fields=_ACTIVITY_CONFIG_REQUIRED_FIELDS,
            )
            activity_config = KinaseActivityConfig(
                enabled=require_bool(
                    activity_payload.get("enabled"),
                    field_name=f"{scope}.activity_config.enabled",
                ),
                method=_parse_activity_method(
                    require_str(
                        activity_payload.get("method"),
                        field_name=f"{scope}.activity_config.method",
                    ),
                    field_name=f"{scope}.activity_config.method",
                ),
                threshold=require_float(
                    activity_payload.get("threshold"),
                    field_name=f"{scope}.activity_config.threshold",
                ),
                min_substrates=require_int(
                    activity_payload.get("min_substrates"),
                    field_name=f"{scope}.activity_config.min_substrates",
                ),
                top_n_substrates=require_int(
                    activity_payload.get("top_n_substrates"),
                    field_name=f"{scope}.activity_config.top_n_substrates",
                ),
                ksea_min_substrates=require_int(
                    activity_payload.get("ksea_min_substrates"),
                    field_name=f"{scope}.activity_config.ksea_min_substrates",
                ),
                ksea_evidence_threshold=(
                    None
                    if activity_payload.get("ksea_evidence_threshold") is None
                    else require_float(
                        activity_payload.get("ksea_evidence_threshold"),
                        field_name=f"{scope}.activity_config.ksea_evidence_threshold",
                    )
                ),
                ksea_p_value_method=_parse_activity_p_value_method(
                    require_str(
                        activity_payload.get("ksea_p_value_method"),
                        field_name=f"{scope}.activity_config.ksea_p_value_method",
                    ),
                    field_name=f"{scope}.activity_config.ksea_p_value_method",
                ),
                ksea_adjust_p_values=require_bool(
                    activity_payload.get("ksea_adjust_p_values"),
                    field_name=f"{scope}.activity_config.ksea_adjust_p_values",
                ),
                ssgsea_min_substrates=require_int(
                    activity_payload.get(
                        "ssgsea_min_substrates",
                        KINASE_ACTIVITY_SSGSEA_DEFAULT_MIN_SUBSTRATES,
                    ),
                    field_name=f"{scope}.activity_config.ssgsea_min_substrates",
                ),
                ssgsea_ranking_direction=_parse_ssgsea_ranking_direction(
                    require_str(
                        activity_payload.get(
                            "ssgsea_ranking_direction",
                            KINASE_ACTIVITY_SSGSEA_RANKING_DIRECTION_DESCENDING,
                        ),
                        field_name=(
                            f"{scope}.activity_config.ssgsea_ranking_direction"
                        ),
                    ),
                    field_name=f"{scope}.activity_config.ssgsea_ranking_direction",
                ),
                ssgsea_permutations=require_int(
                    activity_payload.get(
                        "ssgsea_permutations",
                        KINASE_ACTIVITY_SSGSEA_DEFAULT_PERMUTATIONS,
                    ),
                    field_name=f"{scope}.activity_config.ssgsea_permutations",
                ),
                ssgsea_random_seed=(
                    None
                    if activity_payload.get("ssgsea_random_seed", 0) is None
                    else require_int(
                        activity_payload.get(
                            "ssgsea_random_seed",
                            KINASE_ACTIVITY_SSGSEA_DEFAULT_RANDOM_SEED,
                        ),
                        field_name=f"{scope}.activity_config.ssgsea_random_seed",
                    )
                ),
                ssgsea_adjust_p_values=require_bool(
                    activity_payload.get(
                        "ssgsea_adjust_p_values",
                        KINASE_ACTIVITY_SSGSEA_DEFAULT_ADJUST_P_VALUES,
                    ),
                    field_name=f"{scope}.activity_config.ssgsea_adjust_p_values",
                ),
            )
        reliability_profile = _parse_optional_reliability_profile(
            scoring_payload.get("reliability_profile"),
            field_name=f"{scope}.scoring_config.reliability_profile",
        )
        requested_reliability_profile = _parse_optional_reliability_profile(
            scoring_payload.get("requested_reliability_profile"),
            field_name=f"{scope}.scoring_config.requested_reliability_profile",
        )
        scoring_kwargs: dict[str, object] = {}
        if requested_reliability_profile is not None:
            scoring_kwargs["reliability_profile"] = requested_reliability_profile
        elif (
            reliability_profile is not None and str(reliability_profile) == "production"
        ):
            scoring_kwargs["reliability_profile"] = reliability_profile
        localisation_requirement = _parse_localisation_requirement(
            scoring_payload.get("localisation_requirement"),
            field_name=f"{scope}.scoring_config.localisation_requirement",
        )
        attrition_policy = _parse_attrition_policy(
            scoring_payload.get("attrition_policy"),
            field_name=(f"{scope}.scoring_config.attrition_policy"),
        )
        scoring_config = KinaseScoringConfig(
            min_substrates=require_int(
                scoring_payload.get("min_substrates"),
                field_name=f"{scope}.scoring_config.min_substrates",
            ),
            scoring_mode=_parse_scoring_mode(
                require_str(
                    scoring_payload.get(
                        "scoring_mode",
                        KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED,
                    ),
                    field_name=f"{scope}.scoring_config.scoring_mode",
                ),
                field_name=f"{scope}.scoring_config.scoring_mode",
            ),
            include_diagnostic_scoring_tables=require_bool(
                scoring_payload.get("include_diagnostic_scoring_tables"),
                field_name=(
                    f"{scope}.scoring_config.include_diagnostic_scoring_tables"
                ),
            ),
            include_substrate_contributions=require_bool(
                scoring_payload.get("include_substrate_contributions", False),
                field_name=(f"{scope}.scoring_config.include_substrate_contributions"),
            ),
            profile_missing_value_strategy=(
                _parse_profile_missing_value_strategy(
                    require_str(
                        scoring_payload.get("profile_missing_value_strategy"),
                        field_name=(
                            f"{scope}.scoring_config.profile_missing_value_strategy"
                        ),
                    ),
                    field_name=(
                        f"{scope}.scoring_config.profile_missing_value_strategy"
                    ),
                )
            ),
            profile_self_inclusion_policy=_parse_profile_self_inclusion_policy(
                require_str(
                    scoring_payload.get(
                        "profile_self_inclusion_policy",
                        str(KINASE_PROFILE_SELF_INCLUSION_POLICY_ALLOW),
                    ),
                    field_name=(
                        f"{scope}.scoring_config.profile_self_inclusion_policy"
                    ),
                ),
                field_name=(f"{scope}.scoring_config.profile_self_inclusion_policy"),
            ),
            allow_mixed_total_protein_quantitative_meaning=require_bool(
                scoring_payload.get(
                    "allow_mixed_total_protein_quantitative_meaning", False
                ),
                field_name=(
                    f"{scope}.scoring_config."
                    "allow_mixed_total_protein_quantitative_meaning"
                ),
            ),
            reference_context_compatibility_policy=(
                _parse_reference_context_compatibility_policy(
                    require_str(
                        scoring_payload.get(
                            "reference_context_compatibility_policy",
                            str(
                                ReferenceContextCompatibilityPolicy.REQUIRE_KNOWN_MATCH
                            ),
                        ),
                        field_name=(
                            f"{scope}.scoring_config."
                            "reference_context_compatibility_policy"
                        ),
                    ),
                    field_name=(
                        f"{scope}.scoring_config.reference_context_compatibility_policy"
                    ),
                )
            ),
            attrition_policy=attrition_policy,
            localisation_requirement=localisation_requirement,
            **scoring_kwargs,
        )
        if reliability_profile is not None and (
            scoring_config.effective_reliability_profile is not reliability_profile
        ):
            raise PhosPyInputError(
                f"{scope}.scoring_config.reliability_profile does not match "
                "the effective profile implied by the restored scoring config"
            )
        return cls(
            scoring_config=scoring_config,
            prediction_config=KinasePredictionConfig(
                top_k=require_int(
                    prediction_payload.get("top_k"),
                    field_name=f"{scope}.prediction_config.top_k",
                ),
                deterministic_max_selected_kinases=require_int(
                    prediction_payload.get("deterministic_max_selected_kinases"),
                    field_name=(
                        f"{scope}.prediction_config.deterministic_max_selected_kinases"
                    ),
                ),
                adaptive_ensemble_runs=require_int(
                    prediction_payload.get("adaptive_ensemble_runs"),
                    field_name=f"{scope}.prediction_config.adaptive_ensemble_runs",
                ),
                mode=_parse_prediction_mode(
                    require_str(
                        prediction_payload.get("mode"),
                        field_name=f"{scope}.prediction_config.mode",
                    ),
                    field_name=f"{scope}.prediction_config.mode",
                ),
                adaptive_policy=_parse_adaptive_policy(
                    require_str(
                        prediction_payload.get("adaptive_policy"),
                        field_name=f"{scope}.prediction_config.adaptive_policy",
                    ),
                    field_name=f"{scope}.prediction_config.adaptive_policy",
                ),
                n_iterations=require_int(
                    prediction_payload.get("n_iterations"),
                    field_name=f"{scope}.prediction_config.n_iterations",
                ),
                random_state=(
                    None
                    if prediction_payload.get("random_state") is None
                    else require_int(
                        prediction_payload.get("random_state"),
                        field_name=f"{scope}.prediction_config.random_state",
                    )
                ),
            ),
            activity_config=activity_config,
            _include_scoring_mode="scoring_mode" in scoring_payload,
            _include_attrition_policy="attrition_policy" in scoring_payload,
        )


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


def _parse_profile_missing_value_strategy(
    value: str, *, field_name: str
) -> KinaseProfileMissingValueStrategy:
    if value not in KINASE_PROFILE_MISSING_VALUE_STRATEGIES:
        allowed = ", ".join(sorted(KINASE_PROFILE_MISSING_VALUE_STRATEGIES))
        raise PhosPyInputError(f"{field_name} must be one of: {allowed}")
    return cast(KinaseProfileMissingValueStrategy, value)


def _parse_profile_self_inclusion_policy(
    value: str, *, field_name: str
) -> ProfileSelfInclusionPolicy:
    try:
        return ProfileSelfInclusionPolicy(value)
    except ValueError as exc:
        allowed = ", ".join(
            sorted(str(policy) for policy in KINASE_PROFILE_SELF_INCLUSION_POLICIES)
        )
        raise PhosPyInputError(f"{field_name} must be one of: {allowed}") from exc


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


def _parse_optional_reliability_profile(
    value: object,
    *,
    field_name: str,
) -> KinaseReliabilityProfile | None:
    if value is None:
        return None
    text = require_str(value, field_name=field_name)
    try:
        return KinaseReliabilityProfile(text)
    except ValueError as exc:
        allowed = ", ".join(
            sorted(str(profile) for profile in KINASE_RELIABILITY_PROFILES)
        )
        raise PhosPyInputError(f"{field_name} must be one of: {allowed}") from exc


def _localisation_requirement_to_payload(
    requirement: LocalisationRequirement,
) -> dict[str, object]:
    return {
        "policy": str(requirement.policy),
        "require_present": bool(requirement.require_present),
        "minimum_probability": (
            None
            if requirement.minimum_probability is None
            else float(requirement.minimum_probability)
        ),
    }


def _parse_localisation_requirement(
    value: object,
    *,
    field_name: str,
) -> LocalisationRequirement:
    if value is None:
        return LocalisationRequirement()
    payload = require_mapping(value, field_name=field_name)
    _reject_unsupported_fields(
        payload,
        field_name=field_name,
        allowed_fields=_LOCALISATION_REQUIREMENT_ALLOWED_FIELDS,
    )
    _require_fields(
        payload,
        field_name=field_name,
        required_fields=_LOCALISATION_REQUIREMENT_REQUIRED_FIELDS,
    )
    requirement = LocalisationRequirement(
        require_present=require_bool(
            payload.get("require_present"),
            field_name=f"{field_name}.require_present",
        ),
        minimum_probability=(
            None
            if payload.get("minimum_probability") is None
            else require_float(
                payload.get("minimum_probability"),
                field_name=f"{field_name}.minimum_probability",
            )
        ),
    )
    if "policy" in payload:
        policy = require_str(payload.get("policy"), field_name=f"{field_name}.policy")
        if policy != str(requirement.policy):
            raise PhosPyInputError(
                f"{field_name}.policy must match require_present and "
                "minimum_probability"
            )
    return requirement


def _parse_scoring_mode(value: str, *, field_name: str) -> KinaseScoringMode:
    if value not in KINASE_SCORING_MODES:
        allowed = ", ".join(sorted(KINASE_SCORING_MODES))
        raise PhosPyInputError(f"{field_name} must be one of: {allowed}")
    return cast(KinaseScoringMode, value)


def _parse_attrition_policy(
    value: object,
    *,
    field_name: str,
) -> KinaseAttritionPolicy:
    if value is None:
        return KinaseAttritionPolicy()
    payload = require_mapping(value, field_name=field_name)
    _reject_unsupported_fields(
        payload,
        field_name=field_name,
        allowed_fields=_ATTRITION_POLICY_ALLOWED_FIELDS,
    )
    _require_fields(
        payload,
        field_name=field_name,
        required_fields=_ATTRITION_POLICY_REQUIRED_FIELDS,
    )
    return KinaseAttritionPolicy(
        minimum_reference_overlap_fraction=require_float(
            payload.get("minimum_reference_overlap_fraction"),
            field_name=f"{field_name}.minimum_reference_overlap_fraction",
        ),
        minimum_sequence_supported_fraction=require_float(
            payload.get("minimum_sequence_supported_fraction"),
            field_name=f"{field_name}.minimum_sequence_supported_fraction",
        ),
        minimum_scored_fraction=require_float(
            payload.get("minimum_scored_fraction"),
            field_name=f"{field_name}.minimum_scored_fraction",
        ),
        on_violation=_parse_attrition_violation_mode(
            require_str(
                payload.get("on_violation"),
                field_name=f"{field_name}.on_violation",
            ),
            field_name=f"{field_name}.on_violation",
        ),
    )


def _attrition_policy_to_payload(policy: KinaseAttritionPolicy) -> dict[str, object]:
    return {
        "minimum_reference_overlap_fraction": (
            policy.minimum_reference_overlap_fraction
        ),
        "minimum_sequence_supported_fraction": (
            policy.minimum_sequence_supported_fraction
        ),
        "minimum_scored_fraction": policy.minimum_scored_fraction,
        "on_violation": policy.on_violation,
    }


def _parse_attrition_violation_mode(
    value: str,
    *,
    field_name: str,
) -> KinaseAttritionViolationMode:
    if value not in KINASE_ATTRITION_POLICY_ON_VIOLATION_MODES:
        allowed = ", ".join(sorted(KINASE_ATTRITION_POLICY_ON_VIOLATION_MODES))
        raise PhosPyInputError(f"{field_name} must be one of: {allowed}")
    return cast(KinaseAttritionViolationMode, value)


def _parse_prediction_mode(value: str, *, field_name: str) -> KinasePredictionMode:
    if value not in KINASE_PREDICTION_MODES:
        allowed = ", ".join(sorted(KINASE_PREDICTION_MODES))
        raise PhosPyInputError(f"{field_name} must be one of: {allowed}")
    return cast(KinasePredictionMode, value)


def _parse_adaptive_policy(value: str, *, field_name: str) -> KinaseAdaptivePolicy:
    if value not in KINASE_ADAPTIVE_POLICIES:
        allowed = ", ".join(sorted(KINASE_ADAPTIVE_POLICIES))
        raise PhosPyInputError(f"{field_name} must be one of: {allowed}")
    return cast(KinaseAdaptivePolicy, value)


def _parse_activity_method(value: str, *, field_name: str) -> KinaseActivityMethod:
    if value not in KINASE_ACTIVITY_METHODS:
        allowed = ", ".join(sorted(KINASE_ACTIVITY_METHODS))
        raise PhosPyInputError(f"{field_name} must be one of: {allowed}")
    return cast(KinaseActivityMethod, value)


def _parse_activity_p_value_method(
    value: str, *, field_name: str
) -> KinaseActivityPValueMethod:
    if value not in KINASE_ACTIVITY_KSEA_P_VALUE_METHODS:
        allowed = ", ".join(sorted(KINASE_ACTIVITY_KSEA_P_VALUE_METHODS))
        raise PhosPyInputError(f"{field_name} must be one of: {allowed}")
    return cast(KinaseActivityPValueMethod, value)


def _parse_ssgsea_ranking_direction(
    value: str,
    *,
    field_name: str,
) -> KinaseActivitySsgseaRankingDirection:
    if value not in KINASE_ACTIVITY_SSGSEA_RANKING_DIRECTIONS:
        allowed = ", ".join(sorted(KINASE_ACTIVITY_SSGSEA_RANKING_DIRECTIONS))
        raise PhosPyInputError(f"{field_name} must be one of: {allowed}")
    return cast(KinaseActivitySsgseaRankingDirection, value)
