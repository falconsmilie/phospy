"""Contracts and policy declarations for peptide-evidence dataset resolution."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from enum import Enum
from typing import Any, cast

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.science.evidence.multi_site import (
    MULTI_SITE_POLICY_ERROR,
    MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING,
    MULTI_SITE_POLICY_KEEP_JOINT,
    MULTI_SITE_POLICY_SPLIT_EQUAL_WEIGHT,
    MultiSiteHandlingConfig,
)

DATASET_SITE_RESOLUTION_MODE_SITE_LEVEL_RESOLVED = "site_level_resolved"
DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE = "peptide_evidence"
SUPPORTED_DATASET_SITE_RESOLUTION_MODES: tuple[str, ...] = (
    DATASET_SITE_RESOLUTION_MODE_SITE_LEVEL_RESOLVED,
    DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
)

DATASET_MULTI_SITE_POLICY_REJECT = "reject"
DATASET_MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING = (
    "exclude_from_sequence_scoring"
)
DATASET_MULTI_SITE_POLICY_KEEP_JOINT = "keep_joint"
DATASET_MULTI_SITE_POLICY_SPLIT = "split"
SUPPORTED_DATASET_MULTI_SITE_POLICIES: tuple[str, ...] = (
    DATASET_MULTI_SITE_POLICY_REJECT,
    DATASET_MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING,
    DATASET_MULTI_SITE_POLICY_KEEP_JOINT,
    DATASET_MULTI_SITE_POLICY_SPLIT,
)

_POLICY_TO_MULTI_SITE_HANDLING_POLICY: dict[str, str] = {
    DATASET_MULTI_SITE_POLICY_REJECT: MULTI_SITE_POLICY_ERROR,
    DATASET_MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING: (
        MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING
    ),
    DATASET_MULTI_SITE_POLICY_KEEP_JOINT: MULTI_SITE_POLICY_KEEP_JOINT,
    DATASET_MULTI_SITE_POLICY_SPLIT: MULTI_SITE_POLICY_SPLIT_EQUAL_WEIGHT,
}

DATASET_PEPTIDE_TO_SITE_AGGREGATION_POLICY_MAPPING_WEIGHTED_MEAN = (
    "mapping_weighted_mean"
)
DATASET_PEPTIDE_TO_SITE_AGGREGATION_POLICY_LEGACY_ALIAS = (
    "legacy_alias_for_arithmetic_mean_of_allocated_signals"
)
DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_POLICY_EXPLICIT_OR_DERIVED_EQUAL = (
    "explicit_mapping_weight_when_supplied_else_equal_fraction_per_resolved_site"
)
DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_EXPLICIT = "explicit_mapping_weight"
DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_DERIVED_EQUAL = (
    "derived_equal_weight_per_mapped_site"
)
DATASET_PEPTIDE_MAPPING_WEIGHT_NORMALIZATION_POLICY_SUM_TO_ONE_PER_PEPTIDE_ROW = (
    "sum_to_one_per_peptide_evidence_row"
)
DATASET_PEPTIDE_MAPPING_WEIGHT_NORMALISATION_UNIT_PER_PEPTIDE = (
    DATASET_PEPTIDE_MAPPING_WEIGHT_NORMALIZATION_POLICY_SUM_TO_ONE_PER_PEPTIDE_ROW
)
DATASET_PEPTIDE_SIGNAL_ALLOCATION_POLICY_MULTIPLY_BY_MAPPING_FRACTION = (
    "multiply_peptide_signal_by_mapping_fraction"
)
DATASET_PEPTIDE_SITE_SUMMARISATION_POLICY_ARITHMETIC_MEAN_OF_ALLOCATED_SIGNALS = (
    "arithmetic_mean_of_allocated_signals"
)
DATASET_PEPTIDE_DUPLICATE_EVIDENCE_POLICY_RETAIN_DUPLICATE_ROWS = (
    "retain_duplicate_peptide_evidence_rows_as_separate_observations"
)
DATASET_PEPTIDE_DUPLICATE_POLICY_RETAIN_ALL_ROWS = (
    DATASET_PEPTIDE_DUPLICATE_EVIDENCE_POLICY_RETAIN_DUPLICATE_ROWS
)
DATASET_PEPTIDE_MIXED_AMBIGUITY_POLICY_COMBINE_ALLOCATED_SIGNALS = (
    "combine_ambiguous_and_unambiguous_allocated_signals_in_site_mean"
)
DATASET_PEPTIDE_MIXED_AMBIGUITY_POLICY_SHARED_WEIGHTED_MEAN = (
    DATASET_PEPTIDE_MIXED_AMBIGUITY_POLICY_COMBINE_ALLOCATED_SIGNALS
)
DATASET_PEPTIDE_LOCALISATION_AGGREGATION_POLICY_ARITHMETIC_MEAN_OF_FINITE_VALUES = (
    "arithmetic_mean_of_finite_reported_localisation_values"
)
DATASET_PEPTIDE_SITE_SEQUENCE_POLICY_VALIDATE_WITHOUT_REPAIR = "validate_without_repair"

_LEGACY_MAPPING_WEIGHT_NORMALISATION_UNIT_PER_PEPTIDE = "sum_to_one_per_peptide_row"
_LEGACY_DUPLICATE_PEPTIDE_POLICY_RETAIN_ALL_ROWS = (
    "retain_all_peptide_rows_as_independent_observations"
)
_LEGACY_MIXED_AMBIGUITY_POLICY_SHARED_WEIGHTED_MEAN = (
    "mixed_ambiguous_and_unambiguous_rows_share_same_weighted_mean_aggregation"
)
MAPPING_FRACTION_COLUMN = "mapping_fraction"
MAPPING_WEIGHT_SUM_TOLERANCE = 1e-6
SITE_SEQUENCE_SOURCE_PROVIDED = "provided_site_sequence"
SITE_SEQUENCE_SOURCE_PEPTIDE_CONTEXT = "peptide_sequence_site_string"
SITE_SEQUENCE_SOURCE_MISSING = "missing"


class _MappingWeightSourcePolicy(str, Enum):
    EXPLICIT_OR_DERIVED_EQUAL = (
        DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_POLICY_EXPLICIT_OR_DERIVED_EQUAL
    )


class _MappingWeightNormalizationPolicy(str, Enum):
    SUM_TO_ONE_PER_PEPTIDE_ROW = (
        DATASET_PEPTIDE_MAPPING_WEIGHT_NORMALIZATION_POLICY_SUM_TO_ONE_PER_PEPTIDE_ROW
    )


class _SignalAllocationPolicy(str, Enum):
    MULTIPLY_BY_MAPPING_FRACTION = (
        DATASET_PEPTIDE_SIGNAL_ALLOCATION_POLICY_MULTIPLY_BY_MAPPING_FRACTION
    )


class _SiteSummarisationPolicy(str, Enum):
    ARITHMETIC_MEAN_OF_ALLOCATED_SIGNALS = (
        DATASET_PEPTIDE_SITE_SUMMARISATION_POLICY_ARITHMETIC_MEAN_OF_ALLOCATED_SIGNALS
    )


class _DuplicateEvidencePolicy(str, Enum):
    RETAIN_DUPLICATE_ROWS = (
        DATASET_PEPTIDE_DUPLICATE_EVIDENCE_POLICY_RETAIN_DUPLICATE_ROWS
    )


class _MixedAmbiguityPolicy(str, Enum):
    COMBINE_ALLOCATED_SIGNALS = (
        DATASET_PEPTIDE_MIXED_AMBIGUITY_POLICY_COMBINE_ALLOCATED_SIGNALS
    )


class _LocalisationAggregationPolicy(str, Enum):
    ARITHMETIC_MEAN_OF_FINITE_VALUES = (
        DATASET_PEPTIDE_LOCALISATION_AGGREGATION_POLICY_ARITHMETIC_MEAN_OF_FINITE_VALUES
    )


@dataclass(frozen=True, slots=True)
class _PeptideEvidenceResolutionPolicies:
    mapping_weight_source_policy: _MappingWeightSourcePolicy
    mapping_weight_normalization_policy: _MappingWeightNormalizationPolicy
    signal_allocation_policy: _SignalAllocationPolicy
    site_summarisation_policy: _SiteSummarisationPolicy
    duplicate_evidence_policy: _DuplicateEvidencePolicy
    mixed_ambiguity_policy: _MixedAmbiguityPolicy
    localisation_aggregation_policy: _LocalisationAggregationPolicy
    legacy_aggregation_policy: str
    aggregation_formula: str

    def to_payload(self) -> dict[str, str]:
        return {
            "mapping_weight_source_policy": (self.mapping_weight_source_policy.value),
            "mapping_weight_normalization_policy": (
                self.mapping_weight_normalization_policy.value
            ),
            "signal_allocation_policy": self.signal_allocation_policy.value,
            "site_summarisation_policy": self.site_summarisation_policy.value,
            "duplicate_evidence_policy": self.duplicate_evidence_policy.value,
            "mixed_ambiguity_policy": self.mixed_ambiguity_policy.value,
            "localisation_aggregation_policy": (
                self.localisation_aggregation_policy.value
            ),
            "aggregation_policy": self.legacy_aggregation_policy,
            "aggregation_formula": self.aggregation_formula,
        }


CURRENT_RESOLUTION_POLICIES = _PeptideEvidenceResolutionPolicies(
    mapping_weight_source_policy=(_MappingWeightSourcePolicy.EXPLICIT_OR_DERIVED_EQUAL),
    mapping_weight_normalization_policy=(
        _MappingWeightNormalizationPolicy.SUM_TO_ONE_PER_PEPTIDE_ROW
    ),
    signal_allocation_policy=_SignalAllocationPolicy.MULTIPLY_BY_MAPPING_FRACTION,
    site_summarisation_policy=(
        _SiteSummarisationPolicy.ARITHMETIC_MEAN_OF_ALLOCATED_SIGNALS
    ),
    duplicate_evidence_policy=_DuplicateEvidencePolicy.RETAIN_DUPLICATE_ROWS,
    mixed_ambiguity_policy=_MixedAmbiguityPolicy.COMBINE_ALLOCATED_SIGNALS,
    localisation_aggregation_policy=(
        _LocalisationAggregationPolicy.ARITHMETIC_MEAN_OF_FINITE_VALUES
    ),
    legacy_aggregation_policy=(DATASET_PEPTIDE_TO_SITE_AGGREGATION_POLICY_LEGACY_ALIAS),
    aggregation_formula=(
        "a[p,s,j] = mapping_fraction[p,s] * peptide_signal[p,j]; "
        "site_signal[s,j] = arithmetic_mean(a[p,s,j] for retained peptide rows "
        "p mapped to s)"
    ),
)


@dataclass(frozen=True, slots=True)
class PeptideEvidenceResolutionSummary:
    """Structured summary for peptide-to-site resolution provenance."""

    input_mode: str
    multi_site_policy: str | None
    peptide_observations_received: int
    unique_site_ids_produced: int
    ambiguous_observations: int
    excluded_observations: int
    split_observations: int
    mapping_weight_source_policy: str
    mapping_weight_normalization_policy: str
    signal_allocation_policy: str
    site_summarisation_policy: str
    duplicate_evidence_policy: str
    mixed_ambiguity_policy: str
    localisation_aggregation_policy: str
    aggregation_policy: str
    aggregation_formula: str
    mapping_weight_source: str
    mapping_weight_normalisation: str
    duplicate_peptide_policy: str
    duplicate_peptide_rows: int
    site_sequence_column_present: bool
    provided_site_sequence_count: int
    accepted_site_sequence_count: int
    rejected_site_sequence_count: int
    provided_site_sequence_used_count: int
    peptide_context_derived_site_sequence_count: int
    missing_site_sequence_count: int
    site_sequence_policy: str

    def to_payload(self) -> dict[str, object]:
        return {
            "input_mode": self.input_mode,
            "multi_site_policy": self.multi_site_policy,
            "peptide_observations_received": int(self.peptide_observations_received),
            "unique_site_ids_produced": int(self.unique_site_ids_produced),
            "ambiguous_observations": int(self.ambiguous_observations),
            "excluded_observations": int(self.excluded_observations),
            "split_observations": int(self.split_observations),
            "mapping_weight_source_policy": self.mapping_weight_source_policy,
            "mapping_weight_normalization_policy": (
                self.mapping_weight_normalization_policy
            ),
            "signal_allocation_policy": self.signal_allocation_policy,
            "site_summarisation_policy": self.site_summarisation_policy,
            "duplicate_evidence_policy": self.duplicate_evidence_policy,
            "mixed_ambiguity_policy": self.mixed_ambiguity_policy,
            "localisation_aggregation_policy": self.localisation_aggregation_policy,
            "aggregation_policy": self.aggregation_policy,
            "aggregation_formula": self.aggregation_formula,
            "mapping_weight_source": self.mapping_weight_source,
            "mapping_weight_normalisation": self.mapping_weight_normalisation,
            "duplicate_peptide_policy": self.duplicate_peptide_policy,
            "duplicate_peptide_rows": int(self.duplicate_peptide_rows),
            "site_sequence_column_present": bool(self.site_sequence_column_present),
            "provided_site_sequence_count": int(self.provided_site_sequence_count),
            "accepted_site_sequence_count": int(self.accepted_site_sequence_count),
            "rejected_site_sequence_count": int(self.rejected_site_sequence_count),
            "provided_site_sequence_used_count": int(
                self.provided_site_sequence_used_count
            ),
            "peptide_context_derived_site_sequence_count": int(
                self.peptide_context_derived_site_sequence_count
            ),
            "missing_site_sequence_count": int(self.missing_site_sequence_count),
            "site_sequence_policy": self.site_sequence_policy,
        }

    @classmethod
    def from_payload(
        cls,
        payload: Mapping[str, object],
    ) -> PeptideEvidenceResolutionSummary:
        """Deserialize current or legacy peptide-resolution summary payloads."""

        normalized = normalise_peptide_evidence_resolution_summary_payload(payload)
        return cls(
            input_mode=str(normalized["input_mode"]),
            multi_site_policy=(
                None
                if normalized.get("multi_site_policy") is None
                else str(normalized["multi_site_policy"])
            ),
            peptide_observations_received=_payload_int(
                normalized["peptide_observations_received"]
            ),
            unique_site_ids_produced=_payload_int(
                normalized["unique_site_ids_produced"]
            ),
            ambiguous_observations=_payload_int(normalized["ambiguous_observations"]),
            excluded_observations=_payload_int(normalized["excluded_observations"]),
            split_observations=_payload_int(normalized["split_observations"]),
            mapping_weight_source_policy=str(
                normalized["mapping_weight_source_policy"]
            ),
            mapping_weight_normalization_policy=str(
                normalized["mapping_weight_normalization_policy"]
            ),
            signal_allocation_policy=str(normalized["signal_allocation_policy"]),
            site_summarisation_policy=str(normalized["site_summarisation_policy"]),
            duplicate_evidence_policy=str(normalized["duplicate_evidence_policy"]),
            mixed_ambiguity_policy=str(normalized["mixed_ambiguity_policy"]),
            localisation_aggregation_policy=str(
                normalized["localisation_aggregation_policy"]
            ),
            aggregation_policy=str(normalized["aggregation_policy"]),
            aggregation_formula=str(normalized["aggregation_formula"]),
            mapping_weight_source=str(normalized["mapping_weight_source"]),
            mapping_weight_normalisation=str(
                normalized["mapping_weight_normalisation"]
            ),
            duplicate_peptide_policy=str(normalized["duplicate_peptide_policy"]),
            duplicate_peptide_rows=_payload_int(normalized["duplicate_peptide_rows"]),
            site_sequence_column_present=bool(
                normalized["site_sequence_column_present"]
            ),
            provided_site_sequence_count=_payload_int(
                normalized["provided_site_sequence_count"]
            ),
            accepted_site_sequence_count=_payload_int(
                normalized["accepted_site_sequence_count"]
            ),
            rejected_site_sequence_count=_payload_int(
                normalized["rejected_site_sequence_count"]
            ),
            provided_site_sequence_used_count=_payload_int(
                normalized["provided_site_sequence_used_count"]
            ),
            peptide_context_derived_site_sequence_count=_payload_int(
                normalized["peptide_context_derived_site_sequence_count"]
            ),
            missing_site_sequence_count=_payload_int(
                normalized["missing_site_sequence_count"]
            ),
            site_sequence_policy=str(normalized["site_sequence_policy"]),
        )


def _payload_int(value: object) -> int:
    return int(cast(Any, value))


def normalise_peptide_evidence_resolution_summary_payload(
    payload: Mapping[str, object],
) -> dict[str, object]:
    """Return a current-policy peptide-resolution payload from old or new data."""

    normalized = dict(payload)
    current = CURRENT_RESOLUTION_POLICIES.to_payload()

    normalized["mapping_weight_source_policy"] = str(
        normalized.get("mapping_weight_source_policy")
        or current["mapping_weight_source_policy"]
    )

    legacy_normalization = normalized.get("mapping_weight_normalisation")
    if "mapping_weight_normalization_policy" not in normalized:
        normalized["mapping_weight_normalization_policy"] = (
            current["mapping_weight_normalization_policy"]
            if legacy_normalization
            in (
                None,
                _LEGACY_MAPPING_WEIGHT_NORMALISATION_UNIT_PER_PEPTIDE,
                current["mapping_weight_normalization_policy"],
            )
            else str(legacy_normalization)
        )
    normalized["mapping_weight_normalisation"] = str(
        normalized.get("mapping_weight_normalisation")
        or normalized["mapping_weight_normalization_policy"]
    )
    if (
        normalized["mapping_weight_normalisation"]
        == _LEGACY_MAPPING_WEIGHT_NORMALISATION_UNIT_PER_PEPTIDE
    ):
        normalized["mapping_weight_normalisation"] = str(
            normalized["mapping_weight_normalization_policy"]
        )

    normalized["signal_allocation_policy"] = str(
        normalized.get("signal_allocation_policy")
        or current["signal_allocation_policy"]
    )
    normalized["site_summarisation_policy"] = str(
        normalized.get("site_summarisation_policy")
        or current["site_summarisation_policy"]
    )

    legacy_duplicate_policy = normalized.get("duplicate_peptide_policy")
    if "duplicate_evidence_policy" not in normalized:
        normalized["duplicate_evidence_policy"] = (
            current["duplicate_evidence_policy"]
            if legacy_duplicate_policy
            in (
                None,
                _LEGACY_DUPLICATE_PEPTIDE_POLICY_RETAIN_ALL_ROWS,
                current["duplicate_evidence_policy"],
            )
            else str(legacy_duplicate_policy)
        )
    normalized["duplicate_peptide_policy"] = str(
        normalized.get("duplicate_peptide_policy")
        or normalized["duplicate_evidence_policy"]
    )
    if (
        normalized["duplicate_peptide_policy"]
        == _LEGACY_DUPLICATE_PEPTIDE_POLICY_RETAIN_ALL_ROWS
    ):
        normalized["duplicate_peptide_policy"] = str(
            normalized["duplicate_evidence_policy"]
        )

    mixed_ambiguity_policy = normalized.get("mixed_ambiguity_policy")
    normalized["mixed_ambiguity_policy"] = (
        current["mixed_ambiguity_policy"]
        if mixed_ambiguity_policy
        in (
            None,
            _LEGACY_MIXED_AMBIGUITY_POLICY_SHARED_WEIGHTED_MEAN,
            current["mixed_ambiguity_policy"],
        )
        else str(mixed_ambiguity_policy)
    )
    normalized["localisation_aggregation_policy"] = str(
        normalized.get("localisation_aggregation_policy")
        or current["localisation_aggregation_policy"]
    )

    aggregation_policy = normalized.get("aggregation_policy")
    normalized["aggregation_policy"] = (
        current["aggregation_policy"]
        if aggregation_policy
        in (
            None,
            DATASET_PEPTIDE_TO_SITE_AGGREGATION_POLICY_MAPPING_WEIGHTED_MEAN,
            current["aggregation_policy"],
        )
        else str(aggregation_policy)
    )
    normalized["aggregation_formula"] = str(
        normalized.get("aggregation_formula") or current["aggregation_formula"]
    )
    return normalized


@dataclass(frozen=True, slots=True)
class PeptideEvidenceResolutionResult:
    """Site-level matrices produced from peptide-level evidence."""

    phospho: pd.DataFrame
    site_metadata: pd.DataFrame
    summary: PeptideEvidenceResolutionSummary


@dataclass(frozen=True, slots=True)
class PeptideEvidenceResolutionInputMetrics:
    """Request-level counts collected before scientific resolution stages."""

    peptide_observations_received: int
    ambiguous_observations: int
    excluded_observations: int
    split_observations: int
    duplicate_peptide_rows: int
    site_sequence_column_present: bool
    provided_site_sequence_count: int


def build_multi_site_handling_config_for_dataset_policy(
    *,
    multi_site_policy: str,
) -> MultiSiteHandlingConfig:
    """Translate dataset-builder multi-site policy to evidence config."""

    validate_dataset_multi_site_policy(
        multi_site_policy,
        field_name="dataset build request multi_site_policy",
    )
    resolved_policy = _POLICY_TO_MULTI_SITE_HANDLING_POLICY[multi_site_policy]
    return MultiSiteHandlingConfig(
        statistical_modeling_policy=resolved_policy,
        kinase_sequence_scoring_policy=resolved_policy,
    )


def validate_dataset_multi_site_policy(policy: object, *, field_name: str) -> None:
    if (
        not isinstance(policy, str)
        or policy not in SUPPORTED_DATASET_MULTI_SITE_POLICIES
    ):
        supported = ", ".join(
            repr(value) for value in SUPPORTED_DATASET_MULTI_SITE_POLICIES
        )
        raise PhosPyInputError(f"{field_name} must be one of: {supported}")


__all__ = [
    "CURRENT_RESOLUTION_POLICIES",
    "DATASET_MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING",
    "DATASET_MULTI_SITE_POLICY_KEEP_JOINT",
    "DATASET_MULTI_SITE_POLICY_REJECT",
    "DATASET_MULTI_SITE_POLICY_SPLIT",
    "DATASET_PEPTIDE_DUPLICATE_EVIDENCE_POLICY_RETAIN_DUPLICATE_ROWS",
    "DATASET_PEPTIDE_DUPLICATE_POLICY_RETAIN_ALL_ROWS",
    "DATASET_PEPTIDE_LOCALISATION_AGGREGATION_POLICY_ARITHMETIC_MEAN_OF_FINITE_VALUES",
    "DATASET_PEPTIDE_MAPPING_WEIGHT_NORMALISATION_UNIT_PER_PEPTIDE",
    "DATASET_PEPTIDE_MAPPING_WEIGHT_NORMALIZATION_POLICY_SUM_TO_ONE_PER_PEPTIDE_ROW",
    "DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_DERIVED_EQUAL",
    "DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_EXPLICIT",
    "DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_POLICY_EXPLICIT_OR_DERIVED_EQUAL",
    "DATASET_PEPTIDE_MIXED_AMBIGUITY_POLICY_COMBINE_ALLOCATED_SIGNALS",
    "DATASET_PEPTIDE_MIXED_AMBIGUITY_POLICY_SHARED_WEIGHTED_MEAN",
    "DATASET_PEPTIDE_SIGNAL_ALLOCATION_POLICY_MULTIPLY_BY_MAPPING_FRACTION",
    "DATASET_PEPTIDE_SITE_SEQUENCE_POLICY_VALIDATE_WITHOUT_REPAIR",
    "DATASET_PEPTIDE_SITE_SUMMARISATION_POLICY_ARITHMETIC_MEAN_OF_ALLOCATED_SIGNALS",
    "DATASET_PEPTIDE_TO_SITE_AGGREGATION_POLICY_LEGACY_ALIAS",
    "DATASET_PEPTIDE_TO_SITE_AGGREGATION_POLICY_MAPPING_WEIGHTED_MEAN",
    "DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE",
    "DATASET_SITE_RESOLUTION_MODE_SITE_LEVEL_RESOLVED",
    "MAPPING_FRACTION_COLUMN",
    "MAPPING_WEIGHT_SUM_TOLERANCE",
    "PeptideEvidenceResolutionInputMetrics",
    "PeptideEvidenceResolutionResult",
    "PeptideEvidenceResolutionSummary",
    "SITE_SEQUENCE_SOURCE_MISSING",
    "SITE_SEQUENCE_SOURCE_PEPTIDE_CONTEXT",
    "SITE_SEQUENCE_SOURCE_PROVIDED",
    "SUPPORTED_DATASET_MULTI_SITE_POLICIES",
    "SUPPORTED_DATASET_SITE_RESOLUTION_MODES",
    "build_multi_site_handling_config_for_dataset_policy",
    "normalise_peptide_evidence_resolution_summary_payload",
    "validate_dataset_multi_site_policy",
]
