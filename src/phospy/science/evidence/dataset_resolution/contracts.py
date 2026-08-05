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
from phospy.science.transformations.models import (
    IntensityScaleKind,
    QuantitativeMeaning,
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
DATASET_PEPTIDE_LOCALISATION_SUMMARY_POLICY_DESCRIPTIVE_MEAN = (
    "descriptive_mean_of_finite_reported_localisation_confidence_values"
)
DATASET_PEPTIDE_LOCALISATION_SUMMARY_SEMANTICS = (
    "descriptive_arithmetic_mean_not_calibrated_posterior_probability"
)
DATASET_PEPTIDE_LOCALISATION_SUMMARY_COLUMN = "localisation_confidence_descriptive_mean"
DATASET_PEPTIDE_LOCALISATION_COMPATIBILITY_ALIAS_COLUMN = "localisation_confidence"
DATASET_PEPTIDE_LOCALISATION_SUMMARY_SEMANTICS_COLUMN = (
    "localisation_confidence_summary_semantics"
)
DATASET_PEPTIDE_SITE_SEQUENCE_POLICY_VALIDATE_WITHOUT_REPAIR = "validate_without_repair"
DATASET_PEPTIDE_INPUT_QUANTITATIVE_MEANING_PEPTIDE_ABUNDANCE = "peptide_abundance"
DATASET_PEPTIDE_INPUT_QUANTITATIVE_MEANING_PEPTIDE_LOG2_ABUNDANCE = (
    "peptide_log2_abundance"
)
DATASET_PEPTIDE_ALLOCATION_DOMAIN_LINEAR_ABUNDANCE = "linear_abundance"
DATASET_PEPTIDE_ALLOCATION_DOMAIN_DECLARED_SCALE_UNIT_MAPPING_PASSTHROUGH = (
    "declared_scale_unit_mapping_passthrough"
)
DATASET_PEPTIDE_MAPPING_WEIGHT_SEMANTICS = (
    "unitless_fraction_of_one_peptide_evidence_row_allocated_to_each_resolved_site"
)
DATASET_PEPTIDE_MISSING_VALUE_POLICY_FINITE_MEAN = (
    "mean_finite_allocated_values_per_site_sample_preserve_missing_if_none_finite"
)
DATASET_PEPTIDE_SIGNAL_CONSERVATION_POLICY_NOT_CONSERVED = (
    "not_signal_conserving_after_per_site_arithmetic_mean"
)
DATASET_PEPTIDE_TO_SITE_AGGREGATION_POLICY_LINEAR_ALLOCATED_MEAN_V1 = (
    "peptide_to_site_linear_abundance_fractional_allocation_arithmetic_mean_v1"
)
DATASET_PEPTIDE_TO_SITE_UNCERTAINTY_LIMITATION_DESCRIPTIVE = (
    "no_model_based_uncertainty_or_posterior_localisation_combination"
)
DATASET_PEPTIDE_TO_SITE_UNCERTAINTY_LIMITATION_DEPENDENCE = (
    "peptide_evidence_rows_are_not_modelled_as_independent_replicates"
)

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


class PeptideEvidenceInputQuantitativeMeaning(str, Enum):
    """Supported quantitative meanings for peptide evidence entering site resolution."""

    PEPTIDE_ABUNDANCE = DATASET_PEPTIDE_INPUT_QUANTITATIVE_MEANING_PEPTIDE_ABUNDANCE
    PEPTIDE_LOG2_ABUNDANCE = (
        DATASET_PEPTIDE_INPUT_QUANTITATIVE_MEANING_PEPTIDE_LOG2_ABUNDANCE
    )


class PeptideToSiteAllocationDomain(str, Enum):
    """Numeric domain in which peptide-to-site allocation is defined."""

    LINEAR_ABUNDANCE = DATASET_PEPTIDE_ALLOCATION_DOMAIN_LINEAR_ABUNDANCE
    DECLARED_SCALE_UNIT_MAPPING_PASSTHROUGH = (
        DATASET_PEPTIDE_ALLOCATION_DOMAIN_DECLARED_SCALE_UNIT_MAPPING_PASSTHROUGH
    )


_SUPPORTED_INPUT_SCALES: frozenset[IntensityScaleKind] = frozenset(
    {IntensityScaleKind.LINEAR, IntensityScaleKind.LOG2}
)
_SUPPORTED_INPUT_QUANTITATIVE_MEANINGS: frozenset[
    PeptideEvidenceInputQuantitativeMeaning
] = frozenset(PeptideEvidenceInputQuantitativeMeaning)
_INPUT_MEANING_BY_SCALE: dict[
    IntensityScaleKind, PeptideEvidenceInputQuantitativeMeaning
] = {
    IntensityScaleKind.LINEAR: PeptideEvidenceInputQuantitativeMeaning.PEPTIDE_ABUNDANCE,
    IntensityScaleKind.LOG2: PeptideEvidenceInputQuantitativeMeaning.PEPTIDE_LOG2_ABUNDANCE,
}
_OUTPUT_MEANING_BY_SCALE: dict[IntensityScaleKind, QuantitativeMeaning] = {
    IntensityScaleKind.LINEAR: QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE,
    IntensityScaleKind.LOG2: QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE,
}
_POLICY_UNCERTAINTY_LIMITATIONS: tuple[str, ...] = (
    DATASET_PEPTIDE_TO_SITE_UNCERTAINTY_LIMITATION_DESCRIPTIVE,
    DATASET_PEPTIDE_TO_SITE_UNCERTAINTY_LIMITATION_DEPENDENCE,
)
_LEGACY_NOT_RECORDED = "not_recorded_legacy"


@dataclass(frozen=True, slots=True)
class PeptideToSiteAggregationPolicy:
    """Run-specific peptide-to-site quantitative estimand and allocation contract."""

    policy_id: str
    input_scale: IntensityScaleKind | str
    input_quantitative_meaning: PeptideEvidenceInputQuantitativeMeaning | str
    output_scale: IntensityScaleKind | str
    output_quantitative_meaning: QuantitativeMeaning | str
    allocation_domain: PeptideToSiteAllocationDomain | str
    fractional_mapping_present: bool
    supported_input_scales: frozenset[IntensityScaleKind] = _SUPPORTED_INPUT_SCALES
    supported_input_quantitative_meanings: frozenset[
        PeptideEvidenceInputQuantitativeMeaning
    ] = _SUPPORTED_INPUT_QUANTITATIVE_MEANINGS
    mapping_weight_source_policy: str = (
        DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_POLICY_EXPLICIT_OR_DERIVED_EQUAL
    )
    mapping_weight_normalization_policy: str = (
        DATASET_PEPTIDE_MAPPING_WEIGHT_NORMALIZATION_POLICY_SUM_TO_ONE_PER_PEPTIDE_ROW
    )
    mapping_weight_semantics: str = DATASET_PEPTIDE_MAPPING_WEIGHT_SEMANTICS
    signal_allocation_policy: str = (
        DATASET_PEPTIDE_SIGNAL_ALLOCATION_POLICY_MULTIPLY_BY_MAPPING_FRACTION
    )
    site_summarisation_policy: str = (
        DATASET_PEPTIDE_SITE_SUMMARISATION_POLICY_ARITHMETIC_MEAN_OF_ALLOCATED_SIGNALS
    )
    duplicate_evidence_policy: str = (
        DATASET_PEPTIDE_DUPLICATE_EVIDENCE_POLICY_RETAIN_DUPLICATE_ROWS
    )
    mixed_ambiguity_policy: str = (
        DATASET_PEPTIDE_MIXED_AMBIGUITY_POLICY_COMBINE_ALLOCATED_SIGNALS
    )
    missing_value_policy: str = DATASET_PEPTIDE_MISSING_VALUE_POLICY_FINITE_MEAN
    localisation_aggregation_policy: str = (
        DATASET_PEPTIDE_LOCALISATION_AGGREGATION_POLICY_ARITHMETIC_MEAN_OF_FINITE_VALUES
    )
    localisation_summary_policy: str = (
        DATASET_PEPTIDE_LOCALISATION_SUMMARY_POLICY_DESCRIPTIVE_MEAN
    )
    localisation_summary_semantics: str = DATASET_PEPTIDE_LOCALISATION_SUMMARY_SEMANTICS
    localisation_output_column: str = DATASET_PEPTIDE_LOCALISATION_SUMMARY_COLUMN
    localisation_compatibility_alias_column: str = (
        DATASET_PEPTIDE_LOCALISATION_COMPATIBILITY_ALIAS_COLUMN
    )
    signal_conservation_policy: str = (
        DATASET_PEPTIDE_SIGNAL_CONSERVATION_POLICY_NOT_CONSERVED
    )
    uncertainty_limitations: tuple[str, ...] = _POLICY_UNCERTAINTY_LIMITATIONS

    def __post_init__(self) -> None:
        input_scale = _normalize_intensity_scale_kind(
            self.input_scale,
            field_name="peptide_to_site_aggregation_policy.input_scale",
        )
        output_scale = _normalize_intensity_scale_kind(
            self.output_scale,
            field_name="peptide_to_site_aggregation_policy.output_scale",
        )
        input_meaning = _normalize_input_quantitative_meaning(
            self.input_quantitative_meaning
        )
        output_meaning = _normalize_output_quantitative_meaning(
            self.output_quantitative_meaning
        )
        allocation_domain = _normalize_allocation_domain(self.allocation_domain)
        supported_scales = frozenset(
            _normalize_intensity_scale_kind(
                value,
                field_name="peptide_to_site_aggregation_policy.supported_input_scales",
            )
            for value in self.supported_input_scales
        )
        supported_meanings = frozenset(
            _normalize_input_quantitative_meaning(value)
            for value in self.supported_input_quantitative_meanings
        )
        fractional_mapping_present = _require_bool(
            self.fractional_mapping_present,
            field_name="peptide_to_site_aggregation_policy.fractional_mapping_present",
        )
        policy_id = _canonical_non_empty_string(
            self.policy_id,
            field_name="peptide_to_site_aggregation_policy.policy_id",
        )
        if input_scale not in supported_scales:
            supported = ", ".join(sorted(scale.value for scale in supported_scales))
            raise PhosPyInputError(
                "peptide-to-site aggregation policy does not support input scale "
                f"{input_scale.value!r}; supported_input_scales=[{supported}]"
            )
        if input_meaning not in supported_meanings:
            supported = ", ".join(
                sorted(meaning.value for meaning in supported_meanings)
            )
            raise PhosPyInputError(
                "peptide-to-site aggregation policy does not support input "
                f"quantitative meaning {input_meaning.value!r}; "
                f"supported_input_quantitative_meanings=[{supported}]"
            )
        expected_input_meaning = _INPUT_MEANING_BY_SCALE[input_scale]
        if input_meaning is not expected_input_meaning:
            raise PhosPyInputError(
                "peptide-to-site aggregation policy input scale and quantitative "
                "meaning are incoherent: "
                f"input_intensity_scale={input_scale.value!r}, "
                f"input_quantitative_meaning={input_meaning.value!r}; expected "
                f"{expected_input_meaning.value!r}"
            )
        if output_scale is not input_scale:
            raise PhosPyInputError(
                "the current peptide-to-site aggregation policy preserves the "
                "declared input scale; output_scale must equal input_scale"
            )
        expected_output_meaning = _OUTPUT_MEANING_BY_SCALE[output_scale]
        if output_meaning is not expected_output_meaning:
            raise PhosPyInputError(
                "peptide-to-site aggregation policy output scale and quantitative "
                "meaning are incoherent: "
                f"output_intensity_scale={output_scale.value!r}, "
                f"output_quantitative_meaning={output_meaning.value!r}; expected "
                f"{expected_output_meaning.value!r}"
            )
        if input_scale is IntensityScaleKind.LINEAR:
            if allocation_domain is not PeptideToSiteAllocationDomain.LINEAR_ABUNDANCE:
                raise PhosPyInputError(
                    "linear peptide evidence must use allocation_domain="
                    f"{PeptideToSiteAllocationDomain.LINEAR_ABUNDANCE.value!r}"
                )
        elif allocation_domain is not (
            PeptideToSiteAllocationDomain.DECLARED_SCALE_UNIT_MAPPING_PASSTHROUGH
        ):
            raise PhosPyInputError(
                "log2 peptide evidence can only use allocation_domain="
                f"{PeptideToSiteAllocationDomain.DECLARED_SCALE_UNIT_MAPPING_PASSTHROUGH.value!r} "
                "for unit peptide-to-site mappings"
            )
        if (
            fractional_mapping_present
            and allocation_domain is not PeptideToSiteAllocationDomain.LINEAR_ABUNDANCE
        ):
            raise PhosPyInputError(
                "fractional peptide-to-site mapping is defined only in the "
                "linear_abundance allocation domain"
            )
        object.__setattr__(self, "policy_id", policy_id)
        object.__setattr__(self, "input_scale", input_scale)
        object.__setattr__(self, "output_scale", output_scale)
        object.__setattr__(self, "input_quantitative_meaning", input_meaning)
        object.__setattr__(self, "output_quantitative_meaning", output_meaning)
        object.__setattr__(self, "allocation_domain", allocation_domain)
        object.__setattr__(self, "supported_input_scales", supported_scales)
        object.__setattr__(
            self, "supported_input_quantitative_meanings", supported_meanings
        )
        object.__setattr__(
            self, "fractional_mapping_present", fractional_mapping_present
        )
        object.__setattr__(
            self,
            "uncertainty_limitations",
            tuple(
                _canonical_non_empty_string(
                    item,
                    field_name="peptide_to_site_aggregation_policy.uncertainty_limitations",
                )
                for item in self.uncertainty_limitations
            ),
        )

    @property
    def aggregation_formula(self) -> str:
        """Return the run-specific formula with explicit units and scale."""

        if self.input_scale is IntensityScaleKind.LOG2:
            return (
                "unit mapping only: a[p,s,j] [log2 peptide-abundance units] = "
                "x[p,j] [log2 peptide-abundance units] because w[p,s]=1.0 "
                "[unitless allocation fraction]; y[s,j] [log2 phosphosite-"
                "abundance units] = arithmetic_mean(a[p,s,j] over finite retained "
                "peptide evidence rows p mapped to site s for sample j)"
            )
        return (
            "a[p,s,j] [linear abundance units] = w[p,s] [unitless allocation "
            "fraction] * x[p,j] [linear peptide-abundance units]; y[s,j] "
            "[linear phosphosite-abundance estimate units] = arithmetic_mean("
            "a[p,s,j] over finite retained peptide evidence rows p mapped to "
            "site s for sample j)"
        )

    def validate_allocation_rows(self, rows: pd.DataFrame) -> None:
        """Validate that mapped rows are eligible for this typed allocation policy."""

        if not isinstance(rows, pd.DataFrame):
            raise PhosPyInputError(
                "typed peptide-to-site aggregation policy validation requires a "
                "pandas DataFrame of resolved mapping rows"
            )
        mapping_fractions = _mapping_fractions(rows)
        fractional_mask = _fractional_mapping_mask(mapping_fractions)
        observed_fractional_mapping = bool(fractional_mask.any())
        if observed_fractional_mapping != self.fractional_mapping_present:
            raise PhosPyInputError(
                "typed peptide-to-site aggregation policy fractional_mapping_present "
                "does not match resolved mapping rows"
            )
        if (
            observed_fractional_mapping
            and self.input_scale is not IntensityScaleKind.LINEAR
        ):
            _raise_fractional_allocation_for_non_linear_input(
                rows=rows,
                input_intensity_scale=cast(IntensityScaleKind, self.input_scale),
                fractional_mask=fractional_mask,
            )

    def to_payload(self) -> dict[str, object]:
        """Return JSON-safe policy and estimand provenance for this run."""

        input_scale = cast(IntensityScaleKind, self.input_scale)
        output_scale = cast(IntensityScaleKind, self.output_scale)
        input_meaning = cast(
            PeptideEvidenceInputQuantitativeMeaning,
            self.input_quantitative_meaning,
        )
        output_meaning = cast(QuantitativeMeaning, self.output_quantitative_meaning)
        allocation_domain = cast(PeptideToSiteAllocationDomain, self.allocation_domain)
        return {
            "peptide_to_site_aggregation_policy_id": self.policy_id,
            "supported_input_scales": [
                scale.value for scale in sorted(self.supported_input_scales, key=str)
            ],
            "supported_input_quantitative_meanings": [
                meaning.value
                for meaning in sorted(
                    self.supported_input_quantitative_meanings,
                    key=str,
                )
            ],
            "input_intensity_scale": input_scale.value,
            "input_quantitative_meaning": input_meaning.value,
            "output_intensity_scale": output_scale.value,
            "output_quantitative_meaning": output_meaning.value,
            "allocation_domain": allocation_domain.value,
            "fractional_mapping_present": bool(self.fractional_mapping_present),
            "mapping_weight_source_policy": self.mapping_weight_source_policy,
            "mapping_weight_normalization_policy": (
                self.mapping_weight_normalization_policy
            ),
            "mapping_weight_semantics": self.mapping_weight_semantics,
            "signal_allocation_policy": self.signal_allocation_policy,
            "site_summarisation_policy": self.site_summarisation_policy,
            "duplicate_evidence_policy": self.duplicate_evidence_policy,
            "mixed_ambiguity_policy": self.mixed_ambiguity_policy,
            "missing_value_policy": self.missing_value_policy,
            "localisation_aggregation_policy": self.localisation_aggregation_policy,
            "localisation_summary_policy": self.localisation_summary_policy,
            "localisation_summary_semantics": self.localisation_summary_semantics,
            "localisation_output_column": self.localisation_output_column,
            "localisation_compatibility_alias_column": (
                self.localisation_compatibility_alias_column
            ),
            "signal_conservation_policy": self.signal_conservation_policy,
            "uncertainty_limitations": list(self.uncertainty_limitations),
            "aggregation_policy": DATASET_PEPTIDE_TO_SITE_AGGREGATION_POLICY_LEGACY_ALIAS,
            "aggregation_formula": self.aggregation_formula,
        }


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


def build_peptide_to_site_aggregation_policy(
    *,
    input_intensity_scale: IntensityScaleKind | str,
    input_quantitative_meaning: PeptideEvidenceInputQuantitativeMeaning
    | str
    | None = None,
    mapping_rows: pd.DataFrame,
) -> PeptideToSiteAggregationPolicy:
    """Build the current typed sample-intensity peptide-to-site estimand contract."""

    input_scale = _normalize_intensity_scale_kind(
        input_intensity_scale,
        field_name="peptide_to_site_aggregation_policy.input_intensity_scale",
    )
    resolved_input_meaning = (
        _INPUT_MEANING_BY_SCALE[input_scale]
        if input_quantitative_meaning is None
        else _normalize_input_quantitative_meaning(input_quantitative_meaning)
    )
    mapping_fractions = _mapping_fractions(mapping_rows)
    fractional_mask = _fractional_mapping_mask(mapping_fractions)
    fractional_mapping_present = bool(fractional_mask.any())
    if fractional_mapping_present and input_scale is not IntensityScaleKind.LINEAR:
        _raise_fractional_allocation_for_non_linear_input(
            rows=mapping_rows,
            input_intensity_scale=input_scale,
            fractional_mask=fractional_mask,
        )
    allocation_domain = (
        PeptideToSiteAllocationDomain.LINEAR_ABUNDANCE
        if input_scale is IntensityScaleKind.LINEAR
        else PeptideToSiteAllocationDomain.DECLARED_SCALE_UNIT_MAPPING_PASSTHROUGH
    )
    return PeptideToSiteAggregationPolicy(
        policy_id=DATASET_PEPTIDE_TO_SITE_AGGREGATION_POLICY_LINEAR_ALLOCATED_MEAN_V1,
        input_scale=input_scale,
        input_quantitative_meaning=resolved_input_meaning,
        output_scale=input_scale,
        output_quantitative_meaning=_OUTPUT_MEANING_BY_SCALE[input_scale],
        allocation_domain=allocation_domain,
        fractional_mapping_present=fractional_mapping_present,
    )


@dataclass(frozen=True, slots=True)
class PeptideEvidenceResolutionSummary:
    """Structured summary for peptide-to-site resolution provenance."""

    input_mode: str
    multi_site_policy: str | None
    peptide_to_site_aggregation_policy_id: str
    supported_input_scales: tuple[str, ...]
    supported_input_quantitative_meanings: tuple[str, ...]
    input_intensity_scale: str
    input_quantitative_meaning: str
    output_intensity_scale: str
    output_quantitative_meaning: str
    allocation_domain: str
    fractional_mapping_present: bool
    peptide_observations_received: int
    mapped_peptide_observations: int
    site_mapping_rows: int
    allocated_evidence_rows: int
    unique_site_ids_produced: int
    ambiguous_observations: int
    unambiguous_observations: int
    excluded_observations: int
    split_observations: int
    fractional_mapping_rows: int
    unit_mapping_rows: int
    mapping_weight_source_policy: str
    mapping_weight_normalization_policy: str
    mapping_weight_semantics: str
    signal_allocation_policy: str
    site_summarisation_policy: str
    missing_value_policy: str
    duplicate_evidence_policy: str
    mixed_ambiguity_policy: str
    localisation_aggregation_policy: str
    localisation_summary_policy: str
    localisation_summary_semantics: str
    localisation_output_column: str
    localisation_compatibility_alias_column: str
    signal_conservation_policy: str
    uncertainty_limitations: tuple[str, ...]
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
            "peptide_to_site_aggregation_policy_id": (
                self.peptide_to_site_aggregation_policy_id
            ),
            "supported_input_scales": list(self.supported_input_scales),
            "supported_input_quantitative_meanings": list(
                self.supported_input_quantitative_meanings
            ),
            "input_intensity_scale": self.input_intensity_scale,
            "input_quantitative_meaning": self.input_quantitative_meaning,
            "output_intensity_scale": self.output_intensity_scale,
            "output_quantitative_meaning": self.output_quantitative_meaning,
            "allocation_domain": self.allocation_domain,
            "fractional_mapping_present": bool(self.fractional_mapping_present),
            "peptide_observations_received": int(self.peptide_observations_received),
            "mapped_peptide_observations": int(self.mapped_peptide_observations),
            "site_mapping_rows": int(self.site_mapping_rows),
            "allocated_evidence_rows": int(self.allocated_evidence_rows),
            "unique_site_ids_produced": int(self.unique_site_ids_produced),
            "ambiguous_observations": int(self.ambiguous_observations),
            "unambiguous_observations": int(self.unambiguous_observations),
            "excluded_observations": int(self.excluded_observations),
            "split_observations": int(self.split_observations),
            "fractional_mapping_rows": int(self.fractional_mapping_rows),
            "unit_mapping_rows": int(self.unit_mapping_rows),
            "mapping_weight_source_policy": self.mapping_weight_source_policy,
            "mapping_weight_normalization_policy": (
                self.mapping_weight_normalization_policy
            ),
            "mapping_weight_semantics": self.mapping_weight_semantics,
            "signal_allocation_policy": self.signal_allocation_policy,
            "site_summarisation_policy": self.site_summarisation_policy,
            "missing_value_policy": self.missing_value_policy,
            "duplicate_evidence_policy": self.duplicate_evidence_policy,
            "mixed_ambiguity_policy": self.mixed_ambiguity_policy,
            "localisation_aggregation_policy": self.localisation_aggregation_policy,
            "localisation_summary_policy": self.localisation_summary_policy,
            "localisation_summary_semantics": self.localisation_summary_semantics,
            "localisation_output_column": self.localisation_output_column,
            "localisation_compatibility_alias_column": (
                self.localisation_compatibility_alias_column
            ),
            "signal_conservation_policy": self.signal_conservation_policy,
            "uncertainty_limitations": list(self.uncertainty_limitations),
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
            peptide_to_site_aggregation_policy_id=str(
                normalized["peptide_to_site_aggregation_policy_id"]
            ),
            supported_input_scales=_payload_str_tuple(
                normalized["supported_input_scales"]
            ),
            supported_input_quantitative_meanings=_payload_str_tuple(
                normalized["supported_input_quantitative_meanings"]
            ),
            input_intensity_scale=str(normalized["input_intensity_scale"]),
            input_quantitative_meaning=str(normalized["input_quantitative_meaning"]),
            output_intensity_scale=str(normalized["output_intensity_scale"]),
            output_quantitative_meaning=str(normalized["output_quantitative_meaning"]),
            allocation_domain=str(normalized["allocation_domain"]),
            fractional_mapping_present=bool(normalized["fractional_mapping_present"]),
            peptide_observations_received=_payload_int(
                normalized["peptide_observations_received"]
            ),
            mapped_peptide_observations=_payload_int(
                normalized["mapped_peptide_observations"]
            ),
            site_mapping_rows=_payload_int(normalized["site_mapping_rows"]),
            allocated_evidence_rows=_payload_int(normalized["allocated_evidence_rows"]),
            unique_site_ids_produced=_payload_int(
                normalized["unique_site_ids_produced"]
            ),
            ambiguous_observations=_payload_int(normalized["ambiguous_observations"]),
            unambiguous_observations=_payload_int(
                normalized["unambiguous_observations"]
            ),
            excluded_observations=_payload_int(normalized["excluded_observations"]),
            split_observations=_payload_int(normalized["split_observations"]),
            fractional_mapping_rows=_payload_int(normalized["fractional_mapping_rows"]),
            unit_mapping_rows=_payload_int(normalized["unit_mapping_rows"]),
            mapping_weight_source_policy=str(
                normalized["mapping_weight_source_policy"]
            ),
            mapping_weight_normalization_policy=str(
                normalized["mapping_weight_normalization_policy"]
            ),
            mapping_weight_semantics=str(normalized["mapping_weight_semantics"]),
            signal_allocation_policy=str(normalized["signal_allocation_policy"]),
            site_summarisation_policy=str(normalized["site_summarisation_policy"]),
            missing_value_policy=str(normalized["missing_value_policy"]),
            duplicate_evidence_policy=str(normalized["duplicate_evidence_policy"]),
            mixed_ambiguity_policy=str(normalized["mixed_ambiguity_policy"]),
            localisation_aggregation_policy=str(
                normalized["localisation_aggregation_policy"]
            ),
            localisation_summary_policy=str(normalized["localisation_summary_policy"]),
            localisation_summary_semantics=str(
                normalized["localisation_summary_semantics"]
            ),
            localisation_output_column=str(normalized["localisation_output_column"]),
            localisation_compatibility_alias_column=str(
                normalized["localisation_compatibility_alias_column"]
            ),
            signal_conservation_policy=str(normalized["signal_conservation_policy"]),
            uncertainty_limitations=_payload_str_tuple(
                normalized["uncertainty_limitations"]
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


def _payload_str_tuple(value: object) -> tuple[str, ...]:
    if isinstance(value, str):
        return (value,)
    if isinstance(value, tuple | list):
        return tuple(str(item) for item in value)
    return tuple(str(item) for item in cast(Any, value))


def normalise_peptide_evidence_resolution_summary_payload(
    payload: Mapping[str, object],
) -> dict[str, object]:
    """Return a current-policy peptide-resolution payload from old or new data."""

    normalized = dict(payload)
    current = CURRENT_RESOLUTION_POLICIES.to_payload()
    legacy_policy_defaults = _legacy_policy_payload_defaults()

    for key, value in legacy_policy_defaults.items():
        normalized.setdefault(key, value)

    if "input_scale" in normalized and "input_intensity_scale" not in normalized:
        normalized["input_intensity_scale"] = normalized["input_scale"]
    if "output_scale" in normalized and "output_intensity_scale" not in normalized:
        normalized["output_intensity_scale"] = normalized["output_scale"]

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
    normalized["mapping_weight_semantics"] = str(
        normalized.get("mapping_weight_semantics")
        or DATASET_PEPTIDE_MAPPING_WEIGHT_SEMANTICS
    )
    normalized["missing_value_policy"] = str(
        normalized.get("missing_value_policy")
        or DATASET_PEPTIDE_MISSING_VALUE_POLICY_FINITE_MEAN
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
    normalized["localisation_summary_policy"] = str(
        normalized.get("localisation_summary_policy")
        or DATASET_PEPTIDE_LOCALISATION_SUMMARY_POLICY_DESCRIPTIVE_MEAN
    )
    normalized["localisation_summary_semantics"] = str(
        normalized.get("localisation_summary_semantics")
        or DATASET_PEPTIDE_LOCALISATION_SUMMARY_SEMANTICS
    )
    normalized["localisation_output_column"] = str(
        normalized.get("localisation_output_column")
        or DATASET_PEPTIDE_LOCALISATION_SUMMARY_COLUMN
    )
    normalized["localisation_compatibility_alias_column"] = str(
        normalized.get("localisation_compatibility_alias_column")
        or DATASET_PEPTIDE_LOCALISATION_COMPATIBILITY_ALIAS_COLUMN
    )
    normalized["signal_conservation_policy"] = str(
        normalized.get("signal_conservation_policy")
        or DATASET_PEPTIDE_SIGNAL_CONSERVATION_POLICY_NOT_CONSERVED
    )
    normalized["uncertainty_limitations"] = list(
        _payload_str_tuple(
            normalized.get(
                "uncertainty_limitations",
                _POLICY_UNCERTAINTY_LIMITATIONS,
            )
        )
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
    unambiguous_observations = normalized.get("unambiguous_observations")
    normalized["unambiguous_observations"] = (
        _payload_int(unambiguous_observations)
        if unambiguous_observations is not None
        else (
            _payload_int(normalized["peptide_observations_received"])
            - _payload_int(normalized["ambiguous_observations"])
        )
    )
    mapped_peptide_observations = normalized.get("mapped_peptide_observations")
    normalized["mapped_peptide_observations"] = (
        _payload_int(mapped_peptide_observations)
        if mapped_peptide_observations is not None
        else _payload_int(normalized["peptide_observations_received"])
    )
    site_mapping_rows = normalized.get("site_mapping_rows")
    normalized["site_mapping_rows"] = (
        _payload_int(site_mapping_rows)
        if site_mapping_rows is not None
        else _payload_int(normalized["unique_site_ids_produced"])
    )
    allocated_evidence_rows = normalized.get("allocated_evidence_rows")
    normalized["allocated_evidence_rows"] = (
        _payload_int(allocated_evidence_rows)
        if allocated_evidence_rows is not None
        else _payload_int(normalized["site_mapping_rows"])
    )
    fractional_mapping_rows = normalized.get("fractional_mapping_rows")
    normalized["fractional_mapping_rows"] = (
        _payload_int(fractional_mapping_rows)
        if fractional_mapping_rows is not None
        else 0
    )
    unit_mapping_rows = normalized.get("unit_mapping_rows")
    normalized["unit_mapping_rows"] = (
        _payload_int(unit_mapping_rows)
        if unit_mapping_rows is not None
        else (
            _payload_int(normalized["allocated_evidence_rows"])
            - _payload_int(normalized["fractional_mapping_rows"])
        )
    )
    return normalized


def _legacy_policy_payload_defaults() -> dict[str, object]:
    return {
        "peptide_to_site_aggregation_policy_id": (
            DATASET_PEPTIDE_TO_SITE_AGGREGATION_POLICY_LINEAR_ALLOCATED_MEAN_V1
        ),
        "supported_input_scales": [
            scale.value for scale in sorted(_SUPPORTED_INPUT_SCALES, key=str)
        ],
        "supported_input_quantitative_meanings": [
            meaning.value
            for meaning in sorted(_SUPPORTED_INPUT_QUANTITATIVE_MEANINGS, key=str)
        ],
        "input_intensity_scale": _LEGACY_NOT_RECORDED,
        "input_quantitative_meaning": _LEGACY_NOT_RECORDED,
        "output_intensity_scale": _LEGACY_NOT_RECORDED,
        "output_quantitative_meaning": _LEGACY_NOT_RECORDED,
        "allocation_domain": _LEGACY_NOT_RECORDED,
        "fractional_mapping_present": False,
        "mapping_weight_semantics": DATASET_PEPTIDE_MAPPING_WEIGHT_SEMANTICS,
        "missing_value_policy": DATASET_PEPTIDE_MISSING_VALUE_POLICY_FINITE_MEAN,
        "localisation_summary_policy": (
            DATASET_PEPTIDE_LOCALISATION_SUMMARY_POLICY_DESCRIPTIVE_MEAN
        ),
        "localisation_summary_semantics": (
            DATASET_PEPTIDE_LOCALISATION_SUMMARY_SEMANTICS
        ),
        "localisation_output_column": DATASET_PEPTIDE_LOCALISATION_SUMMARY_COLUMN,
        "localisation_compatibility_alias_column": (
            DATASET_PEPTIDE_LOCALISATION_COMPATIBILITY_ALIAS_COLUMN
        ),
        "signal_conservation_policy": (
            DATASET_PEPTIDE_SIGNAL_CONSERVATION_POLICY_NOT_CONSERVED
        ),
        "uncertainty_limitations": list(_POLICY_UNCERTAINTY_LIMITATIONS),
    }


def _mapping_fractions(rows: pd.DataFrame) -> pd.Series:
    if MAPPING_FRACTION_COLUMN not in rows.columns:
        raise PhosPyInputError(
            "peptide-to-site allocation requires resolved mapping rows with "
            f"{MAPPING_FRACTION_COLUMN!r}; build a typed aggregation policy after "
            "mapping-fraction resolution"
        )
    return pd.to_numeric(rows.loc[:, MAPPING_FRACTION_COLUMN], errors="coerce")


def _fractional_mapping_mask(mapping_fractions: pd.Series) -> pd.Series:
    return mapping_fractions.sub(1.0).abs() > MAPPING_WEIGHT_SUM_TOLERANCE


def _raise_fractional_allocation_for_non_linear_input(
    *,
    rows: pd.DataFrame,
    input_intensity_scale: IntensityScaleKind,
    fractional_mask: pd.Series,
) -> None:
    fractional_rows = rows.loc[fractional_mask, ["peptide_row_id", "site_id"]].head(5)
    preview = ", ".join(
        f"{str(row.peptide_row_id)!r}->{str(row.site_id)!r}"
        for row in fractional_rows.itertuples(index=False)
    )
    suffix = "" if int(fractional_mask.sum()) <= 5 else " ..."
    raise PhosPyInputError(
        "dataset peptide-evidence mode declared input_intensity_scale="
        f"{input_intensity_scale.value!r} cannot use fractional allocation: "
        "mapping_fraction contains non-unit values for peptide-to-site mappings "
        f"({preview}{suffix}). Fractional allocation is defined only for linear "
        "peptide abundance. PhosPy does not invert log2 peptide evidence with "
        "2**x at this boundary because complete pseudocount and transformation "
        "provenance is not available. Supported corrective action: provide "
        "linear peptide evidence before applying split/fractional allocation, "
        "or use only unit/unambiguous peptide-to-site mappings until a validated "
        "scale-aware peptide allocation estimator is implemented."
    )


def _normalize_intensity_scale_kind(
    value: IntensityScaleKind | str,
    *,
    field_name: str,
) -> IntensityScaleKind:
    if isinstance(value, IntensityScaleKind):
        return value
    try:
        return IntensityScaleKind(str(value))
    except ValueError as exc:
        supported = ", ".join(member.value for member in IntensityScaleKind)
        raise PhosPyInputError(f"{field_name} must be one of: {supported}") from exc


def _normalize_input_quantitative_meaning(
    value: PeptideEvidenceInputQuantitativeMeaning | str,
) -> PeptideEvidenceInputQuantitativeMeaning:
    if isinstance(value, PeptideEvidenceInputQuantitativeMeaning):
        return value
    try:
        return PeptideEvidenceInputQuantitativeMeaning(str(value))
    except ValueError as exc:
        supported = ", ".join(
            member.value for member in PeptideEvidenceInputQuantitativeMeaning
        )
        raise PhosPyInputError(
            "peptide-to-site aggregation policy input_quantitative_meaning must "
            f"be one of: {supported}"
        ) from exc


def _normalize_output_quantitative_meaning(
    value: QuantitativeMeaning | str,
) -> QuantitativeMeaning:
    if isinstance(value, QuantitativeMeaning):
        return value
    try:
        return QuantitativeMeaning(str(value))
    except ValueError as exc:
        supported = ", ".join(member.value for member in QuantitativeMeaning)
        raise PhosPyInputError(
            "peptide-to-site aggregation policy output_quantitative_meaning must "
            f"be one of: {supported}"
        ) from exc


def _normalize_allocation_domain(
    value: PeptideToSiteAllocationDomain | str,
) -> PeptideToSiteAllocationDomain:
    if isinstance(value, PeptideToSiteAllocationDomain):
        return value
    try:
        return PeptideToSiteAllocationDomain(str(value))
    except ValueError as exc:
        supported = ", ".join(member.value for member in PeptideToSiteAllocationDomain)
        raise PhosPyInputError(
            "peptide-to-site aggregation policy allocation_domain must be one of: "
            f"{supported}"
        ) from exc


def _canonical_non_empty_string(value: object, *, field_name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise PhosPyInputError(f"{field_name} must be a non-empty string")
    return value.strip()


def _require_bool(value: object, *, field_name: str) -> bool:
    if not isinstance(value, bool):
        raise PhosPyInputError(f"{field_name} must be a bool")
    return value


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
    unambiguous_observations: int
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
    "DATASET_PEPTIDE_ALLOCATION_DOMAIN_DECLARED_SCALE_UNIT_MAPPING_PASSTHROUGH",
    "DATASET_PEPTIDE_ALLOCATION_DOMAIN_LINEAR_ABUNDANCE",
    "DATASET_PEPTIDE_INPUT_QUANTITATIVE_MEANING_PEPTIDE_ABUNDANCE",
    "DATASET_PEPTIDE_INPUT_QUANTITATIVE_MEANING_PEPTIDE_LOG2_ABUNDANCE",
    "DATASET_PEPTIDE_LOCALISATION_AGGREGATION_POLICY_ARITHMETIC_MEAN_OF_FINITE_VALUES",
    "DATASET_PEPTIDE_LOCALISATION_COMPATIBILITY_ALIAS_COLUMN",
    "DATASET_PEPTIDE_LOCALISATION_SUMMARY_COLUMN",
    "DATASET_PEPTIDE_LOCALISATION_SUMMARY_POLICY_DESCRIPTIVE_MEAN",
    "DATASET_PEPTIDE_LOCALISATION_SUMMARY_SEMANTICS",
    "DATASET_PEPTIDE_LOCALISATION_SUMMARY_SEMANTICS_COLUMN",
    "DATASET_PEPTIDE_MAPPING_WEIGHT_NORMALISATION_UNIT_PER_PEPTIDE",
    "DATASET_PEPTIDE_MAPPING_WEIGHT_NORMALIZATION_POLICY_SUM_TO_ONE_PER_PEPTIDE_ROW",
    "DATASET_PEPTIDE_MAPPING_WEIGHT_SEMANTICS",
    "DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_DERIVED_EQUAL",
    "DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_EXPLICIT",
    "DATASET_PEPTIDE_MAPPING_WEIGHT_SOURCE_POLICY_EXPLICIT_OR_DERIVED_EQUAL",
    "DATASET_PEPTIDE_MISSING_VALUE_POLICY_FINITE_MEAN",
    "DATASET_PEPTIDE_MIXED_AMBIGUITY_POLICY_COMBINE_ALLOCATED_SIGNALS",
    "DATASET_PEPTIDE_MIXED_AMBIGUITY_POLICY_SHARED_WEIGHTED_MEAN",
    "DATASET_PEPTIDE_SIGNAL_ALLOCATION_POLICY_MULTIPLY_BY_MAPPING_FRACTION",
    "DATASET_PEPTIDE_SIGNAL_CONSERVATION_POLICY_NOT_CONSERVED",
    "DATASET_PEPTIDE_SITE_SEQUENCE_POLICY_VALIDATE_WITHOUT_REPAIR",
    "DATASET_PEPTIDE_SITE_SUMMARISATION_POLICY_ARITHMETIC_MEAN_OF_ALLOCATED_SIGNALS",
    "DATASET_PEPTIDE_TO_SITE_AGGREGATION_POLICY_LINEAR_ALLOCATED_MEAN_V1",
    "DATASET_PEPTIDE_TO_SITE_AGGREGATION_POLICY_LEGACY_ALIAS",
    "DATASET_PEPTIDE_TO_SITE_AGGREGATION_POLICY_MAPPING_WEIGHTED_MEAN",
    "DATASET_PEPTIDE_TO_SITE_UNCERTAINTY_LIMITATION_DEPENDENCE",
    "DATASET_PEPTIDE_TO_SITE_UNCERTAINTY_LIMITATION_DESCRIPTIVE",
    "DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE",
    "DATASET_SITE_RESOLUTION_MODE_SITE_LEVEL_RESOLVED",
    "MAPPING_FRACTION_COLUMN",
    "MAPPING_WEIGHT_SUM_TOLERANCE",
    "PeptideEvidenceInputQuantitativeMeaning",
    "PeptideEvidenceResolutionInputMetrics",
    "PeptideEvidenceResolutionResult",
    "PeptideEvidenceResolutionSummary",
    "PeptideToSiteAggregationPolicy",
    "PeptideToSiteAllocationDomain",
    "SITE_SEQUENCE_SOURCE_MISSING",
    "SITE_SEQUENCE_SOURCE_PEPTIDE_CONTEXT",
    "SITE_SEQUENCE_SOURCE_PROVIDED",
    "SUPPORTED_DATASET_MULTI_SITE_POLICIES",
    "SUPPORTED_DATASET_SITE_RESOLUTION_MODES",
    "build_multi_site_handling_config_for_dataset_policy",
    "build_peptide_to_site_aggregation_policy",
    "normalise_peptide_evidence_resolution_summary_payload",
    "validate_dataset_multi_site_policy",
]
