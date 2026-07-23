"""Public enrichment workflow configuration models."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import cast

from phospy.errors.validation import ContractValidationError
from phospy.science.configs.enrichment import (
    ENRICHMENT_IDENTIFIER_KIND_DISPLAY_ID,
    ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
    ENRICHMENT_IDENTIFIER_KIND_PHOSPHOSITE,
    ENRICHMENT_IDENTIFIER_KIND_PROTEIN_ID,
    ENRICHMENT_IDENTIFIER_KIND_SITE_KEY,
    ENRICHMENT_METHOD_OVER_REPRESENTATION,
    ENRICHMENT_OUTSIDE_BACKGROUND_POLICY_DROP,
    ENRICHMENT_OUTSIDE_BACKGROUND_POLICY_ERROR,
    GENE_LEVEL_ENRICHMENT_IDENTIFIER_KINDS,
    MULTIPLE_TESTING_CORRECTION_BENJAMINI_HOCHBERG,
    MULTIPLE_TESTING_CORRECTION_BENJAMINI_YEKUTIELI,
    MULTIPLE_TESTING_CORRECTION_BONFERRONI,
    MULTIPLE_TESTING_CORRECTION_HOLM,
    MULTIPLE_TESTING_CORRECTION_NONE,
    PTM_LEVEL_ENRICHMENT_IDENTIFIER_KINDS,
    SUPPORTED_ENRICHMENT_IDENTIFIER_KINDS,
    SUPPORTED_ENRICHMENT_METHODS,
    SUPPORTED_ENRICHMENT_OUTSIDE_BACKGROUND_POLICIES,
    SUPPORTED_MULTIPLE_TESTING_CORRECTIONS,
    EnrichmentIdentifierKind,
    EnrichmentMethod,
    EnrichmentOutsideBackgroundPolicy,
    MultipleTestingCorrection,
)


@dataclass(frozen=True, slots=True)
class EnrichmentConfig:
    """Public configuration for native enrichment.

    Records the enrichment method, multiple-testing correction, and optional
    minimum and maximum set-size filters. It does not compute enrichment
    statistics, load resources, or derive a background universe.
    """

    method: EnrichmentMethod = ENRICHMENT_METHOD_OVER_REPRESENTATION
    multiple_testing_correction: MultipleTestingCorrection = (
        MULTIPLE_TESTING_CORRECTION_BENJAMINI_HOCHBERG
    )
    min_set_size: int | None = None
    max_set_size: int | None = None
    selected_outside_background_policy: EnrichmentOutsideBackgroundPolicy = (
        ENRICHMENT_OUTSIDE_BACKGROUND_POLICY_ERROR
    )
    set_member_outside_background_policy: EnrichmentOutsideBackgroundPolicy = (
        ENRICHMENT_OUTSIDE_BACKGROUND_POLICY_DROP
    )
    minimum_retained_foreground_fraction: float | None = None

    def __post_init__(self) -> None:
        if self.method not in SUPPORTED_ENRICHMENT_METHODS:
            supported = ", ".join(repr(value) for value in SUPPORTED_ENRICHMENT_METHODS)
            raise ContractValidationError(
                f"enrichment.method must be one of: {supported}"
            )
        if (
            self.multiple_testing_correction
            not in SUPPORTED_MULTIPLE_TESTING_CORRECTIONS
        ):
            supported = ", ".join(
                repr(value) for value in SUPPORTED_MULTIPLE_TESTING_CORRECTIONS
            )
            raise ContractValidationError(
                f"enrichment.multiple_testing_correction must be one of: {supported}"
            )
        object.__setattr__(
            self,
            "method",
            cast(EnrichmentMethod, self.method),
        )
        object.__setattr__(
            self,
            "multiple_testing_correction",
            cast(MultipleTestingCorrection, self.multiple_testing_correction),
        )
        min_set_size = _normalise_optional_set_size_threshold(
            self.min_set_size,
            field_name="enrichment.min_set_size",
        )
        max_set_size = _normalise_optional_set_size_threshold(
            self.max_set_size,
            field_name="enrichment.max_set_size",
        )
        if (
            min_set_size is not None
            and max_set_size is not None
            and min_set_size > max_set_size
        ):
            raise ContractValidationError(
                "enrichment.min_set_size must be less than or equal to "
                "enrichment.max_set_size"
            )
        object.__setattr__(self, "min_set_size", min_set_size)
        object.__setattr__(self, "max_set_size", max_set_size)
        object.__setattr__(
            self,
            "selected_outside_background_policy",
            _normalise_outside_background_policy(
                self.selected_outside_background_policy,
                field_name="enrichment.selected_outside_background_policy",
            ),
        )
        object.__setattr__(
            self,
            "set_member_outside_background_policy",
            _normalise_outside_background_policy(
                self.set_member_outside_background_policy,
                field_name="enrichment.set_member_outside_background_policy",
            ),
        )
        object.__setattr__(
            self,
            "minimum_retained_foreground_fraction",
            _normalise_optional_unit_interval(
                self.minimum_retained_foreground_fraction,
                field_name="enrichment.minimum_retained_foreground_fraction",
            ),
        )


def _normalise_optional_set_size_threshold(
    value: object | None,
    *,
    field_name: str,
) -> int | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, int):
        raise ContractValidationError(f"{field_name} must be an int or None")
    if value < 1:
        raise ContractValidationError(
            f"{field_name} must be greater than or equal to 1"
        )
    return value


def _normalise_outside_background_policy(
    value: object,
    *,
    field_name: str,
) -> EnrichmentOutsideBackgroundPolicy:
    if value not in SUPPORTED_ENRICHMENT_OUTSIDE_BACKGROUND_POLICIES:
        supported = ", ".join(
            repr(policy) for policy in SUPPORTED_ENRICHMENT_OUTSIDE_BACKGROUND_POLICIES
        )
        raise ContractValidationError(f"{field_name} must be one of: {supported}")
    return cast(EnrichmentOutsideBackgroundPolicy, value)


def _normalise_optional_unit_interval(
    value: object | None,
    *,
    field_name: str,
) -> float | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise ContractValidationError(f"{field_name} must be numeric or None")
    normalised = float(value)
    if not math.isfinite(normalised) or normalised < 0.0 or normalised > 1.0:
        raise ContractValidationError(f"{field_name} must be within [0.0, 1.0]")
    return normalised


__all__ = [
    "ENRICHMENT_IDENTIFIER_KIND_DISPLAY_ID",
    "ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL",
    "ENRICHMENT_IDENTIFIER_KIND_PHOSPHOSITE",
    "ENRICHMENT_IDENTIFIER_KIND_PROTEIN_ID",
    "ENRICHMENT_IDENTIFIER_KIND_SITE_KEY",
    "ENRICHMENT_METHOD_OVER_REPRESENTATION",
    "ENRICHMENT_OUTSIDE_BACKGROUND_POLICY_DROP",
    "ENRICHMENT_OUTSIDE_BACKGROUND_POLICY_ERROR",
    "GENE_LEVEL_ENRICHMENT_IDENTIFIER_KINDS",
    "MULTIPLE_TESTING_CORRECTION_BENJAMINI_HOCHBERG",
    "MULTIPLE_TESTING_CORRECTION_BENJAMINI_YEKUTIELI",
    "MULTIPLE_TESTING_CORRECTION_BONFERRONI",
    "MULTIPLE_TESTING_CORRECTION_HOLM",
    "MULTIPLE_TESTING_CORRECTION_NONE",
    "PTM_LEVEL_ENRICHMENT_IDENTIFIER_KINDS",
    "SUPPORTED_ENRICHMENT_IDENTIFIER_KINDS",
    "SUPPORTED_ENRICHMENT_METHODS",
    "SUPPORTED_ENRICHMENT_OUTSIDE_BACKGROUND_POLICIES",
    "SUPPORTED_MULTIPLE_TESTING_CORRECTIONS",
    "EnrichmentConfig",
    "EnrichmentIdentifierKind",
    "EnrichmentMethod",
    "EnrichmentOutsideBackgroundPolicy",
    "MultipleTestingCorrection",
]
