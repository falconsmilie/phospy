"""Public enrichment workflow configuration models."""

from __future__ import annotations

from dataclasses import dataclass
from typing import cast

from phospy.errors.validation import WorkflowValidationError
from phospy.science.enrichment.models import (
    ENRICHMENT_IDENTIFIER_KIND_DISPLAY_ID,
    ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
    ENRICHMENT_IDENTIFIER_KIND_PHOSPHOSITE,
    ENRICHMENT_IDENTIFIER_KIND_PROTEIN_ID,
    ENRICHMENT_IDENTIFIER_KIND_SITE_KEY,
    ENRICHMENT_METHOD_OVER_REPRESENTATION,
    GENE_LEVEL_ENRICHMENT_IDENTIFIER_KINDS,
    MULTIPLE_TESTING_CORRECTION_BENJAMINI_HOCHBERG,
    MULTIPLE_TESTING_CORRECTION_BENJAMINI_YEKUTIELI,
    MULTIPLE_TESTING_CORRECTION_BONFERRONI,
    MULTIPLE_TESTING_CORRECTION_HOLM,
    MULTIPLE_TESTING_CORRECTION_NONE,
    PTM_LEVEL_ENRICHMENT_IDENTIFIER_KINDS,
    SUPPORTED_ENRICHMENT_IDENTIFIER_KINDS,
    SUPPORTED_ENRICHMENT_METHODS,
    SUPPORTED_MULTIPLE_TESTING_CORRECTIONS,
    EnrichmentIdentifierKind,
    EnrichmentMethod,
    MultipleTestingCorrection,
)


@dataclass(frozen=True, slots=True)
class EnrichmentConfig:
    """Public configuration for native enrichment.

    The initial contract records method and multiple-testing intent only. It
    does not compute enrichment statistics, load resources, or derive a
    background universe.
    """

    method: EnrichmentMethod = ENRICHMENT_METHOD_OVER_REPRESENTATION
    multiple_testing_correction: MultipleTestingCorrection = (
        MULTIPLE_TESTING_CORRECTION_BENJAMINI_HOCHBERG
    )
    min_set_size: int | None = None
    max_set_size: int | None = None

    def __post_init__(self) -> None:
        if self.method not in SUPPORTED_ENRICHMENT_METHODS:
            supported = ", ".join(repr(value) for value in SUPPORTED_ENRICHMENT_METHODS)
            raise WorkflowValidationError(
                f"enrichment.method must be one of: {supported}"
            )
        if (
            self.multiple_testing_correction
            not in SUPPORTED_MULTIPLE_TESTING_CORRECTIONS
        ):
            supported = ", ".join(
                repr(value) for value in SUPPORTED_MULTIPLE_TESTING_CORRECTIONS
            )
            raise WorkflowValidationError(
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
            raise WorkflowValidationError(
                "enrichment.min_set_size must be less than or equal to "
                "enrichment.max_set_size"
            )
        object.__setattr__(self, "min_set_size", min_set_size)
        object.__setattr__(self, "max_set_size", max_set_size)


def _normalise_optional_set_size_threshold(
    value: object | None,
    *,
    field_name: str,
) -> int | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, int):
        raise WorkflowValidationError(f"{field_name} must be an int or None")
    if value < 1:
        raise WorkflowValidationError(
            f"{field_name} must be greater than or equal to 1"
        )
    return value


__all__ = [
    "ENRICHMENT_IDENTIFIER_KIND_DISPLAY_ID",
    "ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL",
    "ENRICHMENT_IDENTIFIER_KIND_PHOSPHOSITE",
    "ENRICHMENT_IDENTIFIER_KIND_PROTEIN_ID",
    "ENRICHMENT_IDENTIFIER_KIND_SITE_KEY",
    "ENRICHMENT_METHOD_OVER_REPRESENTATION",
    "GENE_LEVEL_ENRICHMENT_IDENTIFIER_KINDS",
    "MULTIPLE_TESTING_CORRECTION_BENJAMINI_HOCHBERG",
    "MULTIPLE_TESTING_CORRECTION_BENJAMINI_YEKUTIELI",
    "MULTIPLE_TESTING_CORRECTION_BONFERRONI",
    "MULTIPLE_TESTING_CORRECTION_HOLM",
    "MULTIPLE_TESTING_CORRECTION_NONE",
    "PTM_LEVEL_ENRICHMENT_IDENTIFIER_KINDS",
    "SUPPORTED_ENRICHMENT_IDENTIFIER_KINDS",
    "SUPPORTED_ENRICHMENT_METHODS",
    "SUPPORTED_MULTIPLE_TESTING_CORRECTIONS",
    "EnrichmentConfig",
    "EnrichmentIdentifierKind",
    "EnrichmentMethod",
    "MultipleTestingCorrection",
]
