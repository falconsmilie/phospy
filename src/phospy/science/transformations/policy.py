"""Transformation-domain quantitative scale policy definitions."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Final, cast

from phospy.errors.transformations import InvalidTransformationStateError

IDENTITY_INTENSITY_SCALE_ESTABLISHER: Final[str] = (
    "phospy.science.transformations.transformers.identity"
)

QUANTITATIVE_MEANING_PROVENANCE_SCHEMA_VERSION_V1: Final[int] = 1

QUANTITATIVE_MEANING_OPERATION_CALLER_DECLARATION: Final[str] = (
    "phospy.dataset_builder.quantitative_meaning.declaration"
)

QUANTITATIVE_MEANING_OPERATION_SCALE_CONTRACT_INFERENCE: Final[str] = (
    "phospy.dataset_builder.quantitative_meaning.infer_from_scale_contract"
)

QUANTITATIVE_MEANING_OPERATION_TOTAL_PROTEIN_SUBTRACT_LOG_TOTAL: Final[str] = (
    "phospy.dataset_preprocessing.total_protein_correction.subtract_log_total"
)

QUANTITATIVE_MEANING_OPERATION_LEGACY_BUNDLE_MIGRATION: Final[str] = (
    "phospy.bundle.legacy_intensity_scale_state_quantitative_meaning_migration"
)

QUANTITATIVE_MEANING_USER_DECLARED_CAVEAT_CODE: Final[str] = (
    "quantitative_meaning_user_declared"
)

QUANTITATIVE_MEANING_LEGACY_UNVERIFIED_CAVEAT_CODE: Final[str] = (
    "quantitative_meaning_legacy_unverified"
)

DATASET_QUANTITATIVE_MEANING_AUTHORITY_SOURCE: Final[str] = (
    "phospy.science.datasets.preprocessing.state_builder"
)

BUNDLE_QUANTITATIVE_MEANING_AUTHORITY_SOURCE: Final[str] = (
    "phospy.io.bundles._shared.intensity_scale_state"
)


class IntensityScaleKind(str, Enum):
    """Supported quantitative-intensity scales."""

    LINEAR = "linear"
    LOG2 = "log2"


class IntensityScaleEstablishmentMode(str, Enum):
    """How intensity-scale state was established."""

    DECLARED = "declared"
    TRANSFORMED = "transformed"
    IDENTITY = "identity"
    DERIVED = "derived"


class DeclaredIntensityScaleDiagnosticPolicy(str, Enum):
    """Policy for suspicious declared input intensity-scale diagnostics."""

    WARN = "warn"
    ERROR = "error"


class IntensityScaleEstablishmentSource(str, Enum):
    """Provenance source for how intensity-scale truth was established."""

    TRANSFORMED_BY_PHOSPY = "transformed_by_phospy"
    DECLARED_BY_USER = "declared_by_user"
    RESTORED_FROM_TRUSTED_PROVENANCE = "restored_from_trusted_provenance"


class IntensityScaleEvidenceLevel(str, Enum):
    """Evidence level supporting an intensity-scale transition event."""

    OBSERVED_TRANSFORMATION = "observed_transformation"
    DECLARED_BY_USER = "declared_by_user"
    INFERRED_FROM_METADATA = "inferred_from_metadata"
    UNKNOWN = "unknown"


class QuantitativeMeaning(str, Enum):
    """Scientific interpretation of phospho matrix values."""

    PHOSPHOSITE_ABUNDANCE = "phosphosite_abundance"
    PHOSPHOSITE_LOG_ABUNDANCE = "phosphosite_log_abundance"
    PHOSPHO_TOTAL_LOG_RATIO = "phospho_total_log_ratio"
    CONTRAST_LOG2_FOLD_CHANGE = "contrast_log2_fold_change"
    DIFFERENTIAL_EFFECT_SIZE = "differential_effect_size"
    ACTIVITY_SCORE = "activity_score"
    MIXED_PHOSPHO_TOTAL_LOG_RATIO_AND_PHOSPHOSITE_LOG_ABUNDANCE = (
        "mixed_phospho_total_log_ratio_and_phosphosite_log_abundance"
    )
    UNKNOWN = "unknown"


class QuantitativeMeaningEvidenceMode(str, Enum):
    """Evidence mode supporting quantitative-meaning state."""

    DERIVED_BY_PHOSPY_OPERATION = "derived_by_phospy_operation"
    DECLARED_BY_CALLER = "declared_by_caller"
    RESTORED_FROM_TRUSTED_SERIALIZED_PROVENANCE = (
        "restored_from_trusted_serialized_provenance"
    )
    INFERRED_FROM_SCALE_CONTRACT = "inferred_from_scale_contract"
    LEGACY_UNVERIFIED = "legacy_unverified"


CALLER_DECLARABLE_QUANTITATIVE_MEANINGS: Final[frozenset[QuantitativeMeaning]] = (
    frozenset(
        {
            QuantitativeMeaning.UNKNOWN,
            QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE,
            QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE,
        }
    )
)


@dataclass(frozen=True, slots=True)
class QuantitativeMeaningScaleRule:
    """Allowed scale and semantic role for one quantitative meaning."""

    meaning: QuantitativeMeaning
    allowed_scales: frozenset[IntensityScaleKind]
    semantic_role: str


_QUANTITATIVE_MEANING_SCALE_RULES: Final[tuple[QuantitativeMeaningScaleRule, ...]] = (
    QuantitativeMeaningScaleRule(
        meaning=QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE,
        allowed_scales=frozenset({IntensityScaleKind.LINEAR}),
        semantic_role="phosphosite_abundance_matrix",
    ),
    QuantitativeMeaningScaleRule(
        meaning=QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE,
        allowed_scales=frozenset({IntensityScaleKind.LOG2}),
        semantic_role="phosphosite_abundance_matrix",
    ),
    QuantitativeMeaningScaleRule(
        meaning=QuantitativeMeaning.PHOSPHO_TOTAL_LOG_RATIO,
        allowed_scales=frozenset({IntensityScaleKind.LOG2}),
        semantic_role="total_corrected_log_ratio_matrix",
    ),
    QuantitativeMeaningScaleRule(
        meaning=QuantitativeMeaning.CONTRAST_LOG2_FOLD_CHANGE,
        allowed_scales=frozenset({IntensityScaleKind.LOG2}),
        semantic_role="contrast_or_effect_matrix",
    ),
    QuantitativeMeaningScaleRule(
        meaning=QuantitativeMeaning.DIFFERENTIAL_EFFECT_SIZE,
        allowed_scales=frozenset({IntensityScaleKind.LOG2}),
        semantic_role="contrast_or_effect_matrix",
    ),
    QuantitativeMeaningScaleRule(
        meaning=QuantitativeMeaning.ACTIVITY_SCORE,
        allowed_scales=frozenset({IntensityScaleKind.LINEAR, IntensityScaleKind.LOG2}),
        semantic_role="activity_score_matrix",
    ),
    QuantitativeMeaningScaleRule(
        meaning=(
            QuantitativeMeaning.MIXED_PHOSPHO_TOTAL_LOG_RATIO_AND_PHOSPHOSITE_LOG_ABUNDANCE
        ),
        allowed_scales=frozenset({IntensityScaleKind.LOG2}),
        semantic_role="mixed_matrix",
    ),
    QuantitativeMeaningScaleRule(
        meaning=QuantitativeMeaning.UNKNOWN,
        allowed_scales=frozenset({IntensityScaleKind.LINEAR, IntensityScaleKind.LOG2}),
        semantic_role="unknown_matrix",
    ),
)

_QUANTITATIVE_MEANING_SCALE_RULE_BY_MEANING: Final[
    dict[QuantitativeMeaning, QuantitativeMeaningScaleRule]
] = {rule.meaning: rule for rule in _QUANTITATIVE_MEANING_SCALE_RULES}
if set(_QUANTITATIVE_MEANING_SCALE_RULE_BY_MEANING) != set(QuantitativeMeaning):
    raise RuntimeError(
        "QuantitativeMeaning scale rules must cover every QuantitativeMeaning member"
    )


def normalize_intensity_scale_evidence_level(
    evidence_level: IntensityScaleEvidenceLevel | str,
) -> IntensityScaleEvidenceLevel:
    raw_evidence_level = cast(object, evidence_level)
    if isinstance(raw_evidence_level, IntensityScaleEvidenceLevel):
        return raw_evidence_level
    try:
        return IntensityScaleEvidenceLevel(str(raw_evidence_level).strip())
    except ValueError as exc:
        supported = ", ".join(member.value for member in IntensityScaleEvidenceLevel)
        raise InvalidTransformationStateError(
            "unsupported intensity-scale evidence level "
            f"{raw_evidence_level!r}; supported: {supported}"
        ) from exc


def normalize_quantitative_meaning(
    quantity: QuantitativeMeaning | str | None,
) -> QuantitativeMeaning | None:
    if quantity is None:
        return None
    if isinstance(quantity, QuantitativeMeaning):
        return quantity
    try:
        return QuantitativeMeaning(str(quantity))
    except ValueError as exc:
        supported = ", ".join(member.value for member in QuantitativeMeaning)
        raise InvalidTransformationStateError(
            "unsupported intensity-scale quantitative meaning "
            f"'{quantity}'; supported: {supported}"
        ) from exc


def normalize_optional_quantitative_meaning(
    quantity: QuantitativeMeaning | str | None,
) -> QuantitativeMeaning | None:
    return normalize_quantitative_meaning(quantity)


def normalize_required_quantitative_meaning(
    quantity: QuantitativeMeaning | str,
    *,
    field_name: str,
) -> QuantitativeMeaning:
    normalized = normalize_quantitative_meaning(quantity)
    if normalized is None:
        raise InvalidTransformationStateError(f"{field_name} must not be None")
    return normalized


def default_quantitative_meaning_for_scale_kind(
    kind: IntensityScaleKind,
) -> QuantitativeMeaning:
    """Return the base quantitative meaning implied by a scale contract."""

    if kind is IntensityScaleKind.LINEAR:
        return QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE
    if kind is IntensityScaleKind.LOG2:
        return QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE
    return QuantitativeMeaning.UNKNOWN


def is_caller_declarable_quantitative_meaning(
    meaning: QuantitativeMeaning | str,
) -> bool:
    """Return whether a public caller may declare this direct input meaning."""

    normalized = normalize_required_quantitative_meaning(
        meaning,
        field_name="quantitative_meaning",
    )
    return normalized in CALLER_DECLARABLE_QUANTITATIVE_MEANINGS


def caller_declarable_quantitative_meaning_values() -> tuple[str, ...]:
    """Return stable public caller-declarable quantitative meaning values."""

    return tuple(
        meaning.value
        for meaning in (
            QuantitativeMeaning.UNKNOWN,
            QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE,
            QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE,
        )
    )


def validate_quantitative_meaning_kind_coherence(
    *,
    quantity: QuantitativeMeaning,
    kind: IntensityScaleKind,
) -> None:
    rule = _QUANTITATIVE_MEANING_SCALE_RULE_BY_MEANING[quantity]
    if kind in rule.allowed_scales:
        return
    allowed = ", ".join(sorted(scale.value for scale in rule.allowed_scales))
    if len(rule.allowed_scales) == 1:
        allowed = f"{allowed} intensity scale"
    else:
        allowed = f"one of these intensity scales: {allowed}"
    if quantity is QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE:
        raise InvalidTransformationStateError(
            "quantitative meaning 'phosphosite_abundance' requires linear "
            "intensity scale"
        )
    raise InvalidTransformationStateError(
        f"quantitative meaning '{quantity.value}' requires {allowed}"
    )


def resolve_establishment_source(
    *,
    authority_source: str,
    establishment_mode: IntensityScaleEstablishmentMode,
) -> IntensityScaleEstablishmentSource:
    if authority_source == "phospy.io.bundles._shared.intensity_scale_state":
        return IntensityScaleEstablishmentSource.RESTORED_FROM_TRUSTED_PROVENANCE
    if establishment_mode is IntensityScaleEstablishmentMode.DECLARED:
        return IntensityScaleEstablishmentSource.DECLARED_BY_USER
    return IntensityScaleEstablishmentSource.TRANSFORMED_BY_PHOSPY


__all__ = [
    "CALLER_DECLARABLE_QUANTITATIVE_MEANINGS",
    "DeclaredIntensityScaleDiagnosticPolicy",
    "IDENTITY_INTENSITY_SCALE_ESTABLISHER",
    "IntensityScaleEstablishmentMode",
    "IntensityScaleEstablishmentSource",
    "IntensityScaleEvidenceLevel",
    "IntensityScaleKind",
    "QUANTITATIVE_MEANING_LEGACY_UNVERIFIED_CAVEAT_CODE",
    "QUANTITATIVE_MEANING_OPERATION_CALLER_DECLARATION",
    "QUANTITATIVE_MEANING_OPERATION_LEGACY_BUNDLE_MIGRATION",
    "QUANTITATIVE_MEANING_OPERATION_SCALE_CONTRACT_INFERENCE",
    "QUANTITATIVE_MEANING_OPERATION_TOTAL_PROTEIN_SUBTRACT_LOG_TOTAL",
    "QUANTITATIVE_MEANING_PROVENANCE_SCHEMA_VERSION_V1",
    "QUANTITATIVE_MEANING_USER_DECLARED_CAVEAT_CODE",
    "QuantitativeMeaning",
    "QuantitativeMeaningEvidenceMode",
    "QuantitativeMeaningScaleRule",
    "caller_declarable_quantitative_meaning_values",
    "default_quantitative_meaning_for_scale_kind",
    "is_caller_declarable_quantitative_meaning",
]
