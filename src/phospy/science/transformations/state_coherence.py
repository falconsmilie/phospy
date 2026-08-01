"""Science-owned intensity-scale-state coherence checks."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import NoReturn, cast

import numpy as np
import pandas as pd

from phospy.errors.validation import TransformationValidationError
from phospy.science.transformations.models import (
    IntensityScaleKind,
    IntensityScaleState,
    QuantitativeMeaning,
)


class ObservedNumericDomain(str, Enum):
    """Observed sign domain for finite numeric matrix values."""

    ZERO_ONLY = "zero_only"
    NON_NEGATIVE = "non_negative"
    SIGNED = "signed"
    NEGATIVE_ONLY = "negative_only"


@dataclass(frozen=True, slots=True)
class NumericDomainObservation:
    """Observed numeric sign-domain summary for one matrix."""

    table_name: str
    observed_domain: ObservedNumericDomain
    value_count: int
    negative_count: int
    zero_count: int
    positive_count: int
    minimum: float
    maximum: float


_NON_NEGATIVE_DOMAINS = frozenset(
    {
        ObservedNumericDomain.ZERO_ONLY,
        ObservedNumericDomain.NON_NEGATIVE,
    }
)
_SIGNED_DOMAINS = frozenset(ObservedNumericDomain)
_PHOSPHO_NUMERIC_DOMAIN_RULES: dict[
    QuantitativeMeaning, frozenset[ObservedNumericDomain]
] = {
    QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE: _NON_NEGATIVE_DOMAINS,
    QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE: _SIGNED_DOMAINS,
    QuantitativeMeaning.PHOSPHO_TOTAL_LOG_RATIO: _SIGNED_DOMAINS,
    QuantitativeMeaning.CONTRAST_LOG2_FOLD_CHANGE: _SIGNED_DOMAINS,
    QuantitativeMeaning.DIFFERENTIAL_EFFECT_SIZE: _SIGNED_DOMAINS,
    QuantitativeMeaning.ACTIVITY_SCORE: _SIGNED_DOMAINS,
    QuantitativeMeaning.MIXED_PHOSPHO_TOTAL_LOG_RATIO_AND_PHOSPHOSITE_LOG_ABUNDANCE: (
        _SIGNED_DOMAINS
    ),
    QuantitativeMeaning.UNKNOWN: frozenset(),
}
if set(_PHOSPHO_NUMERIC_DOMAIN_RULES) != set(QuantitativeMeaning):
    raise RuntimeError(
        "Phospho numeric-domain rules must cover every QuantitativeMeaning member"
    )


def require_intensity_scale_state_coherence(
    intensity_scale_state: IntensityScaleState,
    *,
    has_total_matrix: bool,
    require_established: bool = False,
) -> IntensityScaleState:
    """Return a coherent intensity-scale state or raise a validation error."""

    if not isinstance(cast(object, intensity_scale_state), IntensityScaleState):
        raise TransformationValidationError(
            "dataset.intensity_scale_state must be an IntensityScaleState instance"
        )
    if require_established and not intensity_scale_state.is_established:
        raise TransformationValidationError(
            "dataset.intensity_scale_state must be established through a "
            "supported PhosPy path; use AnalysisReadyDatasetBuilder or a "
            "supported transformer/bundle reconstruction path"
        )
    if require_established and intensity_scale_state.quantity is None:
        raise TransformationValidationError(
            "dataset.intensity_scale_state must have established quantitative "
            "meaning provenance"
        )
    if (
        require_established
        and intensity_scale_state.quantitative_meaning_provenance is None
    ):
        raise TransformationValidationError(
            "dataset.intensity_scale_state must carry quantitative meaning "
            "provenance separate from intensity-scale establishment provenance"
        )
    if has_total_matrix and intensity_scale_state.total is None:
        raise TransformationValidationError(
            "intensity_scale_state.total is required when dataset.total is provided"
        )
    if not has_total_matrix and intensity_scale_state.total is not None:
        raise TransformationValidationError(
            "intensity_scale_state.total must be None when dataset.total is absent"
        )
    if (
        intensity_scale_state.total is not None
        and intensity_scale_state.total.kind is not intensity_scale_state.phospho.kind
    ):
        raise TransformationValidationError(
            "phospho and total matrices must share one intensity scale kind"
        )
    return intensity_scale_state


def require_quantitative_numeric_domain_coherence(
    *,
    phospho: pd.DataFrame,
    total: pd.DataFrame | None,
    intensity_scale_state: IntensityScaleState,
    allow_numeric_semantic_domain_waiver: bool = False,
    error_type: type[Exception] = TransformationValidationError,
) -> tuple[NumericDomainObservation, ...]:
    """Require matrix numeric domains to match established quantitative meaning.

    This is intentionally meaning-aware. It is not a table-schema positivity
    check: signed log effects, centred log quantities, and residual-style
    quantities remain valid when their quantitative meaning permits signed
    values.
    """

    if not isinstance(cast(object, intensity_scale_state), IntensityScaleState):
        _raise(
            error_type,
            "dataset.intensity_scale_state must be an IntensityScaleState instance",
        )
    observations: list[NumericDomainObservation] = []
    phospho_observation = observe_numeric_domain(
        phospho,
        table_name="dataset.phospho",
        error_type=error_type,
    )
    observations.append(phospho_observation)
    _require_phospho_domain_matches_meaning(
        observation=phospho_observation,
        intensity_scale_state=intensity_scale_state,
        allow_numeric_semantic_domain_waiver=allow_numeric_semantic_domain_waiver,
        error_type=error_type,
    )
    if total is not None:
        total_observation = observe_numeric_domain(
            total,
            table_name="dataset.total",
            error_type=error_type,
        )
        observations.append(total_observation)
        _require_total_domain_matches_scale(
            observation=total_observation,
            intensity_scale_state=intensity_scale_state,
            allow_numeric_semantic_domain_waiver=allow_numeric_semantic_domain_waiver,
            error_type=error_type,
        )
    return tuple(observations)


def observe_numeric_domain(
    matrix: pd.DataFrame,
    *,
    table_name: str,
    error_type: type[Exception] = TransformationValidationError,
) -> NumericDomainObservation:
    """Return the observed finite numeric sign domain for ``matrix``."""

    try:
        values = matrix.to_numpy(dtype="float64", copy=False).reshape(-1)
    except (AttributeError, TypeError, ValueError) as exc:
        _raise(
            error_type,
            f"{table_name} must contain numeric values before numeric-domain "
            "coherence can be assessed",
            cause=exc,
        )
    if values.size == 0:
        _raise(
            error_type,
            f"{table_name} must contain at least one value before numeric-domain "
            "coherence can be assessed",
        )
    finite_mask = np.isfinite(values)
    if not bool(finite_mask.all()):
        _raise(
            error_type,
            f"{table_name} must contain only finite numeric values before "
            "numeric-domain coherence can be assessed",
        )
    negative_count = int(np.count_nonzero(values < 0.0))
    zero_count = int(np.count_nonzero(values == 0.0))
    positive_count = int(np.count_nonzero(values > 0.0))
    if negative_count == 0 and positive_count == 0:
        observed_domain = ObservedNumericDomain.ZERO_ONLY
    elif negative_count == 0:
        observed_domain = ObservedNumericDomain.NON_NEGATIVE
    elif positive_count == 0:
        observed_domain = ObservedNumericDomain.NEGATIVE_ONLY
    else:
        observed_domain = ObservedNumericDomain.SIGNED
    return NumericDomainObservation(
        table_name=table_name,
        observed_domain=observed_domain,
        value_count=int(values.size),
        negative_count=negative_count,
        zero_count=zero_count,
        positive_count=positive_count,
        minimum=float(np.min(values)),
        maximum=float(np.max(values)),
    )


def _require_phospho_domain_matches_meaning(
    *,
    observation: NumericDomainObservation,
    intensity_scale_state: IntensityScaleState,
    allow_numeric_semantic_domain_waiver: bool,
    error_type: type[Exception],
) -> None:
    meaning = intensity_scale_state.quantity
    if meaning is None:
        _raise(
            error_type,
            _numeric_domain_error_message(
                observation=observation,
                scale=intensity_scale_state.phospho.kind.value,
                meaning=None,
                expected="an established quantitative meaning",
                reason=(
                    "dataset.phospho cannot be promoted to analysis-ready "
                    "numeric-semantic state without quantitative meaning evidence"
                ),
            ),
        )
    allowed_domains = _PHOSPHO_NUMERIC_DOMAIN_RULES[meaning]
    if observation.observed_domain in allowed_domains:
        return
    if allow_numeric_semantic_domain_waiver:
        return
    if meaning is QuantitativeMeaning.UNKNOWN:
        reason = (
            "unknown quantitative meaning has no numeric-domain contract; "
            "provide evidence for a concrete meaning or an explicit typed "
            "numeric_semantic_domain waiver"
        )
        expected = "a concrete quantitative meaning with a numeric-domain rule"
    elif meaning is QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE:
        reason = "linear phosphosite_abundance must be non-negative"
        expected = "non-negative values"
    else:
        reason = (
            f"quantitative meaning {meaning.value!r} does not allow observed "
            f"numeric domain {observation.observed_domain.value!r}"
        )
        expected = _format_allowed_domains(allowed_domains)
    _raise(
        error_type,
        _numeric_domain_error_message(
            observation=observation,
            scale=intensity_scale_state.phospho.kind.value,
            meaning=meaning.value,
            expected=expected,
            reason=reason,
        ),
    )


def _require_total_domain_matches_scale(
    *,
    observation: NumericDomainObservation,
    intensity_scale_state: IntensityScaleState,
    allow_numeric_semantic_domain_waiver: bool,
    error_type: type[Exception],
) -> None:
    total_state = intensity_scale_state.total
    if total_state is None:
        return
    if total_state.kind is not IntensityScaleKind.LINEAR:
        return
    if observation.observed_domain in _NON_NEGATIVE_DOMAINS:
        return
    if allow_numeric_semantic_domain_waiver:
        return
    _raise(
        error_type,
        _numeric_domain_error_message(
            observation=observation,
            scale=total_state.kind.value,
            meaning="total_protein_abundance",
            expected="non-negative values",
            reason="linear total_protein_abundance must be non-negative",
        ),
    )


def _numeric_domain_error_message(
    *,
    observation: NumericDomainObservation,
    scale: str,
    meaning: str | None,
    expected: str,
    reason: str,
) -> str:
    meaning_label = "<missing>" if meaning is None else meaning
    return (
        "analysis-ready numeric-semantic coherence failed: "
        f"table={observation.table_name!r}, scale={scale!r}, "
        f"meaning={meaning_label!r}, "
        f"observed_numeric_domain={observation.observed_domain.value!r}, "
        f"min={observation.minimum:.12g}, max={observation.maximum:.12g}, "
        f"negative_count={observation.negative_count}, "
        f"zero_count={observation.zero_count}, "
        f"positive_count={observation.positive_count}, "
        f"expected={expected}; {reason}"
    )


def _format_allowed_domains(domains: frozenset[ObservedNumericDomain]) -> str:
    if not domains:
        return "no numeric domains"
    return ", ".join(sorted(domain.value for domain in domains))


def _raise(
    error_type: type[Exception],
    message: str,
    *,
    cause: BaseException | None = None,
) -> NoReturn:
    error = error_type(message)
    if cause is not None:
        raise error from cause
    raise error


__all__ = [
    "NumericDomainObservation",
    "ObservedNumericDomain",
    "observe_numeric_domain",
    "require_intensity_scale_state_coherence",
    "require_quantitative_numeric_domain_coherence",
]
