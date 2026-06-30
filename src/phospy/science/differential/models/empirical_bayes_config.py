"""Empirical-Bayes configuration models for differential analysis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import cast

from phospy.errors.input import PhosPyInputError

EMPIRICAL_BAYES_METHOD_STANDARD = "standard"
EMPIRICAL_BAYES_METHOD_ROBUST = "robust"
SUPPORTED_EMPIRICAL_BAYES_METHODS: tuple[str, ...] = (
    EMPIRICAL_BAYES_METHOD_STANDARD,
    EMPIRICAL_BAYES_METHOD_ROBUST,
)


@dataclass(frozen=True, slots=True)
class EmpiricalBayesConfig:
    """Empirical-Bayes configuration for moderated statistics."""

    method: str = EMPIRICAL_BAYES_METHOD_STANDARD
    trend: bool = False
    winsor_tail_p: tuple[float, float] = (0.05, 0.1)

    def __post_init__(self) -> None:
        if self.method not in SUPPORTED_EMPIRICAL_BAYES_METHODS:
            supported = ", ".join(
                repr(value) for value in SUPPORTED_EMPIRICAL_BAYES_METHODS
            )
            raise PhosPyInputError(
                f"empirical_bayes.method must be one of: {supported}"
            )
        if not isinstance(cast(object, self.trend), bool):
            raise PhosPyInputError("empirical_bayes.trend must be a bool")
        winsor_tail_p = self.winsor_tail_p
        if (
            not isinstance(cast(object, winsor_tail_p), tuple)
            or len(winsor_tail_p) != 2
            or not all(
                isinstance(cast(object, value), int | float) for value in winsor_tail_p
            )
        ):
            raise PhosPyInputError(
                "empirical_bayes.winsor_tail_p must be a tuple of two numeric values"
            )
        left_tail_p = float(winsor_tail_p[0])
        right_tail_p = float(winsor_tail_p[1])
        if not (0.0 <= left_tail_p < 1.0 and 0.0 <= right_tail_p < 1.0):
            raise PhosPyInputError(
                "empirical_bayes.winsor_tail_p values must each be in [0.0, 1.0)"
            )
        if left_tail_p + right_tail_p >= 1.0:
            raise PhosPyInputError(
                "empirical_bayes.winsor_tail_p values must sum to less than 1.0"
            )
        object.__setattr__(self, "winsor_tail_p", (left_tail_p, right_tail_p))


__all__ = [
    "EMPIRICAL_BAYES_METHOD_ROBUST",
    "EMPIRICAL_BAYES_METHOD_STANDARD",
    "SUPPORTED_EMPIRICAL_BAYES_METHODS",
    "EmpiricalBayesConfig",
]
