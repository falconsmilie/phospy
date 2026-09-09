"""Empirical-Bayes configuration models for differential analysis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, cast

from phospy.errors.input import PhosPyInputError

EMPIRICAL_BAYES_METHOD_STANDARD = "standard"
EMPIRICAL_BAYES_METHOD_ROBUST = "robust"
EMPIRICAL_BAYES_TREND_COVARIATE_MEAN_INTENSITY = "mean_intensity"
EMPIRICAL_BAYES_TREND_COVARIATE_QUANTIFICATION_DEPTH = "quantification_depth"
QUANTIFICATION_DEPTH_KIND_PSM_COUNT = "psm_count"
QUANTIFICATION_DEPTH_KIND_PEPTIDE_COUNT = "peptide_count"
EmpiricalBayesMethod = Literal["standard", "robust"]
EmpiricalBayesTrendCovariate = Literal["mean_intensity", "quantification_depth"]
QuantificationDepthKind = Literal["psm_count", "peptide_count"]
SUPPORTED_EMPIRICAL_BAYES_METHODS: tuple[str, ...] = (
    EMPIRICAL_BAYES_METHOD_STANDARD,
    EMPIRICAL_BAYES_METHOD_ROBUST,
)
SUPPORTED_EMPIRICAL_BAYES_TREND_COVARIATES: tuple[str, ...] = (
    EMPIRICAL_BAYES_TREND_COVARIATE_MEAN_INTENSITY,
    EMPIRICAL_BAYES_TREND_COVARIATE_QUANTIFICATION_DEPTH,
)
SUPPORTED_QUANTIFICATION_DEPTH_KINDS: tuple[str, ...] = (
    QUANTIFICATION_DEPTH_KIND_PSM_COUNT,
    QUANTIFICATION_DEPTH_KIND_PEPTIDE_COUNT,
)


@dataclass(frozen=True, slots=True)
class EmpiricalBayesConfig:
    """Empirical-Bayes configuration for moderated statistics."""

    method: EmpiricalBayesMethod = EMPIRICAL_BAYES_METHOD_STANDARD
    trend: bool = False
    trend_covariate: EmpiricalBayesTrendCovariate = (
        EMPIRICAL_BAYES_TREND_COVARIATE_MEAN_INTENSITY
    )
    quantification_depth_kind: QuantificationDepthKind | None = None
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
        if self.trend_covariate not in SUPPORTED_EMPIRICAL_BAYES_TREND_COVARIATES:
            supported = ", ".join(
                repr(value) for value in SUPPORTED_EMPIRICAL_BAYES_TREND_COVARIATES
            )
            raise PhosPyInputError(
                f"empirical_bayes.trend_covariate must be one of: {supported}"
            )
        if self.trend_covariate == EMPIRICAL_BAYES_TREND_COVARIATE_QUANTIFICATION_DEPTH:
            if not self.trend:
                raise PhosPyInputError(
                    "empirical_bayes.trend must be True when "
                    "empirical_bayes.trend_covariate is 'quantification_depth'"
                )
            if (
                self.quantification_depth_kind
                not in SUPPORTED_QUANTIFICATION_DEPTH_KINDS
            ):
                supported = ", ".join(
                    repr(value) for value in SUPPORTED_QUANTIFICATION_DEPTH_KINDS
                )
                raise PhosPyInputError(
                    "empirical_bayes.quantification_depth_kind must be one of "
                    f"{supported} when empirical_bayes.trend_covariate is "
                    "'quantification_depth'"
                )
        elif self.quantification_depth_kind is not None:
            raise PhosPyInputError(
                "empirical_bayes.quantification_depth_kind must be None when "
                "empirical_bayes.trend_covariate is 'mean_intensity'"
            )
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
    "EMPIRICAL_BAYES_TREND_COVARIATE_MEAN_INTENSITY",
    "EMPIRICAL_BAYES_TREND_COVARIATE_QUANTIFICATION_DEPTH",
    "EMPIRICAL_BAYES_METHOD_ROBUST",
    "EMPIRICAL_BAYES_METHOD_STANDARD",
    "EmpiricalBayesMethod",
    "SUPPORTED_EMPIRICAL_BAYES_METHODS",
    "EmpiricalBayesTrendCovariate",
    "EmpiricalBayesConfig",
    "QUANTIFICATION_DEPTH_KIND_PEPTIDE_COUNT",
    "QUANTIFICATION_DEPTH_KIND_PSM_COUNT",
    "QuantificationDepthKind",
    "SUPPORTED_EMPIRICAL_BAYES_TREND_COVARIATES",
    "SUPPORTED_QUANTIFICATION_DEPTH_KINDS",
]
