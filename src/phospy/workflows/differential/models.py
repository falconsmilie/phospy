"""Differential workflow stage-boundary models and contracts."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

from phospy.api.requests import MultipleTestingConfig
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.differential.models import (
    ContrastMatrix,
    DesignMatrix,
    DifferentialAnalysisResult,
    EmpiricalBayesConfig,
)
from phospy.differential.models import (
    DifferentialAnalysisRequest as DifferentialComputationRequest,
)


@dataclass(frozen=True, slots=True)
class ValidatedDifferentialAnalysisRequest:
    """Validated differential request passed to interpretation."""

    dataset: AnalysisReadyPhosphoDataset
    design: DesignMatrix
    contrasts: ContrastMatrix
    empirical_bayes: EmpiricalBayesConfig
    multiple_testing: MultipleTestingConfig


@dataclass(frozen=True, slots=True)
class InterpretedDifferentialAnalysisRequest:
    """Execution-ready differential request produced by interpretation."""

    computation_request: DifferentialComputationRequest
    multiple_testing: MultipleTestingConfig
    design_rank: int
    residual_degrees_of_freedom: float


class DifferentialAnalysisValidatorContract(Protocol):
    """Internal contract for differential workflow validation."""

    def run(self, request: object) -> ValidatedDifferentialAnalysisRequest: ...


class DifferentialAnalysisInterpreterContract(Protocol):
    """Internal contract for differential workflow interpretation."""

    def run(
        self, request: ValidatedDifferentialAnalysisRequest
    ) -> InterpretedDifferentialAnalysisRequest: ...


class DifferentialAnalysisExecutorContract(Protocol):
    """Internal contract for differential workflow execution."""

    def run(
        self, request: InterpretedDifferentialAnalysisRequest
    ) -> DifferentialAnalysisResult: ...


__all__ = [
    "DifferentialAnalysisExecutorContract",
    "DifferentialAnalysisInterpreterContract",
    "DifferentialAnalysisValidatorContract",
    "InterpretedDifferentialAnalysisRequest",
    "ValidatedDifferentialAnalysisRequest",
]
