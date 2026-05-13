"""Differential workflow stage-boundary models and contracts."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Protocol

from phospy.api.configs import DifferentialAnalysisConfig
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.design.models import Contrast, ExperimentalDesign
from phospy.science.differential.models import (
    ContrastMatrix,
    DesignMatrix,
    DifferentialAnalysisResult,
)
from phospy.science.differential.models import (
    DifferentialAnalysisRequest as DifferentialComputationRequest,
)


@dataclass(frozen=True, slots=True)
class ValidatedDifferentialAnalysisRequest:
    """Validated differential request passed to interpretation."""

    dataset: AnalysisReadyPhosphoDataset
    design: ExperimentalDesign
    contrasts: tuple[Contrast, ...]
    analysis_sample_ids: tuple[str, ...]
    design_matrix: DesignMatrix
    contrast_matrix: ContrastMatrix
    config: DifferentialAnalysisConfig
    workflow_provenance: Mapping[str, object] | None = None


@dataclass(frozen=True, slots=True)
class InterpretedDifferentialAnalysisRequest:
    """Execution-ready differential request produced by interpretation."""

    computation_request: DifferentialComputationRequest
    config: DifferentialAnalysisConfig
    design_rank: int
    residual_degrees_of_freedom: float
    workflow_provenance: Mapping[str, object] | None = None


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
