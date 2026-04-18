"""Dataset builder contracts."""

from __future__ import annotations

from dataclasses import dataclass
from os import PathLike
from pathlib import Path
from typing import TYPE_CHECKING, Protocol

import pandas as pd

from phospy.references.models import Organism
from phospy.transformations.models import TransformationState

if TYPE_CHECKING:
    from phospy.api.requests import DatasetBuildRequest
    from phospy.datasets.models import AnalysisReadyPhosphoDataset

DatasetInput = pd.DataFrame | str | Path | PathLike[str]


@dataclass(frozen=True, slots=True)
class InterpretedDatasetBuildRequest:
    """Resolved builder input after request interpretation."""

    phospho: pd.DataFrame
    site_metadata: pd.DataFrame
    sample_metadata: pd.DataFrame | None
    total: pd.DataFrame | None
    organism: Organism | None
    transformation_state: TransformationState | None


class DatasetBuildValidatorContract(Protocol):
    """Internal contract for dataset build request validation."""

    def run(self, request: DatasetBuildRequest) -> DatasetBuildRequest:
        """Validate builder request and return the same request when valid."""


class DatasetBuildInterpreterContract(Protocol):
    """Internal contract for request interpretation into executable inputs."""

    def run(self, request: DatasetBuildRequest) -> InterpretedDatasetBuildRequest:
        """Resolve supported inputs into concrete DataFrame values."""


class DatasetBuildExecutorContract(Protocol):
    """Internal contract for constructing the final dataset."""

    def run(
        self, request: InterpretedDatasetBuildRequest
    ) -> AnalysisReadyPhosphoDataset:
        """Execute builder logic and return an analysis-ready dataset."""


__all__ = [
    "DatasetBuildExecutorContract",
    "DatasetBuildInterpreterContract",
    "DatasetBuildValidatorContract",
    "DatasetInput",
    "InterpretedDatasetBuildRequest",
]
