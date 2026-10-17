"""Validation-owned protocols for batch-correction workflow inputs."""

from __future__ import annotations

from typing import Protocol, runtime_checkable

import pandas as pd

from phospy.preprocessing.batch_correction.contracts.missingness import (
    CorrectionMissingnessPolicy,
)
from phospy.preprocessing.batch_correction.contracts.request import (
    InternalBatchCorrectionRequest,
)
from phospy.preprocessing.control_sites.contracts import ControlSiteSet


@runtime_checkable
class BatchCorrectionWorkflowRequestProtocol(Protocol):
    """Request shape consumed by validation without importing workflow classes."""

    @property
    def phospho(self) -> pd.DataFrame: ...

    @property
    def config(self) -> InternalBatchCorrectionRequest: ...

    @property
    def sample_metadata(self) -> pd.DataFrame | None: ...

    @property
    def control_site_set(self) -> ControlSiteSet | None: ...

    @property
    def missingness_policy(self) -> CorrectionMissingnessPolicy | None: ...

    @property
    def upstream_observation_mask(self) -> pd.DataFrame | None: ...

    @property
    def site_metadata(self) -> pd.DataFrame | None: ...

    @property
    def dataset_organism(self) -> object | None: ...


__all__ = ["BatchCorrectionWorkflowRequestProtocol"]
