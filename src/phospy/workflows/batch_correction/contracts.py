"""Batch-correction workflow boundary contracts."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Protocol

import pandas as pd

from phospy.contracts.configs.preprocessing import (
    CorrectionMissingnessPolicy,
    InternalBatchCorrectionRequest,
)
from phospy.provenance.models import BatchCorrectionProvenance, JsonValue
from phospy.science.datasets.preprocessing.control_sites import (
    ControlSiteMapping,
    ControlSiteSet,
)
from phospy.validation.datasets.batch_correction import ResolvedBatchDesignMetadata
from phospy.workflows.batch_correction.interpreter import ResolvedBatchCorrectionPlan


@dataclass(frozen=True, slots=True)
class BatchCorrectionWorkflowRequest:
    """Inputs required to orchestrate native batch correction."""

    phospho: pd.DataFrame
    config: InternalBatchCorrectionRequest
    sample_metadata: pd.DataFrame | None
    control_site_set: ControlSiteSet | None
    missingness_policy: CorrectionMissingnessPolicy | None = None
    upstream_observation_mask: pd.DataFrame | None = None
    site_metadata: pd.DataFrame | None = None
    dataset_organism: object | None = None


@dataclass(frozen=True, slots=True)
class BatchCorrectionWorkflowResult:
    """Structured output from the batch-correction workflow."""

    corrected_matrix: pd.DataFrame
    diagnostics: Mapping[str, JsonValue]
    warnings: tuple[str, ...]
    provenance: BatchCorrectionProvenance
    corrected_preprocessing_output: object | None = None

    @property
    def corrected(self) -> pd.DataFrame:
        """Return the corrected phosphosite matrix."""

        return self.corrected_matrix


class BatchCorrectionRequestValidatorContract(Protocol):
    """Validate the workflow request shell."""

    def run(self, request: object) -> BatchCorrectionWorkflowRequest: ...


class BatchCorrectionDesignValidatorContract(Protocol):
    """Validate and resolve design metadata before interpretation."""

    def run(
        self, *, request: BatchCorrectionWorkflowRequest
    ) -> ResolvedBatchDesignMetadata: ...


class BatchCorrectionControlSiteValidatorContract(Protocol):
    """Validate and resolve control-site mappings before interpretation."""

    def run(self, *, request: BatchCorrectionWorkflowRequest) -> ControlSiteMapping: ...


class BatchCorrectionStageOrderValidatorContract(Protocol):
    """Validate requested stage-order policy before interpretation."""

    def run(self, *, config: InternalBatchCorrectionRequest) -> None: ...


class BatchCorrectionMissingnessValidatorContract(Protocol):
    """Validate missingness policy before interpretation and execution."""

    def run(
        self, *, request: BatchCorrectionWorkflowRequest
    ) -> CorrectionMissingnessPolicy: ...


class BatchCorrectionInterpreterContract(Protocol):
    """Resolve validated inputs into an execution plan."""

    def run(
        self,
        *,
        config: InternalBatchCorrectionRequest,
        dataset_metadata: ResolvedBatchDesignMetadata,
        control_site_mapping: ControlSiteMapping,
        missingness_policy: CorrectionMissingnessPolicy,
    ) -> ResolvedBatchCorrectionPlan: ...


class BatchCorrectionExecutorDiagnosticsContract(Protocol):
    """Diagnostics emitted by a numerical batch-correction executor."""

    def to_payload(self) -> dict[str, object]: ...


class BatchCorrectionExecutorResultContract(Protocol):
    """Executor result fields consumed by the workflow shell."""

    @property
    def corrected_matrix(self) -> pd.DataFrame: ...

    @property
    def diagnostics(self) -> BatchCorrectionExecutorDiagnosticsContract: ...

    @property
    def warnings(self) -> Sequence[str]: ...

    @property
    def output_observation_mask(self) -> pd.DataFrame: ...

    @property
    def provenance_payload(self) -> Mapping[str, object]: ...

    @property
    def rejected_rows(self) -> Sequence[str]: ...

    @property
    def rejected_cells(self) -> Sequence[tuple[str, str]]: ...

    @property
    def withheld_rows(self) -> Sequence[str]: ...

    @property
    def withheld_cells(self) -> Sequence[tuple[str, str]]: ...


class BatchCorrectionExecutorContract(Protocol):
    """Execute a resolved batch-correction plan."""

    def run(
        self,
        *,
        phospho: pd.DataFrame,
        plan: ResolvedBatchCorrectionPlan,
    ) -> BatchCorrectionExecutorResultContract: ...


class BatchCorrectionProvenanceRecorderContract(Protocol):
    """Assemble workflow provenance after execution."""

    def run(
        self,
        *,
        request: BatchCorrectionWorkflowRequest,
        dataset_metadata: ResolvedBatchDesignMetadata,
        control_site_mapping: ControlSiteMapping,
        missingness_policy: CorrectionMissingnessPolicy,
        plan: ResolvedBatchCorrectionPlan,
        executor_result: BatchCorrectionExecutorResultContract,
    ) -> BatchCorrectionProvenance: ...


__all__ = [
    "BatchCorrectionControlSiteValidatorContract",
    "BatchCorrectionDesignValidatorContract",
    "BatchCorrectionExecutorContract",
    "BatchCorrectionExecutorDiagnosticsContract",
    "BatchCorrectionExecutorResultContract",
    "BatchCorrectionInterpreterContract",
    "BatchCorrectionMissingnessValidatorContract",
    "BatchCorrectionProvenanceRecorderContract",
    "BatchCorrectionRequestValidatorContract",
    "BatchCorrectionStageOrderValidatorContract",
    "BatchCorrectionWorkflowRequest",
    "BatchCorrectionWorkflowResult",
]
