"""Adapters for invoking batch-correction workflows from preprocessing orchestration."""

from __future__ import annotations

import pandas as pd

from phospy.contracts.configs.preprocessing import (
    CorrectionMissingnessPolicy,
    InternalBatchCorrectionRequest,
)
from phospy.science.datasets.preprocessing.control_sites import ControlSiteSet
from phospy.science.datasets.preprocessing.stages.batch_correction import (
    SpsRuvStyleBatchCorrectionResult,
)
from phospy.workflows.batch_correction.contracts import (
    BatchCorrectionWorkflowRequest,
)
from phospy.workflows.batch_correction.workflow import BatchCorrectionWorkflow


class SpsRuvStyleBatchCorrectionWorkflowRunner:
    """Adapter from preprocessing runner protocol to batch-correction workflow."""

    def __init__(self, *, workflow: BatchCorrectionWorkflow | None = None) -> None:
        self._workflow = workflow or BatchCorrectionWorkflow()

    def run(
        self,
        *,
        phospho: pd.DataFrame,
        config: object,
        sample_metadata: pd.DataFrame | None,
        control_site_set: ControlSiteSet,
        missingness_policy: object,
        upstream_observation_mask: pd.DataFrame | None,
        site_metadata: pd.DataFrame,
    ) -> SpsRuvStyleBatchCorrectionResult:
        return self._workflow.run(
            BatchCorrectionWorkflowRequest(
                phospho=phospho,
                config=_require_internal_config(config),
                sample_metadata=sample_metadata,
                control_site_set=control_site_set,
                missingness_policy=_require_missingness_policy(missingness_policy),
                upstream_observation_mask=upstream_observation_mask,
                site_metadata=site_metadata,
            )
        )


def _require_internal_config(value: object) -> InternalBatchCorrectionRequest:
    if isinstance(value, InternalBatchCorrectionRequest):
        return value
    raise TypeError(
        "SPS/RUV-style batch correction requires InternalBatchCorrectionRequest"
    )


def _require_missingness_policy(value: object) -> CorrectionMissingnessPolicy:
    if isinstance(value, CorrectionMissingnessPolicy):
        return value
    raise TypeError(
        "SPS/RUV-style batch correction requires CorrectionMissingnessPolicy"
    )


__all__ = ["SpsRuvStyleBatchCorrectionWorkflowRunner"]
