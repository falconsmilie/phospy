"""Adapters for invoking batch-correction workflows from preprocessing orchestration."""

from __future__ import annotations

import pandas as pd

from phospy.preprocessing.batch_correction.contracts.missingness import (
    CorrectionMissingnessPolicy,
)
from phospy.preprocessing.batch_correction.contracts.request import (
    InternalBatchCorrectionRequest,
)
from phospy.preprocessing.batch_correction.contracts.workflow import (
    BatchCorrectionWorkflowRequest,
    BatchCorrectionWorkflowResult,
)
from phospy.preprocessing.batch_correction.workflow import BatchCorrectionWorkflow
from phospy.preprocessing.control_sites.contracts import ControlSiteSet


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
    ) -> BatchCorrectionWorkflowResult:
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
