"""Request-shell validation for batch-correction workflows."""

from __future__ import annotations

import pandas as pd

from phospy.contracts.configs.preprocessing import InternalBatchCorrectionRequest
from phospy.errors.input import PhosPyInputError
from phospy.validation.configs.preprocessing import (
    reject_unsupported_ruv_iii_style_method,
)
from phospy.validation.workflows.batch_correction.protocols import (
    BatchCorrectionWorkflowRequestProtocol,
)


class BatchCorrectionWorkflowRequestValidator:
    """Validate workflow request containers without scientific interpretation."""

    def run(self, request: object) -> BatchCorrectionWorkflowRequestProtocol:
        if not isinstance(request, BatchCorrectionWorkflowRequestProtocol):
            raise PhosPyInputError(
                "batch-correction workflow request must be "
                "BatchCorrectionWorkflowRequest-compatible"
            )
        if not isinstance(request.config, InternalBatchCorrectionRequest):
            raise PhosPyInputError(
                "batch-correction workflow request.config must be "
                "InternalBatchCorrectionRequest"
            )
        reject_unsupported_ruv_iii_style_method(
            request.config.method,
            field_name="batch-correction workflow request.config.method",
        )
        if not isinstance(request.phospho, pd.DataFrame):
            raise PhosPyInputError(
                "batch-correction workflow request.phospho must be a pandas DataFrame"
            )
        if request.sample_metadata is not None and not isinstance(
            request.sample_metadata,
            pd.DataFrame,
        ):
            raise PhosPyInputError(
                "batch-correction workflow request.sample_metadata must be a pandas "
                "DataFrame when provided"
            )
        if request.site_metadata is not None and not isinstance(
            request.site_metadata,
            pd.DataFrame,
        ):
            raise PhosPyInputError(
                "batch-correction workflow request.site_metadata must be a pandas "
                "DataFrame when provided"
            )
        if request.upstream_observation_mask is not None and not isinstance(
            request.upstream_observation_mask,
            pd.DataFrame,
        ):
            raise PhosPyInputError(
                "batch-correction workflow request.upstream_observation_mask must "
                "be a pandas DataFrame when provided"
            )
        return request


__all__ = ["BatchCorrectionWorkflowRequestValidator"]
