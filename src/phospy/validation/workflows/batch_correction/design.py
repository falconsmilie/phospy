"""Design validation wrapper for batch-correction workflows."""

from __future__ import annotations

from phospy.contracts.configs.preprocessing import InternalBatchCorrectionMethod
from phospy.validation.datasets.batch_correction import (
    BatchCorrectionAdequacyValidator,
    BatchDesignMetadataValidator,
    ResolvedBatchDesignMetadata,
)
from phospy.workflows.batch_correction.contracts import BatchCorrectionWorkflowRequest


class BatchCorrectionWorkflowDesignValidator:
    """Resolve sample metadata and validate batch/condition design adequacy."""

    def __init__(
        self,
        *,
        metadata_validator: BatchDesignMetadataValidator | None = None,
        adequacy_validator: BatchCorrectionAdequacyValidator | None = None,
    ) -> None:
        self._metadata_validator = metadata_validator or BatchDesignMetadataValidator()
        self._adequacy_validator = (
            adequacy_validator or BatchCorrectionAdequacyValidator()
        )

    def run(
        self,
        *,
        request: BatchCorrectionWorkflowRequest,
    ) -> ResolvedBatchDesignMetadata:
        config = request.config
        metadata = self._metadata_validator.run(
            phospho=request.phospho,
            sample_metadata=request.sample_metadata,
            batch_column=config.batch_column,
            condition_columns=config.condition_columns,
            replicate_column=config.replicate_column,
            require_replicate_column=(
                config.method is InternalBatchCorrectionMethod.RUV_III_STYLE
            ),
            context="batch-correction workflow",
        )
        self._adequacy_validator.run(
            batch_by_sample=metadata.batch_by_sample,
            condition_by_sample=metadata.condition_by_sample,
            sample_order=metadata.sample_order,
            preserve_condition_effects=True,
        )
        return metadata


__all__ = ["BatchCorrectionWorkflowDesignValidator"]
