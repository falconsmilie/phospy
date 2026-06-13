"""Batch-correction preprocessing intent configuration."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from phospy.validation.configs.preprocessing import (
    validate_batch_correction_config,
)

DATASET_BATCH_CORRECTION_METHOD_NONE = "none"
DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH = "linear_residualize_batch"
DatasetBatchCorrectionMethod = Literal["none", "linear_residualize_batch"]
DATASET_BATCH_CORRECTION_METHODS = frozenset(
    {
        DATASET_BATCH_CORRECTION_METHOD_NONE,
        DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH,
    }
)


@dataclass(frozen=True, slots=True)
class DatasetBatchCorrectionConfig:
    """Public batch-correction intent for dataset preprocessing.

    - `"none"`: do not request batch correction.
    - `"linear_residualize_batch"`: run fixed-effect residualisation of batch
      terms while preserving condition effects by design during dataset
      preprocessing.

    `"linear_residualize_batch"` is not ComBat, not RUV, and not limma
    `removeBatchEffect` parity. Dataset-build execution resolves sample
    metadata, validates design adequacy, applies correction, and records a
    typed preprocessing report.
    """

    method: DatasetBatchCorrectionMethod = DATASET_BATCH_CORRECTION_METHOD_NONE
    batch_column: str = "batch"
    condition_column: str = "condition"
    preserve_condition_effects: Literal[True] = True

    def __post_init__(self) -> None:
        validate_batch_correction_config(
            method=self.method,
            batch_column=self.batch_column,
            condition_column=self.condition_column,
            preserve_condition_effects=self.preserve_condition_effects,
            supported_methods=DATASET_BATCH_CORRECTION_METHODS,
        )


__all__ = [
    "DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH",
    "DATASET_BATCH_CORRECTION_METHOD_NONE",
    "DATASET_BATCH_CORRECTION_METHODS",
    "DatasetBatchCorrectionConfig",
    "DatasetBatchCorrectionMethod",
]
