"""Batch-correction provenance validation compatibility routes."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from typing import Any, cast

_VALIDATION_DATASET_MODULE_PARTS = ("phospy", "validation", "datasets")


def _validation_provenance_module() -> object:
    # Dynamic import avoids a static science->validation package cycle while keeping
    # applied-output validation implementation in validation.datasets.
    return __import__(
        ".".join((*_VALIDATION_DATASET_MODULE_PARTS, "batch_correction_provenance")),
        fromlist=("validate_applied_native_sps_ruv_correction_provenance",),
    )


def _validation_controls_module() -> object:
    # Dynamic import avoids a static science->validation package cycle while keeping
    # selected-control normalization implementation in validation.datasets.
    return __import__(
        ".".join((*_VALIDATION_DATASET_MODULE_PARTS, "batch_correction_controls")),
        fromlist=("normalize_applied_selected_site_key_rows",),
    )


validate_applied_native_sps_ruv_correction_provenance = cast(
    Callable[..., None],
    cast(
        Any,
        _validation_provenance_module(),
    ).validate_applied_native_sps_ruv_correction_provenance,
)

normalize_applied_selected_site_key_rows = cast(
    Callable[[Sequence[object]], tuple[str, ...]],
    cast(
        Any,
        _validation_controls_module(),
    ).normalize_applied_selected_site_key_rows,
)

__all__ = [
    "normalize_applied_selected_site_key_rows",
    "validate_applied_native_sps_ruv_correction_provenance",
]
