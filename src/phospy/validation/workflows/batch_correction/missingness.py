"""Missingness validation wrapper for batch-correction workflows."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pandas as pd

from phospy.contracts.configs.preprocessing import (
    CorrectionMissingnessCompatibilityValidator,
    CorrectionMissingnessPolicy,
    ObservationMask,
    OriginallyMissingCellTracking,
    RowSampleEligibilityImpact,
)
from phospy.errors.input import PhosPyInputError
from phospy.workflows.batch_correction.contracts import BatchCorrectionWorkflowRequest

_UNSAFE_UPSTREAM_ELIGIBILITY_IMPACTS = frozenset(
    {
        RowSampleEligibilityImpact.EXCLUDE_SAMPLES_WITH_ORIGINALLY_MISSING_VALUES,
        RowSampleEligibilityImpact.REQUIRE_COMPLETE_CASES,
        RowSampleEligibilityImpact.UNSUPPORTED,
    }
)


class BatchCorrectionWorkflowMissingnessValidator:
    """Validate correction missingness policy and provide a complete-data default."""

    def __init__(
        self,
        *,
        compatibility_validator: CorrectionMissingnessCompatibilityValidator
        | None = None,
    ) -> None:
        self._compatibility_validator = (
            compatibility_validator or CorrectionMissingnessCompatibilityValidator()
        )

    def run(
        self,
        *,
        request: BatchCorrectionWorkflowRequest,
    ) -> CorrectionMissingnessPolicy:
        policy = _resolve_upstream_observation_mask_policy(
            request=request,
            policy=request.missingness_policy,
        )
        self._compatibility_validator.run(
            phospho=request.phospho,
            policy=policy,
            allow_complete_matrix_with_observation_mask=(
                request.upstream_observation_mask is not None
            ),
            context="batch-correction workflow missingness validation",
        )
        if policy is None:
            return CorrectionMissingnessPolicy()
        return policy


def _resolve_upstream_observation_mask_policy(
    *,
    request: BatchCorrectionWorkflowRequest,
    policy: CorrectionMissingnessPolicy | None,
) -> CorrectionMissingnessPolicy | None:
    upstream_mask = request.upstream_observation_mask
    if upstream_mask is None:
        return policy

    upstream_observation_mask = _observation_mask_from_upstream_frame(
        upstream_mask,
        phospho=request.phospho,
    )
    if not upstream_observation_mask.originally_missing_cells:
        if (
            policy is not None
            and policy.originally_missing_cells_tracked_by
            is OriginallyMissingCellTracking.OBSERVATION_MASK
            and policy.observation_mask is None
        ):
            return replace(policy, observation_mask=upstream_observation_mask)
        return policy
    if policy is None:
        raise PhosPyInputError(
            "batch-correction workflow missingness validation found "
            "upstream-imputed originally missing cells, but no explicit "
            "missingness policy was provided; upstream-imputed cells cannot be "
            "treated as observed evidence during SPS/RUV-style correction"
        )
    if (
        policy.originally_missing_cells_tracked_by
        is not OriginallyMissingCellTracking.OBSERVATION_MASK
    ):
        raise PhosPyInputError(
            "batch-correction workflow missingness validation found "
            "upstream-imputed originally missing cells, but the requested "
            "missingness policy does not track originally missing cells with an "
            "observation mask; upstream-imputed cells cannot be treated as "
            "observed evidence during SPS/RUV-style correction"
        )
    if policy.row_sample_eligibility_impact in _UNSAFE_UPSTREAM_ELIGIBILITY_IMPACTS:
        raise PhosPyInputError(
            "batch-correction workflow missingness validation found "
            "upstream-imputed originally missing cells, but the requested "
            f"row/sample eligibility policy "
            f"{policy.row_sample_eligibility_impact.value!r} cannot preserve "
            "that distinction during SPS/RUV-style correction"
        )

    resolved_mask = _merge_observation_masks(
        policy_mask=policy.observation_mask,
        upstream_mask=upstream_observation_mask,
    )
    return replace(policy, observation_mask=resolved_mask)


def _observation_mask_from_upstream_frame(
    upstream_mask: pd.DataFrame,
    *,
    phospho: pd.DataFrame,
) -> ObservationMask:
    if not upstream_mask.index.equals(phospho.index):
        raise PhosPyInputError(
            "batch-correction workflow missingness validation upstream "
            "observation mask alignment failed: mask index must match phospho "
            "index so upstream-imputed cells are not treated as observed "
            "evidence during SPS/RUV-style correction"
        )
    if not upstream_mask.columns.equals(phospho.columns):
        raise PhosPyInputError(
            "batch-correction workflow missingness validation upstream "
            "observation mask alignment failed: mask columns must match "
            "phospho columns so upstream-imputed cells are not treated as "
            "observed evidence during SPS/RUV-style correction"
        )
    _require_boolean_mask(upstream_mask)

    observed = upstream_mask.astype(bool)
    feature_ids = tuple(str(value) for value in phospho.index.tolist())
    sample_ids = tuple(str(value) for value in phospho.columns.tolist())
    originally_missing_cells: list[tuple[str, str]] = []
    values = observed.to_numpy(dtype=bool, copy=True)
    for row_position, feature_id in enumerate(feature_ids):
        for column_position, sample_id in enumerate(sample_ids):
            if not bool(values[row_position, column_position]):
                originally_missing_cells.append((feature_id, sample_id))
    try:
        return ObservationMask(
            feature_ids=feature_ids,
            sample_ids=sample_ids,
            originally_missing_cells=tuple(originally_missing_cells),
        )
    except PhosPyInputError as exc:
        raise PhosPyInputError(
            "batch-correction workflow missingness validation could not "
            "serialise the upstream observation mask for provenance; "
            "upstream-imputed cells cannot be treated as observed evidence "
            "during SPS/RUV-style correction"
        ) from exc


def _merge_observation_masks(
    *,
    policy_mask: ObservationMask | None,
    upstream_mask: ObservationMask,
) -> ObservationMask:
    if policy_mask is None:
        return upstream_mask
    if policy_mask.feature_ids != upstream_mask.feature_ids:
        raise PhosPyInputError(
            "batch-correction workflow missingness validation observation mask "
            "alignment failed: policy mask feature_ids must match the upstream "
            "mask and phospho index"
        )
    if policy_mask.sample_ids != upstream_mask.sample_ids:
        raise PhosPyInputError(
            "batch-correction workflow missingness validation observation mask "
            "alignment failed: policy mask sample_ids must match the upstream "
            "mask and phospho columns"
        )
    merged_cells = tuple(
        sorted(
            {
                *policy_mask.originally_missing_cells,
                *upstream_mask.originally_missing_cells,
            }
        )
    )
    return ObservationMask(
        feature_ids=upstream_mask.feature_ids,
        sample_ids=upstream_mask.sample_ids,
        originally_missing_cells=merged_cells,
    )


def _require_boolean_mask(mask: pd.DataFrame) -> None:
    values = mask.to_numpy(dtype="object", copy=True)
    for row_position, row_id in enumerate(mask.index.tolist()):
        for column_position, column_id in enumerate(mask.columns.tolist()):
            value = values[row_position, column_position]
            if isinstance(value, (bool, np.bool_)):
                continue
            raise PhosPyInputError(
                "batch-correction workflow missingness validation upstream "
                "observation mask cannot be serialised as a faithful boolean "
                "mask; invalid value at "
                f"({str(row_id)!r}, {str(column_id)!r})"
            )


__all__ = ["BatchCorrectionWorkflowMissingnessValidator"]
