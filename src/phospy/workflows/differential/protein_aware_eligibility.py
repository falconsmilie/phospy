"""Protein-aware differential post-computation eligibility reconciliation."""

from __future__ import annotations

from typing import NoReturn

import pandas as pd

from phospy.errors.workflows import WorkflowBoundaryError
from phospy.science.differential.models import (
    DIFFERENTIAL_RESULT_STATUS_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_TESTED,
)
from phospy.science.differential.models.protein_aware import (
    ProteinAwareDifferentialComputationResult,
)
from phospy.workflows.differential.models import (
    DifferentialFeatureEligibilityInputs,
    ProteinAwareDifferentialResolvedInputs,
)


def protein_aware_feature_eligibility_after_computation(
    *,
    resolved_inputs: ProteinAwareDifferentialResolvedInputs,
    computation_result: ProteinAwareDifferentialComputationResult,
) -> DifferentialFeatureEligibilityInputs:
    """Merge protein-aware computation failures into full-index eligibility."""

    feature_metadata = pd.DataFrame(
        resolved_inputs.feature_eligibility_inputs.feature_metadata,
        copy=True,
    )
    result_status = pd.Series(
        resolved_inputs.feature_eligibility_inputs.result_status,
        copy=True,
    ).astype(str)
    if DIFFERENTIAL_RESULT_STATUS_COLUMN not in feature_metadata.columns:
        feature_metadata[DIFFERENTIAL_RESULT_STATUS_COLUMN] = result_status.to_numpy()
    if DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN not in feature_metadata.columns:
        feature_metadata[DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN] = ""

    successful_site_ids = set(computation_result.tested_site_ids)
    for site_id in computation_result.tested_site_ids:
        feature_metadata.loc[site_id, DIFFERENTIAL_RESULT_STATUS_COLUMN] = (
            DIFFERENTIAL_RESULT_STATUS_TESTED
        )
        feature_metadata.loc[site_id, DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN] = ""
        result_status.loc[site_id] = DIFFERENTIAL_RESULT_STATUS_TESTED
        if "protein_aware_tested" in feature_metadata.columns:
            feature_metadata.loc[site_id, "protein_aware_tested"] = True

    failure_diagnostics = computation_result.site_failure_diagnostics
    for site_id, row in failure_diagnostics.iterrows():
        site_key = str(site_id)
        if site_key in successful_site_ids:
            _raise_protein_aware_eligibility_error(
                seam="failure_overlap",
                next_action=(
                    "do not report a phosphosite as both successfully tested and "
                    "withheld by protein-aware computation"
                ),
                details={"site_key": site_key},
            )
        status = str(row[DIFFERENTIAL_RESULT_STATUS_COLUMN])
        reason = str(row[DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN])
        feature_metadata.loc[site_key, DIFFERENTIAL_RESULT_STATUS_COLUMN] = status
        feature_metadata.loc[
            site_key,
            DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN,
        ] = reason
        result_status.loc[site_key] = status
        if "protein_aware_tested" in feature_metadata.columns:
            feature_metadata.loc[site_key, "protein_aware_tested"] = False
        if "protein_aware_failure_message" in feature_metadata.columns:
            feature_metadata.loc[site_key, "protein_aware_failure_message"] = str(
                row.get("failure_message", "")
            )

    return DifferentialFeatureEligibilityInputs(
        feature_metadata=feature_metadata,
        result_status=pd.Series(
            feature_metadata[DIFFERENTIAL_RESULT_STATUS_COLUMN].astype(str),
            index=feature_metadata.index.copy(),
            name=DIFFERENTIAL_RESULT_STATUS_COLUMN,
        ),
        testable_feature_ids=tuple(
            str(value) for value in computation_result.tested_site_ids
        ),
        attach_to_result_tables=True,
    )


def _raise_protein_aware_eligibility_error(
    *,
    seam: str,
    next_action: str,
    details: dict[str, object] | None = None,
) -> NoReturn:
    raise WorkflowBoundaryError(
        seam=f"differential.protein_aware_eligibility.{seam}",
        next_action=next_action,
        details=details,
        message_prefix="differential workflow boundary validation failed",
    )


__all__ = ["protein_aware_feature_eligibility_after_computation"]
