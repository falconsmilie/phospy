"""Internal interpreter for differential workflow requests."""

from __future__ import annotations

from typing import cast

import numpy as np
import pandas as pd

from phospy.errors.workflows import WorkflowBoundaryError
from phospy.science.differential.models import (
    ContrastMatrix,
    DesignMatrix,
)
from phospy.science.differential.models import (
    DifferentialAnalysisRequest as DifferentialComputationRequest,
)
from phospy.science.sites.validation import require_site_key_series
from phospy.validation.workflows.differential import (
    ExperimentalDesignContractValidator,
)
from phospy.workflows.differential.models import (
    InterpretedDifferentialAnalysisRequest,
    ValidatedDifferentialAnalysisRequest,
)
from phospy.workflows.differential.provenance import (
    build_differential_policy_provenance,
)
from phospy.workflows.differential.replicates import (
    TechnicalReplicateAggregator,
)


class DifferentialAnalysisInterpreter:
    """Resolve a validated differential request into execution-ready inputs."""

    def __init__(
        self,
        *,
        design_validator: ExperimentalDesignContractValidator | None = None,
        technical_replicate_aggregator: TechnicalReplicateAggregator | None = None,
    ) -> None:
        self._design_validator = (
            design_validator or ExperimentalDesignContractValidator()
        )
        self._technical_replicate_aggregator = (
            technical_replicate_aggregator or TechnicalReplicateAggregator()
        )

    def run(
        self, request: ValidatedDifferentialAnalysisRequest
    ) -> InterpretedDifferentialAnalysisRequest:
        aggregation_plan = request.technical_replicate_aggregation_plan
        resolved_dataset = request.dataset
        resolved_design = request.design
        resolved_contrasts = request.contrasts
        resolved_analysis_sample_ids = request.analysis_sample_ids
        resolved_design_matrix = request.design_matrix
        resolved_contrast_matrix = request.contrast_matrix
        resolved_workflow_provenance = request.workflow_provenance

        if aggregation_plan is not None and aggregation_plan.requires_aggregation:
            technical_replicate_resolution = self._technical_replicate_aggregator.run(
                dataset=request.dataset,
                design=request.design,
                aggregation_plan=aggregation_plan,
            )
            resolved_dataset = technical_replicate_resolution.dataset
            resolved_design = technical_replicate_resolution.design
            resolved_workflow_provenance = (
                technical_replicate_resolution.workflow_provenance
            )
            resolved_design_contract = self._design_validator.run(
                dataset=resolved_dataset,
                design=resolved_design,
                contrasts=request.contrasts,
                allow_design_subset=request.config.allow_design_subset,
                minimum_condition_replicates=request.config.minimum_condition_replicates,
            )
            resolved_contrasts = resolved_design_contract.contrasts
            resolved_analysis_sample_ids = resolved_design_contract.analysis_sample_ids
            resolved_design_matrix = DesignMatrix(resolved_design_contract.design_frame)
            resolved_contrast_matrix = ContrastMatrix(
                resolved_design_contract.contrast_frame
            )

        analysis_sample_ids = resolved_analysis_sample_ids
        matrix = resolved_dataset._borrow_phospho_frame().loc[
            :, list(analysis_sample_ids)
        ]
        matrix = _prefer_site_key_index_for_differential_results(
            matrix=matrix,
            site_metadata=resolved_dataset._borrow_site_metadata_frame(),  # pyright: ignore[reportPrivateUsage] - workflow boundary reads trusted internal dataset snapshots
        )
        design_aligned = resolved_design_matrix.frame
        contrasts_aligned = resolved_contrast_matrix.frame

        matrix_samples = pd.Index(matrix.columns)
        design_samples = pd.Index(design_aligned.index)
        if not matrix_samples.equals(design_samples):
            raise WorkflowBoundaryError(
                seam="differential.interpreter.sample_label_alignment",
                next_action=(
                    "ensure validated design sample IDs exactly match the "
                    "analysis matrix sample order"
                ),
                details={
                    "matrix_samples": matrix_samples.astype(str).tolist(),
                    "design_samples": design_samples.astype(str).tolist(),
                },
                message_prefix="differential workflow boundary validation failed",
            )
        matrix_aligned = cast(pd.DataFrame, matrix.copy(deep=True))
        result_identity_metadata = _build_result_identity_metadata(
            site_metadata=resolved_dataset._borrow_site_metadata_frame(),  # pyright: ignore[reportPrivateUsage] - workflow boundary reads trusted internal dataset snapshots
            expected_index=matrix_aligned.index,
        )

        design_values = design_aligned.to_numpy(dtype=float)
        design_shape = cast(tuple[int, int], design_values.shape)
        sample_count = int(design_shape[0])
        coefficient_count = int(design_shape[1])
        rank = int(np.linalg.matrix_rank(design_values))
        if rank < coefficient_count:
            raise WorkflowBoundaryError(
                seam="differential.interpreter.design_rank",
                next_action=(
                    "remove collinear design terms or simplify the design matrix so "
                    "it is full column rank"
                ),
                details={"rank": rank, "columns": coefficient_count},
                message_prefix="differential workflow boundary validation failed",
            )

        residual_dof = float(sample_count - rank)
        if residual_dof <= 0.0:
            raise WorkflowBoundaryError(
                seam="differential.interpreter.residual_dof",
                next_action=(
                    "increase sample count or reduce design terms so residual "
                    "degrees of freedom stays positive"
                ),
                details={
                    "samples": sample_count,
                    "rank": rank,
                    "residual_dof": residual_dof,
                },
                message_prefix="differential workflow boundary validation failed",
            )

        xtx_inv = np.linalg.pinv(design_values.T @ design_values)
        contrast_values = contrasts_aligned.to_numpy(dtype=float)
        contrast_covariance = contrast_values.T @ xtx_inv @ contrast_values
        contrast_scale = np.sqrt(np.diag(contrast_covariance))
        if np.any(~np.isfinite(contrast_scale)) or np.any(contrast_scale <= 0.0):
            raise WorkflowBoundaryError(
                seam="differential.interpreter.non_estimable_contrast",
                next_action=(
                    "update contrasts so each column is estimable under the "
                    "resolved design matrix"
                ),
                details={
                    "contrast_names": contrasts_aligned.columns.astype(str).tolist(),
                },
                message_prefix="differential workflow boundary validation failed",
            )

        computation_request = DifferentialComputationRequest(
            matrix=cast(pd.DataFrame, matrix_aligned),
            design=DesignMatrix(cast(pd.DataFrame, design_aligned.copy(deep=True))),
            contrasts=ContrastMatrix(
                cast(pd.DataFrame, contrasts_aligned.copy(deep=True))
            ),
            empirical_bayes=request.config.empirical_bayes,
        )
        provenance_request = ValidatedDifferentialAnalysisRequest(
            dataset=resolved_dataset,
            design=resolved_design,
            contrasts=resolved_contrasts,
            analysis_sample_ids=resolved_analysis_sample_ids,
            design_matrix=resolved_design_matrix,
            contrast_matrix=resolved_contrast_matrix,
            config=request.config,
            technical_replicate_aggregation_plan=aggregation_plan,
            workflow_provenance=resolved_workflow_provenance,
            dataset_preprocessing_report=resolved_dataset.preprocessing_report,
        )
        policy_provenance = build_differential_policy_provenance(
            request=provenance_request,
            design_rank=rank,
            residual_degrees_of_freedom=residual_dof,
        )
        return InterpretedDifferentialAnalysisRequest(
            computation_request=computation_request,
            result_identity_metadata=result_identity_metadata,
            config=request.config,
            design_rank=rank,
            residual_degrees_of_freedom=residual_dof,
            policy_provenance=policy_provenance,
            workflow_provenance=resolved_workflow_provenance,
            dataset_preprocessing_report=resolved_dataset.preprocessing_report,
        )


def _build_result_identity_metadata(
    *,
    site_metadata: pd.DataFrame,
    expected_index: pd.Index,
) -> pd.DataFrame:
    required_columns = ("site_key", "display_id", "gene_symbol", "site")
    optional_columns = (
        "protein_id",
        "protein_accession",
        "isoform_id",
    )
    missing = [
        column for column in required_columns if column not in site_metadata.columns
    ]
    if missing:
        joined = ", ".join(missing)
        raise WorkflowBoundaryError(
            seam="differential.interpreter.result_identity_columns",
            next_action=(
                "ensure resolved dataset.site_metadata includes site_key, "
                "display_id, gene_symbol, and site"
            ),
            details={"missing_columns": missing},
            message_prefix=(
                "differential workflow boundary validation failed: "
                f"missing required result identity columns: {joined}"
            ),
        )
    try:
        aligned = cast(pd.DataFrame, site_metadata.loc[expected_index].copy(deep=True))
    except KeyError as exc:
        raise WorkflowBoundaryError(
            seam="differential.interpreter.result_identity_alignment",
            next_action=(
                "ensure resolved dataset.site_metadata is indexed by the exact "
                "site_key labels used by the differential matrix"
            ),
            details={"expected_index_count": int(expected_index.size)},
            message_prefix="differential workflow boundary validation failed",
        ) from exc
    if not aligned.index.equals(expected_index):
        raise WorkflowBoundaryError(
            seam="differential.interpreter.result_identity_alignment",
            next_action=(
                "ensure resolved dataset.site_metadata index order exactly matches "
                "the differential matrix site_key order"
            ),
            details={
                "expected_index_count": int(expected_index.size),
                "actual_index_count": int(aligned.index.size),
            },
            message_prefix="differential workflow boundary validation failed",
        )
    identity = cast(pd.DataFrame, aligned.copy(deep=True))
    site_key_values = identity.loc[:, "site_key"].fillna("").astype(str).str.strip()
    try:
        require_site_key_series(
            site_key_values,
            field_name=("differential workflow request dataset.site_metadata.site_key"),
            error_type=ValueError,
        )
    except ValueError as exc:
        raise WorkflowBoundaryError(
            seam="differential.interpreter.result_identity_site_key",
            next_action=(
                "ensure resolved dataset.site_metadata.site_key contains valid "
                "protein-scoped site keys"
            ),
            details={"error": str(exc)},
            message_prefix="differential workflow boundary validation failed",
        ) from exc
    expected_values = expected_index.astype(str).tolist()
    if site_key_values.tolist() != expected_values:
        raise WorkflowBoundaryError(
            seam="differential.interpreter.result_identity_site_key",
            next_action=(
                "ensure resolved dataset.site_metadata.site_key exactly matches "
                "the differential matrix site_key index"
            ),
            details={"expected_index_count": int(expected_index.size)},
            message_prefix="differential workflow boundary validation failed",
        )
    display_id_values = identity.loc[:, "display_id"].fillna("").astype(str).str.strip()
    if (display_id_values == "").any():
        raise WorkflowBoundaryError(
            seam="differential.interpreter.result_identity_display_id",
            next_action=(
                "ensure resolved dataset.site_metadata.display_id is present for "
                "every differential matrix row"
            ),
            details={"empty_display_id_count": int((display_id_values == "").sum())},
            message_prefix="differential workflow boundary validation failed",
        )
    identity.index = pd.Index(site_key_values.tolist(), name="site_key")
    identity.loc[:, "site_key"] = site_key_values.tolist()
    identity.loc[:, "display_id"] = display_id_values.tolist()
    selected_columns = ("site_key", "display_id", "gene_symbol", "site") + tuple(
        column for column in optional_columns if column in aligned.columns
    )
    return cast(pd.DataFrame, identity.loc[:, list(selected_columns)])


def _prefer_site_key_index_for_differential_results(
    *,
    matrix: pd.DataFrame,
    site_metadata: pd.DataFrame,
) -> pd.DataFrame:
    if "site_key" not in site_metadata.columns:
        raise WorkflowBoundaryError(
            seam="differential.interpreter.result_identity_columns",
            next_action=(
                "ensure resolved dataset.site_metadata includes site_key before "
                "differential execution"
            ),
            message_prefix="differential workflow boundary validation failed",
        )
    aligned_metadata = site_metadata.reindex(matrix.index)
    site_keys = aligned_metadata.loc[:, "site_key"].fillna("").astype(str).str.strip()
    try:
        require_site_key_series(
            site_keys,
            field_name="differential workflow request dataset.site_metadata.site_key",
            error_type=ValueError,
        )
    except ValueError as exc:
        raise WorkflowBoundaryError(
            seam="differential.interpreter.result_identity_site_key",
            next_action=(
                "ensure resolved dataset.site_metadata.site_key contains valid "
                "protein-scoped site keys"
            ),
            details={"error": str(exc)},
            message_prefix="differential workflow boundary validation failed",
        ) from exc
    matrix_index_values = matrix.index.astype(str).tolist()
    if site_keys.tolist() != matrix_index_values:
        raise WorkflowBoundaryError(
            seam="differential.interpreter.result_identity_site_key",
            next_action=(
                "ensure the analysis-ready dataset uses site_key as the phospho "
                "matrix row index"
            ),
            details={"matrix_index_count": int(matrix.index.size)},
            message_prefix="differential workflow boundary validation failed",
        )
    remapped = matrix.copy(deep=True)
    remapped.index = pd.Index(site_keys.tolist(), name="site_key")
    return remapped


__all__ = ["DifferentialAnalysisInterpreter"]
