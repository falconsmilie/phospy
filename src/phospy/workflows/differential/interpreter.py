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
from phospy.workflows.differential.models import (
    InterpretedDifferentialAnalysisRequest,
    ValidatedDifferentialAnalysisRequest,
)
from phospy.workflows.differential.provenance import (
    build_differential_policy_provenance,
)


class DifferentialAnalysisInterpreter:
    """Resolve a validated differential request into execution-ready inputs."""

    def run(
        self, request: ValidatedDifferentialAnalysisRequest
    ) -> InterpretedDifferentialAnalysisRequest:
        analysis_sample_ids = request.analysis_sample_ids
        matrix = request.dataset._borrow_phospho_frame().loc[
            :, list(analysis_sample_ids)
        ]
        design_aligned = request.design_matrix.frame
        contrasts_aligned = request.contrast_matrix.frame

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
        policy_provenance = build_differential_policy_provenance(
            request=request,
            design_rank=rank,
            residual_degrees_of_freedom=residual_dof,
        )
        return InterpretedDifferentialAnalysisRequest(
            computation_request=computation_request,
            config=request.config,
            design_rank=rank,
            residual_degrees_of_freedom=residual_dof,
            policy_provenance=policy_provenance,
            workflow_provenance=request.workflow_provenance,
            dataset_preprocessing_report=request.dataset_preprocessing_report,
        )


__all__ = ["DifferentialAnalysisInterpreter"]
