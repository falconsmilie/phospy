"""Internal validator for signalome workflow requests."""

from __future__ import annotations

import numpy as np
import pandas as pd

from phospy.api.configs import SignalomeConfig
from phospy.api.requests import SignalomeWorkflowRequest
from phospy.api.results import SimpleKinaseWorkflowResult
from phospy.errors.validation import WorkflowValidationError


class SignalomeWorkflowValidator:
    """Validate `SignalomeWorkflowRequest` before interpretation."""

    _MIN_CUTOFF = 0.0
    _MAX_CUTOFF = 1.0

    def run(self, request: SignalomeWorkflowRequest) -> SignalomeWorkflowRequest:
        if not isinstance(request, SignalomeWorkflowRequest):
            raise WorkflowValidationError(
                "signalome workflow input must be a SignalomeWorkflowRequest"
            )
        if not isinstance(request.kinase_result, SimpleKinaseWorkflowResult):
            raise WorkflowValidationError(
                "signalome workflow request kinase_result must be SimpleKinaseWorkflowResult"
            )
        if not isinstance(request.config, SignalomeConfig):
            raise WorkflowValidationError(
                "signalome workflow request config must be SignalomeConfig"
            )
        cutoff = request.config.signalome_cutoff
        if isinstance(cutoff, bool) or not isinstance(cutoff, (int, float)):
            raise WorkflowValidationError(
                "signalome workflow request config.signalome_cutoff must be a float "
                "between 0.0 and 1.0"
            )
        if not self._MIN_CUTOFF <= float(cutoff) <= self._MAX_CUTOFF:
            raise WorkflowValidationError(
                "signalome workflow request config.signalome_cutoff must be between "
                "0.0 and 1.0"
            )

        prediction_matrix = self._validated_numeric_matrix(
            request.kinase_result.prediction_result.pred_mat,
            context="kinase_result.prediction_result.pred_mat",
        )
        if prediction_matrix.shape[1] == 0:
            raise WorkflowValidationError(
                "signalome workflow request kinase_result.prediction_result.pred_mat "
                "must contain at least one kinase column"
            )

        scoring_result = request.kinase_result.scoring_result
        score_matrix_source = scoring_result.combined_scores
        if score_matrix_source is None:
            score_matrix_source = scoring_result.profile_scores
        score_matrix = self._validated_numeric_matrix(
            score_matrix_source,
            context="kinase_result.scoring_result.combined_scores",
        )
        if score_matrix.shape[1] == 0:
            raise WorkflowValidationError(
                "signalome workflow request kinase_result.scoring_result.combined_scores "
                "must contain at least one kinase column"
            )

        dataset_index = pd.Index(
            request.kinase_result.dataset.phospho.index.astype(str),
            name="site_id",
        )
        self._validate_site_alignment(
            expected_index=dataset_index,
            observed_index=prediction_matrix.index,
            context="kinase_result.prediction_result.pred_mat",
        )
        self._validate_site_alignment(
            expected_index=dataset_index,
            observed_index=score_matrix.index,
            context="kinase_result.scoring_result.combined_scores",
        )
        self._validate_kinase_compatibility(
            score_kinases=score_matrix.columns,
            prediction_kinases=prediction_matrix.columns,
        )
        return request

    @staticmethod
    def _validated_numeric_matrix(matrix: object, *, context: str) -> pd.DataFrame:
        if not isinstance(matrix, pd.DataFrame):
            raise WorkflowValidationError(
                f"signalome workflow request {context} must be a pandas DataFrame"
            )
        if matrix.empty:
            raise WorkflowValidationError(
                f"signalome workflow request {context} must not be empty"
            )
        try:
            numeric_matrix = matrix.astype(float)
        except (TypeError, ValueError) as exc:
            raise WorkflowValidationError(
                f"signalome workflow request {context} must contain numeric values"
            ) from exc
        if not np.isfinite(numeric_matrix.to_numpy(dtype=float, copy=False)).all():
            raise WorkflowValidationError(
                f"signalome workflow request {context} must contain finite numeric values"
            )
        return numeric_matrix

    @staticmethod
    def _validate_site_alignment(
        *,
        expected_index: pd.Index,
        observed_index: pd.Index,
        context: str,
    ) -> None:
        observed = pd.Index(observed_index.astype(str), name="site_id")
        if list(expected_index) == list(observed):
            return
        observed_set = set(observed.tolist())
        missing = [site_id for site_id in expected_index if site_id not in observed_set]
        if missing:
            preview = ", ".join(missing[:3])
            suffix = "..." if len(missing) > 3 else ""
            raise WorkflowValidationError(
                f"signalome workflow request {context} is missing dataset sites: "
                f"{preview}{suffix}"
            )
        raise WorkflowValidationError(
            f"signalome workflow request {context} site order must match "
            "kinase_result.dataset.phospho index"
        )

    @staticmethod
    def _validate_kinase_compatibility(
        *,
        score_kinases: pd.Index,
        prediction_kinases: pd.Index,
    ) -> None:
        score_set = set(score_kinases.astype(str).tolist())
        missing = sorted(
            {
                str(kinase)
                for kinase in prediction_kinases.astype(str).tolist()
                if str(kinase) not in score_set
            }
        )
        if not missing:
            return
        preview = ", ".join(missing[:3])
        suffix = "..." if len(missing) > 3 else ""
        raise WorkflowValidationError(
            "signalome workflow request kinase_result.prediction_result.pred_mat "
            "contains kinases absent from scoring outputs: "
            f"{preview}{suffix}"
        )
