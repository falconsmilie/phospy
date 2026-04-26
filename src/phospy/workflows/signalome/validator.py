"""Internal validator for signalome workflow requests."""

from __future__ import annotations

import numpy as np
import pandas as pd

from phospy.api.requests import SignalomeWorkflowRequest
from phospy.api.results import KinaseWorkflowResult
from phospy.errors.validation import WorkflowValidationError
from phospy.prediction.scoring import select_downstream_score_matrix
from phospy.validation.common.dataframes import (
    require_columns,
    require_dataframe,
    require_exact_index_match,
    require_non_empty_string_column,
    require_numeric_dataframe,
    require_unique_columns,
    require_unique_index,
)
from phospy.validation.workflows.configs import SignalomeConfigValidator

SIGNALOME_PROTEIN_IDENTITY_CONTRACT_NOTE = (
    "Supported signalome execution requires explicit non-empty "
    "site_metadata.protein_id for every interpreted site. "
    "Gene-symbol site-ID prefixes encode site identity, not protein identity, "
    "and are not a fallback substitute. "
    "This is an intentional scientific boundary for protein-aware signalome "
    "grouping and module assignment."
)


class SignalomeWorkflowValidator:
    """Validate `SignalomeWorkflowRequest` before interpretation."""

    def __init__(
        self, *, config_validator: SignalomeConfigValidator | None = None
    ) -> None:
        self._config_validator = config_validator or SignalomeConfigValidator()

    def run(self, request: SignalomeWorkflowRequest) -> SignalomeWorkflowRequest:
        if not isinstance(request, SignalomeWorkflowRequest):
            raise WorkflowValidationError(
                "signalome workflow input must be a SignalomeWorkflowRequest"
            )
        if not isinstance(request.kinase_result, KinaseWorkflowResult):
            raise WorkflowValidationError(
                "signalome workflow request kinase_result must be KinaseWorkflowResult"
            )
        self._config_validator.run(request.config)

        prediction_matrix = require_dataframe(
            request.kinase_result.prediction_result.pred_mat,
            field_name=(
                "signalome workflow request kinase_result.prediction_result.pred_mat"
            ),
            allow_empty=False,
            error_type=WorkflowValidationError,
        )
        require_numeric_dataframe(
            prediction_matrix,
            field_name=(
                "signalome workflow request kinase_result.prediction_result.pred_mat"
            ),
            error_type=WorkflowValidationError,
        )
        self._require_no_missing_or_infinite(
            prediction_matrix,
            field_name=(
                "signalome workflow request kinase_result.prediction_result.pred_mat"
            ),
            allow_missing=True,
        )
        prediction_matrix = require_unique_index(
            prediction_matrix,
            field_name=(
                "signalome workflow request kinase_result.prediction_result.pred_mat"
            ),
            error_type=WorkflowValidationError,
        )
        prediction_matrix = require_unique_columns(
            prediction_matrix,
            field_name=(
                "signalome workflow request kinase_result.prediction_result.pred_mat"
            ),
            error_type=WorkflowValidationError,
        )
        if prediction_matrix.shape[1] == 0:
            raise WorkflowValidationError(
                "signalome workflow request kinase_result.prediction_result.pred_mat "
                "must contain at least one kinase column"
            )

        scoring_result = request.kinase_result.scoring_result
        downstream_score_matrix, downstream_score_source = (
            select_downstream_score_matrix(
                profile_scores=scoring_result.profile_scores,
                rank_weighted_fusion_scores=scoring_result.rank_weighted_fusion_scores,
            )
        )
        score_field_name = (
            "signalome workflow request kinase_result.scoring_result."
            f"{downstream_score_source}"
        )
        score_matrix = require_dataframe(
            downstream_score_matrix,
            field_name=score_field_name,
            allow_empty=False,
            error_type=WorkflowValidationError,
        )
        require_numeric_dataframe(
            score_matrix,
            field_name=score_field_name,
            error_type=WorkflowValidationError,
        )
        self._require_no_missing_or_infinite(
            score_matrix,
            field_name=score_field_name,
            # Correlation-based kinase scoring can legitimately emit missing
            # values (for example, zero-variance denominator collapse).
            # Missingness is preconditioned downstream by the interpreter.
            allow_missing=True,
        )
        score_matrix = require_unique_index(
            score_matrix,
            field_name=score_field_name,
            error_type=WorkflowValidationError,
        )
        score_matrix = require_unique_columns(
            score_matrix,
            field_name=score_field_name,
            error_type=WorkflowValidationError,
        )
        if score_matrix.shape[1] == 0:
            raise WorkflowValidationError(
                f"{score_field_name} must contain at least one kinase column"
            )
        self._require_explicit_site_metadata_protein_identity(
            site_metadata=request.kinase_result.dataset.site_metadata,
            dataset_sites=request.kinase_result.dataset.phospho.index,
        )
        return request

    @staticmethod
    def _require_no_missing_or_infinite(
        frame: pd.DataFrame,
        *,
        field_name: str,
        allow_missing: bool,
    ) -> None:
        if not allow_missing and frame.isna().to_numpy().any():
            raise WorkflowValidationError(
                f"{field_name} must not contain missing values"
            )
        if np.isinf(frame.to_numpy(copy=False)).any():
            raise WorkflowValidationError(
                f"{field_name} must contain finite numeric values"
            )

    @staticmethod
    def _require_explicit_site_metadata_protein_identity(
        *,
        site_metadata: object,
        dataset_sites: pd.Index,
    ) -> None:
        field_name = "signalome workflow request kinase_result.dataset.site_metadata"
        site_metadata_frame = require_dataframe(
            site_metadata,
            field_name=field_name,
            allow_empty=False,
            error_type=WorkflowValidationError,
        )
        site_metadata_frame = require_unique_index(
            site_metadata_frame,
            field_name=field_name,
            error_type=WorkflowValidationError,
        )
        require_exact_index_match(
            left=site_metadata_frame.index,
            right=dataset_sites,
            left_name=f"{field_name}.index",
            right_name="signalome workflow request kinase_result.dataset.phospho.index",
            error_type=WorkflowValidationError,
        )
        try:
            require_columns(
                site_metadata_frame,
                field_name=field_name,
                required_columns=("protein_id",),
                error_type=WorkflowValidationError,
            )
        except WorkflowValidationError as exc:
            raise WorkflowValidationError(
                f"{exc}. {SIGNALOME_PROTEIN_IDENTITY_CONTRACT_NOTE}"
            ) from exc
        try:
            require_non_empty_string_column(
                site_metadata_frame,
                field_name=field_name,
                column_name="protein_id",
                error_type=WorkflowValidationError,
            )
        except WorkflowValidationError as exc:
            raise WorkflowValidationError(
                f"{exc}. {SIGNALOME_PROTEIN_IDENTITY_CONTRACT_NOTE}"
            ) from exc
