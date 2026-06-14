"""Internal validator for signalome workflow requests."""

from __future__ import annotations

import pandas as pd

from phospy.contracts.configs.localisation import LocalisationRequirement
from phospy.contracts.requests import SignalomeWorkflowRequest
from phospy.contracts.results import KinaseWorkflowResult
from phospy.errors.validation import WorkflowValidationError
from phospy.science.datasets.internal_view import DatasetInternalView
from phospy.science.prediction.internal_view import (
    KinasePredictionInternalView,
    KinaseScoringInternalView,
)
from phospy.science.scoring.policy_models import DownstreamScoreSource
from phospy.validation.common.dataframes import (
    require_dataframe,
    require_exact_index_match,
    require_finite_numeric_dataframe,
    require_non_empty_dataframe,
    require_numeric_dataframe,
    require_unique_columns,
    require_unique_index,
)
from phospy.validation.datasets.site_metadata import (
    enforce_localisation_requirement,
    enforce_required_non_empty_string_column,
)
from phospy.validation.workflows.configs import (
    SignalomeConfigValidator,
    reject_mixed_total_protein_quantitative_meaning,
)
from phospy.validation.workflows.identity import (
    SIGNALOME_IDENTITY_CONTRACT,
    enforce_workflow_site_identity_contract,
)

SIGNALOME_PROTEIN_GROUPING_METADATA_NOTE = (
    "Signalome uses dataset.site_metadata.protein_id as algorithm-specific "
    "protein grouping metadata. This grouping field is separate from the "
    "dataset-level protein-scoped row identity contract based on site_key, "
    "display_id, organism, protein_namespace, protein_identifier, gene_symbol, "
    "site, and site_sequence where required. "
    "Signalome does not infer protein grouping from gene_symbol or display_id "
    "and does not repair invalid site_key identity."
)

SIGNALOME_SITE_IDENTITY_CONTRACT_NOTE = (
    "Base dataset phosphosite identity must already be valid and site_key-indexed "
    "before signalome validation; signalome does not repair weak dataset identity."
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
        config = self._config_validator.run(request.config)
        dataset = request.kinase_result.dataset
        dataset_view = DatasetInternalView(dataset)
        prediction_result = request.kinase_result.prediction_result
        scoring_result = request.kinase_result.scoring_result
        prediction_view = KinasePredictionInternalView(prediction_result)
        scoring_view = KinaseScoringInternalView(scoring_result)
        reject_mixed_total_protein_quantitative_meaning(
            dataset=dataset,
            allow_mixed=config.validation.allow_mixed_total_protein_quantitative_meaning,
            context="signalome workflow request kinase_result.dataset",
        )

        prediction_matrix = require_dataframe(
            prediction_view.pred_mat,
            field_name=(
                "signalome workflow request kinase_result.prediction_result.pred_mat"
            ),
            allow_empty=True,
            error_type=WorkflowValidationError,
        )
        require_non_empty_dataframe(
            prediction_matrix,
            field_name=(
                "signalome workflow request kinase_result.prediction_result.pred_mat"
            ),
            error_type=WorkflowValidationError,
        )
        require_numeric_dataframe(
            prediction_matrix,
            field_name=(
                "signalome workflow request kinase_result.prediction_result.pred_mat"
            ),
            error_type=WorkflowValidationError,
        )
        require_finite_numeric_dataframe(
            prediction_matrix,
            field_name=(
                "signalome workflow request kinase_result.prediction_result.pred_mat"
            ),
            allow_missing=True,
            error_type=WorkflowValidationError,
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

        downstream_score_source = DownstreamScoreSource.parse(
            scoring_result.score_source,
            field_name="signalome workflow request downstream score source",
        )
        downstream_score_matrix = scoring_view.authoritative_scores
        score_field_name = (
            "signalome workflow request kinase_result.scoring_result."
            f"{downstream_score_source}"
        )
        score_matrix = require_dataframe(
            downstream_score_matrix,
            field_name=score_field_name,
            allow_empty=True,
            error_type=WorkflowValidationError,
        )
        require_non_empty_dataframe(
            score_matrix,
            field_name=score_field_name,
            error_type=WorkflowValidationError,
        )
        require_numeric_dataframe(
            score_matrix,
            field_name=score_field_name,
            error_type=WorkflowValidationError,
        )
        require_finite_numeric_dataframe(
            score_matrix,
            field_name=score_field_name,
            # Correlation-based kinase scoring can legitimately emit missing
            # values (for example, zero-variance denominator collapse).
            # Missingness is preconditioned downstream by the interpreter.
            allow_missing=True,
            error_type=WorkflowValidationError,
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
        self._require_site_identity_and_protein_grouping_metadata(
            site_metadata=dataset_view.site_metadata,
            dataset_sites=dataset_view.phospho.index,
            prediction_sites=prediction_matrix.index,
            score_sites=score_matrix.index,
            score_field_name=score_field_name,
            localisation_requirement=config.validation.localisation_requirement,
            allow_opaque_site_values=dataset.opaque_site_values_allowed,
        )
        return request

    @staticmethod
    def _require_site_identity_and_protein_grouping_metadata(
        *,
        site_metadata: object,
        dataset_sites: pd.Index,
        prediction_sites: pd.Index,
        score_sites: pd.Index,
        score_field_name: str,
        localisation_requirement: LocalisationRequirement,
        allow_opaque_site_values: bool,
    ) -> None:
        field_name = "signalome workflow request kinase_result.dataset.site_metadata"
        site_metadata_frame = require_dataframe(
            site_metadata,
            field_name=field_name,
            allow_empty=False,
            error_type=WorkflowValidationError,
        )
        try:
            enforce_workflow_site_identity_contract(
                site_metadata=site_metadata_frame,
                expected_index=dataset_sites,
                expected_index_field_name=(
                    "signalome workflow request kinase_result.dataset.phospho.index"
                ),
                field_name=field_name,
                contract=SIGNALOME_IDENTITY_CONTRACT,
                error_type=WorkflowValidationError,
                allow_opaque_site_values=allow_opaque_site_values,
            )
        except WorkflowValidationError as exc:
            raise WorkflowValidationError(
                f"{exc}. {SIGNALOME_SITE_IDENTITY_CONTRACT_NOTE}"
            ) from exc
        enforce_localisation_requirement(
            site_metadata=site_metadata_frame,
            field_name=field_name,
            workflow_name="signalome workflow request",
            requirement=localisation_requirement,
            error_type=WorkflowValidationError,
        )
        SignalomeWorkflowValidator._require_signalome_protein_grouping_metadata(
            site_metadata=site_metadata_frame,
            field_name=field_name,
        )
        require_exact_index_match(
            left=prediction_sites,
            right=site_metadata_frame.index,
            left_name=(
                "signalome workflow request "
                "kinase_result.prediction_result.pred_mat.index"
            ),
            right_name=f"{field_name}.index",
            error_type=WorkflowValidationError,
        )
        require_exact_index_match(
            left=score_sites,
            right=site_metadata_frame.index,
            left_name=f"{score_field_name}.index",
            right_name=f"{field_name}.index",
            error_type=WorkflowValidationError,
        )

    @staticmethod
    def _require_signalome_protein_grouping_metadata(
        *,
        site_metadata: pd.DataFrame,
        field_name: str,
    ) -> None:
        try:
            enforce_required_non_empty_string_column(
                site_metadata=site_metadata,
                field_name=field_name,
                workflow_name="signalome protein grouping metadata",
                column_name="protein_id",
                error_type=WorkflowValidationError,
            )
        except WorkflowValidationError as exc:
            raise WorkflowValidationError(
                "Missing signalome protein grouping metadata: protein_id; "
                "signalome protein grouping metadata requirement failed: "
                f"{exc}. {SIGNALOME_PROTEIN_GROUPING_METADATA_NOTE}"
            ) from exc
