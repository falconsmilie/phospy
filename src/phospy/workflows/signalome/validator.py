"""Internal validator for signalome workflow requests."""

from __future__ import annotations

import pandas as pd

from phospy.api.configs.localisation import LocalisationRequirement
from phospy.api.requests import SignalomeWorkflowRequest
from phospy.api.results import KinaseWorkflowResult
from phospy.errors.validation import WorkflowValidationError
from phospy.prediction.scoring import select_downstream_score_matrix
from phospy.sites.identity import build_phosphosite_identity
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

SIGNALOME_PROTEIN_IDENTITY_CONTRACT_NOTE = (
    "Signalome execution requires an explicit site_metadata.protein_id column. "
    "That column must contain explicit non-missing protein_id values for retained sites. "
    "Protein identifiers are resolved for retained signalome sites after "
    "downstream-score preconditioning. "
    "Gene-symbol site-ID prefixes encode site identity, not protein identity."
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
        prediction_result = request.kinase_result.prediction_result
        scoring_result = request.kinase_result.scoring_result
        reject_mixed_total_protein_quantitative_meaning(
            dataset=dataset,
            allow_mixed=config.validation.allow_mixed_total_protein_quantitative_meaning,
            context="signalome workflow request kinase_result.dataset",
        )

        prediction_matrix = require_dataframe(
            prediction_result._borrow_pred_mat_frame(),
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

        downstream_score_matrix, downstream_score_source = (
            select_downstream_score_matrix(
                profile_scores=scoring_result._borrow_profile_scores_frame(),
                rank_weighted_fusion_scores=(
                    scoring_result._borrow_rank_weighted_fusion_scores_frame()
                ),
            )
        )
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
        self._require_explicit_site_metadata_protein_identity(
            site_metadata=dataset._borrow_site_metadata_frame(),
            dataset_sites=dataset._borrow_phospho_frame().index,
            localisation_requirement=config.validation.localisation_requirement,
        )
        return request

    @staticmethod
    def _require_explicit_site_metadata_protein_identity(
        *,
        site_metadata: object,
        dataset_sites: pd.Index,
        localisation_requirement: LocalisationRequirement,
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
            for site_id in site_metadata_frame.index.tolist():
                _ = build_phosphosite_identity(
                    display_id=site_id,
                    gene_symbol=site_metadata_frame.at[site_id, "gene_symbol"],
                    site=site_metadata_frame.at[site_id, "site"],
                    allow_opaque_site_values=True,
                    protein_id=(
                        site_metadata_frame.at[site_id, "protein_id"]
                        if "protein_id" in site_metadata_frame.columns
                        else None
                    ),
                    protein_accession=(
                        site_metadata_frame.at[site_id, "protein_accession"]
                        if "protein_accession" in site_metadata_frame.columns
                        else None
                    ),
                    field_name=f"{field_name}[{site_id!r}]",
                    error_type=WorkflowValidationError,
                )
            enforce_required_non_empty_string_column(
                site_metadata=site_metadata_frame,
                field_name=field_name,
                workflow_name="signalome workflow request",
                column_name="protein_id",
                error_type=WorkflowValidationError,
            )
            enforce_localisation_requirement(
                site_metadata=site_metadata_frame,
                field_name=field_name,
                workflow_name="signalome workflow request",
                requirement=localisation_requirement,
                error_type=WorkflowValidationError,
            )
        except WorkflowValidationError as exc:
            raise WorkflowValidationError(
                f"{exc}. {SIGNALOME_PROTEIN_IDENTITY_CONTRACT_NOTE}"
            ) from exc
