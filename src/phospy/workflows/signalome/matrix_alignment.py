"""Numeric coercion and shared-index alignment for signalome matrices."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from phospy.errors.validation import WorkflowValidationError
from phospy.policy_models import DownstreamScoreSource
from phospy.validation.common.dataframes import (
    require_aligned_dataframe_shape,
    require_no_duplicate_labels,
    require_non_empty_index_intersection,
    require_string_index,
)
from phospy.workflows.signalome.boundary_errors import (
    raise_signalome_boundary_error,
    raise_wrapped_signalome_boundary_error,
)
from phospy.workflows.signalome.constants import (
    SIGNALOME_INTERPRETER_KINASE_OVERLAP_SEAM,
    SIGNALOME_INTERPRETER_SITE_ALIGNMENT_SEAM,
)


@dataclass(frozen=True, slots=True)
class SignalomeMatrixAlignment:
    dataset_site_index: pd.Index
    resolved_prediction_matrix: pd.DataFrame
    resolved_downstream_score_matrix: pd.DataFrame
    aligned_prediction_matrix: pd.DataFrame
    aligned_downstream_score_matrix: pd.DataFrame
    aligned_site_index: pd.Index
    aligned_kinase_index: pd.Index


class SignalomeMatrixAligner:
    """Resolve numeric frames and align sites/kinases across signalome inputs."""

    _SITE_ID_COLUMN = "site_id"
    _KINASE_COLUMN = "kinase"

    def run(
        self,
        *,
        dataset_sites: pd.Index,
        prediction_matrix: pd.DataFrame,
        downstream_score_matrix: pd.DataFrame,
        downstream_score_source: DownstreamScoreSource,
    ) -> SignalomeMatrixAlignment:
        score_field_name = (
            "signalome workflow request kinase_result.scoring_result."
            f"{downstream_score_source}"
        )
        resolved_downstream_score_matrix = self._as_aligned_numeric_frame(
            downstream_score_matrix,
            field_name=score_field_name,
            stage_name="signalome.score_matrix_conversion",
            matrix_label="downstream score matrix",
            seam="signalome.interpreter.downstream_score_matrix_conversion",
            next_action=(
                "ensure kinase scoring outputs contain numeric finite values before "
                "running SignalomeWorkflow"
            ),
            index_name=self._SITE_ID_COLUMN,
            columns_name=self._KINASE_COLUMN,
        )
        resolved_prediction_matrix = self._as_aligned_numeric_frame(
            prediction_matrix,
            field_name="signalome workflow request kinase_result.prediction_result.pred_mat",
            stage_name="signalome.prediction_matrix_conversion",
            matrix_label="prediction matrix",
            seam="signalome.interpreter.prediction_matrix_conversion",
            next_action=(
                "ensure kinase prediction outputs contain numeric finite values "
                "before running SignalomeWorkflow"
            ),
            index_name=self._SITE_ID_COLUMN,
            columns_name=self._KINASE_COLUMN,
        )
        dataset_site_index = pd.Index(
            dataset_sites.astype(str),
            name=self._SITE_ID_COLUMN,
        )
        aligned_site_index = self._resolve_shared_site_index(
            dataset_sites=dataset_site_index,
            prediction_sites=resolved_prediction_matrix.index,
            score_sites=resolved_downstream_score_matrix.index,
        )
        aligned_kinase_index = self._resolve_shared_kinase_index(
            prediction_kinases=resolved_prediction_matrix.columns,
            score_kinases=resolved_downstream_score_matrix.columns,
        )
        try:
            aligned_prediction_matrix = resolved_prediction_matrix.reindex(
                index=aligned_site_index,
                columns=aligned_kinase_index,
            )
            aligned_downstream_score_matrix = resolved_downstream_score_matrix.reindex(
                index=aligned_site_index,
                columns=aligned_kinase_index,
            )
            require_aligned_dataframe_shape(
                left=aligned_prediction_matrix,
                right=aligned_downstream_score_matrix,
                left_name="signalome.aligned_prediction_matrix",
                right_name="signalome.aligned_downstream_score_matrix",
                error_type=WorkflowValidationError,
            )
        except (KeyError, TypeError, ValueError) as exc:
            raise_wrapped_signalome_boundary_error(
                stage_name="signalome.index_alignment",
                seam="signalome.interpreter.aligned_matrix_selection",
                field_name=(
                    "signalome workflow request kinase_result.prediction_result."
                    "pred_mat and downstream score matrix"
                ),
                operation="selecting aligned shared site and kinase labels via .loc",
                next_action=(
                    "ensure prediction and downstream scoring tables expose stable "
                    "canonical site IDs and kinase labels"
                ),
                original_error=exc,
                aligned_sites=int(aligned_site_index.size),
                aligned_kinases=int(aligned_kinase_index.size),
            )
        except WorkflowValidationError as exc:
            raise_wrapped_signalome_boundary_error(
                stage_name="signalome.index_alignment",
                seam="signalome.interpreter.aligned_matrix_selection",
                field_name=(
                    "signalome workflow request kinase_result.prediction_result."
                    "pred_mat and downstream score matrix"
                ),
                operation="validating aligned matrix shape after shared-index selection",
                next_action=(
                    "ensure prediction and downstream scoring tables expose identical "
                    "aligned row/column dimensions for shared site and kinase labels"
                ),
                original_error=exc,
                aligned_sites=int(aligned_site_index.size),
                aligned_kinases=int(aligned_kinase_index.size),
            )
        return SignalomeMatrixAlignment(
            dataset_site_index=dataset_site_index,
            resolved_prediction_matrix=resolved_prediction_matrix,
            resolved_downstream_score_matrix=resolved_downstream_score_matrix,
            aligned_prediction_matrix=aligned_prediction_matrix,
            aligned_downstream_score_matrix=aligned_downstream_score_matrix,
            aligned_site_index=aligned_site_index,
            aligned_kinase_index=aligned_kinase_index,
        )

    @staticmethod
    def _as_aligned_numeric_frame(
        frame: pd.DataFrame,
        *,
        field_name: str,
        stage_name: str,
        matrix_label: str,
        seam: str,
        next_action: str,
        index_name: str,
        columns_name: str,
    ) -> pd.DataFrame:
        try:
            require_string_index(
                frame.index,
                field_name=f"{field_name}.index",
                error_type=WorkflowValidationError,
            )
            require_string_index(
                frame.columns,
                field_name=f"{field_name}.columns",
                error_type=WorkflowValidationError,
            )
            require_no_duplicate_labels(
                frame.index,
                field_name=f"{field_name}.index",
                error_type=WorkflowValidationError,
            )
            require_no_duplicate_labels(
                frame.columns,
                field_name=f"{field_name}.columns",
                error_type=WorkflowValidationError,
            )
            resolved = frame.astype(float)
            resolved.index = pd.Index(frame.index, name=index_name)
            resolved.columns = pd.Index(frame.columns, name=columns_name)
            return resolved
        except (TypeError, ValueError, WorkflowValidationError) as exc:
            raise_wrapped_signalome_boundary_error(
                stage_name=stage_name,
                seam=seam,
                field_name=field_name,
                operation=f"converting {matrix_label} to float",
                next_action=next_action,
                original_error=exc,
            )

    def _resolve_shared_site_index(
        self,
        *,
        dataset_sites: pd.Index,
        prediction_sites: pd.Index,
        score_sites: pd.Index,
    ) -> pd.Index:
        dataset_site_index = pd.Index(
            dataset_sites.astype(str), name=self._SITE_ID_COLUMN
        )
        prediction_site_index = pd.Index(
            prediction_sites.astype(str),
            name=self._SITE_ID_COLUMN,
        )
        score_site_index = pd.Index(score_sites.astype(str), name=self._SITE_ID_COLUMN)
        try:
            shared_dataset_prediction = require_non_empty_index_intersection(
                left=dataset_site_index,
                right=prediction_site_index,
                left_name="kinase_result.dataset.phospho.index",
                right_name="kinase_result.prediction_result.pred_mat.index",
                error_type=WorkflowValidationError,
            )
            shared_sites = require_non_empty_index_intersection(
                left=shared_dataset_prediction,
                right=score_site_index,
                left_name="shared dataset/prediction phosphosite IDs",
                right_name="kinase_result.scoring_result.downstream_score_matrix.index",
                error_type=WorkflowValidationError,
            )
        except WorkflowValidationError:
            raise_signalome_boundary_error(
                seam=SIGNALOME_INTERPRETER_SITE_ALIGNMENT_SEAM,
                next_action=(
                    "ensure prediction and scoring outputs share phosphosite IDs with "
                    "kinase_result.dataset.phospho.index"
                ),
                dataset_sites=int(dataset_site_index.size),
                prediction_sites=int(prediction_site_index.size),
                score_sites=int(score_site_index.size),
                shared_sites=0,
            )
        return pd.Index(shared_sites, name=self._SITE_ID_COLUMN)

    def _resolve_shared_kinase_index(
        self,
        *,
        prediction_kinases: pd.Index,
        score_kinases: pd.Index,
    ) -> pd.Index:
        prediction_kinase_index = pd.Index(
            prediction_kinases.astype(str),
            name=self._KINASE_COLUMN,
        )
        score_kinase_index = pd.Index(
            score_kinases.astype(str), name=self._KINASE_COLUMN
        )
        try:
            shared_kinases = require_non_empty_index_intersection(
                left=prediction_kinase_index,
                right=score_kinase_index,
                left_name="kinase_result.prediction_result.pred_mat.columns",
                right_name="kinase_result.scoring_result.downstream_score_matrix.columns",
                error_type=WorkflowValidationError,
            )
        except WorkflowValidationError:
            raise_signalome_boundary_error(
                seam=SIGNALOME_INTERPRETER_KINASE_OVERLAP_SEAM,
                next_action=(
                    "rerun kinase workflow so scoring_result and "
                    "prediction_result are generated from the same kinase lane"
                ),
                prediction_kinases=int(prediction_kinase_index.size),
                score_kinases=int(score_kinase_index.size),
                shared_kinases=0,
            )
        return pd.Index(shared_kinases, name=self._KINASE_COLUMN)


__all__ = [
    "SignalomeMatrixAligner",
    "SignalomeMatrixAlignment",
]
