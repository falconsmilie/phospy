"""Internal interpreter for signalome workflow requests."""

from __future__ import annotations

import pandas as pd

from phospy.api.requests import SignalomeWorkflowRequest
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.errors.workflows import WorkflowStageError
from phospy.workflows.signalome.contracts import ResolvedSignalomeWorkflowRequest


class SignalomeWorkflowInterpreter:
    """Resolve a validated signalome request into execution-ready inputs."""

    _SITE_ID_COLUMN = "site_id"
    _KINASE_COLUMN = "kinase"
    _PROTEIN_COLUMN = "protein_id"
    _GENE_SYMBOL_COLUMN = "gene_symbol"

    def run(
        self, request: SignalomeWorkflowRequest
    ) -> ResolvedSignalomeWorkflowRequest:
        scoring_result = request.kinase_result.scoring_result
        score_matrix = scoring_result.combined_scores
        if score_matrix is None:
            score_matrix = scoring_result.profile_scores
        resolved_score_matrix = self._as_aligned_numeric_frame(
            score_matrix,
            index_name=self._SITE_ID_COLUMN,
            columns_name=self._KINASE_COLUMN,
        )
        resolved_prediction_matrix = self._as_aligned_numeric_frame(
            request.kinase_result.prediction_result.pred_mat,
            index_name=self._SITE_ID_COLUMN,
            columns_name=self._KINASE_COLUMN,
        )
        site_to_protein = self._resolve_site_to_protein(
            dataset=request.kinase_result.dataset,
            site_index=resolved_prediction_matrix.index,
        )
        return ResolvedSignalomeWorkflowRequest(
            dataset=request.kinase_result.dataset,
            kinase_result=request.kinase_result,
            config=request.config,
            score_matrix=resolved_score_matrix,
            prediction_matrix=resolved_prediction_matrix,
            site_to_protein=site_to_protein,
        )

    @staticmethod
    def _as_aligned_numeric_frame(
        frame: pd.DataFrame,
        *,
        index_name: str,
        columns_name: str,
    ) -> pd.DataFrame:
        resolved = frame.astype(float).copy(deep=True)
        resolved.index = pd.Index(resolved.index.astype(str), name=index_name)
        resolved.columns = pd.Index(resolved.columns.astype(str), name=columns_name)
        return resolved

    def _resolve_site_to_protein(
        self,
        *,
        dataset: AnalysisReadyPhosphoDataset,
        site_index: pd.Index,
    ) -> pd.Series:
        resolved = pd.Series(
            [self._protein_from_site_id(site_id) for site_id in site_index],
            index=site_index.copy(),
            dtype=str,
            name=self._PROTEIN_COLUMN,
        )
        metadata = dataset.site_metadata
        if self._GENE_SYMBOL_COLUMN in metadata.columns:
            metadata_genes = metadata.reindex(site_index).loc[
                :, self._GENE_SYMBOL_COLUMN
            ]
            non_empty_mask = metadata_genes.notna() & (
                metadata_genes.astype(str).str.strip() != ""
            )
            if non_empty_mask.any():
                resolved.loc[non_empty_mask] = (
                    metadata_genes.loc[non_empty_mask].astype(str).str.strip()
                )
        if (resolved.astype(str).str.strip() == "").any():
            raise WorkflowStageError(
                "signalome interpreter could not resolve non-empty protein identifiers "
                "for all prediction sites"
            )
        return resolved

    @staticmethod
    def _protein_from_site_id(site_id: object) -> str:
        raw = str(site_id).strip()
        if raw == "":
            return ""
        return raw.split(";", 1)[0].strip()
