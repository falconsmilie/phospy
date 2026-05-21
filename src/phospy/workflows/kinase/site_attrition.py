"""Compose compact, user-facing kinase site attrition summaries."""

from __future__ import annotations

import pandas as pd

from phospy.contracts.results import (
    KinaseWorkflowPreprocessingAttritionSummary,
    KinaseWorkflowScoringAttritionSummary,
    KinaseWorkflowSiteAttritionSummary,
)
from phospy.errors.validation import WorkflowValidationError
from phospy.science.datasets.models import (
    PreprocessingSiteAttritionSummary,
    SiteSequenceResolutionReport,
)
from phospy.science.prediction.models import KinasePredictionResult
from phospy.science.prediction.sequence_validation import (
    SEQUENCE_VALIDATION_STATUS_INVALID_SITE_ID,
)
from phospy.science.sites.identifiers import parse_canonical_site_identifier
from phospy.workflows.kinase.component_models import KinaseScoringRunResult
from phospy.workflows.kinase.contracts import ResolvedKinaseWorkflowRequest


class _WorkflowIdentifierError(ValueError):
    """Internal exception bridge for identifier parsing checks."""


class KinaseSiteAttritionSummaryComposer:
    """Compose preprocessing and workflow-owned attrition counters."""

    def run(
        self,
        *,
        request: ResolvedKinaseWorkflowRequest,
        scoring_execution: KinaseScoringRunResult,
        prediction_result: KinasePredictionResult,
        activity_enabled: bool,
    ) -> KinaseWorkflowSiteAttritionSummary:
        preprocessing_summary = self._resolve_preprocessing_summary(request)
        scoring_summary = self._resolve_scoring_summary(
            request=request,
            scoring_execution=scoring_execution,
            prediction_result=prediction_result,
            activity_enabled=activity_enabled,
        )
        return KinaseWorkflowSiteAttritionSummary(
            preprocessing=preprocessing_summary,
            scoring=scoring_summary,
        )

    @staticmethod
    def _resolve_preprocessing_summary(
        request: ResolvedKinaseWorkflowRequest,
    ) -> KinaseWorkflowPreprocessingAttritionSummary:
        report = request.dataset.preprocessing_report
        if report is None:
            output_rows = int(request.dataset.phospho.shape[0])
            return KinaseWorkflowPreprocessingAttritionSummary(
                input_rows=output_rows,
                rows_removed_during_preprocessing=0,
                rows_removed_invalid_or_missing_site_identifiers=0,
                duplicate_sites_merged_or_resolved=0,
                output_rows=output_rows,
                sequence_complete_sites=output_rows,
            )
        summary: PreprocessingSiteAttritionSummary = report.site_attrition_summary()
        sequence_summary: SiteSequenceResolutionReport | None = (
            report.site_sequence_resolution_summary()
        )
        sequence_complete_sites = (
            int(sequence_summary.final_sequence_complete_sites)
            if sequence_summary is not None
            else int(summary.output_rows)
        )
        return KinaseWorkflowPreprocessingAttritionSummary(
            input_rows=int(summary.input_rows),
            rows_removed_during_preprocessing=int(
                summary.rows_removed_during_preprocessing
            ),
            rows_removed_invalid_or_missing_site_identifiers=int(
                summary.rows_removed_invalid_or_missing_site_identifiers
            ),
            duplicate_sites_merged_or_resolved=int(
                summary.duplicate_sites_merged_or_resolved
            ),
            output_rows=int(summary.output_rows),
            sequence_complete_sites=sequence_complete_sites,
        )

    def _resolve_scoring_summary(
        self,
        *,
        request: ResolvedKinaseWorkflowRequest,
        scoring_execution: KinaseScoringRunResult,
        prediction_result: KinasePredictionResult,
        activity_enabled: bool,
    ) -> KinaseWorkflowScoringAttritionSummary:
        final_quantitative_sites = int(request.dataset.phospho.shape[0])
        motif_sequence_validation = (
            scoring_execution.scoring_result.motif_sequence_validation
        )
        sequence_supported_sites = int(len(request.scoring_site_index))
        motif_valid_sites = int(
            sequence_supported_sites
            if motif_sequence_validation is None
            else min(
                int(motif_sequence_validation.valid_sequences),
                sequence_supported_sites,
            )
        )
        motif_invalid_site_ids = set()
        if motif_sequence_validation is not None:
            motif_invalid_site_ids = {
                str(row.site_id)
                for row in motif_sequence_validation.rows
                if row.status == SEQUENCE_VALIDATION_STATUS_INVALID_SITE_ID
            }

        malformed_identifier_sites = self._find_malformed_site_identifiers(
            request.dataset.phospho.index
        )
        invalid_or_missing_identifier_sites = set(malformed_identifier_sites)
        invalid_or_missing_identifier_sites.update(motif_invalid_site_ids)

        reference_evidence_sites = {
            str(site_id)
            for site_ids in scoring_execution.quantified_substrates.values()
            for site_id in site_ids
        }
        final_output_sites = self._count_sites_with_non_missing_scores(
            prediction_result.pred_mat
        )
        activity_sites: int | None = None
        if activity_enabled:
            activity_input_scores = prediction_result.pred_mat.loc[
                prediction_result.pred_mat.index.intersection(
                    request.activity_phospho_matrix.index
                )
            ]
            activity_sites = self._count_sites_with_non_missing_scores(
                activity_input_scores
            )

        return KinaseWorkflowScoringAttritionSummary(
            rows_removed_invalid_or_missing_site_identifiers=int(
                len(invalid_or_missing_identifier_sites)
            ),
            final_quantitative_sites_entering_scoring=final_quantitative_sites,
            sites_with_valid_site_sequence=motif_valid_sites,
            sites_without_usable_site_sequence=max(
                final_quantitative_sites - motif_valid_sites,
                0,
            ),
            sites_eligible_for_motif_scoring=motif_valid_sites,
            sites_with_kinase_substrate_reference_profile_evidence=int(
                len(reference_evidence_sites)
            ),
            sites_contributing_to_final_fused_prediction_scoring_output=(
                final_output_sites
            ),
            sites_contributing_to_activity_scoring=activity_sites,
        )

    @staticmethod
    def _count_sites_with_non_missing_scores(score_frame: pd.DataFrame) -> int:
        if score_frame.empty:
            return 0
        return int(score_frame.notna().any(axis=1).sum())

    def _find_malformed_site_identifiers(self, index: pd.Index) -> set[str]:
        malformed: set[str] = set()
        for site_id in index.tolist():
            raw_site_id = str(site_id)
            try:
                parse_canonical_site_identifier(
                    raw_site_id,
                    field_name="kinase.workflow.dataset.phospho.index",
                    error_type=_WorkflowIdentifierError,
                )
            except _WorkflowIdentifierError:
                malformed.add(raw_site_id)
            except Exception as exc:  # pragma: no cover - defensive boundary guard
                raise WorkflowValidationError(
                    "failed to evaluate site identifier validity for attrition "
                    f"summary: {exc}"
                ) from exc
        return malformed


__all__ = ["KinaseSiteAttritionSummaryComposer"]
