"""Alignment diagnostics for interpreted signalome inputs."""

from __future__ import annotations

import pandas as pd

from phospy.science.signalomes.models import (
    SignalomeAlignmentDiagnostics,
    SignalomeAlignmentInputDiagnostics,
)


class SignalomeAlignmentDiagnosticsBuilder:
    """Build alignment diagnostics across sites, kinases, and proteins."""

    _SITE_ID_COLUMN = "site_id"
    _KINASE_COLUMN = "kinase"
    _REASON_MISSING_FROM_DATASET = "missing_from_dataset"
    _REASON_MISSING_FROM_PREDICTION_SCORES = "missing_from_prediction_scores"
    _REASON_MISSING_FROM_DOWNSTREAM_SCORES = "missing_from_downstream_scores"
    _REASON_MISSING_FROM_KINASE_SUPPORT = "missing_kinase_support"
    _REASON_MISSING_PROTEIN_IDENTIFIER = "missing_protein_identifier"
    _REASON_REMOVED_BY_SCORE_PRECONDITIONING = "removed_by_score_preconditioning"
    _REASON_REMOVED_BY_VALIDATION_POLICY = "removed_by_validation_policy"

    def run(
        self,
        *,
        dataset_sites: pd.Index,
        prediction_sites: pd.Index,
        score_sites: pd.Index,
        shared_sites: pd.Index,
        retained_sites: pd.Index,
        prediction_kinases: pd.Index,
        score_kinases: pd.Index,
        shared_kinases: pd.Index,
        interpreted_protein_sites: pd.Index,
        retained_protein_sites: pd.Index,
    ) -> SignalomeAlignmentDiagnostics:
        dataset_site_index = pd.Index(
            dataset_sites.astype(str), name=self._SITE_ID_COLUMN
        )
        prediction_site_index = pd.Index(
            prediction_sites.astype(str),
            name=self._SITE_ID_COLUMN,
        )
        score_site_index = pd.Index(score_sites.astype(str), name=self._SITE_ID_COLUMN)
        retained_site_index = pd.Index(
            retained_sites.astype(str), name=self._SITE_ID_COLUMN
        )
        dataset_sites_diag = self._build_site_alignment_input_diagnostics(
            provided=dataset_site_index,
            retained=retained_site_index,
            missing_from_first=(
                prediction_site_index,
                self._REASON_MISSING_FROM_PREDICTION_SCORES,
            ),
            missing_from_second=(
                score_site_index,
                self._REASON_MISSING_FROM_DOWNSTREAM_SCORES,
            ),
        )
        prediction_sites_diag = self._build_site_alignment_input_diagnostics(
            provided=prediction_site_index,
            retained=retained_site_index,
            missing_from_first=(dataset_site_index, self._REASON_MISSING_FROM_DATASET),
            missing_from_second=(
                score_site_index,
                self._REASON_MISSING_FROM_DOWNSTREAM_SCORES,
            ),
        )
        downstream_sites_diag = self._build_site_alignment_input_diagnostics(
            provided=score_site_index,
            retained=retained_site_index,
            missing_from_first=(dataset_site_index, self._REASON_MISSING_FROM_DATASET),
            missing_from_second=(
                prediction_site_index,
                self._REASON_MISSING_FROM_PREDICTION_SCORES,
            ),
        )

        prediction_kinase_index = pd.Index(
            prediction_kinases.astype(str),
            name=self._KINASE_COLUMN,
        )
        score_kinase_index = pd.Index(
            score_kinases.astype(str), name=self._KINASE_COLUMN
        )
        shared_kinase_index = pd.Index(
            shared_kinases.astype(str),
            name=self._KINASE_COLUMN,
        )
        provided_kinases = prediction_kinase_index.append(
            score_kinase_index.difference(prediction_kinase_index)
        )
        kinase_reasons = {
            self._REASON_MISSING_FROM_PREDICTION_SCORES: 0,
            self._REASON_MISSING_FROM_DOWNSTREAM_SCORES: 0,
            self._REASON_MISSING_FROM_KINASE_SUPPORT: 0,
        }
        retained_kinase_set = set(shared_kinase_index.tolist())
        prediction_kinase_set = set(prediction_kinase_index.tolist())
        score_kinase_set = set(score_kinase_index.tolist())
        for kinase in provided_kinases:
            kinase_id = str(kinase)
            if kinase_id in retained_kinase_set:
                continue
            if kinase_id not in prediction_kinase_set:
                kinase_reasons[self._REASON_MISSING_FROM_PREDICTION_SCORES] += 1
            elif kinase_id not in score_kinase_set:
                kinase_reasons[self._REASON_MISSING_FROM_DOWNSTREAM_SCORES] += 1
            else:
                kinase_reasons[self._REASON_MISSING_FROM_KINASE_SUPPORT] += 1
        kinases_diag = SignalomeAlignmentInputDiagnostics(
            provided_count=int(provided_kinases.size),
            retained_count=int(shared_kinase_index.size),
            dropped_count=int(provided_kinases.size - shared_kinase_index.size),
            dropped_reasons=kinase_reasons,
        )

        interpreted_protein_index = pd.Index(
            interpreted_protein_sites.astype(str),
            name=self._SITE_ID_COLUMN,
        )
        retained_protein_index = pd.Index(
            retained_protein_sites.astype(str),
            name=self._SITE_ID_COLUMN,
        )
        protein_dropped = int(
            interpreted_protein_index.size - retained_protein_index.size
        )
        protein_diag = SignalomeAlignmentInputDiagnostics(
            provided_count=int(interpreted_protein_index.size),
            retained_count=int(retained_protein_index.size),
            dropped_count=protein_dropped,
            dropped_reasons={
                self._REASON_REMOVED_BY_SCORE_PRECONDITIONING: protein_dropped,
                self._REASON_MISSING_PROTEIN_IDENTIFIER: 0,
                self._REASON_REMOVED_BY_VALIDATION_POLICY: 0,
            },
        )

        return SignalomeAlignmentDiagnostics(
            dataset_sites=dataset_sites_diag,
            prediction_score_sites=prediction_sites_diag,
            downstream_score_sites=downstream_sites_diag,
            kinases=kinases_diag,
            protein_identifiers=protein_diag,
        )

    def _build_site_alignment_input_diagnostics(
        self,
        *,
        provided: pd.Index,
        retained: pd.Index,
        missing_from_first: tuple[pd.Index, str],
        missing_from_second: tuple[pd.Index, str],
    ) -> SignalomeAlignmentInputDiagnostics:
        provided_index = pd.Index(provided.astype(str), name=self._SITE_ID_COLUMN)
        retained_index = pd.Index(retained.astype(str), name=self._SITE_ID_COLUMN)
        first_index, first_reason = missing_from_first
        second_index, second_reason = missing_from_second
        first_site_set = set(pd.Index(first_index.astype(str)).tolist())
        second_site_set = set(pd.Index(second_index.astype(str)).tolist())
        retained_site_set = set(retained_index.tolist())
        reasons = {
            first_reason: 0,
            second_reason: 0,
            self._REASON_REMOVED_BY_SCORE_PRECONDITIONING: 0,
            self._REASON_REMOVED_BY_VALIDATION_POLICY: 0,
        }
        for site_id in provided_index:
            site = str(site_id)
            if site in retained_site_set:
                continue
            if site not in first_site_set:
                reasons[first_reason] += 1
            elif site not in second_site_set:
                reasons[second_reason] += 1
            else:
                reasons[self._REASON_REMOVED_BY_SCORE_PRECONDITIONING] += 1
        return SignalomeAlignmentInputDiagnostics(
            provided_count=int(provided_index.size),
            retained_count=int(retained_index.size),
            dropped_count=int(provided_index.size - retained_index.size),
            dropped_reasons=reasons,
        )


__all__ = ["SignalomeAlignmentDiagnosticsBuilder"]
