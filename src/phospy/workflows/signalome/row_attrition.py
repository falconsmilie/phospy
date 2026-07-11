"""Workflow-row attrition provenance for signalome execution."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import pandas as pd

from phospy.provenance.models import RowAttritionRecord, RowAttritionReport
from phospy.workflows.signalome.contracts import ResolvedSignalomeWorkflowRequest

_EXAMPLE_LIMIT = 5
_LOCALISATION_COLUMNS = ("localisation_confidence", "localisation_probability")


@dataclass(frozen=True, slots=True)
class SignalomeRowAttritionProvenance:
    """Signalome site-row attrition provenance."""

    metrics: Mapping[str, object]
    row_attrition: RowAttritionReport | None

    def to_workflow_parameters(self) -> dict[str, object]:
        payload: dict[str, object] = {
            "row_attrition_metrics": dict(self.metrics),
        }
        if self.row_attrition is not None:
            payload["row_attrition"] = self.row_attrition.to_payload()
        return payload


def build_signalome_row_attrition_provenance(
    request: ResolvedSignalomeWorkflowRequest,
) -> SignalomeRowAttritionProvenance:
    """Build standardized row-attrition provenance for signalome execution."""

    site_metadata = request.dataset.site_metadata
    input_site_ids = _index_values(request.dataset.phospho.index)
    retained_site_ids = _index_values(request.downstream_score_matrix.index)
    missing_sequence_ids = _missing_text_ids(site_metadata, "site_sequence")
    missing_localisation_ids = _missing_localisation_probability_ids(site_metadata)
    below_localisation_threshold_ids = _below_localisation_threshold_ids(request)
    missing_protein_ids = _missing_text_ids(site_metadata, "protein_id")
    preconditioning = request.score_preconditioning_diagnostics
    metrics: dict[str, object] = {
        "input_sites": int(len(input_site_ids)),
        "sites_missing_sequence_context": int(len(missing_sequence_ids)),
        "sites_missing_localisation_probability": int(len(missing_localisation_ids)),
        "sites_below_localisation_threshold": int(
            len(below_localisation_threshold_ids)
        ),
        "sites_missing_protein_grouping_metadata": int(len(missing_protein_ids)),
        "sites_removed_by_score_preconditioning": int(
            preconditioning.dropped_all_missing_row_count
        ),
        "sites_retained_for_signalome_scoring_clustering": int(len(retained_site_ids)),
        "site_examples_missing_sequence_context": list(_examples(missing_sequence_ids)),
        "site_examples_missing_localisation_probability": list(
            _examples(missing_localisation_ids)
        ),
        "site_examples_missing_protein_grouping_metadata": list(
            _examples(missing_protein_ids)
        ),
    }
    return SignalomeRowAttritionProvenance(
        metrics=metrics,
        row_attrition=_site_row_attrition_report(
            request=request,
            input_site_ids=input_site_ids,
            retained_site_ids=retained_site_ids,
            missing_sequence_ids=missing_sequence_ids,
            missing_localisation_ids=missing_localisation_ids,
            below_localisation_threshold_ids=below_localisation_threshold_ids,
            missing_protein_ids=missing_protein_ids,
        ),
    )


def _site_row_attrition_report(
    *,
    request: ResolvedSignalomeWorkflowRequest,
    input_site_ids: tuple[str, ...],
    retained_site_ids: tuple[str, ...],
    missing_sequence_ids: tuple[str, ...],
    missing_localisation_ids: tuple[str, ...],
    below_localisation_threshold_ids: tuple[str, ...],
    missing_protein_ids: tuple[str, ...],
) -> RowAttritionReport | None:
    final_site_set = set(retained_site_ids)
    alignment_drops = _site_alignment_drop_ids_by_reason(
        request=request,
        current_site_ids=input_site_ids,
    )
    alignment_drop_site_ids = {
        site_id for _, site_ids in alignment_drops for site_id in site_ids
    }
    preconditioning_drop_ids = _score_preconditioning_drop_ids(
        request=request,
        current_site_ids=tuple(
            site_id
            for site_id in input_site_ids
            if site_id not in alignment_drop_site_ids
        ),
        retained_site_ids=retained_site_ids,
    )
    exact_drop_site_ids = alignment_drop_site_ids.union(preconditioning_drop_ids)
    records: list[RowAttritionRecord] = []
    current_site_ids = input_site_ids
    current_site_ids = _append_site_row_record(
        records,
        stage="signalome_sequence_context",
        reason="sites_missing_sequence_context",
        current_site_ids=current_site_ids,
        removed_site_ids=_metadata_removed_before_final(
            current_site_ids=current_site_ids,
            candidate_site_ids=missing_sequence_ids,
            final_site_ids=final_site_set,
            excluded_site_ids=exact_drop_site_ids,
        ),
    )
    current_site_ids = _append_site_row_record(
        records,
        stage="signalome_localisation_metadata",
        reason="sites_missing_localisation_probability",
        current_site_ids=current_site_ids,
        removed_site_ids=_metadata_removed_before_final(
            current_site_ids=current_site_ids,
            candidate_site_ids=_localisation_missing_record_candidates(
                request=request,
                missing_localisation_ids=missing_localisation_ids,
            ),
            final_site_ids=final_site_set,
            excluded_site_ids=exact_drop_site_ids,
        ),
    )
    current_site_ids = _append_site_row_record(
        records,
        stage="signalome_localisation_metadata",
        reason="sites_below_localisation_threshold",
        current_site_ids=current_site_ids,
        removed_site_ids=_metadata_removed_before_final(
            current_site_ids=current_site_ids,
            candidate_site_ids=below_localisation_threshold_ids,
            final_site_ids=final_site_set,
            excluded_site_ids=exact_drop_site_ids,
        ),
    )
    current_site_ids = _append_site_row_record(
        records,
        stage="signalome_protein_grouping",
        reason="sites_missing_protein_grouping_metadata",
        current_site_ids=current_site_ids,
        removed_site_ids=_metadata_removed_before_final(
            current_site_ids=current_site_ids,
            candidate_site_ids=missing_protein_ids,
            final_site_ids=final_site_set,
            excluded_site_ids=exact_drop_site_ids,
        ),
    )
    for reason, removed_site_ids in alignment_drops:
        current_site_ids = _append_site_row_record(
            records,
            stage="signalome_site_alignment",
            reason=reason,
            current_site_ids=current_site_ids,
            removed_site_ids=removed_site_ids,
        )
    current_site_ids = _append_site_row_record(
        records,
        stage="signalome_score_preconditioning",
        reason="sites_removed_by_score_preconditioning",
        current_site_ids=current_site_ids,
        removed_site_ids=preconditioning_drop_ids,
    )
    current_site_ids = _append_site_row_record(
        records,
        stage="signalome_scoring_clustering_retention",
        reason="not_retained_for_signalome_scoring_clustering",
        current_site_ids=current_site_ids,
        removed_site_ids=tuple(
            site_id for site_id in current_site_ids if site_id not in final_site_set
        ),
    )
    if not records:
        return None
    return RowAttritionReport.from_records(records)


def _append_site_row_record(
    records: list[RowAttritionRecord],
    *,
    stage: str,
    reason: str,
    current_site_ids: tuple[str, ...],
    removed_site_ids: tuple[str, ...],
) -> tuple[str, ...]:
    removed_set = set(removed_site_ids)
    removed_ordered = tuple(
        site_id for site_id in current_site_ids if site_id in removed_set
    )
    if not removed_ordered:
        return current_site_ids
    output_site_ids = tuple(
        site_id for site_id in current_site_ids if site_id not in removed_set
    )
    records.append(
        RowAttritionRecord(
            stage=stage,
            input_rows=int(len(current_site_ids)),
            output_rows=int(len(output_site_ids)),
            removed_rows=int(len(removed_ordered)),
            reason=reason,
            examples=_examples(removed_ordered),
        )
    )
    return output_site_ids


def _metadata_removed_before_final(
    *,
    current_site_ids: tuple[str, ...],
    candidate_site_ids: tuple[str, ...],
    final_site_ids: set[str],
    excluded_site_ids: set[str],
) -> tuple[str, ...]:
    candidates = set(candidate_site_ids)
    return tuple(
        site_id
        for site_id in current_site_ids
        if site_id in candidates
        and site_id not in final_site_ids
        and site_id not in excluded_site_ids
    )


def _localisation_missing_record_candidates(
    *,
    request: ResolvedSignalomeWorkflowRequest,
    missing_localisation_ids: tuple[str, ...],
) -> tuple[str, ...]:
    requirement = request.execution_config.localisation_requirement
    if not requirement.requires_probability_column:
        return ()
    return missing_localisation_ids


def _site_alignment_drop_ids_by_reason(
    *,
    request: ResolvedSignalomeWorkflowRequest,
    current_site_ids: tuple[str, ...],
) -> tuple[tuple[str, tuple[str, ...]], ...]:
    prediction_site_ids = _index_values(
        request.kinase_result.prediction_result.pred_mat.index
    )
    score_site_ids = _index_values(
        request.kinase_result.scoring_result.authoritative_scores.index
    )
    prediction_sites = set(prediction_site_ids)
    score_sites = set(score_site_ids)
    shared_sites = prediction_sites.intersection(score_sites)
    missing_from_prediction: list[str] = []
    missing_from_scores: list[str] = []
    removed_by_validation: list[str] = []
    for site_id in current_site_ids:
        if site_id not in prediction_sites:
            missing_from_prediction.append(site_id)
        elif site_id not in score_sites:
            missing_from_scores.append(site_id)
        elif site_id not in shared_sites:
            removed_by_validation.append(site_id)
    return (
        ("missing_from_prediction_scores", tuple(missing_from_prediction)),
        ("missing_from_downstream_scores", tuple(missing_from_scores)),
        ("removed_by_validation_policy", tuple(removed_by_validation)),
    )


def _score_preconditioning_drop_ids(
    *,
    request: ResolvedSignalomeWorkflowRequest,
    current_site_ids: tuple[str, ...],
    retained_site_ids: tuple[str, ...],
) -> tuple[str, ...]:
    diagnostics = request.score_preconditioning_diagnostics
    expected_drop_count = int(diagnostics.dropped_all_missing_row_count)
    if expected_drop_count <= 0:
        return ()
    retained = set(retained_site_ids)
    not_retained = tuple(
        site_id for site_id in current_site_ids if site_id not in retained
    )
    if not not_retained:
        return ()
    all_missing = _all_missing_score_row_ids(
        request=request,
        candidate_site_ids=not_retained,
    )
    if len(all_missing) == expected_drop_count:
        return all_missing
    if len(not_retained) == expected_drop_count:
        return not_retained
    return all_missing


def _all_missing_score_row_ids(
    *,
    request: ResolvedSignalomeWorkflowRequest,
    candidate_site_ids: tuple[str, ...],
) -> tuple[str, ...]:
    scores = request.kinase_result.scoring_result.authoritative_scores
    columns = request.downstream_score_matrix.columns
    if set(columns.astype(str)).issubset(set(scores.columns.astype(str))):
        scores = scores.copy(deep=False)
        scores.columns = pd.Index(scores.columns.astype(str), name=scores.columns.name)
        matrix = scores.reindex(index=candidate_site_ids, columns=columns.astype(str))
    else:
        matrix = scores.reindex(index=candidate_site_ids)
    all_missing = matrix.isna().all(axis=1)
    return tuple(str(site_id) for site_id in matrix.index[all_missing].tolist())


def _missing_localisation_probability_ids(
    site_metadata: pd.DataFrame,
) -> tuple[str, ...]:
    column_name = next(
        (column for column in _LOCALISATION_COLUMNS if column in site_metadata.columns),
        None,
    )
    if column_name is None:
        return _index_values(site_metadata.index)
    values = site_metadata.loc[:, column_name]
    missing = values.isna() | (values.astype(str).str.strip() == "")
    return tuple(str(site_id) for site_id in site_metadata.index[missing].tolist())


def _below_localisation_threshold_ids(
    request: ResolvedSignalomeWorkflowRequest,
) -> tuple[str, ...]:
    threshold = request.execution_config.localisation_requirement.minimum_probability
    if threshold is None:
        return ()
    site_metadata = request.dataset.site_metadata
    column_name = next(
        (column for column in _LOCALISATION_COLUMNS if column in site_metadata.columns),
        None,
    )
    if column_name is None:
        return ()
    values = pd.to_numeric(site_metadata.loc[:, column_name], errors="coerce")
    below_threshold = values.notna() & (values.astype("float64") < float(threshold))
    return tuple(
        str(site_id) for site_id in site_metadata.index[below_threshold].tolist()
    )


def _missing_text_ids(site_metadata: pd.DataFrame, column_name: str) -> tuple[str, ...]:
    if column_name not in site_metadata.columns:
        return _index_values(site_metadata.index)
    values = site_metadata.loc[:, column_name]
    missing = values.isna() | (values.astype(str).str.strip() == "")
    return tuple(str(site_id) for site_id in site_metadata.index[missing].tolist())


def _index_values(index: pd.Index) -> tuple[str, ...]:
    return tuple(str(value) for value in index.tolist())


def _examples(values: tuple[str, ...]) -> tuple[str, ...]:
    return tuple(str(value) for value in values[:_EXAMPLE_LIMIT])


__all__ = [
    "SignalomeRowAttritionProvenance",
    "build_signalome_row_attrition_provenance",
]
