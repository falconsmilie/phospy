"""Workflow-row attrition provenance for signalome execution."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import pandas as pd

from phospy.errors.workflows import WorkflowStageError
from phospy.provenance.models import RowAttritionRecord, RowAttritionReport
from phospy.workflows._row_attrition import (
    make_row_attrition_record,
    reconcile_row_attrition_report,
)
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
    *,
    final_site_ids: pd.Index | None = None,
) -> SignalomeRowAttritionProvenance:
    """Build standardized row-attrition provenance for signalome execution."""

    site_metadata = request.dataset.site_metadata
    input_site_ids = _index_values(request.dataset.phospho.index)
    final_site_index = (
        request.downstream_score_matrix.index
        if final_site_ids is None
        else final_site_ids
    )
    retained_site_ids = _index_values(final_site_index)
    missing_sequence_ids = _missing_text_ids(site_metadata, "site_sequence")
    missing_localisation_ids = _missing_localisation_probability_ids(site_metadata)
    below_localisation_threshold_ids = _below_localisation_threshold_ids(request)
    missing_protein_ids = _missing_text_ids(site_metadata, "protein_id")
    preconditioning = request.score_preconditioning_diagnostics
    records = _causal_site_row_records(
        request=request,
        final_site_ids=final_site_ids,
    )
    _require_metric_record_agreement(
        record=_record_for_stage(records, "signalome_score_preconditioning"),
        metric_count=int(preconditioning.dropped_all_missing_row_count),
    )
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
        row_attrition=reconcile_row_attrition_report(
            workflow="signalome",
            records=records,
            initial_site_ids=request.dataset.phospho.index,
            final_site_ids=final_site_index,
        ),
    )


def _causal_site_row_records(
    *,
    request: ResolvedSignalomeWorkflowRequest,
    final_site_ids: pd.Index | None,
) -> tuple[RowAttritionRecord, ...]:
    records = list(request.row_attrition_records)
    if final_site_ids is not None:
        retention_record = make_row_attrition_record(
            workflow="signalome",
            stage="signalome_scoring_clustering_retention",
            reason="not_retained_for_signalome_scoring_clustering",
            input_site_ids=request.downstream_score_matrix.index,
            output_site_ids=final_site_ids,
        )
        if retention_record is not None:
            records.append(retention_record)
    return tuple(records)


def _record_for_stage(
    records: tuple[RowAttritionRecord, ...],
    stage: str,
) -> RowAttritionRecord | None:
    for record in records:
        if record.stage == stage:
            return record
    return None


def _require_metric_record_agreement(
    *,
    record: RowAttritionRecord | None,
    metric_count: int,
) -> None:
    if record is None or int(record.removed_rows) == int(metric_count):
        return
    raise WorkflowStageError(
        "signalome row attrition internal consistency error; "
        f"stage={record.stage}; input_count={int(record.input_rows)}; "
        f"output_count={int(record.output_rows)}; "
        "score preconditioning metric count disagrees with causal record; "
        f"metric_count={int(metric_count)}"
    )


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
