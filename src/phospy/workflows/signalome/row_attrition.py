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
    missing_protein_ids = _missing_text_ids(site_metadata, "protein_id")
    preconditioning = request.score_preconditioning_diagnostics
    metrics: dict[str, object] = {
        "input_sites": int(len(input_site_ids)),
        "sites_missing_sequence_context": int(len(missing_sequence_ids)),
        "sites_missing_localisation_probability": int(len(missing_localisation_ids)),
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
        ),
    )


def _site_row_attrition_report(
    *,
    request: ResolvedSignalomeWorkflowRequest,
    input_site_ids: tuple[str, ...],
    retained_site_ids: tuple[str, ...],
) -> RowAttritionReport | None:
    diagnostics = request.alignment_diagnostics.dataset_sites
    dropped_count = int(diagnostics.dropped_count)
    if dropped_count <= 0:
        return None
    dropped_examples = _examples(
        tuple(
            site_id
            for site_id in input_site_ids
            if site_id not in set(retained_site_ids)
        )
    )
    records: list[RowAttritionRecord] = []
    current_rows = int(diagnostics.provided_count)
    accounted = 0
    for reason, count in sorted(diagnostics.dropped_reasons.items()):
        removed_rows = int(count)
        if removed_rows <= 0:
            continue
        current_rows = _append_record(
            records,
            stage="signalome_site_alignment",
            input_rows=current_rows,
            removed_rows=removed_rows,
            reason=str(reason),
            examples=dropped_examples,
        )
        accounted += removed_rows
    remaining = dropped_count - accounted
    if remaining > 0:
        current_rows = _append_record(
            records,
            stage="signalome_site_alignment",
            input_rows=current_rows,
            removed_rows=remaining,
            reason="not_retained_for_signalome_scoring_clustering",
            examples=dropped_examples,
        )
    if not records:
        return None
    _ = current_rows
    return RowAttritionReport.from_records(records)


def _append_record(
    records: list[RowAttritionRecord],
    *,
    stage: str,
    input_rows: int,
    removed_rows: int,
    reason: str,
    examples: tuple[str, ...],
) -> int:
    output_rows = int(input_rows) - int(removed_rows)
    records.append(
        RowAttritionRecord(
            stage=stage,
            input_rows=int(input_rows),
            output_rows=output_rows,
            removed_rows=int(removed_rows),
            reason=reason,
            examples=examples,
        )
    )
    return output_rows


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
