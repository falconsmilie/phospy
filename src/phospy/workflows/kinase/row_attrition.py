"""Workflow-row attrition provenance for kinase execution."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import pandas as pd

from phospy.provenance.models import RowAttritionRecord, RowAttritionReport
from phospy.science.prediction.models import KinaseScoringResult
from phospy.science.tables.kinase import (
    KINASE_PROFILE_SCORE_DIAGNOSTIC_STATUS_UNSCORED,
)
from phospy.workflows._row_attrition import (
    make_row_attrition_record,
    reconcile_row_attrition_report,
)
from phospy.workflows.kinase.contracts import ResolvedKinaseWorkflowRequest
from phospy.workflows.kinase.scoring_mode_contracts import (
    kinase_scoring_mode_input_contract,
)

_LOCALISATION_COLUMNS = ("localisation_confidence", "localisation_probability")
_EXAMPLE_LIMIT = 5


@dataclass(frozen=True, slots=True)
class KinaseRowAttritionProvenance:
    """Kinase site-row and site/kinase-pair attrition provenance."""

    metrics: Mapping[str, object]
    row_attrition: RowAttritionReport | None

    def to_workflow_parameters(self) -> dict[str, object]:
        payload: dict[str, object] = {
            "row_attrition_metrics": dict(self.metrics),
        }
        if self.row_attrition is not None:
            payload["row_attrition"] = self.row_attrition.to_payload()
        return payload


def build_kinase_row_attrition_provenance(
    *,
    request: ResolvedKinaseWorkflowRequest,
    scoring_result: KinaseScoringResult,
) -> KinaseRowAttritionProvenance:
    """Build standardized row-attrition provenance for kinase scoring."""

    input_site_ids = _index_values(request.dataset.phospho.index)
    scoring_site_ids = _index_values(request.scoring_site_index)
    reference_supported_site_ids = _reference_supported_site_ids(request)
    scored_site_ids = _site_ids_with_any_score(scoring_result.authoritative_scores)
    not_present_in_reference_site_ids = tuple(
        site_id
        for site_id in scoring_site_ids
        if site_id not in reference_supported_site_ids
    )
    records = _causal_site_row_records(
        request=request,
        scoring_result=scoring_result,
    )
    sequence_record = _record_for_stage(records, "kinase_sequence_context")
    pair_metrics = _site_kinase_pair_metrics(scoring_result)
    metrics: dict[str, object] = {
        "input_sites": int(len(input_site_ids)),
        "sites_missing_valid_centered_sequence": int(
            0 if sequence_record is None else sequence_record.removed_rows
        ),
        "sites_below_localisation_threshold": int(
            _sites_below_localisation_threshold(request)
        ),
        "sites_not_present_in_reference_resource": int(
            len(not_present_in_reference_site_ids)
        ),
        "sites_with_reference_and_sequence_support": int(
            0
            if request.attrition_metrics is None
            else request.attrition_metrics.scored_sites
        ),
        "sites_scored": int(len(scored_site_ids)),
        "site_examples_missing_valid_centered_sequence": list(
            () if sequence_record is None else sequence_record.examples
        ),
        "site_examples_not_present_in_reference_resource": list(
            _examples(not_present_in_reference_site_ids)
        ),
        **pair_metrics,
    }
    return KinaseRowAttritionProvenance(
        metrics=metrics,
        row_attrition=reconcile_row_attrition_report(
            workflow="kinase",
            records=records,
            initial_site_ids=request.dataset.phospho.index,
            final_site_ids=scoring_result.authoritative_scores.index,
        ),
    )


def _causal_site_row_records(
    *,
    request: ResolvedKinaseWorkflowRequest,
    scoring_result: KinaseScoringResult,
) -> tuple[RowAttritionRecord, ...]:
    records = list(request.row_attrition_records)
    scoring_record = make_row_attrition_record(
        workflow="kinase",
        stage="kinase_scoring_retention",
        reason="sites_removed_by_scoring_retention",
        input_site_ids=request.scoring_site_index,
        output_site_ids=scoring_result.authoritative_scores.index,
    )
    if scoring_record is not None:
        records.append(scoring_record)
    return tuple(records)


def _record_for_stage(
    records: tuple[RowAttritionRecord, ...],
    stage: str,
) -> RowAttritionRecord | None:
    for record in records:
        if record.stage == stage:
            return record
    return None


def _reference_supported_site_ids(
    request: ResolvedKinaseWorkflowRequest,
) -> frozenset[str]:
    mode_contract = kinase_scoring_mode_input_contract(
        request.execution_config.scoring_mode
    )
    if not mode_contract.requires_substrate_reference_overlap:
        return frozenset(_index_values(request.scoring_site_index))
    if "substrate_site" not in request.kinase_substrate_map.columns:
        return frozenset()
    return frozenset(
        str(value)
        for value in request.kinase_substrate_map.loc[:, "substrate_site"].tolist()
    )


def _sites_below_localisation_threshold(
    request: ResolvedKinaseWorkflowRequest,
) -> int:
    threshold = request.execution_config.localisation_requirement.minimum_probability
    if threshold is None:
        return 0
    site_metadata = request.dataset.site_metadata
    column_name = next(
        (column for column in _LOCALISATION_COLUMNS if column in site_metadata.columns),
        None,
    )
    if column_name is None:
        return 0
    values = pd.to_numeric(site_metadata.loc[:, column_name], errors="coerce")
    below_threshold = values.notna() & (values.astype("float64") < float(threshold))
    return int(below_threshold.sum())


def _site_kinase_pair_metrics(
    scoring_result: KinaseScoringResult,
) -> dict[str, object]:
    authoritative_scores = scoring_result.authoritative_scores
    pair_count = int(authoritative_scores.shape[0] * authoritative_scores.shape[1])
    scored_pair_count = int(authoritative_scores.notna().sum().sum())
    unscored_pair_count = int(pair_count - scored_pair_count)
    profile_reason_counts = _profile_unscored_reason_counts(
        scoring_result.profile_score_diagnostics
    )
    return {
        "site_kinase_pairs_considered": pair_count,
        "site_kinase_pairs_scored": scored_pair_count,
        "site_kinase_pairs_unscored_due_to_insufficient_evidence": (
            unscored_pair_count
        ),
        "site_kinase_pair_unscored_reason_counts": profile_reason_counts,
        "site_kinase_pair_examples_unscored_due_to_insufficient_evidence": list(
            _unscored_pair_examples(authoritative_scores)
        ),
    }


def _profile_unscored_reason_counts(
    diagnostics: pd.DataFrame | None,
) -> dict[str, int]:
    if diagnostics is None or diagnostics.empty:
        return {}
    unscored = diagnostics.loc[
        diagnostics.loc[:, "status"].astype(str)
        == KINASE_PROFILE_SCORE_DIAGNOSTIC_STATUS_UNSCORED,
        "reason",
    ]
    return {
        str(reason): int(count)
        for reason, count in unscored.astype(str).value_counts().sort_index().items()
    }


def _site_ids_with_any_score(scores: pd.DataFrame) -> frozenset[str]:
    if scores.empty:
        return frozenset()
    scored = scores.notna().any(axis=1)
    return frozenset(str(site_id) for site_id in scores.index[scored].tolist())


def _unscored_pair_examples(scores: pd.DataFrame) -> tuple[str, ...]:
    examples: list[str] = []
    missing = scores.isna()
    for site_id, row in missing.iterrows():
        for kinase, is_missing in row.items():
            if not bool(is_missing):
                continue
            examples.append(f"{site_id}|{kinase}")
            if len(examples) >= _EXAMPLE_LIMIT:
                return tuple(examples)
    return tuple(examples)


def _index_values(index: pd.Index) -> tuple[str, ...]:
    return tuple(str(value) for value in index.tolist())


def _examples(values: tuple[str, ...]) -> tuple[str, ...]:
    return tuple(str(value) for value in values[:_EXAMPLE_LIMIT])


__all__ = [
    "KinaseRowAttritionProvenance",
    "build_kinase_row_attrition_provenance",
]
