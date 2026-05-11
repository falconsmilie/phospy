"""Localisation-confidence eligibility stage for dataset preprocessing."""

from __future__ import annotations

from dataclasses import replace
from numbers import Real

import pandas as pd

from phospy.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_LOCALISATION,
    PreprocessingPlan,
    PreprocessingStageResult,
    PreprocessingState,
    PreprocessingStateTableKey,
    append_row_audit_records,
)
from phospy.datasets.preprocessing.report_rows import report_rows_from_row_audit_rows
from phospy.datasets.preprocessing.report_schema import PreprocessingRowAuditRow
from phospy.datasets.preprocessing.stage_contract import PreprocessingStageContract
from phospy.errors.input import PhosPyInputError
from phospy.policy_models import LocalisationEligibilityMode
from phospy.validation.datasets.site_metadata import (
    assess_localisation_probability_column,
)

_EXAMPLE_LIMIT = 5


def _coerce_confidence_value(value: object) -> float:
    if isinstance(value, bool):
        raise TypeError("localisation confidence cannot be boolean")
    if isinstance(value, Real):
        return float(value)
    if isinstance(value, str):
        return float(value)
    raise TypeError("localisation confidence must be a real number")


class LocalisationConfidenceStage:
    """Apply site-level localisation confidence policy at preprocessing boundary."""

    stage_key = DATASET_PREPROCESSING_STAGE_LOCALISATION

    def run(self, state: PreprocessingState) -> PreprocessingStageResult:
        mode = state.plan.localisation_mode
        column_name = state.plan.localisation_confidence_column
        threshold = float(state.plan.localisation_min_confidence)
        waiver_reason = state.plan.localisation_waiver_reason

        if mode is LocalisationEligibilityMode.IGNORE:
            return PreprocessingStageResult(
                state=state,
                diagnostics={
                    "dropped_row_ids": (),
                    "dropped_row_count": 0,
                    "imputed_cell_count": 0,
                    "imputed_row_ids": (),
                    "notes": "localisation policy ignored by configuration",
                    "diagnostics": {
                        "mode": mode.value,
                        "confidence_column": column_name,
                        "min_confidence": threshold,
                    },
                },
            )

        site_metadata = state.site_metadata
        if column_name not in site_metadata.columns:
            if mode is LocalisationEligibilityMode.ALLOW_MISSING_WITH_WAIVER:
                return PreprocessingStageResult(
                    state=state,
                    diagnostics={
                        "dropped_row_ids": (),
                        "dropped_row_count": 0,
                        "imputed_cell_count": 0,
                        "imputed_row_ids": (),
                        "notes": "localisation validation waived by configuration",
                        "diagnostics": {
                            "mode": mode.value,
                            "waiver_reason": waiver_reason,
                            "confidence_column": column_name,
                            "min_confidence": threshold,
                            "column_present": False,
                            "missing_count": int(site_metadata.shape[0]),
                            "below_threshold_count": 0,
                        },
                    },
                )
            raise PhosPyInputError(
                "dataset build request preprocessing localisation policy "
                f"mode={mode.value} requires site_metadata.{column_name}; "
                f"affected_rows={int(site_metadata.shape[0])}; "
                f"example_site_ids={_site_id_examples(site_metadata.index)}"
            )

        assessment = assess_localisation_probability_column(
            site_metadata=site_metadata,
            field_name="dataset.site_metadata",
            error_type=PhosPyInputError,
            column_name=column_name,
        )
        if assessment is None:  # pragma: no cover - defensive guard
            return PreprocessingStageResult(state=state)

        if assessment.invalid_count > 0:
            invalid_sites = _site_id_examples(
                site_metadata.index[assessment.invalid_mask.to_numpy()]
            )
            raise PhosPyInputError(
                "dataset build request preprocessing localisation policy "
                f"mode={mode.value} found invalid values in site_metadata.{column_name}; "
                f"affected_rows={assessment.invalid_count}; "
                f"example_site_ids={invalid_sites}; "
                f"example_values={_summarise_examples(list(assessment.invalid_examples), limit=3)}"
            )

        below_threshold = assessment.normalized.notna() & (
            assessment.normalized.astype("float64") < threshold
        )
        below_threshold_count = int(below_threshold.sum())

        if mode is LocalisationEligibilityMode.REQUIRE_THRESHOLD:
            if assessment.missing_count > 0:
                missing_sites = _site_id_examples(
                    site_metadata.index[assessment.missing_mask.to_numpy()]
                )
                raise PhosPyInputError(
                    "dataset build request preprocessing localisation policy "
                    f"mode={mode.value} requires non-missing "
                    f"site_metadata.{column_name}; "
                    f"affected_rows={assessment.missing_count}; "
                    f"example_site_ids={missing_sites}"
                )
            if below_threshold_count > 0:
                below_threshold_sites = _site_id_examples(
                    site_metadata.index[below_threshold.to_numpy()]
                )
                raise PhosPyInputError(
                    "dataset build request preprocessing localisation policy "
                    f"mode={mode.value} requires site_metadata.{column_name} >= "
                    f"{threshold:.3f}; affected_rows={below_threshold_count}; "
                    f"example_site_ids={below_threshold_sites}"
                )
            return PreprocessingStageResult(
                state=state,
                diagnostics={
                    "dropped_row_ids": (),
                    "dropped_row_count": 0,
                    "imputed_cell_count": 0,
                    "imputed_row_ids": (),
                    "notes": "stage executed",
                    "diagnostics": {
                        "mode": mode.value,
                        "waiver_reason": None,
                        "confidence_column": column_name,
                        "min_confidence": threshold,
                        "column_present": True,
                        "missing_count": 0,
                        "below_threshold_count": 0,
                    },
                },
            )

        if mode is not LocalisationEligibilityMode.ALLOW_MISSING_WITH_WAIVER:
            raise PhosPyInputError(
                "dataset preprocessing localisation stage received unsupported mode: "
                f"{mode.value!r}"
            )

        row_audit_records = _build_waiver_row_audit_records(
            site_metadata=site_metadata,
            missing_mask=assessment.missing_mask,
            below_threshold_mask=below_threshold,
            confidence_column=column_name,
            threshold=threshold,
            waiver_reason=waiver_reason or "",
            mode=mode.value,
            stage_order=state.plan.stage_order,
        )
        next_state = append_row_audit_records(state, row_audit_records)
        return PreprocessingStageResult(
            state=replace(next_state),
            report_rows=report_rows_from_row_audit_rows(row_audit_records),
            diagnostics={
                "dropped_row_ids": (),
                "dropped_row_count": 0,
                "imputed_cell_count": 0,
                "imputed_row_ids": (),
                "notes": "localisation validation waived by configuration",
                "diagnostics": {
                    "mode": mode.value,
                    "waiver_reason": waiver_reason,
                    "confidence_column": column_name,
                    "min_confidence": threshold,
                    "column_present": True,
                    "missing_count": int(assessment.missing_count),
                    "below_threshold_count": below_threshold_count,
                },
            },
        )


def _build_waiver_row_audit_records(
    *,
    site_metadata: pd.DataFrame,
    missing_mask: pd.Series,
    below_threshold_mask: pd.Series,
    confidence_column: str,
    threshold: float,
    waiver_reason: str,
    mode: str,
    stage_order: tuple[str, ...],
) -> list[PreprocessingRowAuditRow]:
    records: list[PreprocessingRowAuditRow] = []
    for site_id in site_metadata.index[missing_mask.to_numpy()]:
        site_id_text = str(site_id)
        records.append(
            PreprocessingRowAuditRow(
                stage=DATASET_PREPROCESSING_STAGE_LOCALISATION,
                action="retained",
                reason=(
                    "retained with waived missing localisation confidence "
                    f"({confidence_column})"
                ),
                source_row_id=site_id_text,
                site_id=site_id_text,
                retained=True,
                retained_row_id=site_id_text,
                source_rows=(site_id_text,),
                retained_row=site_id_text,
                parameter_snapshot={
                    "mode": mode,
                    "waiver_reason": waiver_reason,
                    "confidence_column": confidence_column,
                    "min_confidence": float(threshold),
                    "issue": "missing_confidence",
                    "stage_order": [str(stage) for stage in stage_order],
                },
            )
        )
    for site_id in site_metadata.index[below_threshold_mask.to_numpy()]:
        site_id_text = str(site_id)
        confidence_value = _coerce_confidence_value(
            site_metadata.at[site_id, confidence_column]
        )
        records.append(
            PreprocessingRowAuditRow(
                stage=DATASET_PREPROCESSING_STAGE_LOCALISATION,
                action="retained",
                reason=(
                    "retained with waived below-threshold localisation confidence "
                    f"({confidence_column})"
                ),
                source_row_id=site_id_text,
                site_id=site_id_text,
                retained=True,
                retained_row_id=site_id_text,
                source_rows=(site_id_text,),
                retained_row=site_id_text,
                parameter_snapshot={
                    "mode": mode,
                    "waiver_reason": waiver_reason,
                    "confidence_column": confidence_column,
                    "min_confidence": float(threshold),
                    "issue": "below_threshold",
                    "observed_confidence": confidence_value,
                    "stage_order": [str(stage) for stage in stage_order],
                },
            )
        )
    return records


def _site_id_examples(index: pd.Index, *, limit: int = _EXAMPLE_LIMIT) -> str:
    labels = [str(value) for value in index.tolist()]
    if not labels:
        return "(none)"
    preview = ", ".join(repr(label) for label in labels[:limit])
    suffix = "" if len(labels) <= limit else f", +{len(labels) - limit} more"
    return f"[{preview}{suffix}]"


def _summarise_examples(values: list[str], *, limit: int = _EXAMPLE_LIMIT) -> str:
    if not values:
        return "(none)"
    preview = ", ".join(values[:limit])
    suffix = "" if len(values) <= limit else f", +{len(values) - limit} more"
    return f"[{preview}{suffix}]"


def _resolve_operation(plan: PreprocessingPlan) -> str:
    return plan.localisation_mode.value


def _resolve_parameters(plan: PreprocessingPlan) -> dict[str, object]:
    return {
        "localisation_mode": plan.localisation_mode.value,
        "localisation_min_confidence": float(plan.localisation_min_confidence),
        "localisation_confidence_column": plan.localisation_confidence_column,
        "localisation_waiver_reason": plan.localisation_waiver_reason,
    }


LOCALISATION_CONFIDENCE_STAGE_CONTRACT = PreprocessingStageContract(
    stage_key=DATASET_PREPROCESSING_STAGE_LOCALISATION,
    display_label=DATASET_PREPROCESSING_STAGE_LOCALISATION,
    provenance_stage=DATASET_PREPROCESSING_STAGE_LOCALISATION,
    operation_name=_resolve_operation,
    serialize_parameters=_resolve_parameters,
    consumed_input_tables=(PreprocessingStateTableKey.DATASET_SITE_METADATA,),
    produced_output_tables=(
        PreprocessingStateTableKey.DATASET_SITE_METADATA,
        PreprocessingStateTableKey.REPORT_ROW_AUDIT,
    ),
    stage_factory=LocalisationConfidenceStage,
    backend="pandas",
    diagnostics_metadata={
        "known_diagnostics_fields": (
            "mode",
            "waiver_reason",
            "confidence_column",
            "min_confidence",
            "column_present",
            "missing_count",
            "below_threshold_count",
        )
    },
)


__all__ = [
    "LOCALISATION_CONFIDENCE_STAGE_CONTRACT",
    "LocalisationConfidenceStage",
]
