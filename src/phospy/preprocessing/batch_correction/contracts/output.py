"""Passive corrected-output contract for the batch-correction capability."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from typing import Protocol, runtime_checkable

import numpy as np
import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.frames.ownership import (
    export_dataframe,
    own_dataframe,
    own_optional_dataframe,
)
from phospy.preparation_records.batch_correction import (
    BATCH_CORRECTION_CONFOUNDING_PASSED,
    BATCH_CORRECTION_DESIGN_PRESERVATION_PRESERVE_CONDITION_EFFECTS,
    BATCH_CORRECTION_STATUS_APPLIED,
    BatchCorrectionDiagnostics,
    BatchCorrectionPolicy,
    BatchCorrectionReport,
)
from phospy.provenance.models.batch_correction import BatchCorrectionProvenance

_BATCH_CORRECTION_STAGE = "batch_correction"
_DOWNSTREAM_WORKFLOWS_STAGE = "downstream_workflows"
_SUPPORTED_STATUS_VALUES = frozenset({"corrected_observed", "restored_missing"})


@runtime_checkable
class _DiagnosticsLike(Protocol):
    @property
    def method(self) -> str: ...

    @property
    def matrix_shape_before(self) -> tuple[int, int]: ...

    @property
    def matrix_shape_after(self) -> tuple[int, int]: ...

    @property
    def control_site_count(self) -> int: ...

    @property
    def warnings(self) -> Sequence[str]: ...

    def to_payload(self) -> dict[str, object]: ...


@runtime_checkable
class _CorrectionResultLike(Protocol):
    @property
    def corrected_matrix(self) -> pd.DataFrame: ...

    @property
    def output_observation_mask(self) -> pd.DataFrame: ...

    @property
    def corrected_cell_status(self) -> pd.DataFrame: ...

    @property
    def diagnostics(self) -> _DiagnosticsLike: ...

    @property
    def warnings(self) -> Sequence[str]: ...

    @property
    def provenance_payload(self) -> Mapping[str, object]: ...


@dataclass(frozen=True, slots=True)
class CorrectedPreprocessingOutput:
    """Resolved correction output consumed by dataset preprocessing."""

    corrected_matrix: pd.DataFrame
    batch_correction_report: BatchCorrectionReport
    diagnostics: Mapping[str, object]
    output_observation_mask: pd.DataFrame | None = None
    corrected_cell_status: pd.DataFrame | None = None
    provenance: BatchCorrectionProvenance | Mapping[str, object] | None = None
    stage_order: tuple[str, ...] = (
        "missing_data",
        _BATCH_CORRECTION_STAGE,
        _DOWNSTREAM_WORKFLOWS_STAGE,
    )
    consumed_by_downstream: bool = False

    def __post_init__(self) -> None:
        corrected = own_dataframe(
            self.corrected_matrix,
            field_name="corrected preprocessing output.corrected_matrix",
            error_type=PhosPyInputError,
        )
        mask = own_optional_dataframe(
            self.output_observation_mask,
            field_name="corrected preprocessing output.output_observation_mask",
            error_type=PhosPyInputError,
        )
        status = own_optional_dataframe(
            self.corrected_cell_status,
            field_name="corrected preprocessing output.corrected_cell_status",
            error_type=PhosPyInputError,
        )
        _require_numeric_complete_matrix(corrected)
        if mask is not None:
            _require_aligned_frame(
                mask,
                corrected,
                field_name="corrected preprocessing output.output_observation_mask",
            )
            _require_boolean_mask(mask)
        if status is not None:
            _require_aligned_frame(
                status,
                corrected,
                field_name="corrected preprocessing output.corrected_cell_status",
            )
            _require_status_values(status)
        _require_stage_order_precedes_downstream(self.stage_order)
        if not isinstance(self.batch_correction_report, BatchCorrectionReport):
            raise PhosPyInputError(
                "corrected preprocessing output.batch_correction_report must be "
                "BatchCorrectionReport"
            )
        object.__setattr__(self, "corrected_matrix", corrected)
        object.__setattr__(self, "output_observation_mask", mask)
        object.__setattr__(self, "corrected_cell_status", status)
        object.__setattr__(
            self,
            "stage_order",
            tuple(str(stage).strip() for stage in self.stage_order),
        )

    @property
    def corrected(self) -> pd.DataFrame:
        """Return a mutation-isolated corrected matrix snapshot."""

        return export_dataframe(self.corrected_matrix)

    def with_provenance(
        self,
        provenance: BatchCorrectionProvenance | Mapping[str, object],
    ) -> CorrectedPreprocessingOutput:
        """Return a copy with workflow provenance attached."""

        return replace(self, provenance=provenance)

    @classmethod
    def from_sps_ruv_style_result(
        cls,
        result: _CorrectionResultLike,
        *,
        provenance: BatchCorrectionProvenance | Mapping[str, object] | None = None,
        stage_order: Sequence[str] | None = None,
    ) -> CorrectedPreprocessingOutput:
        """Adapt a resolved SPS/RUV-style executor result for preprocessing."""

        diagnostics = result.diagnostics.to_payload()
        seed_data = _sps_ruv_provenance_seed_data(result.provenance_payload)
        design_summary = _mapping_or_empty(diagnostics.get("design_summary"))
        algorithm_description = str(
            diagnostics.get("algorithm_description")
            or "native SPS/RUV-style correction preserves matrix shape and "
            "condition-design terms"
        )
        term_roles = _mapping_or_empty(diagnostics.get("term_roles"))
        batch_term_role = term_roles.get("batch_terms")
        limitations = (
            algorithm_description,
            *(
                (str(batch_term_role),)
                if batch_term_role is not None and str(batch_term_role).strip()
                else ()
            ),
            "native SPS/RUV-style correction preserves matrix shape and "
            "condition-design terms",
        )
        condition_columns = _object_sequence(seed_data.get("condition_columns"))
        normalized_condition_columns = tuple(
            str(column) for column in condition_columns
        )
        batch_report = BatchCorrectionReport(
            status=BATCH_CORRECTION_STATUS_APPLIED,
            policy=BatchCorrectionPolicy(
                method=str(result.diagnostics.method),
                batch_column=_optional_string(seed_data.get("batch_column")),
                condition_column=(
                    None
                    if not normalized_condition_columns
                    else normalized_condition_columns[0]
                ),
                condition_columns=normalized_condition_columns,
                design_preservation_policy=(
                    BATCH_CORRECTION_DESIGN_PRESERVATION_PRESERVE_CONDITION_EFFECTS
                ),
                preserve_condition_effects=True,
            ),
            diagnostics=BatchCorrectionDiagnostics(
                number_of_batches=_optional_int(
                    design_summary.get("number_of_batches")
                ),
                batch_levels=tuple(
                    str(level)
                    for level in _object_sequence(design_summary.get("batch_levels"))
                ),
                condition_levels=tuple(
                    str(level)
                    for level in _object_sequence(
                        design_summary.get("condition_levels")
                    )
                )
                or tuple(
                    str(term)
                    for term in _object_sequence(
                        diagnostics.get("protected_design_terms")
                    )
                    if str(term) != "intercept"
                ),
                confounding_check_status=BATCH_CORRECTION_CONFOUNDING_PASSED,
                matrix_shape_before=_matrix_shape(
                    result.diagnostics.matrix_shape_before
                ),
                matrix_shape_after=_matrix_shape(result.diagnostics.matrix_shape_after),
                warnings=tuple(str(warning) for warning in result.warnings),
                limitations=limitations,
            ),
        )
        return cls(
            corrected_matrix=result.corrected_matrix,
            output_observation_mask=result.output_observation_mask,
            corrected_cell_status=result.corrected_cell_status,
            batch_correction_report=batch_report,
            diagnostics={
                "executor": diagnostics,
                "provenance_payload": dict(result.provenance_payload),
            },
            provenance=provenance,
            stage_order=(
                tuple(str(stage) for stage in stage_order)
                if stage_order is not None
                else (
                    "missing_data",
                    _BATCH_CORRECTION_STAGE,
                    _DOWNSTREAM_WORKFLOWS_STAGE,
                )
            ),
        )


def _require_numeric_complete_matrix(matrix: pd.DataFrame) -> None:
    if matrix.empty:
        raise PhosPyInputError(
            "corrected preprocessing output.corrected_matrix must not be empty"
        )
    non_numeric = [
        str(column)
        for column in matrix.columns.tolist()
        if not pd.api.types.is_numeric_dtype(matrix.loc[:, column])
        or pd.api.types.is_bool_dtype(matrix.loc[:, column])
    ]
    if non_numeric:
        raise PhosPyInputError(
            "corrected preprocessing output.corrected_matrix must contain numeric "
            "non-bool columns. Non-numeric columns: " + ", ".join(non_numeric)
        )
    values = matrix.astype("float64").to_numpy(copy=True)
    if not np.isfinite(values).all():
        raise PhosPyInputError(
            "corrected preprocessing output.corrected_matrix must be finite and "
            "missing-value-free for AnalysisReadyPhosphoDataset construction"
        )


def _object_sequence(value: object) -> tuple[object, ...]:
    if isinstance(value, Sequence) and not isinstance(value, str | bytes | bytearray):
        return tuple(value)
    return ()


def _mapping_or_empty(value: object) -> Mapping[str, object]:
    return value if isinstance(value, Mapping) else {}


def _sps_ruv_provenance_seed_data(
    provenance_payload: Mapping[str, object],
) -> Mapping[str, object]:
    seed = provenance_payload.get("provenance_seed_data")
    if isinstance(seed, Mapping):
        return seed
    plan = provenance_payload.get("resolved_plan")
    if isinstance(plan, Mapping):
        plan_seed = plan.get("provenance_seed_data")
        if isinstance(plan_seed, Mapping):
            return plan_seed
    return {}


def _optional_string(value: object) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return None if text == "" else text


def _optional_int(value: object) -> int | None:
    if isinstance(value, bool) or not isinstance(value, int):
        return None
    return int(value)


def _matrix_shape(value: tuple[int, int]) -> tuple[int, int]:
    return (int(value[0]), int(value[1]))


def _require_aligned_frame(
    frame: pd.DataFrame,
    expected: pd.DataFrame,
    *,
    field_name: str,
) -> None:
    if not frame.index.equals(expected.index):
        raise PhosPyInputError(f"{field_name}.index must match dataset.phospho.index")
    if not frame.columns.equals(expected.columns):
        raise PhosPyInputError(
            f"{field_name}.columns must match dataset.phospho.columns"
        )


def _require_boolean_mask(mask: pd.DataFrame) -> None:
    for row_position, row_id in enumerate(mask.index.tolist()):
        for column_position, column_id in enumerate(mask.columns.tolist()):
            value = mask.iat[row_position, column_position]
            if isinstance(value, (bool, np.bool_)):
                continue
            raise PhosPyInputError(
                "corrected preprocessing output.output_observation_mask must "
                "contain only boolean values; "
                f"invalid_cell=({row_id!r}, {column_id!r})"
            )


def _require_status_values(status: pd.DataFrame) -> None:
    values = status.to_numpy(dtype="object", copy=True)
    for row_position, row_id in enumerate(status.index.tolist()):
        for column_position, column_id in enumerate(status.columns.tolist()):
            value = str(values[row_position, column_position]).strip()
            if value in _SUPPORTED_STATUS_VALUES:
                continue
            raise PhosPyInputError(
                "corrected preprocessing output.corrected_cell_status contains "
                f"unsupported status {value!r} at ({row_id!r}, {column_id!r})"
            )


def _require_stage_order_precedes_downstream(stage_order: Sequence[str]) -> None:
    normalized = tuple(str(stage).strip() for stage in stage_order)
    if _BATCH_CORRECTION_STAGE not in normalized:
        raise PhosPyInputError(
            "corrected preprocessing output.stage_order must include 'batch_correction'"
        )
    if _DOWNSTREAM_WORKFLOWS_STAGE not in normalized:
        return
    correction_index = normalized.index(_BATCH_CORRECTION_STAGE)
    downstream_index = normalized.index(_DOWNSTREAM_WORKFLOWS_STAGE)
    if correction_index > downstream_index:
        raise PhosPyInputError(
            "corrected preprocessing output.stage_order must place "
            "batch_correction before downstream_workflows"
        )


__all__ = ["CorrectedPreprocessingOutput"]
