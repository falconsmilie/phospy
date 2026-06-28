"""Preprocessing boundary integration for resolved correction outputs."""

from __future__ import annotations

import hashlib
import json
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
from phospy.provenance.hashing import (
    fingerprint_matrix,
    fingerprint_optional_table,
    hash_table_tolerance,
)
from phospy.provenance.models import (
    PREPROCESSING_STAGE_DETERMINISM_PURE,
    PREPROCESSING_STAGE_PROVENANCE_SCHEMA_VERSION_V3,
    BatchCorrectionProvenance,
    JsonValue,
    TableFingerprint,
)
from phospy.provenance.serialization import batch_correction_provenance_from_payload
from phospy.science.datasets.preprocessing.batch_correction import (
    BATCH_CORRECTION_CONFOUNDING_PASSED,
    BATCH_CORRECTION_DESIGN_PRESERVATION_PRESERVE_CONDITION_EFFECTS,
    BATCH_CORRECTION_STATUS_APPLIED,
    BatchCorrectionDiagnostics,
    BatchCorrectionPolicy,
    BatchCorrectionReport,
)
from phospy.science.datasets.preprocessing.batch_correction_provenance import (
    build_native_batch_correction_provenance,
)
from phospy.science.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION,
    PreprocessingStageExecution,
    PreprocessingState,
    PreprocessingStateTableKey,
)
from phospy.validation.datasets.batch_correction import (
    normalize_applied_selected_site_key_rows,
    validate_applied_native_sps_ruv_correction_provenance,
)

_DOWNSTREAM_WORKFLOWS_STAGE = "downstream_workflows"
_SUPPORTED_STATUS_VALUES = frozenset({"corrected_observed", "restored_missing"})
_SUPPORTED_EXTERNAL_CORRECTED_OUTPUT_METHODS = frozenset({"sps_ruv_style"})
_EXTERNAL_CORRECTED_OUTPUT_PROVENANCE_REQUIREMENT = (
    "external corrected outputs require applied correction status and a supported "
    "applied correction method with method-specific provenance"
)


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
    """Resolved correction output accepted by dataset preprocessing only."""

    corrected_matrix: pd.DataFrame
    batch_correction_report: BatchCorrectionReport
    diagnostics: Mapping[str, object]
    output_observation_mask: pd.DataFrame | None = None
    corrected_cell_status: pd.DataFrame | None = None
    provenance: BatchCorrectionProvenance | Mapping[str, object] | None = None
    stage_order: tuple[str, ...] = (
        "missing_data",
        DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION,
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
                    DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION,
                    _DOWNSTREAM_WORKFLOWS_STAGE,
                )
            ),
        )


class CorrectedPreprocessingOutputIntegrator:
    """Apply a resolved correction output to preprocessing state."""

    def run(
        self,
        *,
        state: PreprocessingState,
        correction_output: CorrectedPreprocessingOutput,
    ) -> tuple[PreprocessingState, PreprocessingStageExecution]:
        if not isinstance(correction_output, CorrectedPreprocessingOutput):
            raise PhosPyInputError(
                "dataset preprocessing corrected_preprocessing_output must be "
                "CorrectedPreprocessingOutput"
            )
        _validate_corrected_output_report_and_provenance(correction_output)
        if correction_output.consumed_by_downstream:
            raise PhosPyInputError(
                "corrected preprocessing output has already been consumed by a "
                "downstream workflow and cannot be applied to dataset construction"
            )
        _require_aligned_frame(
            correction_output.corrected_matrix,
            state.phospho,
            field_name="corrected preprocessing output.corrected_matrix",
        )
        _validate_selected_controls_in_corrected_matrix(correction_output)
        _validate_applied_provenance_fingerprints(
            state=state,
            correction_output=correction_output,
        )
        _require_stage_order_precedes_downstream(correction_output.stage_order)
        corrected = correction_output.corrected_matrix.copy(deep=True)
        mask = (
            None
            if correction_output.output_observation_mask is None
            else correction_output.output_observation_mask.copy(deep=True)
        )
        previous = state
        current = replace(
            state,
            phospho=corrected,
            imputation_observation_mask=mask,
            batch_correction_report=correction_output.batch_correction_report,
        )
        execution = _build_stage_execution(
            previous=previous,
            current=current,
            correction_output=correction_output,
        )
        return current, execution


def validate_corrected_preprocessing_output(value: object) -> None:
    """Validate a dataset-build request correction output payload."""

    if value is None:
        return
    if not isinstance(value, CorrectedPreprocessingOutput):
        raise PhosPyInputError(
            "dataset build request corrected_preprocessing_output must be "
            "CorrectedPreprocessingOutput"
        )
    _validate_corrected_output_report_and_provenance(value)


def _validate_corrected_output_report_and_provenance(
    correction_output: CorrectedPreprocessingOutput,
) -> None:
    status = str(correction_output.batch_correction_report.status).strip()
    if status != BATCH_CORRECTION_STATUS_APPLIED:
        raise PhosPyInputError(
            "corrected_preprocessing_output "
            f"{_EXTERNAL_CORRECTED_OUTPUT_PROVENANCE_REQUIREMENT}; "
            f"observed status={status!r}"
        )

    method = _normalize_external_correction_method(
        correction_output.batch_correction_report.method
    )
    if method == "none":
        raise PhosPyInputError(
            "corrected_preprocessing_output "
            f"{_EXTERNAL_CORRECTED_OUTPUT_PROVENANCE_REQUIREMENT}; "
            "method='none' is reserved for uncorrected outputs and cannot "
            "describe an external corrected matrix"
        )
    if method not in _SUPPORTED_EXTERNAL_CORRECTED_OUTPUT_METHODS:
        raise PhosPyInputError(
            "corrected_preprocessing_output "
            f"{_EXTERNAL_CORRECTED_OUTPUT_PROVENANCE_REQUIREMENT}; "
            f"observed method={method!r} has no registered strict "
            "method-specific provenance validator"
        )
    validate_applied_native_sps_ruv_correction_provenance(
        method=method,
        status=status,
        provenance=correction_output.provenance,
    )


def _normalize_external_correction_method(method: object) -> str:
    normalized = str(method).strip().lower()
    if normalized == "":
        raise PhosPyInputError(
            "corrected_preprocessing_output "
            f"{_EXTERNAL_CORRECTED_OUTPUT_PROVENANCE_REQUIREMENT}; "
            "observed method is missing or empty"
        )
    return normalized


def _validate_selected_controls_in_corrected_matrix(
    correction_output: CorrectedPreprocessingOutput,
) -> None:
    provenance = correction_output.provenance
    if not isinstance(provenance, BatchCorrectionProvenance):
        return
    if str(correction_output.batch_correction_report.status).strip() != "applied":
        return
    if not provenance.selected_site_key_rows:
        return

    selected_rows = normalize_applied_selected_site_key_rows(
        provenance.selected_site_key_rows
    )
    corrected_index = {
        str(row).strip() for row in correction_output.corrected_matrix.index.tolist()
    }
    absent_rows = tuple(row for row in selected_rows if row not in corrected_index)
    if absent_rows:
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance "
            "selected_site_key_rows must all be present in the corrected matrix "
            f"index; absent rows: {_format_labels(absent_rows)}"
        )


def _validate_applied_provenance_fingerprints(
    *,
    state: PreprocessingState,
    correction_output: CorrectedPreprocessingOutput,
) -> None:
    provenance = correction_output.provenance
    if not isinstance(provenance, BatchCorrectionProvenance):
        return
    if str(correction_output.batch_correction_report.status).strip() != "applied":
        return
    _require_matching_table_fingerprint(
        actual=state.phospho,
        recorded=provenance.input_matrix_fingerprint,
        field_name="BatchCorrectionProvenance.input_matrix_fingerprint",
        actual_description="current pre-correction dataset.phospho matrix",
    )
    if provenance.output_matrix_fingerprint is not None:
        _require_matching_table_fingerprint(
            actual=correction_output.corrected_matrix,
            recorded=provenance.output_matrix_fingerprint,
            field_name="BatchCorrectionProvenance.output_matrix_fingerprint",
            actual_description="corrected preprocessing output.corrected_matrix",
        )
    for recorded_mask in provenance.observation_masks:
        _require_matching_observation_mask_fingerprint(
            correction_output=correction_output,
            recorded=recorded_mask,
        )


def _require_matching_observation_mask_fingerprint(
    *,
    correction_output: CorrectedPreprocessingOutput,
    recorded: TableFingerprint,
) -> None:
    normalized_name = recorded.name.strip().lower()
    if "corrected_cell_status" in normalized_name:
        actual = correction_output.corrected_cell_status
        if actual is None:
            raise PhosPyInputError(
                "corrected_preprocessing_output BatchCorrectionProvenance includes "
                f"observation mask fingerprint {recorded.name!r}, but "
                "corrected preprocessing output.corrected_cell_status is missing "
                "and cannot be verified"
            )
        _require_matching_table_fingerprint(
            actual=actual,
            recorded=recorded,
            field_name="BatchCorrectionProvenance.observation_masks",
            actual_description="corrected preprocessing output.corrected_cell_status",
        )
        return
    if "observation_mask" not in normalized_name:
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance includes "
            f"unsupported observation mask fingerprint {recorded.name!r}; "
            "expected an observation_mask or corrected_cell_status fingerprint"
        )
    actual_mask = correction_output.output_observation_mask
    if actual_mask is None:
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance includes "
            f"observation mask fingerprint {recorded.name!r}, but corrected "
            "preprocessing output.output_observation_mask is missing and cannot "
            "be verified"
        )
    _require_matching_table_fingerprint(
        actual=actual_mask.astype("int8"),
        recorded=recorded,
        field_name="BatchCorrectionProvenance.observation_masks",
        actual_description="corrected preprocessing output.output_observation_mask",
    )


def _require_matching_table_fingerprint(
    *,
    actual: pd.DataFrame,
    recorded: TableFingerprint,
    field_name: str,
    actual_description: str,
) -> None:
    observed = fingerprint_matrix(actual, name=recorded.name)
    if observed == recorded:
        return
    raise PhosPyInputError(
        "corrected_preprocessing_output fingerprint mismatch: "
        f"{field_name} does not match the actual {actual_description}; "
        f"fingerprint_name={recorded.name!r}; "
        f"recorded_exact_hash={recorded.exact_hash_value!r}; "
        f"actual_exact_hash={observed.exact_hash_value!r}"
    )


def _build_stage_execution(
    *,
    previous: PreprocessingState,
    current: PreprocessingState,
    correction_output: CorrectedPreprocessingOutput,
) -> PreprocessingStageExecution:
    consumed_input_tables = _collect_fingerprints(
        previous,
        (PreprocessingStateTableKey.DATASET_PHOSPHO,),
    )
    produced_output_tables = _collect_fingerprints(
        current,
        (
            PreprocessingStateTableKey.DATASET_PHOSPHO,
            PreprocessingStateTableKey.DATASET_IMPUTATION_OBSERVATION_MASK,
        ),
    )
    diagnostics = _json_mapping(
        {
            **dict(correction_output.diagnostics),
            "status": correction_output.batch_correction_report.status,
            "method": correction_output.batch_correction_report.method,
            "condition_column": (
                correction_output.batch_correction_report.condition_column
            ),
            "condition_columns": list(
                correction_output.batch_correction_report.condition_columns
            ),
            "stage_order": list(correction_output.stage_order),
            "has_output_observation_mask": (
                correction_output.output_observation_mask is not None
            ),
            "has_corrected_cell_status": (
                correction_output.corrected_cell_status is not None
            ),
        }
    )
    return PreprocessingStageExecution(
        stage=DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION,
        operation=str(correction_output.batch_correction_report.method),
        parameters={
            "source": "resolved_correction_result",
            "stage_order": list(correction_output.stage_order),
        },
        input_shape=(int(previous.phospho.shape[0]), int(previous.phospho.shape[1])),
        output_shape=(int(current.phospho.shape[0]), int(current.phospho.shape[1])),
        input_hash=_hash_stage_table_fingerprints(
            direction="input",
            table_fingerprints=consumed_input_tables,
        ),
        output_hash=_hash_stage_table_fingerprints(
            direction="output",
            table_fingerprints=produced_output_tables,
        ),
        phospho_input_hash=hash_table_tolerance(
            previous.phospho,
            name="batch_correction.input.phospho",
        ),
        phospho_output_hash=hash_table_tolerance(
            current.phospho,
            name="batch_correction.output.phospho",
        ),
        dropped_row_ids=(),
        dropped_row_count=0,
        schema_version=PREPROCESSING_STAGE_PROVENANCE_SCHEMA_VERSION_V3,
        consumed_input_tables=consumed_input_tables,
        produced_output_tables=produced_output_tables,
        backend="resolved_correction_result",
        random_seed=None,
        determinism=PREPROCESSING_STAGE_DETERMINISM_PURE,
        is_deterministic=True,
        imputed_cell_count=0,
        imputed_row_ids=(),
        notes="resolved correction output integrated before dataset boundary",
        diagnostics=diagnostics,
        batch_correction_provenance=_resolve_batch_correction_provenance(
            previous=previous,
            current=current,
            correction_output=correction_output,
        ),
    )


def _resolve_batch_correction_provenance(
    *,
    previous: PreprocessingState,
    current: PreprocessingState,
    correction_output: CorrectedPreprocessingOutput,
) -> BatchCorrectionProvenance:
    provided = correction_output.provenance
    if isinstance(provided, BatchCorrectionProvenance):
        return provided
    if isinstance(provided, Mapping) and _looks_like_batch_provenance_payload(provided):
        return batch_correction_provenance_from_payload(provided)
    diagnostics: dict[str, object] = dict(correction_output.diagnostics)
    if isinstance(provided, Mapping):
        diagnostics["provided_untyped_provenance"] = dict(provided)
    return build_native_batch_correction_provenance(
        input_matrix=previous.phospho,
        output_matrix=current.phospho,
        plan=previous.plan,
        report=correction_output.batch_correction_report,
        metadata=None,
        diagnostics=diagnostics,
        warnings=correction_output.batch_correction_report.warnings,
        observation_mask=current.imputation_observation_mask,
        corrected_cell_status=correction_output.corrected_cell_status,
        control_site_source={
            "source_type": "not_provided",
            "reason": (
                "resolved correction output did not provide typed control-site "
                "provenance"
            ),
        },
        selected_site_key_rows=(),
        source="resolved_correction_output",
    )


def _looks_like_batch_provenance_payload(value: Mapping[str, object]) -> bool:
    return "schema_version" in value and "requested_method" in value


def _collect_fingerprints(
    state: PreprocessingState,
    table_names: tuple[PreprocessingStateTableKey, ...],
) -> tuple[TableFingerprint, ...]:
    fingerprints: list[TableFingerprint] = []
    for table_name in table_names:
        table = _resolve_state_table(state=state, table_name=table_name)
        fingerprint = fingerprint_optional_table(table, name=table_name.value)
        if fingerprint is not None:
            fingerprints.append(fingerprint)
    return tuple(fingerprints)


def _resolve_state_table(
    *,
    state: PreprocessingState,
    table_name: PreprocessingStateTableKey,
) -> pd.DataFrame | None:
    if table_name is PreprocessingStateTableKey.DATASET_PHOSPHO:
        return state.phospho
    if table_name is PreprocessingStateTableKey.DATASET_IMPUTATION_OBSERVATION_MASK:
        return state.imputation_observation_mask
    return None


def _hash_stage_table_fingerprints(
    *,
    direction: str,
    table_fingerprints: tuple[TableFingerprint, ...],
) -> str:
    payload = {
        "stage": DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION,
        "direction": direction,
        "tables": [
            {
                "name": item.name,
                "rows": int(item.rows),
                "columns": int(item.columns),
                "exact_hash_algorithm": item.exact_hash_algorithm,
                "exact_hash_value": item.exact_hash_value,
                "tolerance_hash_algorithm": item.tolerance_hash_algorithm,
                "tolerance_hash_value": item.tolerance_hash_value,
            }
            for item in table_fingerprints
        ],
    }
    encoded = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _require_numeric_complete_matrix(matrix: pd.DataFrame) -> None:
    if matrix.empty:
        raise PhosPyInputError(
            "corrected preprocessing output.corrected_matrix must not be empty"
        )
    non_numeric = [
        str(column)
        for column in matrix.columns.tolist()
        if (
            not pd.api.types.is_numeric_dtype(matrix.loc[:, column])
            or pd.api.types.is_bool_dtype(matrix.loc[:, column])
        )
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
    if isinstance(value, Mapping):
        return value
    return {}


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


def _format_labels(labels: Sequence[str]) -> str:
    preview = ", ".join(repr(value) for value in tuple(labels)[:5])
    suffix = "" if len(labels) <= 5 else " ..."
    return f"{preview}{suffix}"


def _require_stage_order_precedes_downstream(stage_order: Sequence[str]) -> None:
    normalized = tuple(str(stage).strip() for stage in stage_order)
    if DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION not in normalized:
        raise PhosPyInputError(
            "corrected preprocessing output.stage_order must include 'batch_correction'"
        )
    if _DOWNSTREAM_WORKFLOWS_STAGE not in normalized:
        return
    correction_index = normalized.index(DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION)
    downstream_index = normalized.index(_DOWNSTREAM_WORKFLOWS_STAGE)
    if correction_index > downstream_index:
        raise PhosPyInputError(
            "corrected preprocessing output.stage_order must place "
            "batch_correction before downstream_workflows"
        )


def _json_mapping(value: Mapping[str, object]) -> dict[str, JsonValue]:
    return {str(key): _json_value(item) for key, item in value.items()}


def _json_value(value: object) -> JsonValue:
    if value is None or isinstance(value, str | int | float | bool):
        return value
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, Sequence) and not isinstance(value, str | bytes | bytearray):
        return [_json_value(item) for item in value]
    if isinstance(value, TableFingerprint):
        return _json_value(
            {
                "name": value.name,
                "rows": value.rows,
                "columns": value.columns,
                "index_name": value.index_name,
                "column_names": list(value.column_names),
                "dtypes": list(value.dtypes),
                "exact_hash_algorithm": value.exact_hash_algorithm,
                "exact_hash_value": value.exact_hash_value,
                "tolerance_hash_algorithm": value.tolerance_hash_algorithm,
                "tolerance_hash_value": value.tolerance_hash_value,
            }
        )
    return str(value)


__all__ = [
    "CorrectedPreprocessingOutput",
    "CorrectedPreprocessingOutputIntegrator",
    "validate_corrected_preprocessing_output",
]
