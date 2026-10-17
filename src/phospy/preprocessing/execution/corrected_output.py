"""Preprocessing boundary integration for resolved correction outputs."""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping, Sequence
from dataclasses import replace

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.preprocessing.batch_correction.contracts.output import (
    CorrectedPreprocessingOutput as _CorrectedPreprocessingOutput,
)
from phospy.preprocessing.batch_correction.provenance.assembly import (
    build_native_batch_correction_provenance as _build_native_batch_correction_provenance,
)
from phospy.preprocessing.batch_correction.validation.output import (
    validate_corrected_output_report_and_provenance as _validate_corrected_output_report_and_provenance,
)
from phospy.preprocessing.batch_correction.validation.output import (
    validate_selected_controls_in_corrected_matrix as _validate_selected_controls_in_corrected_matrix,
)
from phospy.preprocessing.contracts import (
    DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION,
    PreprocessingStageExecution,
    PreprocessingState,
    PreprocessingStateTableKey,
)
from phospy.preprocessing.provenance.models import (
    PREPROCESSING_EXTERNAL_NONDETERMINISM_CAVEAT_CODE,
    PREPROCESSING_STAGE_PROVENANCE_SCHEMA_VERSION_V3,
    DeterminismKind,
    ReproducibilityCaveat,
)
from phospy.provenance.hashing import (
    fingerprint_matrix,
    fingerprint_optional_table,
    hash_table_tolerance,
)
from phospy.provenance.models import (
    BatchCorrectionProvenance,
    JsonValue,
    TableFingerprint,
)
from phospy.provenance.serialization.batch_correction import (
    batch_correction_provenance_from_payload,
)

_DOWNSTREAM_WORKFLOWS_STAGE = "downstream_workflows"


class CorrectedPreprocessingOutputIntegrator:
    """Apply a resolved correction output to preprocessing state."""

    def run(
        self,
        *,
        state: PreprocessingState,
        correction_output: _CorrectedPreprocessingOutput,
    ) -> tuple[PreprocessingState, PreprocessingStageExecution]:
        if not isinstance(correction_output, _CorrectedPreprocessingOutput):
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


def _validate_applied_provenance_fingerprints(
    *,
    state: PreprocessingState,
    correction_output: _CorrectedPreprocessingOutput,
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
    correction_output: _CorrectedPreprocessingOutput,
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
    correction_output: _CorrectedPreprocessingOutput,
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
        determinism=DeterminismKind.EXTERNALLY_NONDETERMINISTIC,
        reproducibility_caveats=_external_correction_reproducibility_caveats(
            method=str(correction_output.batch_correction_report.method),
        ),
        is_deterministic=False,
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


def _external_correction_reproducibility_caveats(
    *,
    method: str,
) -> tuple[ReproducibilityCaveat, ...]:
    return (
        ReproducibilityCaveat(
            code=PREPROCESSING_EXTERNAL_NONDETERMINISM_CAVEAT_CODE,
            severity="warning",
            message=(
                "Corrected preprocessing output was produced outside the dataset "
                "preprocessing pipeline; exact reproduction requires the external "
                "execution environment, runtime state, and inputs used to produce "
                "the recorded matrix."
            ),
            details={
                "stage": DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION,
                "operation": method,
                "backend": "resolved_correction_result",
                "determinism_kind": (DeterminismKind.EXTERNALLY_NONDETERMINISTIC.value),
            },
        ),
    )


def _resolve_batch_correction_provenance(
    *,
    previous: PreprocessingState,
    current: PreprocessingState,
    correction_output: _CorrectedPreprocessingOutput,
) -> BatchCorrectionProvenance:
    provided = correction_output.provenance
    if isinstance(provided, BatchCorrectionProvenance):
        return provided
    if isinstance(provided, Mapping) and _looks_like_batch_provenance_payload(provided):
        return batch_correction_provenance_from_payload(provided)
    diagnostics: dict[str, object] = dict(correction_output.diagnostics)
    if isinstance(provided, Mapping):
        diagnostics["provided_untyped_provenance"] = dict(provided)
    return _build_native_batch_correction_provenance(
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


__all__ = ["CorrectedPreprocessingOutputIntegrator"]
