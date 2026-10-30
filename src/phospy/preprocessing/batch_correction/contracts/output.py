"""Corrected-output contract and supported external-output factory."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from typing import Literal, Protocol, cast, runtime_checkable

import numpy as np
import pandas as pd

from phospy.core.errors.input import PhosPyInputError
from phospy.core.frames.ownership import (
    export_dataframe,
    own_dataframe,
    own_optional_dataframe,
)
from phospy.core.preparation_records.batch_correction import (
    BATCH_CORRECTION_CONFOUNDING_NOT_CHECKED,
    BATCH_CORRECTION_CONFOUNDING_PASSED,
    BATCH_CORRECTION_DESIGN_PRESERVATION_PRESERVE_CONDITION_EFFECTS,
    BATCH_CORRECTION_STATUS_APPLIED,
    BatchCorrectionDiagnostics,
    BatchCorrectionPolicy,
    BatchCorrectionReport,
)
from phospy.core.provenance.hashing import fingerprint_matrix
from phospy.core.provenance.models.batch_correction import BatchCorrectionProvenance
from phospy.core.provenance.scalar_validation import JsonValue
from phospy.preprocessing.control_sites.contracts import ControlSiteSet
from phospy.preprocessing.control_sites.validation import (
    ControlSiteEligibilityValidator,
)

_BATCH_CORRECTION_STAGE = "batch_correction"
_DOWNSTREAM_WORKFLOWS_STAGE = "downstream_workflows"
_SUPPORTED_STATUS_VALUES = frozenset({"corrected_observed", "restored_missing"})

_CALLER_ATTESTED_EXTERNAL_CORRECTION_ORIGIN = "caller_supplied_external_correction"
_CALLER_ATTESTED_EXTERNAL_CORRECTION_VERIFICATION = (
    "caller_attested_not_verified_by_phospy"
)
_CALLER_ATTESTED_EXTERNAL_CONTROL_ATTESTATION = "caller_attested_external_controls"
_CALLER_ATTESTED_EXTERNAL_ENVIRONMENT = "external_execution_not_observed_by_phospy"
_SUPPORTED_CALLER_ATTESTED_EXTERNAL_METHODS = frozenset(
    {"sps_ruv_style", "ruv_iii_style"}
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


def create_external_corrected_preprocessing_output(
    *,
    pre_correction_matrix: pd.DataFrame,
    corrected_matrix: pd.DataFrame,
    method: str,
    source_software: str,
    source_version: str,
    control_site_set: ControlSiteSet,
    n_unwanted_factors: int,
    batch_column: str,
    condition_columns: Sequence[str],
    replicate_column: str | None = None,
    preserve_condition_effects: Literal[True] = True,
    parameters: Mapping[str, object] | None = None,
    diagnostics: Mapping[str, object] | None = None,
    output_observation_mask: pd.DataFrame | None = None,
    warnings: Sequence[str] = (),
) -> CorrectedPreprocessingOutput:
    """Create a validated caller-attested output from an external correction.

    The factory records what the caller declares and binds the result to the
    supplied pre-correction and corrected matrices. It does not represent the
    external operation as a PhosPy-executed or PhosPy-verified correction.
    """

    before = own_dataframe(
        pre_correction_matrix,
        field_name="pre_correction_matrix",
        error_type=PhosPyInputError,
    )
    corrected = own_dataframe(
        corrected_matrix,
        field_name="corrected_matrix",
        error_type=PhosPyInputError,
    )
    _require_numeric_complete_matrix(before, field_name="pre_correction_matrix")
    _require_numeric_complete_matrix(corrected, field_name="corrected_matrix")
    _require_unique_axes(before, field_name="pre_correction_matrix")
    _require_unique_axes(corrected, field_name="corrected_matrix")
    _require_aligned_frame(
        corrected,
        before,
        field_name="corrected_matrix",
        expected_description="pre_correction_matrix",
    )

    normalized_method = _required_text(method, field_name="method").lower()
    if normalized_method == "none":
        raise PhosPyInputError(
            "method='none' cannot describe an externally corrected matrix"
        )
    if normalized_method not in _SUPPORTED_CALLER_ATTESTED_EXTERNAL_METHODS:
        supported = ", ".join(sorted(_SUPPORTED_CALLER_ATTESTED_EXTERNAL_METHODS))
        raise PhosPyInputError(
            "method has no registered strict external-correction validator; "
            f"supported methods: {supported}"
        )
    normalized_source_software = _required_text(
        source_software,
        field_name="source_software",
    )
    normalized_source_version = _required_text(
        source_version,
        field_name="source_version",
    )
    normalized_batch_column = _required_text(
        batch_column,
        field_name="batch_column",
    )
    normalized_condition_columns = _string_sequence(
        condition_columns,
        field_name="condition_columns",
    )
    if not normalized_condition_columns:
        raise PhosPyInputError(
            "condition_columns must contain at least one protected condition column"
        )
    if isinstance(n_unwanted_factors, bool) or not isinstance(
        n_unwanted_factors,
        int,
    ):
        raise PhosPyInputError("n_unwanted_factors must be a positive integer")
    if n_unwanted_factors < 1:
        raise PhosPyInputError("n_unwanted_factors must be at least 1")
    normalized_replicate_column = _optional_required_text(
        replicate_column,
        field_name="replicate_column",
    )
    if normalized_method == "ruv_iii_style" and normalized_replicate_column is None:
        raise PhosPyInputError(
            "replicate_column is required for external ruv_iii_style correction"
        )
    if preserve_condition_effects is not True:
        raise PhosPyInputError(
            "preserve_condition_effects must be True for externally corrected "
            "dataset input"
        )
    normalized_warnings = _string_sequence(warnings, field_name="warnings")
    normalized_parameters = _json_mapping_copy(parameters, field_name="parameters")
    normalized_diagnostics = _json_mapping_copy(
        diagnostics,
        field_name="diagnostics",
    )

    required_control_count = (
        max(2, n_unwanted_factors)
        if normalized_method == "ruv_iii_style"
        else n_unwanted_factors + 1
    )
    if not isinstance(control_site_set, ControlSiteSet):
        raise PhosPyInputError("control_site_set must be a ControlSiteSet")
    control_mapping = ControlSiteEligibilityValidator().run(
        control_set=control_site_set,
        method=normalized_method,
        min_eligible_controls=required_control_count,
        site_keys=corrected.index.tolist(),
        n_unwanted_factors=n_unwanted_factors,
        expected_identifier_namespace="site_key",
        control_site_source_type=control_site_set.source_metadata.source_type,
    )
    selected_control_site_keys = tuple(
        str(row.site_key)
        for row in control_mapping.row_eligibility
        if row.is_control and row.site_key is not None
    )

    if output_observation_mask is None:
        mask = pd.DataFrame(True, index=corrected.index, columns=corrected.columns)
        mask_source = "all_cells_observed_declared_by_caller"
    else:
        mask = own_dataframe(
            output_observation_mask,
            field_name="output_observation_mask",
            error_type=PhosPyInputError,
        )
        _require_aligned_frame(
            mask,
            corrected,
            field_name="output_observation_mask",
            expected_description="corrected_matrix",
        )
        _require_boolean_mask(mask, field_name="output_observation_mask")
        mask_source = "caller_supplied"
    status = pd.DataFrame(
        np.where(
            mask.to_numpy(dtype=bool, copy=True),
            "corrected_observed",
            "restored_missing",
        ),
        index=corrected.index.copy(),
        columns=corrected.columns.copy(),
    )

    policy = BatchCorrectionPolicy(
        method=normalized_method,
        batch_column=normalized_batch_column,
        condition_column=None,
        condition_columns=normalized_condition_columns,
        design_preservation_policy=(
            BATCH_CORRECTION_DESIGN_PRESERVATION_PRESERVE_CONDITION_EFFECTS
        ),
        preserve_condition_effects=True,
    )
    shape = (int(corrected.shape[0]), int(corrected.shape[1]))
    limitations = (
        "correction was executed outside PhosPy and is caller-attested",
        "PhosPy validates the supplied state, axes, and fingerprints but does "
        "not verify the external scientific method or execution environment",
    )
    report = BatchCorrectionReport(
        status=BATCH_CORRECTION_STATUS_APPLIED,
        policy=policy,
        diagnostics=BatchCorrectionDiagnostics(
            confounding_check_status=BATCH_CORRECTION_CONFOUNDING_NOT_CHECKED,
            matrix_shape_before=shape,
            matrix_shape_after=shape,
            warnings=normalized_warnings,
            limitations=limitations,
        ),
    )
    provenance = BatchCorrectionProvenance(
        requested_method=normalized_method,
        resolved_parameters=cast(
            Mapping[str, JsonValue],
            {
                "execution_origin": _CALLER_ATTESTED_EXTERNAL_CORRECTION_ORIGIN,
                "execution_verification": (
                    _CALLER_ATTESTED_EXTERNAL_CORRECTION_VERIFICATION
                ),
                "source_software": normalized_source_software,
                "source_version": normalized_source_version,
                "n_unwanted_factors": n_unwanted_factors,
                "external_method_parameters": normalized_parameters,
            },
        ),
        preprocessing_stage_order=(
            "missing_data",
            _BATCH_CORRECTION_STAGE,
            _DOWNSTREAM_WORKFLOWS_STAGE,
        ),
        control_site_source=cast(
            Mapping[str, JsonValue],
            {
                **control_site_set.source_metadata.to_payload(),
                "selection_performed_by_phospy": False,
                "external_control_attestation": (
                    _CALLER_ATTESTED_EXTERNAL_CONTROL_ATTESTATION
                ),
            },
        ),
        selected_site_key_rows=selected_control_site_keys,
        batch_metadata={
            "batch_column": normalized_batch_column,
            "source": "caller_attested",
        },
        replicate_metadata=(
            None
            if normalized_replicate_column is None
            else {
                "replicate_column": normalized_replicate_column,
                "used_for_numerical_factor_estimation": (
                    normalized_method == "ruv_iii_style"
                ),
                "ruv_iii_semantics_enabled": (normalized_method == "ruv_iii_style"),
                "source": "caller_attested",
            }
        ),
        design_metadata={
            "condition_columns": list(policy.condition_columns or ()),
            "design_preservation_policy": policy.design_preservation_policy,
            "preserve_condition_effects": True,
            "source": "caller_attested",
        },
        missing_value_policy={
            "output_observation_mask_source": mask_source,
            "corrected_matrix_requires_complete_finite_values": True,
        },
        observation_masks=(
            fingerprint_matrix(
                mask.astype("int8"),
                name="external_correction.output_observation_mask",
            ),
            fingerprint_matrix(
                status,
                name="external_correction.corrected_cell_status",
            ),
        ),
        input_matrix_fingerprint=fingerprint_matrix(
            before,
            name="external_correction.pre_correction_matrix",
        ),
        output_matrix_fingerprint=fingerprint_matrix(
            corrected,
            name="external_correction.corrected_matrix",
        ),
        diagnostics=cast(Mapping[str, JsonValue], normalized_diagnostics),
        warnings=normalized_warnings,
        phospy_version=_CALLER_ATTESTED_EXTERNAL_ENVIRONMENT,
        python_version=_CALLER_ATTESTED_EXTERNAL_ENVIRONMENT,
        dependency_versions={normalized_source_software: normalized_source_version},
        imputation_policy={
            "source": "caller_attested",
            "output_observation_mask_source": mask_source,
        },
    )
    return CorrectedPreprocessingOutput(
        corrected_matrix=corrected,
        batch_correction_report=report,
        diagnostics={
            "execution_origin": _CALLER_ATTESTED_EXTERNAL_CORRECTION_ORIGIN,
            "execution_verification": (
                _CALLER_ATTESTED_EXTERNAL_CORRECTION_VERIFICATION
            ),
            "source_software": normalized_source_software,
            "source_version": normalized_source_version,
            "external_diagnostics": normalized_diagnostics,
        },
        output_observation_mask=mask,
        corrected_cell_status=status,
        provenance=provenance,
    )


def _require_numeric_complete_matrix(
    matrix: pd.DataFrame,
    *,
    field_name: str = "corrected preprocessing output.corrected_matrix",
) -> None:
    if matrix.empty:
        raise PhosPyInputError(f"{field_name} must not be empty")
    non_numeric = [
        str(column)
        for column in matrix.columns.tolist()
        if not pd.api.types.is_numeric_dtype(matrix.loc[:, column])
        or pd.api.types.is_bool_dtype(matrix.loc[:, column])
    ]
    if non_numeric:
        raise PhosPyInputError(
            f"{field_name} must contain numeric non-bool columns. "
            "Non-numeric columns: " + ", ".join(non_numeric)
        )
    values = matrix.astype("float64").to_numpy(copy=True)
    if not np.isfinite(values).all():
        raise PhosPyInputError(
            f"{field_name} must be finite and missing-value-free for "
            "AnalysisReadyPhosphoDataset construction"
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
    expected_description: str = "dataset.phospho",
) -> None:
    if not frame.index.equals(expected.index):
        raise PhosPyInputError(
            f"{field_name}.index must match {expected_description}.index"
        )
    if not frame.columns.equals(expected.columns):
        raise PhosPyInputError(
            f"{field_name}.columns must match {expected_description}.columns"
        )


def _require_unique_axes(matrix: pd.DataFrame, *, field_name: str) -> None:
    if not matrix.index.is_unique:
        raise PhosPyInputError(f"{field_name}.index must contain unique labels")
    if not matrix.columns.is_unique:
        raise PhosPyInputError(f"{field_name}.columns must contain unique labels")


def _require_boolean_mask(
    mask: pd.DataFrame,
    *,
    field_name: str = "corrected preprocessing output.output_observation_mask",
) -> None:
    for row_position, row_id in enumerate(mask.index.tolist()):
        for column_position, column_id in enumerate(mask.columns.tolist()):
            value = mask.iat[row_position, column_position]
            if isinstance(value, (bool, np.bool_)):
                continue
            raise PhosPyInputError(
                f"{field_name} must contain only boolean values; "
                f"invalid_cell=({row_id!r}, {column_id!r})"
            )


def _required_text(value: object, *, field_name: str) -> str:
    if not isinstance(value, str):
        raise PhosPyInputError(f"{field_name} must be a string")
    normalized = value.strip()
    if normalized == "":
        raise PhosPyInputError(f"{field_name} must be a non-empty string")
    return normalized


def _optional_required_text(value: object, *, field_name: str) -> str | None:
    if value is None:
        return None
    return _required_text(value, field_name=field_name)


def _string_sequence(value: object, *, field_name: str) -> tuple[str, ...]:
    if isinstance(value, str | bytes | bytearray) or not isinstance(value, Sequence):
        raise PhosPyInputError(f"{field_name} must be a sequence of strings")
    normalized = tuple(_required_text(item, field_name=field_name) for item in value)
    return normalized


def _json_mapping_copy(
    value: object,
    *,
    field_name: str,
) -> dict[str, object]:
    if value is None:
        return {}
    if not isinstance(value, Mapping):
        raise PhosPyInputError(f"{field_name} must be a mapping")
    copied: dict[str, object] = {}
    for key, item in value.items():
        if not isinstance(key, str) or key.strip() == "":
            raise PhosPyInputError(
                f"{field_name} must contain only non-empty string keys"
            )
        copied[key] = item
    return copied


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


__all__ = [
    "CorrectedPreprocessingOutput",
    "create_external_corrected_preprocessing_output",
]
