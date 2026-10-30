"""Scientific validation for passive corrected preprocessing outputs."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

from phospy.core.errors.input import PhosPyInputError
from phospy.core.preparation_records.batch_correction import (
    BATCH_CORRECTION_STATUS_APPLIED,
)
from phospy.core.provenance.models.batch_correction import BatchCorrectionProvenance
from phospy.preprocessing.batch_correction.contracts.output import (
    _CALLER_ATTESTED_EXTERNAL_CONTROL_ATTESTATION,
    _CALLER_ATTESTED_EXTERNAL_CORRECTION_ORIGIN,
    _CALLER_ATTESTED_EXTERNAL_CORRECTION_VERIFICATION,
    _CALLER_ATTESTED_EXTERNAL_ENVIRONMENT,
    CorrectedPreprocessingOutput,
)
from phospy.preprocessing.batch_correction.provenance.validation import (
    validate_applied_native_sps_ruv_correction_provenance,
)
from phospy.preprocessing.batch_correction.validation.controls import (
    normalize_applied_selected_site_key_rows,
    require_control_site_source_metadata,
    require_unique_selected_control_site_rows,
)

_SUPPORTED_EXTERNAL_CORRECTED_OUTPUT_METHODS = frozenset(
    {"sps_ruv_style", "ruv_iii_style"}
)
_EXTERNAL_CORRECTED_OUTPUT_PROVENANCE_REQUIREMENT = (
    "external corrected outputs require applied correction status and a supported "
    "applied correction method with method-specific provenance"
)


def validate_corrected_preprocessing_output(value: object) -> None:
    """Validate a dataset-build request correction output payload."""

    if value is None:
        return
    if not isinstance(value, CorrectedPreprocessingOutput):
        raise PhosPyInputError(
            "dataset build request corrected_preprocessing_output must be "
            "CorrectedPreprocessingOutput"
        )
    validate_corrected_output_report_and_provenance(value)


def validate_corrected_output_report_and_provenance(
    correction_output: CorrectedPreprocessingOutput,
) -> None:
    """Require applied, method-specific typed provenance."""

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
    if _is_caller_attested_external_correction(correction_output.provenance):
        _validate_caller_attested_external_correction(
            correction_output=correction_output,
            method=method,
        )
        return
    validate_applied_native_sps_ruv_correction_provenance(
        method=method,
        status=status,
        provenance=correction_output.provenance,
    )


def _is_caller_attested_external_correction(provenance: object) -> bool:
    if not isinstance(provenance, BatchCorrectionProvenance):
        return False
    return (
        provenance.resolved_parameters.get("execution_origin")
        == _CALLER_ATTESTED_EXTERNAL_CORRECTION_ORIGIN
    )


def _validate_caller_attested_external_correction(
    *,
    correction_output: CorrectedPreprocessingOutput,
    method: str,
) -> None:
    provenance = correction_output.provenance
    if not isinstance(provenance, BatchCorrectionProvenance):
        raise PhosPyInputError(
            "caller-attested external corrected output requires typed "
            "BatchCorrectionProvenance"
        )
    parameters = provenance.resolved_parameters
    if (
        parameters.get("execution_verification")
        != _CALLER_ATTESTED_EXTERNAL_CORRECTION_VERIFICATION
    ):
        raise PhosPyInputError(
            "caller-attested external corrected output provenance must state "
            "that PhosPy did not verify the external execution"
        )
    if str(provenance.requested_method).strip().lower() != method:
        raise PhosPyInputError(
            "caller-attested external corrected output provenance method must "
            "match the batch-correction report method"
        )
    for field_name in ("source_software", "source_version"):
        value = parameters.get(field_name)
        if not isinstance(value, str) or value.strip() == "":
            raise PhosPyInputError(
                "caller-attested external corrected output provenance requires "
                f"a non-empty {field_name}"
            )
    report = correction_output.batch_correction_report
    if report.policy.preserve_condition_effects is not True:
        raise PhosPyInputError(
            "caller-attested external corrected output must preserve condition effects"
        )
    if provenance.design_metadata.get("preserve_condition_effects") is not True:
        raise PhosPyInputError(
            "caller-attested external corrected output provenance must record "
            "preserve_condition_effects=True"
        )
    selected_control_rows = normalize_applied_selected_site_key_rows(
        provenance.selected_site_key_rows
    )
    require_unique_selected_control_site_rows(selected_control_rows)
    n_unwanted_factors = parameters.get("n_unwanted_factors")
    if isinstance(n_unwanted_factors, bool) or not isinstance(
        n_unwanted_factors,
        int,
    ):
        raise PhosPyInputError(
            "caller-attested external corrected output provenance requires a "
            "positive integer n_unwanted_factors"
        )
    if n_unwanted_factors < 1:
        raise PhosPyInputError(
            "caller-attested external corrected output provenance requires "
            "n_unwanted_factors >= 1"
        )
    required_control_count = (
        max(2, n_unwanted_factors)
        if method == "ruv_iii_style"
        else n_unwanted_factors + 1
    )
    if len(selected_control_rows) < required_control_count:
        raise PhosPyInputError(
            "caller-attested external corrected output provenance records too "
            "few controls for n_unwanted_factors"
        )
    require_control_site_source_metadata(
        provenance.control_site_source,
        selected_site_key_rows=selected_control_rows,
    )
    if (
        provenance.control_site_source.get("external_control_attestation")
        != _CALLER_ATTESTED_EXTERNAL_CONTROL_ATTESTATION
        or provenance.control_site_source.get("selection_performed_by_phospy")
        is not False
    ):
        raise PhosPyInputError(
            "caller-attested external corrected output provenance must not claim "
            "PhosPy control-site selection"
        )
    if provenance.batch_metadata.get("batch_column") != report.batch_column:
        raise PhosPyInputError(
            "caller-attested external corrected output provenance batch column "
            "must match the batch-correction report"
        )
    recorded_condition_columns = provenance.design_metadata.get("condition_columns")
    if isinstance(
        recorded_condition_columns,
        str | bytes | bytearray,
    ) or not isinstance(recorded_condition_columns, Sequence):
        raise PhosPyInputError(
            "caller-attested external corrected output provenance condition "
            "columns must be a sequence"
        )
    if tuple(recorded_condition_columns) != tuple(report.condition_columns):
        raise PhosPyInputError(
            "caller-attested external corrected output provenance condition "
            "columns must match the batch-correction report"
        )
    if (
        provenance.design_metadata.get("design_preservation_policy")
        != report.design_preservation_policy
    ):
        raise PhosPyInputError(
            "caller-attested external corrected output provenance design policy "
            "must match the batch-correction report"
        )
    if method == "ruv_iii_style":
        replicate_metadata = provenance.replicate_metadata
        if (
            not isinstance(replicate_metadata, Mapping)
            or not isinstance(replicate_metadata.get("replicate_column"), str)
            or replicate_metadata.get("used_for_numerical_factor_estimation")
            is not True
            or replicate_metadata.get("ruv_iii_semantics_enabled") is not True
        ):
            raise PhosPyInputError(
                "caller-attested external RUV-III provenance requires active "
                "replicate metadata"
            )
    expected_stage_order = (
        "missing_data",
        "batch_correction",
        "downstream_workflows",
    )
    if provenance.preprocessing_stage_order != expected_stage_order:
        raise PhosPyInputError(
            "caller-attested external corrected output provenance must use the "
            "supported preprocessing stage order"
        )
    if tuple(correction_output.stage_order) != expected_stage_order:
        raise PhosPyInputError(
            "caller-attested external corrected output must use the supported "
            "preprocessing stage order"
        )
    if provenance.output_matrix_fingerprint is None:
        raise PhosPyInputError(
            "caller-attested external corrected output provenance requires input "
            "and output matrix fingerprints"
        )
    observed_shape = (
        int(correction_output.corrected_matrix.shape[0]),
        int(correction_output.corrected_matrix.shape[1]),
    )
    if report.matrix_shape_before != observed_shape:
        raise PhosPyInputError(
            "caller-attested external corrected output report input shape must "
            "match corrected matrix shape"
        )
    if report.matrix_shape_after != observed_shape:
        raise PhosPyInputError(
            "caller-attested external corrected output report output shape must "
            "match corrected matrix shape"
        )
    if (
        provenance.input_matrix_fingerprint.rows,
        provenance.input_matrix_fingerprint.columns,
    ) != observed_shape or (
        provenance.output_matrix_fingerprint.rows,
        provenance.output_matrix_fingerprint.columns,
    ) != observed_shape:
        raise PhosPyInputError(
            "caller-attested external corrected output provenance fingerprint "
            "shapes must match the corrected matrix"
        )
    if tuple(provenance.warnings) != tuple(report.warnings):
        raise PhosPyInputError(
            "caller-attested external corrected output provenance warnings must "
            "match the batch-correction report"
        )
    mask_names = tuple(
        fingerprint.name.strip().lower() for fingerprint in provenance.observation_masks
    )
    if not any("observation_mask" in name for name in mask_names):
        raise PhosPyInputError(
            "caller-attested external corrected output provenance requires an "
            "observation-mask fingerprint"
        )
    if not any("corrected_cell_status" in name for name in mask_names):
        raise PhosPyInputError(
            "caller-attested external corrected output provenance requires a "
            "corrected-cell-status fingerprint"
        )
    if (
        provenance.phospy_version != _CALLER_ATTESTED_EXTERNAL_ENVIRONMENT
        or provenance.python_version != _CALLER_ATTESTED_EXTERNAL_ENVIRONMENT
    ):
        raise PhosPyInputError(
            "caller-attested external corrected output provenance must not record "
            "the external environment as a PhosPy-observed execution environment"
        )
    source_software = parameters["source_software"]
    source_version = parameters["source_version"]
    if provenance.dependency_versions.get(str(source_software)) != source_version:
        raise PhosPyInputError(
            "caller-attested external corrected output provenance software "
            "version must match dependency_versions"
        )


def validate_selected_controls_in_corrected_matrix(
    correction_output: CorrectedPreprocessingOutput,
) -> None:
    """Require every selected control row in an applied output matrix."""

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


def _normalize_external_correction_method(method: object) -> str:
    normalized = str(method).strip().lower()
    if normalized == "":
        raise PhosPyInputError(
            "corrected_preprocessing_output "
            f"{_EXTERNAL_CORRECTED_OUTPUT_PROVENANCE_REQUIREMENT}; "
            "observed method is missing or empty"
        )
    return normalized


def _format_labels(labels: tuple[str, ...]) -> str:
    preview = ", ".join(repr(value) for value in labels[:5])
    suffix = "" if len(labels) <= 5 else " ..."
    return f"{preview}{suffix}"


__all__ = [
    "validate_corrected_output_report_and_provenance",
    "validate_corrected_preprocessing_output",
    "validate_selected_controls_in_corrected_matrix",
]
