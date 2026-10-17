"""Scientific validation for passive corrected preprocessing outputs."""

from __future__ import annotations

from phospy.errors.input import PhosPyInputError
from phospy.preparation_records.batch_correction import (
    BATCH_CORRECTION_STATUS_APPLIED,
)
from phospy.preprocessing.batch_correction.contracts.output import (
    CorrectedPreprocessingOutput,
)
from phospy.preprocessing.batch_correction.provenance.validation import (
    validate_applied_native_sps_ruv_correction_provenance,
)
from phospy.preprocessing.batch_correction.validation.controls import (
    normalize_applied_selected_site_key_rows,
)
from phospy.provenance.models.batch_correction import BatchCorrectionProvenance

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
    validate_applied_native_sps_ruv_correction_provenance(
        method=method,
        status=status,
        provenance=correction_output.provenance,
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
