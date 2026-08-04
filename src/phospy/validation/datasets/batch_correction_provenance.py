"""Applied SPS/RUV-style provenance validation for batch correction."""
# pyright: reportUnnecessaryIsInstance=false
# Runtime boundary guards are intentionally retained for untyped external callers.

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import cast

from phospy.contracts.configs.preprocessing.internal_batch_correction import (
    SUPPORTED_INTERNAL_BATCH_CORRECTION_EXECUTED_STAGE_ORDER,
)
from phospy.errors.input import PhosPyInputError
from phospy.provenance.environment import BATCH_CORRECTION_ENVIRONMENT_DEPENDENCIES
from phospy.provenance.models import BatchCorrectionProvenance
from phospy.validation.datasets._batch_correction_helpers import (
    MISSING,
    duplicates_in_order,
    format_labels,
    is_missing_environment_text,
    reject_not_provided_required_mapping,
    require_non_empty_mapping,
    unique_in_order,
)
from phospy.validation.datasets.batch_correction_controls import (
    normalize_applied_selected_site_key_rows,
    require_control_site_source_metadata,
    require_unique_selected_control_site_rows,
)

_APPLIED_NATIVE_SPS_RUV_METHODS = frozenset({"sps_ruv_style"})
_UNSUPPORTED_SPS_RUV_METHODS = frozenset({"ruv_iii_style"})
_MISSING_PROVENANCE_MESSAGE = (
    "corrected_preprocessing_output with applied native SPS/RUV-style correction "
    "requires typed BatchCorrectionProvenance"
)


def validate_applied_native_sps_ruv_correction_provenance(
    *,
    method: object,
    status: object,
    provenance: object,
) -> None:
    """Require complete typed provenance for applied SPS/RUV-style outputs."""

    if str(status).strip() != "applied":
        return
    normalized_method = _normalise_method(method)
    if normalized_method in _UNSUPPORTED_SPS_RUV_METHODS:
        raise PhosPyInputError(
            "corrected_preprocessing_output declares unsupported SPS/RUV-style "
            f"batch correction method {normalized_method!r}; regenerate with a "
            "supported native method and complete BatchCorrectionProvenance"
        )
    if not _is_sps_ruv_style_label(normalized_method):
        return
    if normalized_method not in _APPLIED_NATIVE_SPS_RUV_METHODS:
        raise PhosPyInputError(
            "corrected_preprocessing_output declares ambiguous or unsupported "
            f"SPS/RUV-style batch correction method {normalized_method!r}; "
            "applied corrected outputs require a supported native method and "
            "typed BatchCorrectionProvenance"
        )
    if provenance is None:
        raise PhosPyInputError(_MISSING_PROVENANCE_MESSAGE)
    if not isinstance(provenance, BatchCorrectionProvenance):
        raise PhosPyInputError(
            _MISSING_PROVENANCE_MESSAGE
            + "; untyped provenance payloads are not accepted for applied "
            "SPS/RUV-style corrected outputs"
        )
    _validate_complete_sps_ruv_provenance(
        provenance,
        expected_method=normalized_method,
    )


def _validate_complete_sps_ruv_provenance(
    provenance: BatchCorrectionProvenance,
    *,
    expected_method: str,
) -> None:
    requested_method = _normalise_method(provenance.requested_method)
    if requested_method != expected_method:
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance method must "
            f"match the applied correction method; expected {expected_method!r}, "
            f"observed {requested_method!r}"
        )
    if requested_method in _UNSUPPORTED_SPS_RUV_METHODS:
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance declares "
            f"unsupported method {requested_method!r}"
        )
    _require_environment_provenance(provenance)
    selected_site_key_rows = normalize_applied_selected_site_key_rows(
        provenance.selected_site_key_rows
    )
    n_unwanted_factors = _extract_sps_ruv_n_unwanted_factors(provenance)
    _require_selected_control_count_for_unwanted_factors(
        selected_site_key_rows=selected_site_key_rows,
        n_unwanted_factors=n_unwanted_factors,
    )
    require_unique_selected_control_site_rows(selected_site_key_rows)
    require_control_site_source_metadata(
        provenance.control_site_source,
        selected_site_key_rows=selected_site_key_rows,
    )
    require_non_empty_mapping(
        provenance.batch_metadata,
        field_name="BatchCorrectionProvenance.batch_metadata",
    )
    require_non_empty_mapping(
        provenance.design_metadata,
        field_name="BatchCorrectionProvenance.design_metadata",
    )
    require_non_empty_mapping(
        provenance.missing_value_policy,
        field_name="BatchCorrectionProvenance.missing_value_policy",
    )
    if not provenance.observation_masks:
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance must include "
            "observation mask fingerprints for SPS/RUV-style missingness provenance"
        )
    input_matrix_fingerprint = getattr(provenance, "input_matrix_fingerprint", None)
    if input_matrix_fingerprint is None:
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance is missing "
            "input/output matrix fingerprints: input_matrix_fingerprint is required"
        )
    if provenance.output_matrix_fingerprint is None:
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance is missing "
            "input/output matrix fingerprints: output_matrix_fingerprint is required"
        )
    _require_supported_stage_order(provenance.preprocessing_stage_order)
    require_non_empty_mapping(
        provenance.diagnostics,
        field_name="BatchCorrectionProvenance.diagnostics",
    )
    reject_not_provided_required_mapping(
        provenance.resolved_parameters,
        field_name="BatchCorrectionProvenance.resolved_parameters",
    )
    reject_not_provided_required_mapping(
        provenance.batch_metadata,
        field_name="BatchCorrectionProvenance.batch_metadata",
    )
    reject_not_provided_required_mapping(
        provenance.design_metadata,
        field_name="BatchCorrectionProvenance.design_metadata",
    )
    reject_not_provided_required_mapping(
        provenance.missing_value_policy,
        field_name="BatchCorrectionProvenance.missing_value_policy",
    )


def _require_environment_provenance(provenance: BatchCorrectionProvenance) -> None:
    if is_missing_environment_text(provenance.phospy_version):
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance must include "
            "non-empty phospy_version for applied SPS/RUV-style corrected output"
        )
    if is_missing_environment_text(provenance.python_version):
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance must include "
            "non-empty python_version for applied SPS/RUV-style corrected output"
        )
    dependency_versions = provenance.dependency_versions
    if not isinstance(dependency_versions, Mapping) or not dependency_versions:
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance must include "
            "non-empty dependency_versions for applied SPS/RUV-style corrected output"
        )

    missing_dependencies = tuple(
        dependency
        for dependency in BATCH_CORRECTION_ENVIRONMENT_DEPENDENCIES
        if is_missing_environment_text(dependency_versions.get(dependency))
    )
    if missing_dependencies:
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance "
            "dependency_versions must include versions for native "
            "batch-correction dependencies: "
            f"{format_labels(missing_dependencies)}"
        )


def _extract_sps_ruv_n_unwanted_factors(
    provenance: BatchCorrectionProvenance,
) -> int:
    resolved_parameters = provenance.resolved_parameters
    value = resolved_parameters.get("n_unwanted_factors", MISSING)
    if value is MISSING:
        config = resolved_parameters.get("config")
        if isinstance(config, Mapping):
            value = cast(Mapping[str, object], config).get(
                "n_unwanted_factors",
                MISSING,
            )
    if value is MISSING:
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance selected "
            "controls cannot be validated against unwanted-factor count because "
            "resolved_parameters.n_unwanted_factors is missing"
        )
    if isinstance(value, bool) or not isinstance(value, int):
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance selected "
            "controls require a positive integer unwanted-factor count; "
            f"observed n_unwanted_factors={value!r}"
        )
    if value < 1:
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance selected "
            "controls require unwanted-factor count n_unwanted_factors >= 1; "
            f"observed n_unwanted_factors={value}"
        )
    return value


def _require_selected_control_count_for_unwanted_factors(
    *,
    selected_site_key_rows: Sequence[str],
    n_unwanted_factors: int,
) -> None:
    duplicates = duplicates_in_order(selected_site_key_rows)
    selected_count = len(unique_in_order(selected_site_key_rows))
    required_count = n_unwanted_factors + 1
    if selected_count < required_count:
        duplicate_detail = (
            ""
            if not duplicates
            else "; selected_site_key_rows duplicate identifiers: "
            f"{format_labels(duplicates)}"
        )
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance selected "
            "controls are too few for unwanted-factor count; "
            f"unique_selected_controls={selected_count}, "
            f"n_unwanted_factors={n_unwanted_factors}, "
            f"required_selected_controls={required_count}"
            f"{duplicate_detail}"
        )


def _normalise_method(method: object) -> str:
    normalized = str(method).strip().lower()
    if not normalized:
        raise PhosPyInputError(
            "corrected_preprocessing_output applied batch correction method is "
            "missing or empty"
        )
    return normalized


def _is_sps_ruv_style_label(method: str) -> bool:
    return "ruv" in method or method.startswith("sps_") or "_sps_" in method


def _require_supported_stage_order(stage_order: Sequence[str]) -> None:
    normalized = tuple(str(stage).strip() for stage in tuple(stage_order))
    if not normalized:
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance stage order "
            "is missing"
        )
    if normalized != SUPPORTED_INTERNAL_BATCH_CORRECTION_EXECUTED_STAGE_ORDER:
        supported = " -> ".join(
            SUPPORTED_INTERNAL_BATCH_CORRECTION_EXECUTED_STAGE_ORDER
        )
        observed = " -> ".join(normalized)
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance stage order "
            f"is unsupported; observed {observed!r}; supported stage order is "
            f"{supported}"
        )


__all__ = ["validate_applied_native_sps_ruv_correction_provenance"]
