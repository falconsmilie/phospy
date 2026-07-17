"""Preprocessing-domain checks for applied batch-correction provenance."""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from typing import cast

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.provenance import fingerprint_matrix
from phospy.provenance.environment import (
    BATCH_CORRECTION_ENVIRONMENT_DEPENDENCIES,
    collect_batch_correction_environment_provenance,
)
from phospy.provenance.models import (
    BATCH_CORRECTION_SELECTED_SITE_KEY_ROW_SENTINELS,
    BatchCorrectionProvenance,
    JsonValue,
)
from phospy.science.configs.preprocessing.internal_batch_correction import (
    SUPPORTED_INTERNAL_BATCH_CORRECTION_EXECUTED_STAGE_ORDER,
)

_APPLIED_NATIVE_SPS_RUV_METHODS = frozenset({"sps_ruv_style"})
_UNSUPPORTED_SPS_RUV_METHODS = frozenset({"ruv_iii_style"})
_MISSING_PROVENANCE_MESSAGE = (
    "corrected_preprocessing_output with applied native SPS/RUV-style correction "
    "requires typed BatchCorrectionProvenance"
)
_NOT_PROVIDED_VALUES = frozenset({"not_provided", "not provided"})
_MISSING_ENVIRONMENT_VALUES = _NOT_PROVIDED_VALUES | frozenset({"unknown"})
_SELECTED_SITE_KEY_ROW_SENTINELS = (
    BATCH_CORRECTION_SELECTED_SITE_KEY_ROW_SENTINELS | _NOT_PROVIDED_VALUES
)
_STRICT_CONTROL_SOURCE_TYPE_MARKERS = frozenset({"packaged", "reference", "external"})
_CALLER_CONTROL_SOURCE_AUDIT_FIELDS = (
    "organism",
    "identifier_namespace",
    "source_version",
    "license",
    "redistribution",
)
_STRICT_CONTROL_SOURCE_REQUIRED_FIELDS = (
    "organism",
    "identifier_namespace",
    "source_name",
    "source_version",
    "license",
    "redistribution",
    "selection_method",
)
_MISSING = object()


def build_native_batch_correction_provenance(
    *,
    input_matrix: pd.DataFrame,
    output_matrix: pd.DataFrame,
    plan: object,
    report: object,
    metadata: object | None,
    diagnostics: Mapping[str, object] | object,
    warnings: Sequence[object] = (),
    observation_mask: pd.DataFrame | None = None,
    corrected_cell_status: pd.DataFrame | None = None,
    control_site_source: Mapping[str, object] | None = None,
    selected_site_key_rows: Sequence[object] | None = None,
    source: str,
) -> BatchCorrectionProvenance:
    """Build typed provenance for native preprocessing batch correction."""

    environment = collect_batch_correction_environment_provenance()
    method = str(
        _getattr_or(report, "method", _getattr_or(plan, "batch_correction_method", ""))
    ).strip()
    stage_order = tuple(
        str(stage).strip()
        for stage in _object_sequence(
            _getattr_or(
                plan,
                "stage_order",
                SUPPORTED_INTERNAL_BATCH_CORRECTION_EXECUTED_STAGE_ORDER,
            )
        )
    )
    selected_rows = _resolve_selected_site_key_rows(
        explicit=selected_site_key_rows,
        plan=plan,
    )
    observation_masks = ()
    if observation_mask is None:
        method = str(
            _getattr_or(
                report, "method", _getattr_or(plan, "batch_correction_method", "")
            )
        ).strip()
        if method and not _is_sps_ruv_style_label(method):
            observation_mask = input_matrix.notna()
    if observation_mask is not None:
        observation_masks = (
            fingerprint_matrix(
                observation_mask.astype("int8"),
                name="batch_correction.native.observation_mask",
            ),
        )
    if corrected_cell_status is not None:
        observation_masks = (
            *observation_masks,
            fingerprint_matrix(
                corrected_cell_status.astype("string"),
                name="batch_correction.native.corrected_cell_status",
            ),
        )
    control_source = (
        _json_mapping(control_site_source)
        if control_site_source is not None
        else _control_site_source_payload(plan)
    )
    batch_metadata = _metadata_payload(metadata, plan=plan)
    diagnostics_payload = _diagnostics_payload(diagnostics)
    return BatchCorrectionProvenance(
        requested_method=method,
        resolved_parameters=_json_mapping(
            {
                "source": source,
                "batch_column": _getattr_or(
                    plan, "batch_correction_batch_column", None
                ),
                "condition_column": _getattr_or(
                    plan,
                    "batch_correction_condition_column",
                    None,
                ),
                "condition_columns": list(
                    _object_sequence(
                        _getattr_or(plan, "batch_correction_condition_columns", ())
                    )
                ),
                "replicate_column": _getattr_or(
                    plan,
                    "batch_correction_replicate_column",
                    None,
                ),
                "n_unwanted_factors": _extract_plan_unwanted_factors(plan),
            }
        ),
        preprocessing_stage_order=stage_order,
        control_site_source=control_source,
        selected_site_key_rows=selected_rows,
        batch_metadata=batch_metadata,
        replicate_metadata=_replicate_metadata_payload(plan),
        design_metadata=_json_mapping(
            {
                "condition_column": _getattr_or(
                    plan,
                    "batch_correction_condition_column",
                    None,
                ),
                "condition_columns": list(
                    _object_sequence(
                        _getattr_or(plan, "batch_correction_condition_columns", ())
                    )
                ),
            }
        ),
        missing_value_policy=_missingness_payload(plan),
        imputation_policy=_imputation_payload(plan),
        observation_masks=observation_masks,
        input_matrix_fingerprint=fingerprint_matrix(
            input_matrix,
            name="batch_correction.native.input",
        ),
        output_matrix_fingerprint=fingerprint_matrix(
            output_matrix,
            name="batch_correction.native.corrected",
        ),
        diagnostics=diagnostics_payload,
        warnings=tuple(str(warning) for warning in warnings),
        phospy_version=environment.package_version,
        python_version=environment.python_version,
        dependency_versions=environment.dependency_versions,
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


def normalize_applied_selected_site_key_rows(rows: Sequence[object]) -> tuple[str, ...]:
    """Normalize and validate applied selected control row identifiers."""

    if not rows:
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance selected "
            "controls/control provenance must include selected_site_key_rows"
        )

    normalized_rows: list[str] = []
    missing_rows: list[int] = []
    blank_rows: list[int] = []
    sentinel_rows: list[int] = []
    for position, row in enumerate(tuple(rows)):
        if _is_missing_value(row):
            missing_rows.append(position)
            continue
        normalized = str(row).strip()
        if normalized == "":
            blank_rows.append(position)
            continue
        if _is_selected_site_key_row_sentinel(normalized):
            sentinel_rows.append(position)
            continue
        normalized_rows.append(normalized)

    if missing_rows:
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance selected "
            "controls/control provenance selected_site_key_rows contains missing "
            f"site_key rows at positions {_format_positions(missing_rows)}"
        )
    if blank_rows:
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance selected "
            "controls/control provenance selected_site_key_rows contains blank "
            f"site_key rows at positions {_format_positions(blank_rows)}"
        )
    if sentinel_rows:
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance selected "
            "controls/control provenance selected_site_key_rows contains sentinel "
            f"site_key rows at positions {_format_positions(sentinel_rows)}"
        )

    return tuple(normalized_rows)


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
    _require_unique_selected_control_site_rows(selected_site_key_rows)
    _require_control_site_source_metadata(
        provenance.control_site_source,
        selected_site_key_rows=selected_site_key_rows,
    )
    _require_non_empty_mapping(
        provenance.batch_metadata,
        field_name="BatchCorrectionProvenance.batch_metadata",
    )
    _require_non_empty_mapping(
        provenance.design_metadata,
        field_name="BatchCorrectionProvenance.design_metadata",
    )
    _require_non_empty_mapping(
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
    _require_non_empty_mapping(
        provenance.diagnostics,
        field_name="BatchCorrectionProvenance.diagnostics",
    )
    _reject_not_provided_required_mapping(
        provenance.resolved_parameters,
        field_name="BatchCorrectionProvenance.resolved_parameters",
    )
    _reject_not_provided_required_mapping(
        provenance.batch_metadata,
        field_name="BatchCorrectionProvenance.batch_metadata",
    )
    _reject_not_provided_required_mapping(
        provenance.design_metadata,
        field_name="BatchCorrectionProvenance.design_metadata",
    )
    _reject_not_provided_required_mapping(
        provenance.missing_value_policy,
        field_name="BatchCorrectionProvenance.missing_value_policy",
    )


def _require_environment_provenance(provenance: BatchCorrectionProvenance) -> None:
    if _is_missing_environment_text(provenance.phospy_version):
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance must include "
            "non-empty phospy_version for applied SPS/RUV-style corrected output"
        )
    if _is_missing_environment_text(provenance.python_version):
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
        if _is_missing_environment_text(dependency_versions.get(dependency))
    )
    if missing_dependencies:
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance "
            "dependency_versions must include versions for native "
            "batch-correction dependencies: "
            f"{_format_labels(missing_dependencies)}"
        )


def _extract_sps_ruv_n_unwanted_factors(
    provenance: BatchCorrectionProvenance,
) -> int:
    resolved_parameters = provenance.resolved_parameters
    value = resolved_parameters.get("n_unwanted_factors", _MISSING)
    if value is _MISSING:
        config = resolved_parameters.get("config")
        if isinstance(config, Mapping):
            value = cast(Mapping[str, object], config).get(
                "n_unwanted_factors",
                _MISSING,
            )
    if value is _MISSING:
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
    duplicates = _duplicates_in_order(selected_site_key_rows)
    selected_count = len(_unique_in_order(selected_site_key_rows))
    required_count = n_unwanted_factors + 1
    if selected_count < required_count:
        duplicate_detail = (
            ""
            if not duplicates
            else "; selected_site_key_rows duplicate identifiers: "
            f"{_format_labels(duplicates)}"
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


def _is_selected_site_key_row_sentinel(value: object) -> bool:
    return str(value).strip().lower() in _SELECTED_SITE_KEY_ROW_SENTINELS


def _require_unique_selected_control_site_rows(rows: Sequence[str]) -> None:
    duplicates = _duplicates_in_order(rows)
    if duplicates:
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance "
            "selected_site_key_rows contains duplicate selected control row "
            f"identifiers: {_format_labels(duplicates)}"
        )


def _require_control_site_source_metadata(
    source: Mapping[str, object],
    *,
    selected_site_key_rows: Sequence[str],
) -> None:
    _require_non_empty_mapping(
        source,
        field_name="BatchCorrectionProvenance.control_site_source",
    )
    _reject_not_provided_required_mapping(
        source,
        field_name="BatchCorrectionProvenance.control_site_source",
    )
    source_type = _source_type(source)
    if source_type is None or _is_not_provided(source_type):
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance selected "
            "controls/control provenance must include control source metadata"
        )

    if _has_strict_control_source_type(source):
        missing = tuple(
            field_name
            for field_name in _STRICT_CONTROL_SOURCE_REQUIRED_FIELDS
            if not _has_non_missing_text(source.get(field_name))
        )
        if missing:
            raise PhosPyInputError(
                "corrected_preprocessing_output BatchCorrectionProvenance "
                "packaged/reference/external control-source metadata is "
                f"incomplete; missing {_format_labels(missing)}"
            )
        return

    missing_without_reason = tuple(
        field_name
        for field_name in _CALLER_CONTROL_SOURCE_AUDIT_FIELDS
        if not _has_non_missing_text(source.get(field_name))
        and not _has_metadata_missing_reason(
            source,
            field_name,
            selected_site_key_rows=selected_site_key_rows,
        )
    )
    if missing_without_reason:
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance control "
            "source is missing "
            f"{_format_labels(missing_without_reason)} without explicit rationale"
        )

    has_source_name = _has_non_missing_text(source.get("source_name"))
    has_source_version = _has_non_missing_text(source.get("source_version"))
    has_unavailable_reason = _has_non_missing_text(
        source.get("source_version_unavailable_reason")
    )
    has_missing_reason = _has_metadata_missing_reason(
        source,
        "source_version",
        selected_site_key_rows=selected_site_key_rows,
    )
    if has_source_name and not (
        has_source_version or has_unavailable_reason or has_missing_reason
    ):
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance control "
            "source declares source_name without source_version or explicit "
            "source_version_unavailable_reason"
        )
    if source_type == "caller_supplied" and not (
        has_source_version or has_unavailable_reason or has_missing_reason
    ):
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance "
            "caller_supplied control source must record source_version or "
            "source_version_unavailable_reason"
        )


def _source_type(source: Mapping[str, object]) -> str | None:
    for key in ("source_type", "source"):
        value = source.get(key)
        if _has_non_missing_text(value):
            return str(value).strip().lower()
    return None


def _has_strict_control_source_type(source: Mapping[str, object]) -> bool:
    for key in (
        "source_type",
        "source",
        "control_site_set_source_type",
        "source_name",
    ):
        value = source.get(key)
        if _has_non_missing_text(value) and _is_strict_control_source_type(
            str(value).strip().lower()
        ):
            return True
    return False


def _is_strict_control_source_type(source_type: str | None) -> bool:
    if source_type is None:
        return False
    tokens = frozenset(source_type.replace("-", "_").split("_"))
    return bool(tokens & _STRICT_CONTROL_SOURCE_TYPE_MARKERS)


def _has_metadata_missing_reason(
    source: Mapping[str, object],
    field_name: str,
    *,
    selected_site_key_rows: Sequence[str],
) -> bool:
    reasons = source.get("metadata_missing_reason")
    if isinstance(reasons, Mapping) and _has_non_missing_text(
        cast(Mapping[str, object], reasons).get(field_name)
    ):
        return True
    if _has_non_missing_text(source.get(f"{field_name}_missing_reason")):
        return True
    if field_name == "source_version" and _has_non_missing_text(
        source.get("source_version_unavailable_reason")
    ):
        return True
    return _has_metadata_missing_reason_by_site_key(
        source,
        field_name,
        selected_site_key_rows=selected_site_key_rows,
    )


def _has_metadata_missing_reason_by_site_key(
    source: Mapping[str, object],
    field_name: str,
    *,
    selected_site_key_rows: Sequence[str],
) -> bool:
    reasons_by_site_key = source.get("metadata_missing_reason_by_site_key")
    if not isinstance(reasons_by_site_key, Mapping):
        return False
    selected = tuple(str(site_key) for site_key in selected_site_key_rows)
    if not selected:
        return False
    by_site_key = cast(Mapping[str, object], reasons_by_site_key)
    for site_key in selected:
        site_reasons = by_site_key.get(site_key)
        if not isinstance(site_reasons, Mapping):
            return False
        if not _has_non_missing_text(
            cast(Mapping[str, object], site_reasons).get(field_name)
        ):
            return False
    return True


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


def _require_non_empty_mapping(
    value: Mapping[str, object],
    *,
    field_name: str,
) -> None:
    if not isinstance(value, Mapping) or not value:
        raise PhosPyInputError(f"{field_name} must be a non-empty object")


def _reject_not_provided_required_mapping(
    value: Mapping[str, object],
    *,
    field_name: str,
) -> None:
    for key, item in value.items():
        if _is_not_provided(item):
            raise PhosPyInputError(
                f"{field_name}.{key} must not be recorded as not_provided for "
                "applied SPS/RUV-style corrected output"
            )


def _has_non_missing_text(value: object) -> bool:
    return not _is_missing_required_text(value) and not _is_not_provided(value)


def _is_missing_environment_text(value: object) -> bool:
    return (
        _is_missing_required_text(value)
        or str(value).strip().lower() in _MISSING_ENVIRONMENT_VALUES
    )


def _is_missing_required_text(value: object) -> bool:
    return value is None or str(value).strip() == ""


def _is_not_provided(value: object) -> bool:
    return str(value).strip().lower() in _NOT_PROVIDED_VALUES


def _levels_in_order(labels: Sequence[str]) -> tuple[str, ...]:
    levels: list[str] = []
    seen: set[str] = set()
    for label in labels:
        if label in seen:
            continue
        seen.add(label)
        levels.append(label)
    return tuple(levels)


def _unique_in_order(labels: Sequence[str]) -> tuple[str, ...]:
    return _levels_in_order(labels)


def _duplicates_in_order(labels: Sequence[str]) -> tuple[str, ...]:
    counts: dict[str, int] = {}
    for label in labels:
        counts[label] = counts.get(label, 0) + 1
    return tuple(label for label in _levels_in_order(labels) if counts[label] > 1)


def _is_missing_value(value: object) -> bool:
    return bool(pd.Series((value,), dtype="object").isna().iat[0])


def _format_positions(positions: Sequence[int]) -> str:
    preview = ", ".join(str(position) for position in tuple(positions)[:8])
    suffix = "" if len(positions) <= 8 else ", ..."
    return f"[{preview}{suffix}]"


def _format_labels(labels: Sequence[str]) -> str:
    preview = ", ".join(repr(value) for value in tuple(labels)[:5])
    suffix = "" if len(labels) <= 5 else " ..."
    return f"{preview}{suffix}"


__all__ = [
    "build_native_batch_correction_provenance",
    "normalize_applied_selected_site_key_rows",
    "validate_applied_native_sps_ruv_correction_provenance",
]


def _resolve_selected_site_key_rows(
    *,
    explicit: Sequence[object] | None,
    plan: object,
) -> tuple[str, ...]:
    if explicit is not None:
        return tuple(str(row).strip() for row in explicit if str(row).strip())
    control_site_set = _getattr_or(plan, "batch_correction_control_site_set", None)
    rows = _getattr_or(control_site_set, "selected_site_key_rows", None)
    if rows is not None:
        return tuple(
            str(row).strip() for row in _object_sequence(rows) if str(row).strip()
        )
    controls = _getattr_or(control_site_set, "controls", None)
    if controls is not None:
        resolved: list[str] = []
        for item in _object_sequence(controls):
            site_key = _getattr_or(item, "site_key", item)
            if str(site_key).strip():
                resolved.append(str(site_key).strip())
        return tuple(resolved)
    return ()


def _control_site_source_payload(plan: object) -> Mapping[str, JsonValue]:
    control_site_set = _getattr_or(plan, "batch_correction_control_site_set", None)
    if control_site_set is None:
        method = str(_getattr_or(plan, "batch_correction_method", "")).strip()
        if method and not _is_sps_ruv_style_label(method):
            return _json_mapping(
                {
                    "source_type": "not_applicable",
                    "reason": (
                        "control-site selection is not used by this native "
                        "batch-correction method"
                    ),
                }
            )
        return _json_mapping(
            {
                "source_type": "not_provided",
                "reason": "batch-correction plan did not carry a control-site set",
            }
        )
    if hasattr(control_site_set, "to_payload"):
        payload = _payload(control_site_set)
        if isinstance(payload, Mapping):
            return _json_mapping(payload)
    return _json_mapping(
        {
            "source_type": _getattr_or(control_site_set, "source_type", "not_provided"),
            "source_name": _getattr_or(control_site_set, "source_name", "not_provided"),
            "source_version": _getattr_or(
                control_site_set,
                "source_version",
                "not_provided",
            ),
            "license": _getattr_or(control_site_set, "license", "not_provided"),
            "redistribution": _getattr_or(
                control_site_set,
                "redistribution_status",
                "not_provided",
            ),
            "selection_method": _getattr_or(
                control_site_set,
                "selection_method",
                "not_provided",
            ),
        }
    )


def _metadata_payload(
    metadata: object | None,
    *,
    plan: object,
) -> Mapping[str, JsonValue]:
    if metadata is None:
        return _json_mapping({"source": "not_provided"})
    batch_column = _getattr_or(plan, "batch_correction_batch_column", None)
    return _json_mapping(
        {
            "column": batch_column,
            "sample_order": list(
                _object_sequence(_getattr_or(metadata, "sample_order", ()))
            ),
            "batch_by_sample": _object_mapping(
                _getattr_or(metadata, "batch_by_sample", {})
            ),
            "condition_by_sample": _object_mapping(
                _getattr_or(metadata, "condition_by_sample", {})
            ),
        }
    )


def _missingness_payload(plan: object) -> Mapping[str, JsonValue]:
    policy = _getattr_or(plan, "batch_correction_missingness_policy", None)
    if policy is None:
        method = str(_getattr_or(plan, "batch_correction_method", "")).strip()
        if method and not _is_sps_ruv_style_label(method):
            return _json_mapping({"policy": "reject_missing_at_batch_correction"})
        return _json_mapping({"policy": "not_provided"})
    return _json_mapping(_payload(policy))


def _imputation_payload(plan: object) -> Mapping[str, JsonValue]:
    policy = _getattr_or(plan, "batch_correction_imputation_policy", None)
    if policy is None:
        method = str(_getattr_or(plan, "batch_correction_method", "")).strip()
        if method and not _is_sps_ruv_style_label(method):
            return _json_mapping({"policy": "forbid"})
        return _json_mapping({})
    return _json_mapping(_payload(policy))


def _replicate_metadata_payload(plan: object) -> Mapping[str, JsonValue] | None:
    replicate_column = _getattr_or(plan, "batch_correction_replicate_column", None)
    if replicate_column is None:
        return None
    return _json_mapping({"column": str(replicate_column).strip()})


def _diagnostics_payload(
    diagnostics: Mapping[str, object] | object,
) -> Mapping[str, JsonValue]:
    payload = _json_mapping(_payload(diagnostics))
    if "executor" in payload or "stage_diagnostics" in payload:
        return payload
    return _json_mapping({**payload, "stage_diagnostics": dict(payload)})


def _extract_plan_unwanted_factors(plan: object) -> object:
    request = _getattr_or(plan, "batch_correction_internal_request", None)
    config = _getattr_or(request, "config", request)
    return _getattr_or(config, "n_unwanted_factors", None)


def _payload(value: object) -> object:
    to_payload = getattr(value, "to_payload", None)
    if callable(to_payload):
        return cast(Callable[[], object], to_payload)()
    if isinstance(value, Mapping):
        return dict(value)
    return value


def _json_mapping(value: Mapping[str, object] | object) -> Mapping[str, JsonValue]:
    if not isinstance(value, Mapping):
        return {"value": _json_value(value)}
    return {str(key): _json_value(item) for key, item in value.items()}


def _json_value(value: object) -> JsonValue:
    if value is None or isinstance(value, (str, int, float, bool)):
        return cast(JsonValue, value)
    if isinstance(value, Mapping):
        return cast(
            JsonValue,
            {str(key): _json_value(item) for key, item in value.items()},
        )
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return cast(JsonValue, [_json_value(item) for item in value])
    enum_value = getattr(value, "value", _MISSING)
    if enum_value is not _MISSING:
        return _json_value(enum_value)
    to_payload = getattr(value, "to_payload", None)
    if callable(to_payload):
        return _json_value(cast(Callable[[], object], to_payload)())
    return str(value)


def _getattr_or(value: object, name: str, default: object) -> object:
    return getattr(value, name, default)


def _object_sequence(value: object) -> tuple[object, ...]:
    if value is None or isinstance(value, (str, bytes, bytearray)):
        return ()
    if isinstance(value, Sequence):
        return tuple(value)
    return ()


def _object_mapping(value: object) -> dict[str, object]:
    if not isinstance(value, Mapping):
        return {}
    return {str(key): item for key, item in value.items()}
