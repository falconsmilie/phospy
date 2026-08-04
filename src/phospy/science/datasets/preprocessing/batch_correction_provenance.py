"""Preprocessing provenance construction and legacy validation import routes."""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from typing import Any, cast

import pandas as pd

from phospy.provenance import fingerprint_matrix
from phospy.provenance.environment import (
    collect_batch_correction_environment_provenance,
)
from phospy.provenance.models import (
    BatchCorrectionProvenance,
    JsonValue,
)
from phospy.science.configs.preprocessing.internal_batch_correction import (
    SUPPORTED_INTERNAL_BATCH_CORRECTION_EXECUTED_STAGE_ORDER,
)

_VALIDATION_DATASET_MODULE_PARTS = ("phospy", "validation", "datasets")
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


def _validation_provenance_module() -> object:
    # Dynamic import avoids a static science->validation package cycle while keeping
    # applied-output validation implementation in validation.datasets.
    return __import__(
        ".".join((*_VALIDATION_DATASET_MODULE_PARTS, "batch_correction_provenance")),
        fromlist=("validate_applied_native_sps_ruv_correction_provenance",),
    )


def _validation_controls_module() -> object:
    # Dynamic import avoids a static science->validation package cycle while keeping
    # selected-control normalization implementation in validation.datasets.
    return __import__(
        ".".join((*_VALIDATION_DATASET_MODULE_PARTS, "batch_correction_controls")),
        fromlist=("normalize_applied_selected_site_key_rows",),
    )


def _is_sps_ruv_style_label(method: str) -> bool:
    return "ruv" in method or method.startswith("sps_") or "_sps_" in method


validate_applied_native_sps_ruv_correction_provenance = cast(
    Callable[..., None],
    cast(
        Any,
        _validation_provenance_module(),
    ).validate_applied_native_sps_ruv_correction_provenance,
)
normalize_applied_selected_site_key_rows = cast(
    Callable[[Sequence[object]], tuple[str, ...]],
    cast(
        Any,
        _validation_controls_module(),
    ).normalize_applied_selected_site_key_rows,
)

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
