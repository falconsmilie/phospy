"""Preprocessing batch-correction provenance construction."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import pandas as pd

from phospy.provenance import fingerprint_matrix
from phospy.provenance.environment import (
    collect_batch_correction_environment_provenance,
)
from phospy.provenance.models import BatchCorrectionProvenance
from phospy.science.configs.preprocessing.internal_batch_correction import (
    SUPPORTED_INTERNAL_BATCH_CORRECTION_EXECUTED_STAGE_ORDER,
)
from phospy.science.datasets.preprocessing.batch_correction_provenance_payloads import (
    _control_site_source_payload,
    _diagnostics_payload,
    _extract_plan_unwanted_factors,
    _getattr_or,
    _imputation_payload,
    _is_sps_ruv_style_label,
    _json_mapping,
    _metadata_payload,
    _missingness_payload,
    _object_sequence,
    _replicate_metadata_payload,
    _resolve_selected_site_key_rows,
)
from phospy.science.datasets.preprocessing.batch_correction_provenance_validation import (
    normalize_applied_selected_site_key_rows,
    validate_applied_native_sps_ruv_correction_provenance,
)


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


__all__ = [
    "build_native_batch_correction_provenance",
    "normalize_applied_selected_site_key_rows",
    "validate_applied_native_sps_ruv_correction_provenance",
]
